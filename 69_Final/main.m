% IN4393-16: Final Assignment
% 
% Jesse Hagenaars & Michiel Mollema - 04.06.2018

clc; clear; close all

run('C:/Users/jesse/Documents/MATLAB/vlfeat/toolbox/vl_setup')
% run('/home/michiel/Programs/MATLAB/vlfeat/toolbox/vl_setup')

%% Settings
match_threshold = 3;
dist_threshold = 50;


%% Load data

% Images
image_files = dir('model_castle/*.JPG');

for i = 1:length(image_files)
    
    current_image = imread([image_files(i).folder '/' image_files(i).name]);
    
    if i == 1
        images = uint8(zeros([size(current_image) length(image_files)]));
    end
    
    images(:, :, :, i) = uint8(current_image);
   
end

% Convert to greyscale
images_grey = uint8(mean(images, 3));

% SIFT features, obtained using Harris-affine & Hessian-affine SIFT from
%   http://www.robots.ox.ac.uk/~vgg/research/affine/detectors.html
feature_files = dir('features/*.png.harhes.sift');
sift_oxford = {};

for i = 1:length(feature_files)
    
    current_features = dlmread([feature_files(i).folder '/' feature_files(i).name], ' ', 2, 0);
    
    sift_oxford{1, i} = current_features(:, 1:5)';
    sift_oxford{2, i} = current_features(:, 6:end)';
   
end

 
%% Feature point detection & extraction of SIFT descriptors (4 pts)

% Only do if file doesn't exist already
if ~exist('sift_vlfeat.mat', 'file')
  
    % Using VLFeat's SIFT
    sift_vlfeat = {};

    for i = 1:size(images, 4)

        [sift_vlfeat{1, i}, sift_vlfeat{2, i}] = vl_sift(single(images_grey(:, :, i)));

    end
    
    % Save
    save('sift_vlfeat', 'sift_vlfeat')
    
else
    
    % Load already detected features
    load('sift_vlfeat', 'sift_vlfeat')
    
end

% Using Harris-affine & Hessian-affine SIFT from
%   http://www.robots.ox.ac.uk/~vgg/research/affine/detectors.html
%   --> loaded in above section


%%  Normalized 8-point RANSAC to find best matches (4 pts)

% Only do if file doesn't exist already
if ~exist('matches_8pt_RANSAC.mat', 'file')
    
    % Cell array of matches per frame pair
    matches_8pt_RANSAC = {};

    for i = 1:size(sift_vlfeat, 2)  

        disp(['i: ' num2str(i)])

        [F_ransac_denorm, inliers_1, inliers_2, inliers_idx] = do_eightpoint(sift_vlfeat, match_threshold, dist_threshold, i);

        matches_8pt_RANSAC{1, i} = inliers_idx;

        if i == 1

            plot_eightpoint(inliers_1, inliers_2, F_ransac_denorm);

        end

    end

    save('matches_8pt_RANSAC', 'matches_8pt_RANSAC')
    
else
    
    load('matches_8pt_RANSAC', 'matches_8pt_RANSAC')
    
end

%% Chaining (8 pts)

point_view_matrix = chaining(matches_8pt_RANSAC);


%% Stitching (12 pts)

% Cells to store 3D point set for each set of frames
S = {};
M = {};
colours = {};
tform = {};

% Use 4 consecutive frames each time
for f = 0:size(point_view_matrix, 1) - 1
    
    % Shift (cell) array circularly
    pv_matrix_circ = circshift(point_view_matrix, -f, 1);
    sift_vlfeat_circ = circshift(sift_vlfeat, -f, 2);
    
    pv_matrix_circ_next = circshift(point_view_matrix, -f-1, 1);
    sift_vlfeat_circ_next = circshift(sift_vlfeat, -f-1, 2);
    
    % Get x, y for each SIFT descriptor
    points = get_points(sift_vlfeat_circ(1, 1:3), pv_matrix_circ(1:3, :));
    
    points_next = get_points(sift_vlfeat_circ_next(1, 1:3), pv_matrix_circ_next(1:3, :));
    
    points_plot = get_points(sift_vlfeat_circ(1, 1:2), pv_matrix_circ(1:2, :));
    
    showMatchedFeatures(images(:, :, :, f+1), images(:, :, :, f+2), points_plot(1:2, :)', points_plot(3:4, :)')
    
    [~, tform{2, f+1}, tform{1, f+1}] = intersect(points_next(1:end-2, :)', points(3:end, :)', 'rows');
    
    colour = [images(sub2ind(size(images), uint16(points(2, :)), uint16(points(1, :)), ones([1, size(points, 2)]), f+1 * ones([1, size(points, 2)]))); ...
              images(sub2ind(size(images), uint16(points(2, :)), uint16(points(1, :)), 2*ones([1, size(points, 2)]), f+1 * ones([1, size(points, 2)]))); ...
              images(sub2ind(size(images), uint16(points(2, :)), uint16(points(1, :)), 3*ones([1, size(points, 2)]), f+1 * ones([1, size(points, 2)])))];
    
    if size(points, 2) > 2
    
        % Perform structure-from-motion and solve for affine ambiguity
        [S{1, f+1}, M{1, f+1}] = SfM_affine(points);
        colours{1, f+1} = colour;
        
    end
    
end

% Extend this with transformation
S_final = S{1,1};
M_final = [];

counter = 0;

final_colors = colours{1, 1};

scene = pointCloud(S{1, 1}', 'Color', uint8(colours{1, 1}'));
merge_size = 1;

% Stitch various 3D point sets together
for s = 0:size(S, 2) - 1
    
    % Shift cell array circularly
    S_circ = circshift(S, -s, 2);
    M_circ = circshift(M, -s, 2);
    Col_circ = circshift(colours, -s, 2);
    
    % Minimum number of rows (equal rows needed for procrustes)
    min_cols_S = min(cellfun('size', S_circ(1, 1:2), 2));
    min_rows_M = min(cellfun('size', M_circ(1, 1:2), 1));
    
    max_cols_S = max(cellfun('size', S_circ(1, 1:2), 2));
    
    disp(max_cols_S - min_cols_S)
    
    counter = counter + (max_cols_S - min_cols_S);
    
    if min_cols_S > 3 && min_rows_M > 0 % --> cells next to zero also deleted
        
        final_colors = [final_colors Col_circ{1, 2}];
        
        % Get transformation between 3D point sets
        % Transpose since x, y, z have to be columns
        [~, ~, tr] = procrustes(S_circ{1, 1}(:, tform{1,s+1})', S_circ{1, 2}(:, tform{2,s+1})');
        [~, X] = procrustes(M_circ{1, 1}(1:min_rows_M, :), M_circ{1, 2}(1:min_rows_M, :));
        
        % All rows equal
        tr.c = tr.c(1, :);
        
        Z = tr.b * S_circ{1, 2}' * tr.T + tr.c;
        
%         figure; scatter3(S{1, 1}(1, :), S{1, 1}(2, :), S{1, 1}(3, :), 'b')
%         hold on; scatter3(S{1, 1}(1, tform{1, 1}), S{1, 1}(2, tform{1, 1}), S{1, 1}(3, tform{1, 1}), 'g')
%         title('Set 1 and used points')
%         figure; scatter3(S{1, 2}(1, :), S{1, 2}(2, :), S{1, 2}(3, :), 'r')
%         hold on; scatter3(S{1, 2}(1, tform{2, 1}), S{1, 2}(2, tform{2, 1}), S{1, 2}(3, tform{2, 1}), 'g')
%         title('Set 2 and used points')
%         figure; scatter3(S{1, 1}(1, :), S{1, 1}(2, :), S{1, 1}(3, :), 'b')
%         hold on; scatter3(Z(:, 1), Z(:, 2), Z(:, 3), 'r')
%         title('Set 1 and set 2 transformed')
%         
%         figure;
%         title('Used points from set 1 and their transform')
%         scatter3(S{1, 1}(1, tform{1, 1}), S{1, 1}(2, tform{1, 1}), S{1, 1}(3, tform{1, 1}), 'b')
%         hold on; scatter3(Z(tform{2, 1}, 1), Z(tform{2, 1}, 2), Z(tform{2, 1}, 3), 'r')
        
        new_cloud = pointCloud(Z, 'Color', Col_circ{1, 2}');

        % Extend final shape matrix --> check dimensions
        S_final = [S_final Z'];
        
        % Extend final motion matrix
        M_final = [M_final; X];
        
        scene = pcmerge(scene, new_cloud, merge_size);
        
    end
    
end

% Check for close points in some way?


%% Bundle adjustment (4 pts)




%% Eliminate affine ambiguity (4 pts)

% S_final = solve_affine_ambiguity(S_final, M_final);


%% 3D model plotting (4 pts)

% scene = pointCloud(S_final', 'Color', uint8(final_colors'));
pcshow(scene, 'MarkerSize', 100)

