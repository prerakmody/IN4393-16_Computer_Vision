% Point cloud stitching

reference = pointCloud(S{1, 1}');
current = pointCloud(S{1, 2}');

fixed = reference;
moving = current;

tform = pcregistericp(moving, fixed, 'Metric', 'PointToPlane', 'Extrapolate', true);
aligned = pctransform(current, tform);

merge_size = 0.015;
scene = pcmerge(reference, aligned, merge_size);

accum_tform = tform;

pcshow(scene)

for i = 3:4
    
    if size(S{1, i}, 2) > 5
        
        disp(i)
        
        current = pointCloud(S{1, i}');

        fixed = moving;
        moving = current;

        tform = pcregistericp(moving, fixed, 'Metric', 'PointToPlane', 'Extrapolate', true);
        accum_tform = affine3d(tform.T * accum_tform.T);

        aligned = pctransform(current, accum_tform);

        scene = pcmerge(scene, aligned, merge_size);
        
    end
    
end

% pcshow(scene, 'VerticalAxis','Y', 'VerticalAxisDir', 'Down')
% title('Updated world scene')
% xlabel('X (m)')
% ylabel('Y (m)')
% zlabel('Z (m)')
