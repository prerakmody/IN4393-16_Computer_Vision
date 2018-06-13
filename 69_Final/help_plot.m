imshow(images(:, :, :, 1))
hold on
scatter(points(1, 1:6), points(2, 1:6), 'c')

figure
imshow(images(:, :, :, 2))
hold on
scatter(points(3, 1:6), points(4, 1:6), 'r')

figure
imshow(images(:, :, :, 3))
hold on
scatter(points(5, 1:6), points(6, 1:6), 'b')

figure
imshow(images(:, :, :, 4))
hold on
scatter(points(7, 1:6), points(8, 1:6), 'k')

figure;
scatter3(S{1, 1}(1, 1:6), S{1, 1}(2, 1:6), S{1, 1}(3, 1:6))

% scatter3(S(1, 1:4), S(2, 1:4), S(3, 1:4))