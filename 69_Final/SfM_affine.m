function [S_hat, M_hat] = SfM_affine(points)
% SFM_AFFINE Performs affine structure-from-motion.
%
% Inputs:
% - points: matrix of x & y coordinates, shape (2 * N_frames, N_points),
%   rows of x & y coordinates alternate (so 1st 2 rows belong to 1st frame)
%
% Outputs:
% - S_hat: estimated shape matrix (with affine ambiguity)
% - M_hat: estimated motion matrix (with affine ambiguity)
%
% Jesse Hagenaars & Michiel Mollema - 11.06.2018

% Sizes
[N_frames, N_points] = size(points);
N_frames = N_frames / 2;

% Center points
points_center = points - repmat(sum(points, 2) / N_points, 1, N_points);

% SVD
[U, W, V] = svd(points_center);

% Only top 3 singular values
U_acc = U(:, 1:3);
W_acc = W(1:3, 1:3);
V_acc_T = V(:, 1:3)';

% Decompose into measurements M and shape S
M_hat = U_acc * sqrt(W_acc);
S_hat = sqrt(W_acc) * V_acc_T;

% Denormalize -> how?
% S_hat = S_hat + repmat(sum(points, 2) / N_points, 1, N_points)

end