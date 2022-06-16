function [x1, x2, A, V] = epa(correspondences, K)
    % Depending on whether a calibrating matrix 'K' is given,
    % this function calculates either the essential or the fundamental matrix
    % with the eight-point algorithm.
    x1 = correspondences(1:2, :); % euclidean
    x2 = correspondences(3:4, :);
    x1 = [x1; ones(size(correspondences,2))]; %homogenous
    x2 = [x2; ones(size(correspondences,2))];
    
    if exist('K', 'var')
        x1 = inv(K)*x1;
        x2 = inv(K)*x2;
    end

    
end