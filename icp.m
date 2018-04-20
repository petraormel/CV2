function [R, t] = icp(source, target)
    R = eye(3);

    t = 0;

    % Threshold for change in RMS values    
    threshold = 1e-5; 
    
    % Initialization
    prev_RMS = inf ;
    RMS = 0 ; 
    count = 0;
    diff_RMS = inf ;
    
    fprintf("Performing ICP...\n")
    
%     while diff_RMS > threshold
        for i =1:100
        % Find closest points using brute-force
        closest_idx = find_closest_points(source, target);
        
        % Instantiate corresponding points in target set
        target_corr = target(:, closest_idx) ;
        
        % Refine R and t using SVD
        [R, t] = refine_R_t(source, target_corr) ;
    
        % Transform source/base points using R, t
        source_t = zeros(size(source)) ; 
        for i = 1:size(source, 2)
            source_t(:, i) = R * source(:, i) ;
        end
        source = source_t + t ;
        
        % Calculate RMS to determine convergence
        RMS = calculate_RMS(source, target_corr) ;
        
        prev_RMS = RMS ;
        diff_RMS = abs(RMS - prev_RMS) ; 
        
        count = count + 1 ;
        disp("RMS at iteration : %f\n", RMS) ;
    end
    
    fprintf("Completed ICP in %d iterations.\n", count) ;
end

function mins = find_closest_points(source, target)
% Finds indices of closest points in target set corresponding points in
% source set.
% INPUTS:    
%       source : (D x N) array with D dimensionality and N num data points
%       target : (D x N) array with D dimensionality and N num data points
% OUTPUTS:
%       mins   : (N x 1) array containing indices of closest points in
%                       target array

    % Make sure input dimensions are correct
    if size(source,1) > size(source, 2)
        source = source.';
        target = target.';
        fprintf("Transposed source and target matrices for compatability")
    end

    mins = zeros(size(source,2),1);
    minj = 0;
    
    for i=1:size(source,2)
        mindist = inf;
        for j=1:size(target,2)
            dist = norm(source(:,i)-target(:,j));
            if dist < mindist
                
                mindist = dist;
                minj = j;
               
            end
        end
        mins(i) = minj;    
    end
    
end
    

function [R, t] = refine_R_t(source, target)
    
    N = size(source, 2) ;
    D = size(source, 1) ;
   
    % Equal weights for all points
    weights = ones(1, size(source, 2)) ;
    
    % Initialize weighted centroid values
    wc_source = zeros(D, 1) ;
    wc_target = zeros(D, 1) ;
    
    % Find weighted centroids of source and target set
    for i=1:N
        wc_source = wc_source + weights(i) * source(:, i) ;
        wc_target = wc_target + weights(i) * target(:, i) ;
    end
    wc_source = wc_source / sum(weights) ; 
    wc_target = wc_target / sum(weights) ;
    
    % Check dimensions
    if not(size(wc_source, 1) == size(source, 1))
        disp("Weighted centroid vectors do not have correct dimensions\n")
    end
    
    % Calculate centered vectors
    source_c = source - wc_source ;
    target_c = target - wc_target ;
    
    % Check dimensions
    if not(size(source_c) == size(source))
        disp("Centered vectors do not have correct dimensions\n")
    end
    
    % Calculate covariance matrix
    cov_matrix = source_c * diag(weights) * target_c.' ;
    
    % Perform SVD on covariance matrix
    [U,S,V] = svd(cov_matrix) ; 
    
    % Calculate optimal R matrix
    diag_matrix = eye(size(S)) ; 
    diag_matrix(end, end) = det(V*U.') ;
    R = V * diag_matrix * U.' ;
    
    % Check dimensions
    if not(isequal(size(R), [D, D]))
        disp("Rotation matrix does not have the correct dimensions: ")
        disp(size(R))
    end
    
    % Calculate optimal t vector
    t = target_c - R * source_c ;
    
    disp(size(t))
    % Check dimensions
    if not(isequal(size(t), [D, N]))
        disp("Translation vector does not have the correct dimensions: ")
        disp(size(t))
    end
end

function RMS = calculate_RMS(source, target_corr)
    
    % Calculate RMS
    RMS = 1/size(source, 2) * (source - target_corr).^2 ; 
end