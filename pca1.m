function [denoised]=pca1(P, C)
    % Inputs:
    %   P : noisy patch matrix
    %   C : reduced number of dimensions for PCA
    % Outputs:
    %   denoised : denoised patch matrix


    meanpatch = mean(P,2);
    P_wrt_mean = P - meanpatch;                                             % Patches with respect to mean value

    standard_dev = sqrt(sum(P_wrt_mean.^2, 2) / size(P_wrt_mean, 2)); 
    standard_dev_inv = 1 ./ standard_dev;
    normalised = P_wrt_mean .* standard_dev_inv;                            % Normalizing each patch

    [Vals, obt_score, ~] = pca(normalised');                                % Applying PCA

    dim_red = obt_score(:,1:C) * ((Vals(:,1:C))');                          % Patch Data dim_red to C dimensions

    P_wrt_mean_op = dim_red' .* standard_dev;                               %Multiplying with Std Dev

    denoised = P_wrt_mean_op + meanpatch;                                   % Adding mean val to get orignal data                      
end