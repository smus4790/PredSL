function [sum_ssim] = thresholding2(threshold_snow, threshold_lake, dist_k, count, data_k, data_mod, mat_k)
    %threshold_lake = 0.96;
    
    % Update the thresholding logic
    mat_k(:, :, count + 1) = (dist_k > threshold_snow) | ((data_k == 2) & (dist_k > threshold_lake));
    
    % Display the image
    imshow(mat_k(:, :, count + 1));

    % Calculate SSIM
    ssim_accuracy(count + 1) = ssim(1 - mat_k(:, :, count + 1), data_mod(:, :, count + 1));

    % Display SSIM accuracy for the current count
    disp(['SSIM accuracy for count ' num2str(count) ': ' num2str(1-ssim_accuracy(count + 1))]);

    % Sum of SSIM accuracies
    sum_ssim = sum(ssim_accuracy);
end