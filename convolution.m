function output_post_conv = convolution(output_pre_convolution, dimdata)
    % Extract the raster dimensions from input
    [r, c, ~] = size(output_pre_convolution);
    
    % Initialize outputs dynamically
    output_post_conv = zeros(r, c, dimdata);
    nhood = zeros(r, c, dimdata);

    for i = 1:dimdata
        nhood(:, :, i) = output_pre_convolution(:, :, i);
    end

    % Define the convolution kernel
    k_size = 3;
    kernel = ones(k_size, k_size);
    kernel(ceil(k_size / 2), ceil(k_size / 2)) = 0;

    % Apply convolution and thresholding
    for idx = 1:size(nhood, 3)
        currentNhood = nhood(:, :, idx);
        convolutionResult = convn(currentNhood, kernel, 'same');
        
        % Threshold value
        thresholdValue = 4;
        currentNhood(convolutionResult > thresholdValue) = 1;
        currentNhood(convolutionResult <= thresholdValue) = 0;
        
        output_post_conv(:, :, idx) = currentNhood;
    end
end
