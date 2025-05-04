function output_post_conv = convolution(output_pre_convolution, dimdata)
output_post_conv = zeros(4279,6826,dimdata);
nhood = zeros(4279,6826,dimdata);
for i=1:dimdata
    nhood(:,:,i) = output_pre_convolution(:,:,i);
end

k_size = 3;
kernel = ones(k_size, k_size);
kernel((k_size + 1) / 2, (k_size + 1) / 2) = 0;

% Iterate over nhood matrices
for idx = 1:size(nhood, 3)
    currentNhood = nhood(:, :, idx);

    % Perform convolution
    convolutionResult = convn(currentNhood, kernel, 'same');

    % Thresholding based on the convolution result
    thresholdValue = 4;
    currentNhood(convolutionResult > thresholdValue) = 1;
    currentNhood(convolutionResult <= thresholdValue) = 0;

    % Update the nhood matrix
    output_post_conv(:, :, idx) = currentNhood;
end
end
  