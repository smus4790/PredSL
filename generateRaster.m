function output_k_array = generateRaster(best_threshold_snow, best_threshold_lake, dist_k, count, data_k, output_k_array)
    % Assuming you have code to generate the raster using the best thresholds
    % Modify this function based on your specific raster generation logic
    
    % Update the thresholding logic
    output_k_array(:, :, count + 1) = (dist_k > best_threshold_snow) | ((data_k == 2) & (dist_k > best_threshold_lake));
    
    % Display the images using subplots
    subplot(1, 2, 1);
    imshow(output_k_array(:, :, count + 1));
    title(['Predicted - Timestamp ' num2str(count)]);
    
    subplot(1, 2, 2);
    imshow(data_k);
    title(['Actual - Timestamp ' num2str(count)]);
    
    % Pause to allow time for display (optional)
    pause(0.5);
end