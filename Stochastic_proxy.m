%function Stochastic Cellular Automata
rho = xin(1);
alpha = xin(2);
beta = xin(3);
gamma = xin(4);
pp = xin(5);
qq = xin(6);
rr = xin(7);

global mask mask2 incidence elevation n m cells nbors data prop frac means ffunc ssim_accuracy

[r, c, ~] = size(data); % Dynamically get raster dimensions

ssim_accuracy = zeros(1, dimdata);
s = ones(m, n) .* mask;
[x, y, z] = find(s == 1);

b = zeros(m, n) ./ nbors; 
output = zeros(m, n, length(prop));

%-----Automated assignment of values in governing equation----%
steps = 0; 
datastep = 1; 
snowcells = cells; 
measure = ones(2 * length(prop), 1); 
carryon = true;

% Initialize variables using dynamic sizes
scaling = zeros(1, dimdata);
fixed = cell(1, dimdata);
dist_2 = zeros(r, c, dimdata);
data_mod = zeros(r, c, dimdata);

for loop = 1:dimdata
    % Calculate scaling
    scaling(loop) = 1.0 / ((1.0 + alpha * (means(loop)^pp)) * (1.0 + beta * (means(6)^qq)));

    % Access incidence (assumed predefined as incidence1, incidence2, etc.)
    incidence = eval(['incidence' num2str(loop)]);
    incidence(incidence < 0) = 0;

    elevation(elevation < 0) = 0;

    % Compute fixed array
    fixed{loop} = scaling(loop) * ((1 + alpha * (power(incidence, pp))) .* ...
                                    (1 + beta * (power(elevation, qq))));

    % Generate random sequence for current timestamp
    [xx, yy, ii] = find(data(:, :, loop) > 0);
    seq = randperm(length(xx));

    for j = 1:length(seq)
        x = xx(seq(j));
        y = yy(seq(j));
        ffunc = 1 / (fixed{loop}(x, y) * (1 + gamma * (b(x, y)^rr)));
        dist_2(x, y, loop) = ffunc;
    end
end

% Modify data to convert lake class from 2 to 1, preserving other values
for loop = 1:dimdata
    for z = 1:r
        for w = 1:c
            if (data(z, w, loop) == 2)
                data_mod(z, w, loop) = 1;
            else
                data_mod(z, w, loop) = data(z, w, loop);
            end
        end
    end
end


%-------new automation for generating threshold values------------%
% Generate potential threshold values
potential_thresholds_snow = linspace(1, 1.6, 7);  % Adjust the range and granularity as needed
potential_thresholds_lake = linspace(0.75, 0.98, 7);

best_avg = inf;
best_threshold_snow = 1;
best_threshold_lake = 0;

for i = 1:size(potential_thresholds_snow, 2)
    for j = 1:size(potential_thresholds_lake, 2)
        sum = 0;  % Initialize average value for the current thresholds
        snow_threshold = potential_thresholds_snow(i);
        lake_threshold = potential_thresholds_lake(j);
        disp(['Snow Threshold:' num2str(snow_threshold)]);
        disp(['Lake Threshold:' num2str(lake_threshold)]);
        
        % Calculate averages for each timestamp
        for k = 1:(dimdata - 1)
            data_k = data(:, :, k);
            dist_k = dist_2(:, :, k);
            mat_k = zeros(size(data_k));  % Initialize mat_k appropriately
            
            % Assuming thresholding2 is defined somewhere in your code
            sum = sum + thresholding2(snow_threshold, lake_threshold, dist_k, k, data_k, data_mod, mat_k);
            disp(['Sum is:' num2str(sum)]);
        end
        
        % Update best values if the current combination yields a higher average
        avg = sum / (dimdata - 1);  % Calculate the average over the timestamps
        disp(['Average:' num2str(avg)]);
        if avg < best_avg
            best_avg = avg;
            disp(['Best Avg:' num2str(best_avg)]);
            best_threshold_snow = snow_threshold;
            best_threshold_lake = lake_threshold;
        end
    end
end

disp(['Best Threshold Snow: ' num2str(best_threshold_snow)]);
disp(['Best Threshold Lake: ' num2str(best_threshold_lake)]);
disp(['Best Average: ' num2str(best_avg)]);


% Initialize a 3D array to store the outputs
output_k_array = zeros(size(data, 1), size(data, 2), dimdata - 1);

for k = 1:(dimdata - 1)
    % Extract data and dist for the current timestamp
    data_k = data(:, :, k);
    dist_k = dist_2(:, :, k);
    
    % Call generateRaster with the existing output_k_array
    output_k_array = generateRaster(best_threshold_snow, best_threshold_lake, dist_k, k, data_k, output_k_array);
end
copy_output_k_array = output_k_array;
output_post_conv = convolution(copy_output_k_array, dimdata);   

% % -------------------Subplot Predicted vs actual-------------------------%
subplot(1,2,1)
imshow(output_post_conv(:,:,2))
title('Predicted')
subplot(1,2,2)
imshow(output_k_array(:,:,2))
title('Actual')
% % end



