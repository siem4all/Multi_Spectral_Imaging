classdef MultiSpectralImaging
    properties
        wavelengths; % Should be a vector of wavelengths
        total_bands; % Number of bands
        roi; 
        img_bands = {}; 
    end

    methods
        function obj = MultiSpectralImaging(wavelengths, total_bands, roi) 
            obj.wavelengths = wavelengths;
            obj.total_bands = total_bands;
            obj.roi = roi;
            if roi == false
                obj.img_bands{1} = zeros(512);
                obj.img_bands{2} = zeros(512);
                obj.img_bands{3} = zeros(512);
                obj.img_bands{4} = zeros(512);
            else
                obj.img_bands{1} = zeros(234, 272);
                obj.img_bands{2} = zeros(234, 272);
                obj.img_bands{3} = zeros(234, 272);
                obj.img_bands{4} = zeros(234, 272);
            end 
        end
        
        function data = initializeData(obj)
            % Create a structure to hold dynamic variable names
            data = struct();
        end

        function obj = runMultiSpectral(obj, imageDir, numImages, dir_ind)
            % Create a structure to hold dynamic variable names
            data = obj.initializeData();
            % Get a list of all PNG files in the directory
            imageFiles = dir(fullfile(imageDir, '*.png'));
            % Check if any images were found
            if isempty(imageFiles)
                disp('No PNG images found in the specified directory.');
                return; % Exit the function early
            end
        
            % Initialize data arrays
            data.norm_avg = zeros(1, obj.total_bands);
            data.max_avg = zeros(1, obj.total_bands);
            data.min_avg = zeros(1, obj.total_bands);
            data.norm_sigma = zeros(1, obj.total_bands);
            data.max_sigma = zeros(1, obj.total_bands);
            data.min_sigma = zeros(1, obj.total_bands);

            % Process images in groups of four
            for i = 1:numImages
                % Construct the full file name
                filename = fullfile(imageDir, imageFiles(i).name);
                
                % Check if the file exists
                if isfile(filename)
                    % Read the image
                    img = imread(filename);
                    
                    % Convert to grayscale if the image is RGB
                    if size(img, 3) == 3
                        img = rgb2gray(img);
                    end
                    
                    % Call the function imBandCat 
                    obj = obj.imBandCat(img); 
                else
                    disp(['Image not found: ', filename]);
                end
            end
            
            for j = 1:obj.total_bands
                if ~isempty(obj.img_bands{j})  % Check if the band is not empty
                    std_error = std(obj.img_bands{j}(:)) / sqrt(numImages);
                    data.norm_avg(j) = mean(obj.img_bands{j}(:));  % Use (:)
                    data.max_avg(j) = data.norm_avg(j) + std_error;
                    data.min_avg(j) = data.norm_avg(j) - std_error;

                    sigma = obj.cal_absop_to_scatt_ratio(obj.img_bands{j});
                    std_sigma_error = std(sigma(:)) / sqrt(numImages);
                    data.norm_sigma(j) = mean(sigma(:));
                    data.max_sigma(j) = data.norm_sigma(j) + std_sigma_error;
                    data.min_sigma(j) = data.norm_sigma(j) - std_sigma_error;

                    if j == 4  % For the fourth band, display images
                        subplot(1, 2, 1);
                        imshow(obj.img_bands{4} / numImages);
                        title('Average image of the fourth Band');
                        subplot(1, 2, 2);
                        imshow(sigma / numImages);
                        title('Average image of the fourth Sigma(K/S)');
                    end
                end
            end
            
            obj.exportData(data.max_avg, data.norm_avg, data.min_avg, data.max_sigma, data.norm_sigma, data.min_sigma, dir_ind);
        end

        function obj = imBandCat(obj, img)
            [x, y] = size(img);  % Get the dimensions of the image
           
            if obj.roi == true
                % Define regions of interest (ROI)
                row_ban1_init = 150;
                row_ban1_final = 383;
                row_ban3_init = 662;
                row_ban3_final = 895;
                col_ban1_init = 52;
                col_ban1_final = 323;
                col_ban2_init = 564;
                col_ban2_final = 835;
            else
                % Define full image
                row_ban1_init = 1;
                row_ban1_final = x / 2;
                row_ban3_init = x / 2 + 1;
                row_ban3_final = x;
                col_ban1_init = 1;
                col_ban1_final = y / 2;
                col_ban2_init = y / 2 + 1;
                col_ban2_final = y;
            end
            
            % Check if the dimensions are even
            if mod(x, 2) ~= 0 || mod(y, 2) ~= 0
                error('Image dimensions must be even.');
            end
            
            img = im2double(img);

            % Define the row and column indices for each band
           row_indices = [row_ban1_init, row_ban1_final; row_ban1_init, row_ban1_final; row_ban3_init, row_ban3_final; row_ban3_init, row_ban3_final];
           col_indices = [col_ban1_init, col_ban1_final; col_ban2_init, col_ban2_final; col_ban1_init, col_ban1_final; col_ban2_init, col_ban2_final];

          % Update img_bands using one line
          for i = 1:obj.total_bands
                  obj.img_bands{i} = obj.img_bands{i} + img(row_indices(i, 1):row_indices(i, 2), col_indices(i, 1):col_indices(i, 2));
          end
        end

        function sigma = cal_absop_to_scatt_ratio(obj, img)
            % Check if the input image is valid
            if isempty(img) || ~isnumeric(img)
                error('Input image must be a non-empty numeric array.');
            end

            img = im2double(img); % Ensure the image is double
            % Avoid division by zero
            img(img == 0) = 1e-10; % Small value to prevent log(0)

            % Calculate the absorption to scattering ratio for each pixel
            sigma = ((1 - img).^2) ./ (2 * img);
        end

        function obj = exportData(obj, max_avg, norm_avg, min_avg, max_sigma, norm_sigma, min_sigma, dir_ind)
            % Create dynamic filenames using the provided directory index
            avg_sigma = sprintf('Avg_sigma_%d', dir_ind);  % Filename for average sigma data
            avg_Ref = sprintf('Avg_reflection_%d', dir_ind); % Filename for average reflection data
            sigma_filename = strcat(avg_sigma, '.mat');      % Complete filename for sigma data
            avg_filename = strcat(avg_Ref, '.mat');          % Complete filename for reflection data

            % Check if the average reflection file exists and remove it if it does
            if exist(avg_filename, 'file')
                delete(avg_filename);  % Delete the existing file
            end

            % Check if the sigma file exists and remove it if it does
            if exist(sigma_filename, 'file')
                delete(sigma_filename);  % Delete the existing file
            end

            % Save data to the dynamically generated .mat file for average reflection
            save(avg_filename, 'max_avg', 'norm_avg', 'min_avg');
            % Save data to the dynamically generated .mat file for sigma values
            save(sigma_filename, 'max_sigma', 'norm_sigma', 'min_sigma');
        end
    end
end
