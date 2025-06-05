classdef DataVisualize
    properties
        ImageFolder % Folder containing images
        ResultsPath % Path for results
        snv_mean
    end

    methods
        function obj = DataVisualize(imageFolder, resultsPath)
            % Constructor to initialize properties
            obj.ImageFolder = imageFolder;
            obj.ResultsPath = resultsPath;
            obj.snv_mean=[0,0];
        end
        function graph_coloring(obj,i, si_images, wavelengths, exp_result)
            % Step 1: Define the Laplacian kernel
            num_bins = 30;
            if i<5
                value1=histcounts(si_images{i}, num_bins);
                value2=histcounts(si_images{i+4}, num_bins);
                text_pos=min(max(value1), max(value2));
                min_p_value=min(min(si_images{i}(:)),min(si_images{i+4}(:)));
                max_p_value=max(max(si_images{i}(:)),max(si_images{i+4}(:)));
                x_max=max(mean(si_images{i}(:)), mean(si_images{i+4}(:)));
            else
                value1=histcounts(si_images{i+4}, num_bins);
                value2=histcounts(si_images{i+6}, num_bins);
                text_pos=min(max(value1), max(value2));
                min_p_value=min(min(si_images{i+4}(:)),min(si_images{i+6}(:)));
                max_p_value=max(max(si_images{i+4}(:)),max(si_images{i+6}(:)));
                x_max=max(mean(si_images{i+4}(:)), mean(si_images{i+6}(:)));
            end
            band_num   = @(i) mod(i, 4) + (mod(i, 4) == 0) * 4; % Define the lambda expression
            cnt_pos=1;
            params=[];
            %for the FDR
            findex=1;
            fdr_mu=[0, 0];
            fdr_std=[0,0];
            maxfguassian=0;%max guassian value for text of FDR
            figure;
            for j=0:length(wavelengths):5
                band_index=i+j;
                if i>4 && j==0
                    band_index=i+length(wavelengths);
                end
                if i>4 && j~=0
                    band_index=i+length(wavelengths)+2;
                end
                if j==0
                    exp_result_index=1;
                    color=[0, 0.4470, 0.7410];
                else
                    exp_result_index=2;
                    color=[0.8500, 0.3250, 0.0980];
                end
                % Power Spectral Density
                psd=mean2(obj.powerSpectralDensity(si_images{band_index}));
                subplot(2,2,cnt_pos);
                imagesc(si_images{band_index}); % Display the data as an image
                axis square; % Correct the axis orientation
                colorbar; % Add a color bar for reference
                % Step 2: Set color limits for the color map
                clim([min_p_value, max_p_value]); % Set color limits between 30 and 100
                % Step 3: Customize the colormap
                colormap(jet); % Change to the 'jet' colormap
                % Add a title using sprintf
                title(sprintf("%s %s Color Graph", exp_result{exp_result_index}, exp_result{3})); % Dynamic title
                xlabel('X'); % Label the X-axis
                ylabel('Y'); % Label the Y-axis
                hold on;
                % Calculate mean and standard deviation
                % Compute the histogram
                counts = histcounts(si_images{band_index}, num_bins);
                % Total number of pixels
                totalPixels = sum(counts);

                % Calculate probabilities
                probabilities = counts / totalPixels;

                % Remove zero probabilities (to avoid log(0))
                probabilities = probabilities(probabilities > 0);
                [x,y]=size(counts);
                zerovectror=zeros(x, y/2);
                countszeroappended=[counts, zerovectror];
                % Calculate entropy
                entropy = -sum(probabilities .* log2(probabilities));
                % Create a subplot for the histogram
                subplot(2, 2, 2);
                bar(countszeroappended);
                xlim([0, length(counts)+num_bins]);  % Adjust x-axis limits for padding
                ylim([0, text_pos+mean(counts)]+2*num_bins);  % Adjust y-axis limits for padding
                mu=mean2(si_images{band_index});
                %To save snv mean of fused image
                if i==5 && j==0
                    obj.snv_mean(1)=mu;
                end
                if i==5 && j~=0
                    obj.snv_mean(2)=mu;
                end
                sigma=std2(si_images{band_index});
                fdr_mu(findex)=mu;
                fdr_std(findex)=sigma;
                findex=findex+1;% incerement fdr index for mean and std
                skew=skewness(counts);
                kurto=kurtosis(counts);
                energy=mean2((si_images{band_index}).*2);
                % % Add mean, std, skewness, and kurtosis as text labels
                if j==0
                    text('String', sprintf('Mean:%.2f |', mu), ...
                        'Position', [num_bins, text_pos], ...
                        'Color', color, 'FontSize', 6, 'FontWeight', 'bold');
                    text('String', sprintf('Std:%.2f |', sigma), ...
                        'Position', [num_bins, text_pos*0.9], ...
                        'Color', color, 'FontSize', 6, 'FontWeight', 'bold');

                    text('String', sprintf('Skewness:%.2f |', skew), ...
                        'Position', [num_bins, text_pos*0.8], ...
                        'Color', color, 'FontSize', 6, 'FontWeight', 'bold');

                    text('String', sprintf('Kurtosis:%.2f |', kurto), ...
                        'Position', [num_bins, text_pos*0.7], ...
                        'Color', color, 'FontSize', 6, 'FontWeight', 'bold');
                    text('String', sprintf('Entropy:%.2f |', entropy), ...
                        'Position', [num_bins, text_pos*0.6], ...
                        'Color', color, 'FontSize', 6, 'FontWeight', 'bold');
                    text('String', sprintf('Energy:%.2f |', energy), ...
                        'Position', [num_bins, text_pos*0.5], ...
                        'Color', color, 'FontSize', 6, 'FontWeight', 'bold');
                    text('String', sprintf('PSD:%.2f |', psd), ...
                        'Position', [num_bins, text_pos*0.4], ...
                        'Color', color, 'FontSize', 6, 'FontWeight', 'bold');
                else
                    text('String', sprintf('%.2f', mu), ...
                        'Position', [num_bins+20, text_pos], ...
                        'Color', color, 'FontSize', 6, 'FontWeight', 'bold');
                    text('String', sprintf('%.2f', sigma), ...
                        'Position', [num_bins+20,text_pos*0.9], ...
                        'Color', color, 'FontSize', 6, 'FontWeight', 'bold');

                    text('String', sprintf(' %.2f', skew), ...
                        'Position', [num_bins+20, text_pos*0.8], ...
                        'Color', color, 'FontSize', 6, 'FontWeight', 'bold');

                    text('String', sprintf('%.2f', kurto), ...
                        'Position', [num_bins+20, text_pos*0.7], ...
                        'Color', color, 'FontSize', 6, 'FontWeight', 'bold');
                    text('String', sprintf('%.2f', entropy), ...
                        'Position', [num_bins+20, text_pos*0.6], ...
                        'Color', color, 'FontSize', 6, 'FontWeight', 'bold');
                    text('String', sprintf('%.2f', energy), ...
                        'Position', [num_bins+20, text_pos*0.5], ...
                        'Color', color, 'FontSize', 6, 'FontWeight', 'bold');
                    text('String', sprintf('%.2f', psd), ...
                        'Position', [num_bins+20, text_pos*0.4], ...
                        'Color', color, 'FontSize', 6, 'FontWeight', 'bold');
                end
                cnt_pos=cnt_pos+2;
                hold on;
                max_value=2*x_max+12;
                %  Create the new x-axis from 0 to 300, with 500 points
                x = linspace(0, max_value, 500);  % Define a higher range and finer resolution for Gaussian
                min_value=0;
                if i==6
                    max_value=15;
                    min_value=-10;
                    x=linspace(min_value, max_value, 500);
                end
                gaussian=obj.gaussian_distribution(x,sigma, mu);
                maxfguassian=max(max(gaussian * sum(counts)), maxfguassian);
                %  Plot the histogram and Gaussian distribution
                % Plot the histogram
                subplot(2, 2, 4); % Histogram in the first subplot
                plot(x, gaussian * sum(counts), 'color',color, 'LineWidth', 2); % Scale the Gaussian to match the histogram
                xlim([min_value, max_value]);  % Set the x-axis from 0 to 300 for the Gaussian plot
                hold on;
                %The parametrs
                temp_params=[band_index, mu, sigma, skew, kurto, entropy, energy, psd];
                params=[params;temp_params];
            end
            subplot(2, 2, 2);
            title('Histogram');
            xlabel('Pixel Intensity');
            ylabel('Counts');
            hLegend = legend(exp_result{1}, exp_result{2}, 'Location', 'southeast');
            set(hLegend, 'FontSize', 4, 'EdgeColor', 'none'); % Adjust font size and remove edge color
            subplot(2,2,4);
            title('Gaussian Dist.');
            text('String', sprintf('FDR = %.3f', (fdr_mu(1)-fdr_mu(2)).^2./(fdr_std(1)^2-fdr_std(2)^2)), ...
                        'Position', [max(fdr_mu)+3*max(fdr_std)/2, maxfguassian/2+10], ...
                         'FontSize', 8, 'FontWeight', 'bold');
            xlabel('Pixel Intensity');
            ylabel('Frequency');
            hLegend = legend(exp_result{1}, exp_result{2});
            set(hLegend, 'FontSize', 5, 'EdgeColor', 'none'); % Adjust font size and remove edge color
            saveas(gcf, fullfile(obj.ResultsPath, sprintf('color_graph_of_%d_nm_%d.png',wavelengths(band_num(i)),i)));
            hold off;
            %Export The parametrs
            obj.exportparam2excel(wavelengths, params, exp_result)
            % if i<5
            %        value1=histcounts(conv2(si_images{i}, laplacianKernel, 'same'), num_bins);
            %        value2=histcounts(conv2(si_images{i+4}, laplacianKernel, 'same'), num_bins);
            %        text_pos=max(max(value1), max(value2))-61;
            % else
            %       value1=histcounts(conv2(si_images{i+4}, laplacianKernel, 'same'), num_bins);
            %       value2=histcounts(conv2(si_images{i+4+1}, laplacianKernel, 'same'), num_bins);
            %       text_pos=max(max(value1), max(value2))-61;
            % end
            % figure;
            % cnt_pos=1;
            % for j=0:4:5
            %      band_index=i;
            %      exp_result_index=1;
            %      color=[0, 0.4470, 0.7410];
            %     if j~=0 && i<5
            %         band_index=i+j;
            %     end
            %     if j==0 && i>=5
            %        band_index=i+4;
            %     end
            %     if j~=0 && i>=5
            %        band_index=i+4+1;
            %     end
            %
            %     if j~=0
            %         exp_result_index=2;
            %         color=[0.8500, 0.3250, 0.0980];
            %    end
            %
            %    % Apply convolution
            %    img = conv2(si_images{band_index}, laplacianKernel, 'same');
            % %    %Plot histogram and Gaussian curve
            %   subplot(2,2,cnt_pos);
            %   imagesc(img); % Display the data as an image
            %   axis square; % Correct the axis orientation
            %   colorbar; % Add a color bar for reference
            %   % Step 2: Set color limits for the color map
            %  clim([4, 16]); % Set color limits between 30 and 100
            %  % Step 3: Customize the colormap
            %  colormap(jet); % Change to the 'jet' colormap
            %  % Add a title using sprintf
            %  title(sprintf("%s %s 2nd Derivative", exp_result{exp_result_index}, exp_result{3}));
            %  xlabel('X'); % Label the X-axis
            %  ylabel('Y'); % Label the Y-axis
            % hold on;
            % % Apply convolution
            % counts=histcounts(img(:), num_bins);
            % % Total number of pixels
            % totalPixels = sum(counts);
            %
            % % Calculate probabilities
            % probabilities = counts / totalPixels;
            %
            %  % Remove zero probabilities (to avoid log(0))
            %  probabilities = probabilities(probabilities > 0);
            %  [x,y]=size(counts);
            %  zerovectror=zeros(x, y/2);
            %  zeroappended=[counts, zerovectror];
            %  % Calculate entropy
            %  H = -sum(probabilities .* log2(probabilities));
            % subplot(2, 2, 2);
            % bar(zeroappended);
            % xlim([0, length(counts)+num_bins+20]);  % Adjust x-axis limits for padding
            % ylim([0, text_pos+num_bins]);  % Adjust x-axis limits for padding
            % % Add mean and std as text labels
            % if j==0
            %     text('String', sprintf('Mean:%.2f | ', mean2(img)), ...
            %      'Position', [num_bins, text_pos*0.9], ...
            %      'Color', color, 'FontSize', 8, 'FontWeight', 'bold');
            %
            %     text('String', sprintf('Std:%.2f | ', std2(img)), ...
            %      'Position', [num_bins, text_pos*0.8], ...
            %      'Color', color, 'FontSize', 8, 'FontWeight', 'bold');
            %     text('String', sprintf('skewns:%.2f | ', skewness(counts)), ...
            %      'Position', [num_bins, text_pos*0.7], ...
            %      'Color', color, 'FontSize', 8, 'FontWeight', 'bold');
            %     text('String', sprintf('kurtss:  %.2f | ', kurtosis(counts)), ...
            %      'Position', [num_bins, text_pos*0.6], ...
            %      'Color', color, 'FontSize', 8, 'FontWeight', 'bold');
            %     text('String', sprintf('Entropy: %.2f | ', H), ...
            %     'Position', [num_bins, text_pos*0.5], ...
            %     'Color', color, 'FontSize', 8, 'FontWeight', 'bold');
            %    text('String', sprintf('Energy:%.2f | ', sum(sum((img).^2))/numel(img)), ...
            %     'Position', [num_bins, text_pos*0.4], ...
            %     'Color', color, 'FontSize', 8, 'FontWeight', 'bold');
            %      cnt_pos=cnt_pos+2;
            % else
            %    text('String', sprintf(' %.2f', mean2(img)), ...
            %      'Position', [num_bins+30, text_pos*0.9], ...
            %      'Color', color, 'FontSize', 8, 'FontWeight', 'bold');
            %
            %     text('String', sprintf(' %.2f', std2(img)), ...
            %      'Position', [num_bins+30, text_pos*0.8], ...
            %      'Color', color, 'FontSize', 8, 'FontWeight', 'bold');
            %     text('String', sprintf(' %.2f', skewness(counts)), ...
            %      'Position', [num_bins+30, text_pos*0.7], ...
            %      'Color', color, 'FontSize', 8, 'FontWeight', 'bold');
            %     text('String', sprintf(' %.2f', kurtosis(counts)), ...
            %      'Position', [num_bins+30, text_pos*0.6], ...
            %      'Color', color, 'FontSize', 8, 'FontWeight', 'bold');
            %     text('String', sprintf(' %.2f', H), ...
            %     'Position', [num_bins+30, text_pos*0.5], ...
            %     'Color', color, 'FontSize', 8, 'FontWeight', 'bold');
            %    text('String', sprintf(' %.2f', sum(sum((img).^2))/numel(img)), ...
            %     'Position', [num_bins+30, text_pos*0.4], ...
            %     'Color', color, 'FontSize', 8, 'FontWeight', 'bold');
            % end
            % hold on;
            % end
            % subplot(2, 2, 2);
            % title('Histogram');
            % xlabel('Pixel Intensity');
            % ylabel('Counts');
            % hLegend = legend(exp_result{1}, exp_result{2}, 'Location', 'southeast');
            % set(hLegend, 'FontSize', 4, 'EdgeColor', 'none'); % Adjust font size and remove edge color
            % saveas(gcf, fullfile(obj.ResultsPath, sprintf('histogram_of_%d_nm_%d.png',wavelengths(band_num(i)),i)));
            % hold off;
        end
        function gaussian=gaussian_distribution(~, x,sigma, mu)
            %  Calculate the Gaussian normal distribution based on the mean and std
            gaussian = (1 / (sigma * sqrt(2 * pi))) * exp(-0.5 * ((x - mu) / sigma).^2);
        end
        function makePlot(obj, exp_result, si_images)
            % Define arrays for line widths, colors, and styles
            lineWidths      = [1, 1.5, 2, 2.5];
            lineColors      = {'blue', 'red', 'blue', 'red'};
            lineStyles      = {'-', '--', ':', '-.'};
            wavelengths     = [735, 800, 865, 930];
            names           = {'Patient', ' ', 'Patient'};

            % Initialize a figure
            figure;
            prefixes        = {'Avg reflection', 'Avg si'};
            avgs            = {'si_std_error','avg_si_refl'};
            data_before     = load(sprintf('%s %d.mat', prefixes{2}, 1));
            data_after      = load(sprintf('%s %d.mat', prefixes{2}, 2));
            % Access the correct field from the loaded data
            % Assuming the structure of loaded data is consistent
            avg_before      = data_before.(avgs{2});
            se_before       =data_before.(avgs{1});
            avg_after       = data_after.(avgs{2});
            se_after        = data_after.(avgs{1});
            ymax_value      = max(max(avg_before+se_before), max(avg_after+se_after))+35;
            obj.plotLinearFit(wavelengths, avg_before, se_before, avg_after, se_after, si_images, ymax_value);
            % Set labels and title
            xlabel('Wavelength (nm)');
            xlim([700,1000]);
            ylim([15, ymax_value]);
            ylabel('Avg Si (K/S)');
            set(gca, 'XTick', wavelengths, 'XTickLabel', string(wavelengths));
            grid on;
            title(sprintf('Wavelength vs. Si(K/S) of %s', names{2}));
            legend(exp_result{1}, exp_result{2});

            % Save the figure
            saveas(gcf, fullfile(obj.ResultsPath, 'average_and_si.png'));
            hold off; % Release the plot
            obj.polynomialfit(wavelengths, avg_before, se_before, avg_after, se_after, exp_result)
        end

        function changeOfY(obj, roi, isNormalized)
            % This function calculates the percentage decrease of avg si before relative to after
            results                     = {}; % Initialize a cell array to store results
            for i = 1:1
                name                    = 'Patiant';
                % Load the result of before pressure
                before                  = load(sprintf('Avg si %d.mat', i));
                % Load the result of after pressure
                after                   = load(sprintf('Avg si %d.mat', i + 1));
                % Store results in a cell array
                results{end + 1, 1}     = name; % Store the name
                results{end, 2}         = string(roi); % Store the roi
                results{end, 3}         = string(isNormalized); % Store the roi
                for j = 1:length(after.avg_si_refl)
                    % Calculate the difference and store it at index j
                    results{end, j + 3} = abs(after.avg_si_refl(j) - before.avg_si_refl(j)) / after.avg_si_refl(j) * 100;
                end
            end

            % Convert results to a table
            resultsTable                = cell2table(results, 'VariableNames', {'Name', 'ROI','isNormalized', 'Band1', 'Band2', 'Band3', 'Band4'});

            % Define the output file path
            outputFilePath              = fullfile(obj.ResultsPath, 'percentage_decrease_results.csv');

            % Check if the file exists
            if isfile(outputFilePath)
                % Append to the existing file
                writetable(resultsTable, outputFilePath, 'WriteMode', 'append', 'WriteVariableNames', false);
            else
                % Write the table to a new CSV file
                writetable(resultsTable, outputFilePath);
            end
        end

        function [position, x, y] = regionOfInterest(obj, roi, img_proc)
            % Get all image files from the specified folder
            imageFiles        = dir(fullfile(obj.ImageFolder, '*.png')); % Adjust file extension as needed
            if isempty(imageFiles)
                disp('No images found in the specified folder.');
                return; % Exit if no images are found
            end

            % Select an image from the folder
            selectedImageFile = imageFiles(1).name;
            firstImagePath    = fullfile(obj.ImageFolder, selectedImageFile); % Full path to the selected image

            % Load the selected image
            img               = imread(firstImagePath);
            % Convert the image to grayscale
            gray_image        = im2gray(img);
            img_flt=im2double(gray_image);

            % Get the size of the grayscale image
            [x, y]            = size(img_flt);
            if ~roi
                %position      = [0,0,0,0];
            else
                % Extract the top-left quarter of the image
                %img_quarter =  medfilt2(img_flt(floor(x/2):x, floor(y/2):y),[3 3]);
                img_quarter =  img_flt(1:floor(x/2+1), floor(y/2):y);
                %image processing to enhance the contrast
                %img_enhanced=medfilt2(imlocalbrighten(obj.image_processing(img_quarter)), [3 3], "symmetric");
                if img_proc==true
                    img_enhanced = imadjust(img_quarter);
                    img_enhanced = imlocalbrighten(img_enhanced, 0.3, "AlphaBlend",true);
                    img_enhanced=imbilatfilt(img_enhanced);
                else
                    img_enhanced = imlocalbrighten(img_quarter);
                end
                %img_enhanced(img_enhanced > 220) = img_enhanced(img_enhanced > 220) - 30;
                imshow(img_enhanced+50, []); % Display the original image
                imwrite(img_enhanced,"img_enhanced.png")
                title('Draw the ROIs '); % Title for clarity
                text('String', sprintf('Size: %s', mat2str(size(img_enhanced)-1)), ...
                    'Position', [5, 10], ...
                    'Color', 'yellow', 'FontSize', 8, 'FontWeight', 'bold');
                hold on; % Hold on to overlay the ROIs

                % Create a freehand ROI
                h = drawrectangle();

                % Wait for the user to finalize the ROI
                position     = round(h.Position); % Get the coordinates of the ROI
                %position =[125,   227,   195,   142];
                % Extract coordinates from the position
                % Extract the ROI from the original image
                roiImage = img(position(2):position(2)+position(4)-1, position(1):position(1)+position(3)-1);
                text('String', sprintf('Size: %s', mat2str(size(roiImage))), ...
                    'Position', [position(1)+15, position(2)+7], ...
                    'Color', 'yellow', 'FontSize', 8, 'FontWeight', 'bold');
                saveas(gcf, fullfile(obj.ResultsPath, 'roiImage.png'));
            end
        end
        function enhanced_image=image_processing(~, image)
            % Add the path to the directory containing enhance_image function
            addpath('C:\Users\DELL\Documents\MATLAB\infrared-enhancement-main');
            % Call the enhance_image function
            enhanced_image = enhance_image(image);
        end
        function exportparam2excel(~, wavelengths, params, exp_result)
            % This function exports statistical results to a CSV file.
            temp_params=zeros(1,length(params(1,:)));
            temp_params(1)=params(2,1);
            for i=2:length(params(1,:))
                % relative change is abs(after-before)/before*100
                temp_params(i)= abs((params(2,i)-params(1,i))/params(1,i))*100;
            end
            params=[params;temp_params];
            for j=1:length(params(:,1))
                status=exp_result{1};
                if params(j,1)>4
                    status=exp_result{2};
                end
                band_num   = @(i) mod(i, 4) + (mod(i, 4) == 0) * 4; % Define the lambda expression
                % Initialize a cell array to store results
                results = {};
                % Store results in a cell array
                results{end + 1, 1} = string(wavelengths(band_num(params(j,1))));
                results{end, 2} = status;
                results{end, 3} = string(params(j,2));         % Store the mean
                results{end, 4}     = string(params(j,3));   % Store the standard deviation
                results{end, 5}     = string(params(j,4)); % Store the skewness
                results{end, 6}     = string(params(j,5)); % Store the kurtosis
                results{end, 7}     = string(params(j,6));  % Store the entropy
                results{end, 8}     = string(params(j,7));   % Store the energy
                results{end, 9}     = string(params(j,8));   % Store the energy

                % Convert results to a table
                resultsTable = cell2table(results, 'VariableNames', {'Band','Status','Mean', 'Std', 'Skewness', 'Kurtosis', 'Entropy', 'Energy','PSD'});
                % Define the output file path
                outputFilePath = fullfile('results', 'all_parametrs.csv');

                % Check if the file exists
                if isfile(outputFilePath)
                    % Append to the existing file
                    writetable(resultsTable, outputFilePath, 'WriteMode', 'append', 'WriteVariableNames', false);
                else
                    % Write the table to a new CSV file
                    writetable(resultsTable, outputFilePath);
                end
            end
        end
        function psd_uint8=powerSpectralDensity(~, img)
            %  Perform the 2D Fourier Transform
            F = fft2(img);

            %  Shift the zero-frequency component to the center
            Fshift = fftshift(F);

            % Compute the magnitude and then the power spectrum
            magnitude = abs(Fshift);  % Magnitude of the FFT
            psd = magnitude.^2;  % Power Spectrum

            %  Normalize (optional)
            psd = log(1 + psd);  % Apply log scaling for better visualization (optional)
            % Normalize the PSD to the range [0, 255] for proper visualization
            psd_normalized = mat2gray(psd) * 255;

            % Convert to uint8 for saving as an image
            psd_uint8 = mean2(uint8(psd_normalized));

        end
        function plotLinearFit(~, wavelengths, before, se_before, after, se_after, si_images, ymax_value)
            % PLOTLINEARFIT - This function fits a linear model to the provided data
            % Fit a linear model to the 'before' data
            md_before = fitlm(wavelengths', before'); % Fit model with wavelengths as predictor and before as response
            md_after = fitlm(wavelengths', after'); % Fit model with wavelengths as predictor and before as response
            x_test=linspace(720, 945, 100);
            predict_before=predict(md_before, x_test');
            predict_after=predict(md_after, x_test');
            % Scatter plot for 'before' values
            %scatter(wavelengths, before, 'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'Before');
            errorbar(wavelengths, before, se_before, 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'DisplayName', 'Before');
            hold on; % Hold on to plot 'after' values in the same figure
            % Scatter plot for 'after' values
            % scatter(wavelengths, after, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'After');
            errorbar(wavelengths, after, se_after, '*', 'Color', 'r', 'MarkerFaceColor', 'r',  'DisplayName', 'After');
            % Display MSE values as text
            hold on; % Hold on to plot 'after' values in the same figure
            % Plot the predicted values for 'before'
            % errorbar(x_test,predict_before,uint8(255*se_before),'o-','Color', 'b', 'LineWidth', 1.5);
            % errorbar(x_test, predict_after ,uint8(255*se_after),'*-', 'Color', 'r','LineWidth', 1.5);
            plot(x_test,predict_before, 'b-', 'LineWidth', 2);
            plot(x_test, predict_after ,'r-', 'LineWidth', 2);
            ylim([min(min(before), min(after))-10, max(max(before),max(after))+20]);
            % Prepare the model equation text
            int_bfr = md_before.Coefficients.Estimate(1);
            int_afr = md_after.Coefficients.Estimate(1);
            slope_bfr = md_before.Coefficients.Estimate(2);
            slope_afr = md_after.Coefficients.Estimate(2);
            before_eq = sprintf('y = %.2fx + %.2f', slope_bfr, int_bfr);
            after_eq = sprintf('y = %.2fx + %.2f', slope_afr, int_afr);
            % Calculate the y-coordinate for the text above the predicted line
            y_text_before = slope_afr * 735 + int_afr +7; % Adjust the offset as needed
            y_text_after = slope_afr * 735 + int_afr +5; % Adjust the offset as needed
            angle=atan(max(slope_bfr,slope_afr))* (180/pi);
            % Set y-axis limits to ensure text is visible
            text(720, ymax_value-31, before_eq, 'Color', 'blue', 'FontSize', 8,  'Rotation', 0);
            text(720, ymax_value-38, after_eq, 'Color', 'red', 'FontSize', 8,  'Rotation', 0);
            text('String', '\underline{Fused Image Without SNV}', ...
                'Position', [720, ymax_value-6], ...
                'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
            text('String', sprintf('Mean $\\pm$ Std: %.2f $\\pm$ %.2f', mean2(si_images{9}), std2(si_images{9})), ...
                'Position', [720, ymax_value-17], ...
                'Color', 'blue', 'Interpreter', 'latex', 'FontSize', 10, 'FontWeight', 'bold');
            text('String', sprintf('Mean $\\pm$ Std: %.2f $\\pm$ %.2f', mean2(si_images{11}), std2(si_images{11})), ...
                'Position', [720, ymax_value-24], ...
                'Color', 'red', 'Interpreter', 'latex', 'FontSize', 10, 'FontWeight', 'bold');

            hold on;
        end
        function polynomialfit(obj, wavelengths, before, se_before, after, se_after, exp_result)
            % Choose the degree of the polynomial
            degree = 3; % Change this value for different polynomial degrees

            % Fit a polynomial to the data
            before_coef = polyfit(wavelengths, before, degree);
            after_coef = polyfit(wavelengths, after, degree);

            % Generate x values for plotting the polynomial curve
            x_fit = linspace(730, 935, 100);
            before_fit = polyval(before_coef, x_fit);
            after_fit = polyval(after_coef, x_fit);

            % Plot the original data points
            figure;
            errorbar(wavelengths, before, se_before, 'o', 'Color', 'b', 'MarkerFaceColor', 'b',  'DisplayName', 'Before');
            hold on;
            errorbar(wavelengths, after, se_after, '*', 'Color', 'r', 'MarkerFaceColor', 'r',  'DisplayName', 'After');
            hold on;
            % Plot the fitted polynomial curve
            plot(x_fit, before_fit, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Polynomial Fit (degree %d)', degree));
            hold on;
            % Plot the fitted polynomial curve
            plot(x_fit, after_fit, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Polynomial Fit (degree %d)', degree));

            % Add labels and title
            xlabel('Wavelength (nm)');
            ylabel('Avg Si (K/S)');
            set(gca, 'XTick', wavelengths, 'XTickLabel', string(wavelengths));
            title('Polynomial Fitting');
            legend(exp_result{1}, exp_result{2});
            grid on;
            % Save the figure
            saveas(gcf, fullfile(obj.ResultsPath, 'polyfit.png'));
            hold off;


        end
        
    end
end
