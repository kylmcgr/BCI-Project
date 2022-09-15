%% PCA+LDA decoder of various windows to beamform IQ data
% Vary circle/square beamforming, beamforming area sizes, center pixel
% coords. All decoding with 3x3 square aorund center pixel
% Author(s) - Kyle McGraw, Lydia Lin, Whitney Griggs, Sumner Norman
% Written on August 2022
% 
mpiprofile on
%% Define user specified parameters
close all;
clc;

% ADD YOUR PATH TO THE BIN DATA HERE
data_path = 'Y:\Raw Data\20180308_S2R1';

%% load data and convert to PD

% load data
memory_bin_files = load(fullfile(data_path,'memory_bin_files.mat'));
bin_inds_to_use = memory_bin_files.bin_index; 
labels = memory_bin_files.RInds; 

% get # channels & ensemble length from experiment parameters
n_channels = 128;
params = load(fullfile(data_path,'UF.mat'));
ens_length = params.UF.numFrames;

% find bin files in this directory
bin_files = dir(fullfile(data_path,'*.bin'));

% Sort them numerically
[~, sort_ind] = sort([bin_files.datenum]);
bin_files_sorted = bin_files(sort_ind);
bin_index_1D = sort(bin_inds_to_use(:));

% Load searchlight pixels
load('Y:\Raw Data\20180308_S2R1\S2R1_searchlight_results.mat');

% Fibonacci sequence for window sizes
fib = [1,2,3,5,8,13,20,30,40,50,60,70];

% Find all pixels of interest
clear pixels all_pixels
pix_num = 1;
for cz_pixels=1:size(pvalue_mask,1)
    for cx_pixels=1:size(pvalue_mask,2)
        if pvalue_mask(cz_pixels,cx_pixels)
            % If pixel is at edge of screen, use pixel next to it
            cz = cz_pixels;
            cx = cx_pixels;
            if cz_pixels == 1
                cz = cz_pixels + 1;
            elseif cz_pixels == size(pvalue_mask,1)
                cz = cz_pixels - 1;
            elseif cx_pixels == 1
                cx = cx_pixels + 1;
            elseif cx_pixels == size(pvalue_mask,2)
                cx = cx_pixels - 1;
            end
            all_pixels(pix_num,:) = [cz,cx];
            pix_num = pix_num + 1;
        end
    end
end

pixels = all_pixels(randsample(length(all_pixels),100),:);

% define upper limits of loops outside parfor
pix_len = length(pixels);
win_len = length(fib);

% convert the IQ data into power doppler data
% possible to parallelize this using a parfor loop if memory allows
% f = waitbar(0,'');
clear PD_data_array
parfor i_bin = 1:length(bin_index_1D)
    bin_filename = sprintf('fUS_block_%.3d.bin', bin_index_1D(i_bin));

    % convert bin -> IQ (note that we will overwrite this variable on every
    % loop because storing all the IQ data could overload memory)
    IQ_data = bin_to_IQ(fullfile(data_path, bin_filename), n_channels, ens_length); % 3002 total, 2965 IQ
    
    % Changed from 0:1 to 1:2 to make parfor work
    for circle=1%:2
        for pix=1:pix_len
            for window_size = 1:win_len%1:max_window
                % convert IQ -> Power Doppler data
                PD_data_array(:,i_bin,circle,pix,window_size) = IQ_to_PD(IQ_data, pixels(pix,1), pixels(pix,2), fib(window_size), circle-1); 
            end
        end
    end
%     waitbar(i_bin/length(bin_index_1D), f,  ...
%         sprintf('converting bin %i of %i',i_bin, length(bin_index_1D)));
end 
 
accuracy = zeros(1,pix_len,win_len);
% First run with square beamforming
for circle=0%:1
    sq_circ = 'square';
    if circle
        sq_circ = 'circle';
    end
    % Loop through all pixels
    for pix=1:pix_len
        % Loop through window sizes 1 (3x3 square or diameter 3
        % circle (cross)) to full image
        for window_size = 1:win_len
            % Skips window 1 for circle becuase its smaller than
            % 3x3 decoding size
            if circle && fib(window_size) == 1
                fprintf('\n\twindow not large enough: %s window %d\n', sq_circ, fib(window_size));
                continue 
            end

            clear PD_data
            PD_data = PD_data_array(:,:,circle+1,pix,window_size); 

            %% Reshape into 4D structure
%                     iDop = reshape(PD_data, size(PD_data, 1), size(PD_data, 2), size(bin_inds_to_use, 2), []);
            iDop = reshape(PD_data, size(PD_data, 1), size(bin_inds_to_use, 2), []);

            %% Pre-process data into appropriate shape

            iDopP = iDop;
%                     [zPix, xPix, nWindows, nTrials] = size(iDopP);
            [nPixels, nWindows, nTrials] = size(iDopP);

            % run cross validation for each sample in the memory epoch
            % this stacks the images in the epoch side by side
%                     dop3D = reshape(iDopP, [zPix, 3*xPix, nTrials]);
            dop3D = reshape(iDopP, [3*nPixels, nTrials]);

%                     nPixels = size(dop3D,1)*size(dop3D,2);
%                     data = zeros(nTrials,nPixels);
%                     for i = 1:size(dop3D,3)
%                         currentImage = dop3D(:,:,i);
%                         data(i,:) = reshape(currentImage,1,nPixels);
%                     end
            data = permute(dop3D,[2 1]);

            %% cross validated performance is calculated here

            % normalize to z score 
            data = zscore(data);

            % set up class performance
            cp = classperf(labels');     % initializing class performance var
            indices = crossvalind('HoldOut', nTrials, 0.2);

            % create indices for training set & test set
            test = (indices == 0);
            train = ~test;

            % set up training & testing data & labels
            trainData = data(train,:);
            trainLabels = labels(train)';
            testData = data(test,:);

            %% run the BCI

            predictedClass = BCI(trainData, trainLabels, testData);

            % Check performance
            classperf(cp, predictedClass, test);

            %% display results

            % Confusion matrix
%             confusion = cp.CountingMatrix;
            percentCorrect = 100*cp.CorrectRate;
%             f1 = figure();
%             set(f1,'Position',[800 500 700 500])
%             confusionchart(confusion(1:2,:),'Normalization','row-normalized')
%             fprintf('\n\taccuracy: %2.1f%%\n', percentCorrect)

            % Save accuracy and confusion chart for window size
            accuracy(circle+1,pix,window_size)=percentCorrect;
%             exportgraphics(f1,['C:\Users\User\OneDrive - California Institute of Technology\Caltech\Junior\SURF\BCI\confusioncharts_' sq_circ '\window' num2str(window_size) '.png'])
        end
    end
end
save(['C:\Users\User\OneDrive - California Institute of Technology\Caltech\Junior\SURF\BCI\square_accuracy_fib.mat'],'accuracy');
mpiprofile viewer

%% function to decode data
function predictedClass = BCI(trainData, trainLabels, testData)

% PCA
[pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(trainData);
explainedVarianceToKeepAsFraction = 0.95;

numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
trainPredictors = pcaScores(:,1:numComponentsToKeep);   % (PCA transformed trainData)

% LDA
MdlLinear = fitcdiscr(trainPredictors,trainLabels);

% Apply PCA+LDA to the test data
testDataPCA = (testData - pcaCenters) * pcaCoefficients;
predictedClass = predict(MdlLinear, testDataPCA);

end

%% Bin to IQ
function IQ_data = bin_to_IQ(filename, n_channels, ens_length)

% open the .bin file, read the IQ data, close the .bin file
fid = fopen(filename, 'r');
bin_data = fread(fid, 'double');                    % vector <double>
fclose(fid);

% reshape the IQ data to 3D (n_z, n_x*2, ens_length) <double>
IQ_3D = reshape(bin_data, [],n_channels*2, ens_length);

% reshape the IQ data to 3D (2D image + complex conjugate)
% split into real & imaginary components (n_z, n_x, ens_length) <complex>
IQ_data = IQ_3D(:,1:n_channels,:) + ...             % real
    1i * IQ_3D(:, n_channels+1 : 2*n_channels, :);  % imag

end


%% IQ to Power Doppler
function PD_data = IQ_to_PD(IQ_data, cz, cx, window_size, circle)
% Save the 3x3 data so we can find it in the reshaped data
IQ_data_3x3 = IQ_data(cz-1:cz+1,cx-1:cx+1,:);

[nz, nx, nt] = size(IQ_data);
% Crops data to square/circle centered at pixel
[x, z, ~] = meshgrid(1:nx, 1:nz, 1:nt);
if circle
    IQ_data = reshape(IQ_data((x - cx).^2 + (z - cz).^2 <= window_size.^2), [], nt);
else
    IQ_data = reshape(IQ_data(max(cz-window_size,1):min(cz+window_size,end),max(cx-window_size,1):min(cx+window_size,end),:), [], nt);
end

% Find the indices of the 3x3 data
% Get indices of all 3x3 places, ignore the time index
[indices,~] = find(ismember(IQ_data,IQ_data_3x3));
all_indices = reshape(indices, [], nt);
% Make sure we didn't grab a duplicate value's index
assert(range(mean(all_indices)) == 0);
IQ_indices_3x3 = all_indices(:,1);

% flatten IQ_data
% [nz, nx, nt] = size(IQ_data);
% IQ_data = reshape(IQ_data, [nz*nx, nt]);

% compute covariance matrix
cov_matrix = IQ_data'*IQ_data;

% get eigenvalues
[Eig_vect, Eig_val]= eig(cov_matrix);
Eig_vect=fliplr(Eig_vect);
Eig_val=rot90(Eig_val,2); %#ok<NASGU>
M_A = IQ_data*Eig_vect;

% separate tissue motion 
skipped_eigs =  1 : round(0.15 * nt);
IQF_tissu =     M_A(:,skipped_eigs)*Eig_vect(:,skipped_eigs)';
% IQF_tissu =     reshape(IQF_tissu, [nz, nx, nt]);
% IQ_data =       reshape(IQ_data, [nz, nx, nt]);
IQF_corrected = IQ_data-IQF_tissu;

% compute doppler power of corrected (tissue filtered) IQ
PD_data = mean(abs(IQF_corrected).^2,2);

% grab the 3x3 region of the PD using the previous calculate indices
PD_data = PD_data(IQ_indices_3x3);
% make sure we have the correct shape
assert(length(PD_data)==9);

end

