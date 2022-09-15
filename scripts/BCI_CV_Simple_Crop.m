%% PCA+LDA decoder of various cropped IQ data
% Adapted from BCI_CV_Simple.m
% Crops IQData to window size before beamforming
% Author(s) - Kyle McGraw, Lydia Lin, Whitney Griggs, Sumner Norman
% Written on July 2022

%% Define user specified parameters
close all;
clc;

% ADD YOUR PATH TO THE BIN DATA HERE
data_path = 'Z:\Acq_161502';

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

% Set coords for center most important pixel
cz = 70; 
cx = 80;

% Save accuracies
accuracy = zeros(1,81);

% Loop through window sizes 0 (single pixel) to 80 (full image)
% Windows are all square with radius of window_size around the center pixel
for window_size = 0:80
    % Load accuracies in case of error mid-run
    accuracy=load('C:\Users\User\OneDrive - California Institute of Technology\Caltech\Junior\SURF\BCI\accuracy.mat').accuracy;
    % Skips previously ran window sizes
    if accuracy(window_size+1)
        fprintf('\n\tskipping: %d\n', window_size);
        continue 
    end
    fprintf('\n\twindow size: %d\n', window_size);
    fprintf('\n\tz pixels: %d to %d\n', max(cz-window_size,1),min(cz+window_size,103));
    fprintf('\n\tx pixels: %d to %d\n', max(cx-window_size,1),min(cx+window_size,128));
    % convert the IQ data into power doppler data
    % possible to parallelize this using a parfor loop if memory allows
    % f = waitbar(0,'');
    clear PD_data
    parfor i_bin = 1:length(bin_index_1D)
        bin_filename = sprintf('fUS_block_%.3d.bin', bin_index_1D(i_bin));

        % convert bin -> IQ (note that we will overwrite this variable on every
        % loop because storing all the IQ data could overload memory)
        IQ_data = bin_to_IQ(fullfile(data_path, bin_filename), n_channels, ens_length, cz, cx, window_size);   

        % convert IQ -> Power Doppler data
        PD_data(:,:,i_bin) = IQ_to_PD(IQ_data); 
    %     waitbar(i_bin/length(bin_index_1D), f,  ...
    %         sprintf('converting bin %i of %i',i_bin, length(bin_index_1D)));
    end 

    %% Reshape into 4D structure
    iDop = reshape(PD_data, size(PD_data, 1), size(PD_data, 2), size(bin_inds_to_use, 2), []);

    %% Pre-process data into appropriate shape

    iDopP = iDop;
    [zPix, xPix, nWindows, nTrials] = size(iDopP);

    % run cross validation for each sample in the memory epoch
    % this stacks the images in the epoch side by side
    dop3D = reshape(iDopP, [zPix, 3*xPix, nTrials]);

    nPixels = size(dop3D,1)*size(dop3D,2);
    data = zeros(nTrials,nPixels);
    for i = 1:size(dop3D,3)
        currentImage = dop3D(:,:,i);
        data(i,:) = reshape(currentImage,1,nPixels);
    end

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
    confusion = cp.CountingMatrix;
    percentCorrect = 100*cp.CorrectRate;
    f1 = figure();
    set(f1,'Position',[800 500 700 500])
    confusionchart(confusion(1:2,:),'Normalization','row-normalized')
    fprintf('\n\taccuracy: %2.1f%%\n', percentCorrect)
    
    % Save accuracy and confusion chart for window size
    accuracy(window_size+1)=percentCorrect;
    save('C:\Users\User\OneDrive - California Institute of Technology\Caltech\Junior\SURF\BCI\accuracy.mat','accuracy');
    exportgraphics(f1,['C:\Users\User\OneDrive - California Institute of Technology\Caltech\Junior\SURF\BCI\confusioncharts\window' num2str(window_size) '.png'])
end

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
function IQ_data = bin_to_IQ(filename, n_channels, ens_length, cz, cx, window_size)

% open the .bin file, read the IQ data, close the .bin file
fid = fopen(filename, 'r');
bin_data = fread(fid, 'double');                    % vector <double>
fclose(fid);

% reshape the IQ data to 3D (n_z, n_x*2, ens_length) <double>
IQ_3D = reshape(bin_data, [],n_channels*2, ens_length);

% reshape the IQ data to 3D (2D image + complex conjugate)
% split into real & imaginary components (n_z, n_x, ens_length) <complex>
IQ_data_uncropped = IQ_3D(:,1:n_channels,:) + ...             % real
    1i * IQ_3D(:, n_channels+1 : 2*n_channels, :);  % imag

% Crops data to square centered at pixel
IQ_data = IQ_data_uncropped(max(cz-window_size,1):min(cz+window_size,end),max(cx-window_size,1):min(cx+window_size,end),:);

end


%% IQ to Power Doppler
function PD_data = IQ_to_PD(IQ_data)

% flatten IQ_data
[nz, nx, nt] = size(IQ_data);
IQ_data = reshape(IQ_data, [nz*nx, nt]);

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
IQF_tissu =     reshape(IQF_tissu, [nz, nx, nt]);
IQ_data =       reshape(IQ_data, [nz, nx, nt]);
IQF_corrected = IQ_data-IQF_tissu;

% compute doppler power of corrected (tissue filtered) IQ
PD_data = mean(abs(IQF_corrected).^2,3);
end

