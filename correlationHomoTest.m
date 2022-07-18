% this is basically similar to correlationVarFunc02 for the input
% variables and datastructure, meant to run on different chunks of
% the correlation matrix 

% this function gets the mean and variance of the correlation
% values between the given index of genes from the correlation
% matrix of the datasets for one tissue. I had to write a function
% for it, limiting the input matrix since loading all the datasets
% on the ram at a time was not possible and also I can do it in
% parallel using the function
%
% function[] = correlationVarFunc02(ib, ie, jb, je, tissue)
%
% where ib, ie, jb and je define the sub-matrix of the correlation
% matrix and tissue is a string, name of the tissue(which is name
% of the folder for tissue)
% Notes:
% 1. this function is faster than the original
% correlationVarFunc.m, since it uses MATLAB 3dim matrix 
% 2. this function uses the rank normalization. 

% for testing the function
% ib = 1 
% ie = 1000
% jb = 1
% je = 1000
% tissue = 'blood'

function[] = correlationHomoTest(ib, ie, jb, je, tissue)

dataFolder = '~/data/array/';
saveFolder = '/space/grp/marjan/data/corrHomoTest/'
gCount = 18494;

% ib = str2num(ib);
% ie = str2num(ie);
% jb = str2num(jb);
% je = str2num(je);


fileList = dir([dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/*.mat'])
%wholeData = zeros(length(fileList), (ie-ib+1), (je-jb+1));
a = (ie-ib+1)
b = (je-jb+1)
corrKai = zeros(a, b);
corrAvg = zeros(a, b);
dsCount = length(fileList);
sampleCount = zeros(1, dsCount);
       allCorr = zeros(dsCount,a, b);
        allZ = zeros(dsCount, a, b);
        allWZ = zeros(dsCount,a, b);
        
        size(allZ)
        size(corrKai)
       size(allCorr)
 

%reading files and getting the desired correlation
'loading datasets'
    for i = 1:dsCount  % i, j, k are USED. pick l 
        i
        Name = fileList(i).name;
        load([dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/' ...
              fileList(i).name])
        gCount = size((dataSet.mat), 1);
        sampleCount(i) = size(dataSet.mat, 2);
        
        % for each correlation we need the z and zw values
        bigCorr = corr(dataSet.mat'); % pearson correlation
        
        sCorr = bigCorr(ib:ie, jb:je);
        tempZ = .5.*log((1+sCorr)./(1-sCorr));
        tempWZ = (sampleCount(i)-3) .* tempZ;
        
        allCorr(i, :, :) = sCorr;
        allZ(i, :, :) = tempZ;
        allWZ(i, :, :) = tempWZ;
       
        clear sCorr;
        clear tempZ; 
       clear tempWZ;
    end 
    
    'datasets loaded'
    
    for i = 1:(ie-ib+1)
        for j = 1:(je-jb+1)
            temp = 0;
             tCorrSum = sum(allWZ(:, i, j)); 
            tSumSample = sum(sampleCount - 3);
            tAvgCorr = tCorrSum/tSumSample;
            temp = (sampleCount-3)'.*((allZ(:, i, j) - tAvgCorr).^2);
            corrKai(i,j) = sum(temp);
            corrAvg(i,j) = tAvgCorr;
        end
    end

    fileName = sprintf('%scorreHomoAvg_%s_%d_%d_%d_%d.mat', saveFolder,tissue, ib, ie, ...
                       jb, je)
        save(fileName, 'corrAvg', '-v7.3')

    fileName = sprintf('%scorreHomochi2_%s_%d_%d_%d_%d.mat', saveFolder,tissue, ib, ie, ...
                       jb, je);
    save(fileName, 'corrKai', '-v7.3')
    
end
