% this function gets the mean and variance of the correlation
% values between the given index of genes from the correlation
% matrix of the datasets for one tissue. I had to write a function
% for it, limiting the input matrix since loading all the datasets
% on the ram at a time was not possible and also I can do it in
% parallel using the function
%
% function[] = correlationVarFunc01(ib, ie, jb, je, tissue)
%
% where ib, ie, jb and je define the sub-matrix of the correlation
% matrix and tissue is a string, name of the tissue(which is name
% of the folder for tissue)
% Notes:
% 1. this function is faster than the original
% correlationVarFunc.m, since it uses MATLAB 3dim matrix 
% 2. for normalization, this function doesn't use the ranking but
% transfers the correlation values between [0,1], this is faster
% than rank normalization but might not be the best normalization

function[] = correlationVarFunc01(ib, ie, jb, je, tissue)

dataFolder = '~/data/array/';
saveFolder = '/space/grp/marjan/data/corrVar/'
gCount = 18377;

fileList = dir([dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/outlierRemoved/*.mat'])
wholeData = zeros(length(fileList)-2, (ie-ib+1), (je-jb+1));

%reading files and getting the desired correlation
'loading datasets'
    for i = 1:length(fileList)  % i, j, k are USED. pick l 
        i
        Name = fileList(i).name;
        load([dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/outlierRemoved/' ...
              fileList(i).name])
        gCount = size((dataSet.mat), 1);
        
        dataMat = 2.^ dataSet.mat;
        
        % 1.1 network 
        tempMat = corr(dataMat'); % pearson correlation
           
        tmin = min(min(tempMat));
        tmax = max(max(tempMat));
        tempMat = ((tempMat - tmin) / (tmax - tmin))*(2) -1;
        wholeData(i, :, :) = tempMat(ib:ie, jb:je);
        clear tempMat
        clear dataMat
        clear dataSet
        
    end 
    
    'datasets loaded'
    %getting the variance
    varmat = var(wholeData);
    size(varmat)
    
    fmat = zeros(gCount);
    fmat(ib:ie, jb:je) = varmat(1, :, :);
    
    svarmat = sparse(fmat);
    fileName = sprintf('%scorrVar_%s_%d_%d_%d_%d.mat', saveFolder,tissue, ib, ie, ...
                       jb, je);
    save(fileName, 'svarmat', '-v7.3')
end