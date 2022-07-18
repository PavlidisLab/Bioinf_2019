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

function[] = correlationVarFunc02(ib, ie, jb, je, tissue)


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
        'ranking...'
        tic
        tempMat = reshape(tiedrank(tempMat(:)), gCount, gCount);
        
        toc
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
