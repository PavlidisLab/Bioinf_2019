% this function is to get sum of the ranked correlation values for
% all genes in a tissue, given the tissue
%
% the code is part of the main code in foldernetwork.m

function[] = rankSum(tissue)


'inf'
dataFolder = ['~/data/array/' tissue '/matFiles/geneExprMatGemmaMapBER/' ...
              'outlierRemoved/'] 
'datafolder'

fileList = dir([dataFolder '*.mat'])

load([dataFolder fileList(1).name]);
gCount = size(dataSet.mat, 1);

Sum = zeros(gCount, gCount);

for i = 1: length(fileList)
    i
    load([dataFolder fileList(i).name]);
    dataMat = dataSet.mat.^2;
    
    tic
    sib = corr(dataMat');
    tempMat = reshape(tiedrank(sib(:)), gCount, gCount);
    toc
    
    Sum = Sum + tempMat;
    clear tempMat;
    
end

'saving file...'

fileName = sprintf('/space/grp/marjan/data/corrVar/rankSum_%s.mat', ...
                   tissue);
save(fileName, 'Sum', '-v7.3');

'file saved.'

end

