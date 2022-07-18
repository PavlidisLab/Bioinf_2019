
%given the address of a dataset, this file gets the correlation
%network, having the correlation method

function[sparseSingleNet] = getNetwork01(file, correlation, logt, q)

load(file);
gCount = size(dataSet.mat, 1);
sampleCount = size(dataSet.mat, 2);
if(logt)
    dataSet.mat = 2.^dataSet.mat;
end
if(strcmp(correlation, 'SP'))
    sib = corr(dataSet.mat', 'type', 'Spearman')
else
    sib = corr(dataSet.mat');
end

sib = sib - eye(gCount);

upperSingle = triu(sib, 1);
QSingle(i-2) = quantile(upperSingle(:), (1 - q))
singleNet = sib > QSingle(i-2);
sparseSingleNet = sparse(singleNet);

end