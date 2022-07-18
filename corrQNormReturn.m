% dsName: name of the dataset
% corrArray: the list of correlation values
% for the given dataset, this function returns the normalized
% ranked values of the correlation values in the corrArray

function[Qarray] = corrQNormReturn(dsName, corrArray)

load('~/data/general/correlationNorm/dsQNormInf.mat')
load('~/data/general/correlationNorm/correlationQuantiles.mat')
load('~/data/general/correlationNorm/representedValues.mat')

DSlist = {dsQNormInf(:).name};
dsIndCell = strfind(DSlist, dsName);
dsInd = find(~cellfun(@isempty, dsIndCell));
Q = qsamples(:, dsInd);

Qarray = zeros(1, length(corrArray));

for i = 1:length(corrArray)
    if (corrArray(i) < Q(1000))
        index = min(find(Q > corrArray(i)));
        Qarray(i) = repValues(index);
    else
        Qarray(i) = repValues(1000);
    end
end
end
