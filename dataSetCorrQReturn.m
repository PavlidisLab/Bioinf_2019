% this is a function to return the quantile vector of a dataset

function[Q] = dataSetCorrQReturn(dsName)

load('~/data/general/corrBins_expStd.mat')

DSlist = {dsCorrQInf(:).name};
dsIndCell = strfind(DSlist, dsName);
dsInd = find(~cellfun(@isempty, dsIndCell));
Q = dsCorrQInf(dsInd).totalQ;

end


