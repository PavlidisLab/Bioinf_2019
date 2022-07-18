
function[Qarray] = corrQuantileReturn(dsName, corrArray)

load('~/data/general/corrBins_expStd.mat')

DSlist = {dsCorrQInf(:).name};
dsIndCell = strfind(DSlist, dsName);
dsInd = find(~cellfun(@isempty, dsIndCell));
Q = dsCorrQInf(dsInd).totalQ;

Qarray = zeros(1, length(corrArray));

for i = 1:length(corrArray)
    if (corrArray(i) < Q(1000))
        Qarray(i) = min(find(Q > corrArray(i)));
    else
        Qarray(i) = 1000;
    end
end

end
