% this function is a modification of the function:
% continuousTSLinksFunction() 
% this is for random testing, so it takes a list of datasets as the
% "tissue" (which are supposed to be a mixture of all other
% tissues) and the rest of datasets are random. The dataSetList
% should be a list of dataSet IDs. 
% s is the starting point (cause I get in limited number of edges
% each time) - n is the count of edges , I have ~57e+6 edges, for
% brain, for example. 1e+6 takes about half an hour to run. I give
% it 10mil at a time, I guess. 10e+6 takes 5 hours, and probably
% 10G of RAM. 

function[] = continuousTSLinksFunction_randomRep(dataSetList, s, n, ...
                                                 r, rs)

% get the expressed genes from the matrix. 
    load('~/data/general/tissueExpGenes/allDSExprInfo.mat')
    load('~/data/general/linkExprInfo/wholeGeneExpr.mat')
    load('~/data/general/linkExprInfo/dataSetProbeInf.mat')

    totalDSCount = length(dataSetProbeInf)
    % getting the expressed genes
    exprMat = zeros(18494, length(dataSetList));
    [a, dsExpr] = ismember(dataSetList, {allDSExprInfo(:).GSEID});
    for i = 1:length(dataSetList)
        exprMat(:, i) = allDSExprInfo(dsExpr(i)).exprGenes;
    end

    dsThr = ceil(.8 * length(dataSetList))
    sib = sum(exprMat');
    expGenesInd = sib >= dsThr;

    % building the finalmat, regardless of which datasets 
    finalMat = zeros(n, totalDSCount);
    for i = 1:totalDSCount
        i
        dsName = dataSetProbeInf(i).name;
        dsMat = wholeGeneExpr(:, dataSetProbeInf(i).sampleInd(1): ...
                              dataSetProbeInf(i).sampleInd(2));
        dsMat = dsMat(expGenesInd, :);
        %sib = corr(dsMat(1:100, :)');
        sib = corr(dsMat');
        upOnes = sib(logical(triu(ones(size(sib)), 1)));
        clear sib
        arr = upOnes(s:s+n-1);
        clear upOnes
        finalMat(:, i) = corrQuantileReturn(dataSetProbeInf(i).name, ...
                                           arr);

    end

    sib = finalMat < 500;
    finalMat(sib) = 500;
    clear sib

    %getting the IDs for the "tisssue" and "others"
    AllIDs = {dataSetProbeInf(:).name};
    [a, allDSInd] = ismember(dataSetList, AllIDs);
    totalDSCount = length(dataSetProbeInf);


    tic
    allInd = [1:totalDSCount];
    space = zeros(n, 1);
    tissueInd = allDSInd;
    otherInd = allInd(~ismember(allInd, tissueInd));

    for i = 1:length(tissueInd)
        for j = 1:length(otherInd)
            sib = finalMat(:, tissueInd(i)) - finalMat(:, otherInd(j));
            sib(sib < 0) = 0;
            space = space + sib;
        end
    end
    toc
    fileName = sprintf('~/networks/continuousTSLinksResults/%d_random%d_%d_%d_NegAdjusted.mat', ...
                        rs, s, n, r)
    save(fileName, 'space');
end


%     tc = length(tissueInd);
%     oc = length(otherInd);
%     ec = 1000 - totalDSCount;

%         for p = 1: n

%             small = finalMat(p, tissueInd);

%             big = finalMat(p, otherInd);

%             total = ec + tc + oc;

%             % for each step, I should count the space under the curve. 
%             allseqBig = zeros(1, total);
%             allseqSM = zeros(1, total);
%             for i = 1:length(big)
%                 allseqBig(big(i)) = allseqBig(big(i)) + 1; 
%             end

%             for i = 1:length(small)
%                  allseqSM(small(i)) = allseqSM(small(i)) + 1; 
%             end

%             emptyUnit = 1/ec;
%             ts = 0;
%             tsUnit = 1/tc; % 1/tc
%             nontsUnit = 1/oc;
%             ss = 0;
%             steps = 0;
%             tsEmptySurface = 0;
%             while steps < total
%                 currentBig = allseqBig(1000 - steps);
%                 currentSM = allseqSM(1000 - steps);
%                 if(currentBig + currentSM)  == 0 
%                     tsEmptySurface = tsEmptySurface + emptyUnit * ts;
%                     steps = steps + 1;
%                 else
%                     for j = 1:currentBig
%                         ss = ss + (tsEmptySurface * nontsUnit);
%                         steps = steps + 1;
%                     end
%                     for j = 1:currentSM
%                         ts = ts + tsUnit;
%                         steps = steps + 1;
%                     end
%                 end
%             end
%             space(p)= ss;
%         end

%     fileName = sprintf('~/networks/continuousTSLinksResults/random%d_%d_%d.mat', ...
%                         s, n, r)
%     save(fileName, 'space');
%     toc
% end