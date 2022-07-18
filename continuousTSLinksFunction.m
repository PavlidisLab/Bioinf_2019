% this is the main function for running the code in parallel and
% saving the data. 
% where tissue is the tissue and s is the starting point for the
% correlation matrix and n is the number of edges to be checked
% (cause I am running in parallel, so I choose how long each run takes)
%function[] = continuousTSLinksFunction(tissue, s, n)
% s is the starting point (cause I get in limited number of edges
% each time) - n is the count of edges , I have ~57e+6 edges, for
% brain, for example. 1e+6 takes about half an hour to run. I give
% it 10mil at a time, I guess. 10e+6 takes 5 hours, and probably
% 10G of RAM. 

function[] = continuousTSLinksFunction(tissue, s, n)

    load(['~/data/general/' tissue 'ExpGenes.mat']);
    load('~/data/general/linkExprInfo/wholeGeneExpr.mat')
    load('~/data/general/linkExprInfo/dataSetProbeInf.mat')


    finalMat = zeros(n, 53);
    for i = 1:53
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
    % find the tissueInd

    tissueCell = {dataSetProbeInf(:).tissue};
    dsIndCell = strfind(tissueCell, tissue);
    allDSInd = find(~cellfun(@isempty, dsIndCell));

    tic
    allTissues = [1:53];
    space = zeros(n, 1);
    tissueInd = allDSInd;
    otherInd = allTissues(~ismember(allTissues, tissueInd));
    tc = length(tissueInd);
    oc = length(otherInd);
    ec = 947;

        for p = 1: n

            small = finalMat(p, tissueInd);

            big = finalMat(p, otherInd);

            total = ec + tc + oc;

            % for each step, I should count the space under the curve. 
            allseqBig = zeros(1, total);
            allseqSM = zeros(1, total);
            for i = 1:length(big)
                allseqBig(big(i)) = allseqBig(big(i)) + 1; 
            end

            for i = 1:length(small)
                 allseqSM(small(i)) = allseqSM(small(i)) + 1; 
            end

            emptyUnit = 1/ec;
            ts = 0;
            tsUnit = 1/tc; % 1/tc
            nontsUnit = 1/oc;
            ss = 0;
            steps = 0;
            tsEmptySurface = 0;
            while steps < total
                currentBig = allseqBig(1000 - steps);
                currentSM = allseqSM(1000 - steps);
                if(currentBig + currentSM)  == 0 
                    tsEmptySurface = tsEmptySurface + emptyUnit * ts;
                    steps = steps + 1;
                else
                    for j = 1:currentBig
                        ss = ss + (tsEmptySurface * nontsUnit);
                        steps = steps + 1;
                    end
                    for j = 1:currentSM
                        ts = ts + tsUnit;
                        steps = steps + 1;
                    end
                end
            end
            space(p)= ss;
        end

    fileName = sprintf('~/networks/continuousTSLinksResults/%s_%d_%d.mat', ...
                       tissue, s, n)
    save(fileName, 'space');

    toc
end