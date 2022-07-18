% for each dataset, build the sparse coexpression networ with the
% -2% of edges and 2% of edges for correlation. Keep the average
% gene expression level for each gene as well as its STD among the
% samples, as well as its expression rank which is a value between
% 1:6
%
% The goal after this process is to explore the characteristics of
% the TS links VS the common links. 
% 1. getting the netStr file for each of the data sets
% 2. checking the percent of edges in every coexpression category for
% some of the networks
% 3. recording the correlation quantiles
% 4. saving the qExpr (the expression level separately )

% 1. getting the netStr file for each of the data sets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

dataFolder = '~/data/affyArray/tissues/'
folderName = dir(dataFolder)
counter = 0
gCount = 18494;
for t = 2:length(folderName)

    tissue = folderName(t).name;
    exprFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/']
    fileList = dir([exprFolder '*.mat'])


    % getting the single binary networks: 
    % count of gene quantiles
    qCount = 6;
    for i = 1: length(fileList)
        counter = counter + 1
        tic
        %loading new file
        '1'
        load([exprFolder fileList(i).name]);
        
        % getting the gene ranks and average
        % (save the avgExpr and stdExpr for the network)
        avgExpr = mean(dataSet.mat');
        stdExpr = std(dataSet.mat');
        
        kado = [0:1/6:1]
        Qs = quantile(avgExpr, kado)
        qSize= floor(gCount/qCount)
        qExpr = 2 .* ones(1, gCount);
        addVal = [1, 2, 2, 2, 2];
        for j = 2:qCount
            tempInd = avgExpr >= Qs(j);
            qExpr(tempInd) = qExpr(tempInd) + addVal(j-1);
        end
        
        sib = corr(dataSet.mat');
        %sib = sib - eye(gCount);
        
        '2'
        %getting the sparse raw values network
        upperSingle = sib(logical(triu(ones(size(sib)), 1)));
        Qplus = quantile(upperSingle, (1 - .1))
        Qminus = quantile(upperSingle, .06)
        sNet1 = sib > Qplus;
        sNet1(1:gCount+1:gCount*gCount) = 0; % putting the upper
                                             % triangle zero
        net = sib .* sNet1;
        sNet2 = sib < Qminus;
        net = net + sib .* sNet2;

        [r, c] = find(net);

        % getting the sparse corr p-values? or sth similar? 
        % for the .1 upper quantile and .6 lower quantile, we
        % divide it into 100. 
        q = 100;
        bookp = upperSingle(upperSingle > Qplus);
        qbookp = quantile(bookp, q);
        bookm = upperSingle(upperSingle < Qminus);
        qbookm = quantile(bookm, q);

        % find the quantile for each of the net values in the final
        % net.
        
        rankNet = zeros(gCount, gCount);
        qs = [qbookm , qbookp];
        tic
        for ii = 1:length(r)
            temp = find(qs < net(r(ii), c(ii)));
            % if (net(r(i), c(i))) < 0
            %     '<0'
            % end
            if length(temp) < 100 % the corr is negative
                rankNet(r(ii), c(ii)) = -(100 - length(temp));
            else % the corr is positive
                rankNet(r(ii), c(ii)) = temp(end) - 100;
            end
        end
        toc

        % getting the sparse expression values
        sparseExpNet = zeros(gCount, gCount);
        for ii=1:length(r)
            g1 = r(ii);
            g2 = c(ii);
            if(qExpr(g1) < qExpr(g2))
                sparseExpNet(g1,g2) = -(qExpr(g1) * qExpr(g2));            
            else
                sparseExpNet(g1,g2) = qExpr(g1) * qExpr(g2);                        
            end
        end
        
        netStr.expNet = sparse(sparseExpNet .* tril(ones(gCount,gCount),-1));
        netStr.coNt = sparse(net .* tril(ones(gCount,gCount),-1));
        netStr.rankNet = rankNet; 
        netStr.qExpr = qExpr;
        netStr.avgExpr = avgExpr;
        netStr.stdExpr = stdExpr;
        
        '3'
        %saving the single Network
        [a, b] = regexp(fileList(i).name ,'GSE[0123456789]+');
        Name = fileList(i).name(a:b)
        save(sprintf('~/networks/tissues/%s/singleNet/netStr_%s.mat', ...
                     tissue, Name), 'netStr', '-v7.3'); % i is included so I can follow up which experiments were in the i'th aggNet
        
        '4'
        %adding it to the aggregated ranked matrix we had so far and
        %get the results - tied rank takes 10m ?
        % aggMat = aggMat + reshape(tiedrank(sib(:)), gCount, gCount);
        % upperAgg = triu(aggMat, 1);
        % QAgg = quantile(upperAgg(:), (1 - 0.005));
        % aggNet = aggMat > QAgg;df
        % sparseAggNet = sparse(aggNet);
        
        % '5'
        % %saving the AggNetwork 
        % save(sprintf('%s%s/AggNet/%dagg.mat', dataFolder, tissue, (i)), 'aggNet');
        
        % clear dataSet
        toc
    end
end

% 2. checking the percent of edges in every coexpression category for
% some of the networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
fileList = dir('~/networks/liver/singleNet/netStr*.mat')
load(['~/networks/liver/singleNet/' fileList(2).name])
kado = full(netStr.expNet);
kado = netStr.expNet;
kado = abs(netStr.expNet);

sib = [ 1 2 3 5 7 9 4 6 10 14 18 9 15 21 27 25 35 45 49 63 81]
%sib = [1 2 3 4 6 9]
counts = zeros(1, length(sib));
meanCorr = zeros(1, length(sib));
for i = 1: length(sib)
    counts(i) = sum(sum(kado == sib(i)));
    % getting all the edges which have the current value
    tempNet = netStr.coNt(kado == sib(i));
    meanCorr(i) = mean(mean(tempNet(tempNet > 0)));
end

h = figure;
bar(1:length(sib), counts)
set(gca, 'Xtick', 1:length(sib), 'XTickLabel', sib)

% 3. recording the correlation quantiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

addpath '~/codes/MATLAB/myCodes/general/'
dataFolder = '~/data/affyArray/tissues/'
folderName = dir(dataFolder)
counter = 0
gCount = 18494;
for t = 3:length(folderName)

    tissue = folderName(t).name;
    exprFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/']
    fileList = dir([exprFolder '*.mat'])

    % getting the correlation
    for i = 1: length(fileList)
        counter = counter + 1
        tic
        %loading new file
        '1'
        load([exprFolder fileList(i).name]);
        
        % correlation
        stdExpr = std(dataSet.mat');
        
        sib = corr(dataSet.mat');
        %sib = sib - eye(gCount);
        
        '2'
        %getting the sparse raw values network
        upperSingle = sib(logical(triu(ones(size(sib)), 1)));
        Qs = quantile(upperSingle, 1000);
        
        qthds = upperSingle > Qs(1000);
        qthdsCorrs = upperSingle(qthds);
        qthdsQs = quantile(qthdsCorrs, 1000);
        
        % saving the file 
        min(qthdsQs)
        max(Qs)

        totalQ = [Qs qthdsQs];

        dsCorrQInf(counter).name = getGSEIDfromStr(fileList(i).name);
        dsCorrQInf(counter).tissue = tissue;
        dsCorrQInf(counter).totalQ = totalQ;
        dsCorrQInf(counter).stdExpr = stdExpr;
    end
end

save('~/data/general/corrBins_expStd.m', 'dsCorrQInf')

% 3.5 recording the correlation counts in each bin 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

addpath '~/codes/MATLAB/myCodes/general/'
dataFolder = '~/data/affyArray/tissues/'
folderName = dir(dataFolder)
counter = 0
gCount = 18494;
for t = 3:length(folderName)

    tissue = folderName(t).name;
    exprFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/']
    fileList = dir([exprFolder '*.mat'])


    % getting the correlation

    for i = 1: length(fileList)
        counter = counter + 1
        tic
        %loading new file
        '1'
        load([exprFolder fileList(i).name]);
                        
        sib = corr(dataSet.mat');
        sib = sib - eye(gCount);

        maxCorr = max(max(sib));
        minCorr = min(min(sib));
        
        '2'
        %getting the sparse raw values network
        upperSingle = sib(logical(triu(ones(size(sib)), 1)));
        Qs = quantile(upperSingle, 1000);
        
        qthds = upperSingle > Qs(1000);
        qthdsCorrs = upperSingle(qthds);
        qthdsQs = quantile(qthdsCorrs, 1000);
        
        % saving the file 
        min(qthdsQs)
        max(Qs)

        totalQ = [Qs qthdsQs];

        dsCorrQInf(counter).name = getGSEIDfromStr(fileList(i).name);
        dsCorrQInf(counter).tissue = tissue;
        dsCorrQInf(counter).totalQ = totalQ;
        dsCorrQInf(counter).stdExpr = stdExpr;
    end
end

save('~/data/general/corrBins_expStd.m', 'dsCorrQInf')

% 4. saving the qExpr (the expression level separately )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% saving the qExpr 

ts{1} = 'blood'
ts{2} = 'liver'
ts{3} = 'lung'
ts{4} = 'skeletalMuscle'
ts{5} = 'brain'

for t = 1: 5
    tissue = ts{t};
    fileList = dir(['~/networks/tissues/' tissue ['/singleNet/' ...
                        'netStr*.mat']])
    for f = 1:length(fileList)
f
        load(['~/networks/tissues/' tissue '/singleNet/' ...
              fileList(f).name]);
        expLevels = netStr.qExpr;
        save(['~/data/affyArray/tissues/' tissue ['matFiles/' ...
                            'rankedExpLevels/qExp'] fileList(f).name ...
                   ], 'expLevels')
    end
end


% getting the expression levels 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('~/codes/MATLAB/myCodes/general/')

ts{1} = 'blood'
ts{2} = 'liver'
ts{3} = 'lung'
ts{4} = 'skeletalMuscle'
ts{5} = 'brain'


c = 1;
for t = 1: 5
    tissue = ts{t};
    fileList = dir(['~/networks/tissues/' tissue ['/singleNet/' ...
                        'netStr*.mat']])
    for f = 1:length(fileList)
        f
        load(['~/networks/tissues/' tissue '/singleNet/' ...
              fileList(f).name]);

        name = getGSEIDfromStr(fileList(f).name)
        subNetStr(c).qExp = netStr.qExpr;
        subNetStr(c).avgExpr = netStr.avgExpr;
        subNetStr(c).stdExpr = netStr.stdExpr;
        subNetStr(c).tissue = tissue;
        subNetStr(c).name = name;        
        c = c + 1
    end
end

save('~/data/general/subNetStr.mat', 'subNetStr')

