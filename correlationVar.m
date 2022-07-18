% this is the main file for getting the correlation variation among
% different datasets for a tissue. The functions for correlation
% variation are correlationVarFunc.m and a newer version
% correlationVarFunc01.m which is much faster. 

% TODO : 
% back to this project: 
% for the having files, do the following: 
% merging the sd and getting the man for each tissue
% checking the overlap of edges with each of the selected edges
% design this for the random values. which means having 10
% correlation matrix of 10 random data. 


% general variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
dataFolder = '/space/grp/marjan/data/';

tissueList = [{'blood'}, {'skeletalMuscle'}, {'liver'}, {'lung'}];

GOtermsSubFolder = 'GOSimNet/'
gCount = 18377

% merging all the correlationVariances for one tissue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for l = 1:length(tissueList)
    l = 1
    corrVarMat = zeros(gCount);
    tissue = tissueList{l};
    fileList = dir([dataFolder 'corrVar/*' tissue '*.mat'])
    
    for i = 1:length(fileList)
        i
        fileList(i).name
        load([dataFolder 'corrVar/' fileList(i).name])
        corrVarMat = corrVarMat + full(svarmat);
        clear svarmat
    end
    
    fileName = sprintf('/space/grp/marjan/data/corrVar/corrVar_%s.mat', tissue);
    save(fileName, 'corrVarMat', '-v7.3')
    clear corrVarMat
end

%NOTICE : the variance matrix "might be" just the upper triangle
% getting the histogram of the variance 
tissueList = [{'blood'}, {'skeletalMuscle'}, {'liver'}, {'lung'}];
tissue = tissueList{4}
 load(sprintf('/space/grp/marjan/data/corrVar/corrVar_%s.mat', tissue));

% getting the upper triangle
corrTriUp = triu(corrVarMat, 1); 
corrArray = corrTriUp(:);
corrArray = corrArray(corrArray > 0);
corrArray = sqrt(corrArray);
size(corrArray)
hist(sqrt(corrArray), 40)

% having the i and j, getting the index of our desired value in the
% array 
a = ((j -1)*gCount + i) - sum([gCount-j+2:1:gCount]);

% getting the meanCorr for each tissue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for l = 2:length(tissueList)
    l = 4
    sumCorrMat = zeros(gCount);
    tissue = tissueList{l};
    fileList = dir([dataFolder tissue '/geneExprMat/gemmaProbeMap/clean/'])
    
    for i = 3:length(fileList)
        i
        fileList(i).name 
        
        load([dataFolder tissue '/geneExprMat/gemmaProbeMap/clean/' ...
              fileList(i).name])

        
        dataMat = 2.^ dataSet.mat;
        gCount = size(dataMat, 1);
        
        % 1.1 network 
        %tempMat = corr(dataMat'); % pearson correlation
        tempMat = corr(dataMat');
        
        sumCorrMat = sumCorrMat + tempMat;
        
        clear dataSet
        clear tempMat

    end
    
    fileName = sprintf('corrVar_%s.mat', tissue);
    save(fileName, 'sumCorrMat', '-v7.3')
end


meanCorrMat = sumCorrMat/(length(fileList)-2);
%NOTICE work with just the upper triangle
meanTriUp = triu(meanCorrMat, 1);
meanArray = meanTriUp(:);
meanArray(meanArray < 0) = 0.0001;
meanArray = meanArray(meanArray>0);
size(meanArray)
hist(sqrt(meanArray), 40)

p = 100000;
plot(meanArray(1:p), corrArray(1:p), '.', 'color', 'red')
xlabel('Correlation Mean')
ylabel('Correlation standard deviation')


% getting the aggregated network for each tissue to check the mean
% and variances 

%loading the aggNet
load('/space/grp/marjan/data/lung/AggNet/17agg.mat')

% the three methods of picking pairs 
% what I do not want is a standard deviation bigger than the mean -
% it must be half of the mean, considering I am getting the minimum
% mean at 0.3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. getting the index of the "aggNet" in the array 
aggNetCorrMeanTriUp = triu(aggNet, 1);
clear aggNetCorrMean
aggNetCorrMeanArray = aggNetCorrMeanTriUp(:);
clear aggNetCorrMeanTriUp
aggNetBin = aggNetCorrMeanArray;
clear aggNetCorrMeanArray
A = find(aggNetBin>0);

% 2. this is the below diagonal 
meanstd = meanArray./corrArray;
q = quantile(meanstd, (1 - 0.01))
newInPairs = find(meanstd>q);
belowDiagInd = newInPairs;
B = belowDiagInd;

plot(meanArray(newInPiars), corrArray(newInPiars), '.', 'color', 'red')
xlabel('Correlation Mean')
ylabel('Correlation standard deviation')

% 3. this is the low right most square. 
meanBin = (meanArray>0.4);
corrBin = (corrArray<0.2);
meanCorrIn = meanBin .* corrBin;
lowRightBin = meanCorrIn;
C = find(lowRightBin>0);
hold on
plot(meanArray(logical(lowRightBin)), corrArray(logical(lowRightBin)), ...
     '.', 'color', 'green')

%%% plotting the Venn Diagram
lA = length(A);
lB = length(B);
lC = length(C);
lAB = intersect(A, B);
lABC = length(intersect(lAB, C));
lAB = length(lAB);
lAC = length(intersect(A, C));
lBC = length(intersect(B, C));

vennX(floor([lA, lAB, lB, lBC, lC, lAC, lABC]/100), 0.1)

%TODO : save this diagram. 

%%%%%%%%%%%%%%%%%5
% getting the mean corr value for my selected pairs
aggNetCorrMean = aggNet .* meanCorrMat; 
aggNetCorrMeanTriUp = triu(aggNetCorrMean, 1);
clear aggNetCorrMean
aggNetCorrMeanArray = aggNetCorrMeanTriUp(:);
clear aggNetCorrMeanTriUp
aggNetCorrMeanArray(aggNetCorrMeanArray < 0) = 0.0001;
aggNetCorrMeanArray = aggNetCorrMeanArray(aggNetCorrMeanArray>0);
length(aggNetCorrMeanArray)

% getting the var value for my selected pairs
aggNetCorrVar = aggNet .* corrVarMat; 
aggNetCorrVarTriUp = triu(aggNetCorrVar, 1);
clear aggNetCorrVar
aggNetCorrVarArray = aggNetCorrVarTriUp(:);
clear aggNetCorrVarTriUp
aggNetCorrVarArray = aggNetCorrVarArray(aggNetCorrVarArray>0);
aggNetCorrVarArray = sqrt(aggNetCorrVarArray);
length(aggNetCorrVarArray)

pp = randi(length(aggNetCorrMeanArray), 1, 100000);
h = figure;
plot(aggNetCorrMeanArray(1:pp), aggNetCorrVarArray(1:pp), '.', ...
     'color', [80, 20, 20]./255)
xlabel('Correlation Mean')
ylabel('Correlation standard deviation')
title(['tissue : ' tissue])
%HERE : get the plot from other tissues

fileName = ['corrMeanVar_' tissue]
figFolder = [dataFolder 'general/']
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);

% so lets get the 0.005 quantile of maxmean-minstd pairs, and see
% how do they overlap with my simple maxmean pairs. 


% how much my three method of picking "top" pairs are matching with
% each other. 

% the other thing would be that, how much my "right down square "
% is sharing with the others 

% SCRATCH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for each tissue (the l iteration to be added)
for l = 1:length(tissueList)
    tissue = tissueList{l};
    fileList = dir([dataFolder tissue '/geneExprMat/gemmaProbeMap/clean/'])
    tissue
    %    corCell = cell(1, length(fileList)-2)
    wholeData = zeros(length(fileList)-2,gCount, gCount);
    ib = 1
    ie = 4500;
    jb = 1
    je = 4500;
    % getting the correlation list
    for i = 3: 6%length(fileList)  % i, j, k are USED. pick l 
        i
        Name = fileList(i).name;
        load([dataFolder tissue '/geneExprMat/gemmaProbeMap/clean/' ...
              fileList(i).name])
        
        dataMat = 2.^ dataSet.mat;
        gCount = size(dataMat, 1);
        
        % 1.1 network 
        %tempMat = corr(dataMat'); % pearson correlation
        tempMat = corr(dataMat');
        
        tmin = min(min(tempMat));
        tmax = max(max(tempMat));
        tempMat = ((tempMat - tmin) / (tmax - tmin))*(2) -1;
        %corCell{i - 2} = tempMat(ib:ie, jb:je);
        wholeData(i-2, :, :) = tempMat;
        clear tempMat
        clear dataMat
        clear dataSet
        
    end

    % getting the variance
        varmat = zeros(gCount);
       
        array = zeros(1:length(corCell));
        for i = 1:(ie-ib+1)%gCount
tic
            for j = 1+1:(je-jb+1)%gCount
                for k = 1:1:length(corCell)
                    array(k) = corCell{k}(i,j);
                    varmat(i,j) = var(array);
                end
            end
toc
        end
        
        % another way
        
        genesIn = 10000;
        inData = wholeData(1:4,1:genesIn, 1:genesIn);
        size(inData)
        
        tic
                 holu = var(inData);
        toc
        
        tic
        correlationVarFunc01(1, 10000, 1, 10000, 'lung')
        toc
        
        
        tic
        correlationVarFunc01(1, 10000, 10001, 18377, 'lung')
        toc
        
        
        tic
        correlationVarFunc01(10001, 18377, 1, 10000, 'lung')
        toc
        
        
        tic
        correlationVarFunc01(10001, 18377, 10001, 18377, 'lung')
        toc
        
        
        tic
        correlationVarFunc01(1, 10000, 1, 10000, 'skeletalMuscle')
        toc
        
        
        tic
        correlationVarFunc01(1, 10000, 10001, 18377, 'skeletalMuscle')
        toc
        
        
        tic
        correlationVarFunc01(10001, 18377, 1, 10000, 'skeletalMuscle')
        toc
        
        
        tic
        correlationVarFunc01(10001, 18377, 10001, 18377, 'skeletalMuscle')
        toc
        
        
        
%getting some big subplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        
h = figure;
        
aggNetFiles{1} = [dataFolder 'blood/AggNet/10agg.mat'] 
aggNetFiles{2} = [dataFolder 'skeletalMuscle/AggNet/9agg.mat']
aggNetFiles{3} = [dataFolder 'liver/AggNet/12agg.mat']
aggNetFiles{4} = [dataFolder 'lung/AggNet/17agg.mat']

for l = 1:length(tissueList)
    sumCorrMat = zeros(gCount);
    tissue = tissueList{l};
    fileList = dir([dataFolder tissue '/geneExprMat/gemmaProbeMap/clean/'])
    l
    %getting the meanCorr
    for i = 3:length(fileList)
        i
        fileList(i).name 
        
        load([dataFolder tissue '/geneExprMat/gemmaProbeMap/clean/' ...
              fileList(i).name])

        
        dataMat = 2.^ dataSet.mat;
        gCount = size(dataMat, 1);
        
        % 1.1 network 
        %tempMat = corr(dataMat'); % pearson correlation
        tempMat = corr(dataMat');
        
        sumCorrMat = sumCorrMat + tempMat;
        
        clear dataSet
        clear tempMat

    end
    meanCorrMat = sumCorrMat/(length(fileList)-2);
    
    %loading the correlation variation
    load(sprintf('/space/grp/marjan/data/corrVar/corrVar_%s.mat', ...
                 tissue));
    
    %loading the network
    load(aggNetFiles{l})
    'all data loaded'
    % getting the desired values
    aggNetCorrMean = aggNet .* meanCorrMat; 
    aggNetCorrMeanTriUp = triu(aggNetCorrMean, 1);
    clear aggNetCorrMean
    aggNetCorrMeanArray = aggNetCorrMeanTriUp(:);
    clear aggNetCorrMeanTriUp
    aggNetCorrMeanArray(aggNetCorrMeanArray < 0) = 0.0001;
    aggNetCorrMeanArray = aggNetCorrMeanArray(aggNetCorrMeanArray>0);
    length(aggNetCorrMeanArray)

    % getting the var value for my selected pairs
    aggNetCorrVar = aggNet .* corrVarMat; 
    aggNetCorrVarTriUp = triu(aggNetCorrVar, 1);
    clear aggNetCorrVar
    aggNetCorrVarArray = aggNetCorrVarTriUp(:);
    clear aggNetCorrVarTriUp
    aggNetCorrVarArray = aggNetCorrVarArray(aggNetCorrVarArray>0);
    aggNetCorrVarArray = sqrt(aggNetCorrVarArray);
    length(aggNetCorrVarArray)

    pp = randi(length(aggNetCorrMeanArray), 1, 100000);
    subplot(2,2, l)
    plot(aggNetCorrMeanArray(1:pp), aggNetCorrVarArray(1:pp), '.', ...
         'color', [80, 20, 20]./255)
    xlabel('Correlation Mean')
    ylabel('Correlation standard deviation')
    title(['tissue : ' tissue])
end


fileName = ['corrMeanVarAll']
figFolder = [dataFolder 'general/']
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);


%NOTICE work with just the upper triangle
