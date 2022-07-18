% in this file I work on getting some statistics of my aggregated
% tissue networks obtained from the homotest and the average rank
% comparison.

% 1. the homotest networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
dataFolder = '/space/grp/marjan/data/corrHomoTest/'
tissueList = {'blood', 'liver', 'lung', 'skeletalMuscle', 'brain'}
tissueNum = [2, 8, 4, 32, 16]
tissueEdgeCount = zeros(size(tissueNum));
gCount = 18494
totalEdgeCount = (gCount * (gCount-1))/2;
homoNetCorrDensity = struct;

% 1.1. what percentage of the edges are present for each tissue
% 1.2. the distribution of the correlation values for each tissue
% get the values and plot them in the local, it is nicer. put the
% 1.1 and 1.2 in one plot 

%loading the common network
load(['/space/grp/marjan/data/corrHomoTest/commonNetchi005.mat'])
tissueEdgeInfo = struct

for i = 1:length(tissueList)
    tissue = tissueList{i}
    file = dir([dataFolder tissue '/homoNetchi005*Rep.mat'])

    load([dataFolder tissue '/' file.name])

    [a, b, c] = find(net);
    tissueEdgeCount(i) = length(a);
    e = length(a)/ totalEdgeCount;

    % getting the all edge density
    X = randi(length(a), [1, floor(length(a)/10)]);

    tic
    [f, xi] = ksdensity(c(X));
    toc
    allEdgeDensity = [xi; f];
    
    sib = commonNet == tissueNum(i);
    sib2 = net(sib);
    X = randi(length(sib2), [1, floor(length(sib2)/10)]);
    [f, xi] = ksdensity(sib2(X));
    specificEdgeDensity = [xi; f];
    
    
    sib = commonNet == 62;
    sib2 = net(sib);
    X = randi(length(sib2), [1, floor(length(sib2)/10)]);
    [f, xi] = ksdensity(sib2(X));
    commonEdgeDensity = [xi; f];
    
    homoNetCorrDensity(i).tissue = tissue;
    homoNetCorrDensity(i).edgePer = e;    
    homoNetCorrDensity(i).allEdgeDensity = allEdgeDensity;
    homoNetCorrDensity(i).specificEdegeDensity = specificEdgeDensity;
    homoNetCorrDensity(i).commonEdgeDensity = commonEdgeDensity;
    
    %filling the tissue struct
    tissueEdgeInfo(i).tissue = tissueList{i};
    tissueEdgeInfo(i).specificEdge = sum(sum(commonNet == ...
                                         tissueNum(i)));
    tissueEdgeInfo(i).edgeCount = tissueEdgeCount(i);
    tissueEdgeInfo(i).commonEdge = sum(sum(commonNet == 62));
    
end

for i =1:5
    plot(homoNetCorrDensity(i).xi, homoNetCorrDensity(i).f )
end
file = ['~/data/array/homoNetCorrDensity.mat']
save(file, 'homoNetCorrDensity')
file = ['~/data/array/homoNetTissuEdgeInfo.mat']
save(file, 'tissueEdgeInfo')

% 1.3 the bar plot of the precentage of edges in common, the total
% count of edges and the specific edges. For each tissue, get the
% distribution of the correlation for  the two groups. 
load(['/space/grp/marjan/data/corrHomoTest/commonNetchi005.mat'])

tissueEdgeInfo = struct

% the tissue numbers are in the file : correHomoResultMergin 
for i = 1:length(tissueList)
    tissueEdgeInfo(i).tissue = tissueList{i};
    tissueEdgeInfo(i).specificEdge = sum(sum(commonNet == ...
                                         tissueNum(i)));
    tissueEdgeInfo(i).edgeCount = tissueEdgeCount(i);
    tissueEdgeInfo(i).commonEdge = sum(sum(commonNet == 62));
    
    % load the network and get the distribution of the specific and
    % common edges and other edges 
end

% 2. comparison between the two methods 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2.1 what percentage of the average rank are in the homotest
% what is their average correlation

% get the top 005 * total from the selected corrs I have 
clear
dataFolder = '/space/grp/marjan/data/corrHomoTest/'
tissueList = {'blood', 'liver', 'lung', 'skeletalMuscle', 'brain'}
gCount = 18494
totalEdgeCount = (gCount * (gCount-1))/2;
edgeCount = floor(totalEdgeCount * .005);
minTissueCorr = zeros(1, 5)
for i = 1:length(tissueList)
    i
    homoBinMat = (zeros(gCount));
    tissue = tissueList{i}
    file = dir([dataFolder tissue '/homoNetchi005*Rep.mat'])

    load([dataFolder tissue '/' file.name])
    
    [a, b, c] = find(net);
    tic
    [d, e] = sort(c, 'descend');
    toc
    newA = a(e(1:edgeCount));
    newB = b(e(1:edgeCount));
    minTissueCorr(i) = c(e(edgeCount));
    tic
    for j = 1:length(newA)
        homoBinMat(newA(j), newB(j)) = 1;
    end
    toc
    shomoBinMat = sparse(homoBinMat);
    file = [dataFolder tissue '/homoBinMat005_' tissue '.mat']
    save(file, 'shomoBinMat', '-v7.3')
end

% for each tissue, load the networks and count the number of edges
clear
dataFolder = '/space/grp/marjan/data/corrHomoTest/'
tissueList = {'blood', 'liver', 'lung', 'skeletalMuscle', 'brain'}
gCount = 18494
totalEdgeCount = (gCount * (gCount-1))/2;
overlap = zeros(1, 5);
overlapPer = zeros(1, 5);

for i = 1:length(tissueList)
    tissue = tissueList{i};

    % the average net
    file = dir(['~/networks/' tissue '/sparse*.mat'])
    load(['~/networks/' tissue '/' file.name])

    % the homo net
    file = dir([dataFolder tissue '/homoBin*01*.mat'])
    load([dataFolder tissue '/' file.name])
    tot = shomoBinMat + sbinAggNet;
    overlap(i) = sum(sum(tot == 2));
end
overlapPer = overlap ./ edgeCount

%%%%%%%%%%%%%%%%%%
clear
tissueList = {'blood', 'liver', 'lung', 'skeletalMuscle', 'brain'}
tissueNum = [2, 8, 4, 32, 16];
gCount = 18494
totalEdgeCount = (gCount * (gCount-1))/2;
commonNet = zeros(gCount, gCount);

for i = 1:length(tissueList)
    tissue = tissueList{i}
    file = dir(['~/networks/' tissue '/sparse*.mat'])
    load(['~/networks/' tissue '/' file.name])
    sbinAggNet = sbinAggNet .* tissueNum(i);
    commonNet = commonNet + sbinAggNet;
    
end

scommonNet = sparse(commonNet);
% save the commonNet 
file = ['~/networks/rankedAvgCommonNet.mat']
save(file, 'scommonNet', '-v7.3')

%%%%%
% loading the common net

clear
tissueList = {'blood', 'liver', 'lung', 'skeletalMuscle', 'brain'}
tissueNum = [2, 8, 4, 32, 16];
gCount = 18494
totalEdgeCount = (gCount * (gCount-1))/2;
commonNet = zeros(gCount, gCount);
load('~/networks/rankedAvgCommonNet.mat')
tissueEdgeInfo = struct
edgeCount = floor(totalEdgeCount * .005);
% the tissue numbers are in the file : correHomoResultMergin 
for i = 1:length(tissueList)
    tissueEdgeInfo(i).tissue = tissueList{i};
    tissueEdgeInfo(i).specificEdge = sum(sum(scommonNet == ...
                                         tissueNum(i)));
        tissueEdgeInfo(i).commonEdge = sum(sum(scommonNet == 62));
    tissueEdgeInfo(i).edgeCount = edgeCount - tissueEdgeInfo(i).commonEdge ...
        - tissueEdgeInfo(i).specificEdge;
end


file = ['~/data/array/rankedAvgNetTissuEdgeInfo.mat']
save(file, 'tissueEdgeInfo')

% getting the same two values for the homoNet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
clear
dataFolder = '/space/grp/marjan/data/corrHomoTest/'
tissueList = {'blood', 'liver', 'lung', 'skeletalMuscle', 'brain'}
tissueNum = [2, 8, 4, 32, 16];
gCount = 18494
totalEdgeCount = (gCount * (gCount-1))/2;
commonNet = zeros(gCount, gCount);

for i = 1:length(tissueList)
    tissue = tissueList{i}
    file = dir([dataFolder tissue '/homoBinMat005_*.mat'])
    load([dataFolder tissue '/' file.name])
    sum(sum(homoBinMat))
    shomoBinMat = shomoBinMat .* tissueNum(i);
    commonNet = commonNet + shomoBinMat;
end

scommonNet = sparse(commonNet);
% save the commonNet 
file = ['~/networks/homoNetTop005CommonNet.mat']
save(file, 'scommonNet', '-v7.3')

%%%%%
% loading the common net

clear
tissueList = {'blood', 'liver', 'lung', 'skeletalMuscle', 'brain'}
tissueNum = [2, 8, 4, 32, 16];
gCount = 18494
totalEdgeCount = (gCount * (gCount-1))/2;
commonNet = zeros(gCount, gCount);
load('~/networks/homoNetTop005CommonNet.mat')
tissueEdgeInfo = struct
edgeCount = floor(totalEdgeCount * .005);
% the tissue numbers are in the file : correHomoResultMergin 
for i = 1:length(tissueList)
    tissueEdgeInfo(i).tissue = tissueList{i};
    tissueEdgeInfo(i).specificEdge = sum(sum(scommonNet == ...
                                         tissueNum(i)));
        tissueEdgeInfo(i).commonEdge = sum(sum(scommonNet == 62));
    tissueEdgeInfo(i).edgeCount = edgeCount - tissueEdgeInfo(i).commonEdge ...
        - tissueEdgeInfo(i).specificEdge;
end


file = ['~/data/array/homoNetTop005TissuEdgeInfo.mat']
save(file, 'tissueEdgeInfo')
