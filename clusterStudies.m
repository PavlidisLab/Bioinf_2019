% in this file I study clusterings of the genes. How dense each
% cluster is and how connected the clusterings are with each
% other. My inputs are the clusterings and the networks (tissue
% networks). Questions I want to find the answer to are: 
% 1. How dense is the cluster? 
% 2. How is it "clustered"? ND distribution and cluster coefficient
% maybe? Does it have like, sub clusters or genes are smoothly
% clustered?
% 3. how is each cluster connected to the other clusters? 

% loading the clusterings for the tissues 
clear
load('~/data/general/GOdata_GPL570_04.mat')

addpath('~/codes/MATLAB/myCodes/general/')
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}

% loading the R clustering 
minModuleSize = 100
deepSplit = 2
linum = 2
% fid = fopen(['~/data/RDataFiles/' ...
%              'bd_bn_lr_lg_sm_clusterMat_minMod50_DS2.csv'])

% file = sprintf(['~/data/RDataFiles/' ...
%          'bd_bn_lr_lg_sm_clusterMat_MinMod%d_DS%d_Exp0.8_Ind0' ...
%                 '.05_4QDS.csv'], ...
%                minModuleSize, deepSplit)

file = sprintf(['~/data/RDataFiles/finalClusterings/' ...
                'bd_bn_lr_lg_sm_clusterMat_MinMod%d_DS2_Exp0' ...
                '.8_Ind0.10_FDR5e-5.csv'], minModuleSize)

fid = fopen(file)

book = textscan(fid, [repmat('%s', 1, 1)], ...
                'Headerlines', linum - 1, 'Delimiter', '\t')

gCount = length(book{1})
clusteringMat = zeros(gCount, 5);
for i = 1:gCount
    sib = book{1}(i);
    parts = strsplit(sib{1});
    tCell = parts(2:end);
    for j = 1:5
        clusteringMat(i, j) = str2num(tCell{j});
    end
end

% select the tissue
tissueID = 1;
tissue = tissues{tissueID}

% load the tissue network
% load(['~/networks/tissues/' tissue ['/' ...
%                     'binaryNet_0.5DSCountCorr_0.8Expr.mat'])
% net = tempNet;

% load(['~/networks/tissues/' tissue ['/' ...
%                     'binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat']])

load(['~/networks/tissues/' tissue ...
      '/binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
net = binNet;

% network density

(sum(sum(net)) * 2) / ((size(net, 1) * size(net, 1)) - ...
                       size(net, 1))

clustering = clusteringMat(:, tissueID);

% for each cluster
clusterEdgesTotalCount = 0;
clusterEdgeCount = zeros(1, max(clustering));
clusterMemberCount = zeros(1, max(clustering));
for i = 1:max(clustering)
    geneList = clustering == i;
    interNet = net(geneList, :);
    miniNet = interNet(:, geneList);
    clusterEdgesTotalCount = clusterEdgesTotalCount + sum(sum(miniNet));
    clear interNet;
    sum(sum(miniNet)) * 2/ ((size(miniNet, 1) * size(miniNet, 1)) - ...
                            size(miniNet, 1))
    clusterEdgeCount(i) = sum(sum(miniNet));
    clusterMemberCount(i) = sum(geneList);
end

clusterEdgesTotalCount / sum(sum(net))
% for each pair of clusters, how close or far a way from each other

%%%% TODO : FIX THE SUM FOR THE FIRST AND SECOND HERE. 
connectedFirstGenes = zeros(max(clustering), max(clustering));
connectedSecondGenes = zeros(max(clustering), max(clustering));
linkCounts = zeros(max(clustering), max(clustering));
externalLinkVSInternal = zeros(max(clustering), max(clustering));
for i = 1:max(clustering)
    i
    firstSet = clustering == i;
    for j = 1:max(clustering)
        secondSet = clustering == j;
        interNet = net(firstSet, :);
        miniNet = interNet(:, secondSet);
        % sum for the firstSet
        connectedFirstGenes(i, j) = sum(sum(miniNet') > 0);
        % sum for the secondNet
        connectedSecondGenes(i, j) = sum(sum(miniNet) > 0);
        linkCounts(i, j) = sum(sum(miniNet));
        externalLinkVSInternal(i, j) = linkCounts(i, j) / ...
            clusterEdgeCount(i);
        externalLinkVSInternal(j, i) = linkCounts(i, j) / clusterEdgeCount(j);
    end
end

logLinkCounts = log10(linkCounts + 1);
sib = externalLinkVSInternal;
sib(sib > 1) = 1;
heatmap(sib, [], [], [], 'Colorbar', true, 'Colormap', bone)

sum(sum(linkCounts)) / sum(sum(net))

% the hub study: how many do we have, are they local or general,
% where do they belong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IT IS A VERY FULL NETWORK, I DONT LIKE IT! 

nd = sum(net)  + sum(net');

sum(nd > 400)

hist(nd(nd > 0), 200)

sum(nd(nd > 0) < 100)

q = quantile(nd(nd > 0), 100);

ndThr = 90;

q(95)

hubs = nd > q(95);
sum(hubs)

% for each hub, which cluster does it belong to and what is its
% local hubness? 
hubLocal = zeros(1, length(clustering));
for i = 1:max(clustering)
    geneList = clustering == i;
    book = (hubs' + geneList) == 2;
    
    interNet = net(book, :);
    miniNet1 = interNet(:, geneList);
    
    interNet = net(:, book);
    miniNet2 = interNet(geneList, :);
    
    sumMiniNet = miniNet1 + miniNet2';
    
    sib = sum(sumMiniNet');
    length(sib)
    hubLocal(book) = sib;
end

sib = hubLocal ./ nd;
hist(sib(sib > 0))

%%%%%%
% comparison with GTEx and Illumina is remaining. Also, some sort
% of cluster studies. 
