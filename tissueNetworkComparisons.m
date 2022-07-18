% here is the file to compare the tissue networks, this is where I
% do the clustering. For clustering, I use WGCNA in R, here I
% prepare the data for WGCNA. 
% 0. Just loading the networks
% 1. getting the expressed gens 
% 2. the ribosome link check
% 3. clustering
% 4. getting the genes for the R clustering

% just testing I guess!
clear
% loading the binary networks for each tissue 
load('~/networks/tissues/blood/binaryNet_0.5DSCountCorr_0.8Expr.mat')
bloodNet = tempNet;
load('~/networks/tissues/brain/binaryNet_0.5DSCountCorr_0.8Expr.mat')
brainNet = tempNet;
load('~/networks/tissues/lung/binaryNet_0.5DSCountCorr_0.8Expr.mat')
lungNet = tempNet;
load('~/networks/tissues/liver/binaryNet_0.5DSCountCorr_0.8Expr.mat')
liverNet = tempNet;
load('~/networks/tissues/skeletalMuscle/binaryNet_0.5DSCountCorr_0.8Expr.mat')
smNet = tempNet;

clear
% loading the binary networks for each tissue 
load('~/networks/tissues/blood/binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat')
bloodNet = binNet;
load('~/networks/tissues/brain/binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat')
brainNet = binNet;
load('~/networks/tissues/lung/binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat')
lungNet = binNet;
load('~/networks/tissues/liver/binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat')
liverNet = binNet;
load('~/networks/tissues/skeletalMuscle/binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat')
smNet = binNet;


clear sib
sib = ...
    bloodNet + ...
      brainNet + ...
      lungNet + ...
      liverNet + ...
      smNet;

pg = sum(sib) + sum(sib');
sum(pg > 0)
[a, b, c] = find(sib);
sum(c == 5)/length(c)
sum(c == 4)
sum(c == 3)

myNet = sparse(sib);

sibHist = hist(c, [1:1:5])%/length(c)
sibHist = sibHist / sum(sibHist)
% the five networks share so little with each other, but how don
% they do in more general structures? ribosome links, clustering
% and node degree distribution. We already know how many genes are
% shared between them. 
% 1. getting the gene counts, network density, count of edges and
% genes in the network, nd dist, connectivity of the network.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tissue = 'blood'
load(['~/data/general/tissueExpGenes/' tissue 'ExpGenes0.8.mat'])

net = bloodNet;
k = sum(net) + sum(net');

sum(k > 0)

[f, xi] = ksdensity(k(k > 0));

plot(xi, f)

% 2. the ribosome links
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% performance of the ribosome links in each tissue

load('~/data/general/ribosomeIDs.mat')
% how many of them are expressed in all tissues?  41 out of 48,
% total of 820 links. 
load('~/data/general/tissueExpGenes/allExp0.8.mat')
ribTemplate = sparse(zeros(18494, 18494));
ribTemplate(ribIDs, ribIDs) = 1;
ribNet = sib .*ribTemplate;

ribNet = smNet .* ribTemplate;

sum(sum(ribNet == 5))
[a, b, c] = find(ribNet);
ribHist = hist(c, [1:1:5])
ribHist = ribHist / sum(ribHist)
log10(ribHist)

netOverlap.ribosome = ribHist;
netOverlap.all = sibHist;

save('~/resultsAndFigures/tissueNetworkOverlap.mat', 'netOverlap')

% the overlap of the ribosome correlation is much higher between
% the tissues. How do I preenet it? 

ribHist(5)/820 % 38.66% of the ribosome links are presenet in all
               % the tissues networks

% 3. clustering and functional enrichment of the clusters TOP 
%%%%%%%%%%%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% before that, getting the genes in each tissue. 
%%% NOTICE: I am using three thresholds for networks. MODIFY THE
%%% INPUT FILE

clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
genesMat = zeros(18494, 5);
for t = 1:length(tissues)
    tissue = tissues{t};
    %    load(['~/networks/tissues/' tissue '/binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
     load(['~/networks/tissues/' tissue '/topSumNet_5_0.8Expr_Ind0.10.mat'])
    % gettign the edge count for the current net
    binNetT = binNet';
    fullCurrentNet = binNetT + binNet;
    genesIn = sum(fullCurrentNet) > 0;
    genesMat(:, t) = genesIn + 0;
end
%count of neighbours

% load('~/data/general/tissueExpGenes/bloodExpGenes0.8.mat')
% bloodSubNet = bloodNet;
% gCount = 18494;

clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
dsCount = [9, 12, 10, 15, 7]

% this function can modified to work with different networks, load
% the files properly! 
% loading the finaltable for the TS networks
%load('~/resultsAndFigures/TSlinks/finalTable.mat')

load('~/networks/GTEx/allNets.mat')

load('~/resultsAndFigures/TSlinks/finalTables/finalTable_CG23_FC2_FDR0012.mat')

for t = 3:3%length(GTExNet)%length(tissues)
    tissue = tissues{t}
    %    load(['~/networks/tissues/' tissue '/binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat'])
     % load(['~/networks/tissues/' tissue ...
     % '/binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
    
    tissue = finalTable(t).tissue
    binNet = finalTable(t).wholeNet;
    
    % load(['~/networks/tissues/' tissue '/' ...
    %                     'topSumNet_5_0.8Expr_Ind0.10.mat'])
    
    % GTEX
    %%%%
    % binNet = GTExNet(t).net;
    % tissue = ['GTEx_' GTExNet(t).tissue];
    % binNet = finalTable(t).wholeNet;
    % tissue = ['TSlinks_' finalTable(t).tissue];
    %%% End of GTEx and TS
    
    
    % tempMax = max(max(binNet));
    % binNet(binNet > 0) = binNet(binNet > 0) + (dsCount(t) - tempMax);
    % binNet = binNet .^2;
    
    % gettign the edge count for the current net (node degree)
    
    %    binNet = binNet ./ max(max(binNet));
    a1 = full(sum(binNet));
    a2 = full(sum(binNet'));
    K = a1 + a2;

    binNetT = binNet';
    fullCurrentNet = binNetT + binNet;
    
    % this is the Lij
    
    tic
    Lij = fullCurrentNet * fullCurrentNet;
    toc
    
    Lij(1:10, 1:10)

    % b(1) has a(1) neighbours. 
    gCount = 18494
    [a, b] = sort(K, 'descend');
    minK = zeros(gCount, gCount);
    tic
    for i = 1:gCount
        minK(b(i), :) = a(i);
        minK(:, b(i)) = a(i);
    end
    toc

    top1 =  Lij + fullCurrentNet;
    top2 = minK - fullCurrentNet + 1;
    tic
    top = top1 ./ top2;
    toc
    
    sum(sum(top > .5))
    % do the clustering of edges
    % genes with links:

    genesIn = sum(fullCurrentNet) > 0;

    subTop = top(genesIn, :);

    subTop = subTop(:, genesIn);

    sum(sum(subTop > .5))

    % [a, b, c] = find(subTop);

    % subTop = triu(subTop, -1);

    % save(['~/data/general/' tissue 'SubTop.mat'], 'subTop')


    % just saving it to a file. 
    addpath('~/codes/MATLAB/myCodes/general/')

    subTop = full(subTop);

    % it can not write it as a whole for whatever reason, so I am gonna
    % make it to 2k * 2k chunks and write it. boo. boo. boo. 

    subTopSize = size(subTop, 1)

    step = 2000;
    rs = 1
    cs = 1
    re = 0
    ce = 0

    while re < subTopSize
        if (re + step) > subTopSize 
            re = subTopSize
        else 
            re = re + step
        end
        while (ce) < subTopSize
            if (ce + step) > subTopSize
                ce = subTopSize
            else
                ce = ce + step
            end
            subTopMini = subTop(rs:re, cs:ce);
            file = sprintf('~/data/RDataFiles/temp/subTop%s_%d_%d_%d_%d_FDR5e-5.RData', tissue, rs, ...
                           re, cs, ce)
            
            % file = sprintf('~/data/RDataFiles/temp/TSnet_extended_subTop%s_%d_%d_%d_%d_FDR5e-5.RData', tissue, rs, ...
            %                re, cs, ce)
            
                   % file = sprintf('~/data/RDataFiles/temp/TSnet_subTop%s_%d_%d_%d_%d_FDR5e-5.RData', tissue, rs, ...
            %                re, cs, ce)


            % file = sprintf('~/data/RDataFiles/temp/subTop%s_%d_%d_%d_%d_FDR5e-5_ContinuousExtendedNet_power2.RData', tissue, rs, ...
            %                 re, cs, ce)
                            saveR(file, 'subTopMini')
            cs = ce  + 1

        end
        rs = re + 1
        cs = 1
        ce = 0;
    end
end

% 4. getting the genes for the R clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
dsCount = [9, 12, 10, 15, 7]

% this function can modified to work with different networks, load
% the files properly! 
% loading the finaltable for the TS networks
load('~/resultsAndFigures/TSlinks/finalTable.mat')

load('~/networks/GTEx/allNets.mat')
tissues = { GTExNet(14:17).tissue GTExNet(5:7).tissue}
tInd = [14:17, 5:7]

load('~/resultsAndFigures/TSlinks/finalTables/finalTable_CG23_FC2_FDR0012.mat')

genesIn = zeros(18494, 7);
for t = 1:length(tInd)%length(GTExNet)%length(tissues)
           %      tissue = tissues{t}
     % load(['~/networks/tissues/' tissue '/' ...
     %                     'topSumNet_5_0.8Expr_Ind0.10.mat'])

    %  load(['~/networks/tissues/' tissue '/binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat'])
     % load(['~/networks/tissues/' tissue ...
     % '/binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
    
    % GTEX
    %%%%
     binNet = GTExNet(tInd(t)).net;
    % tissue = ['GTEx_' GTExNet(t).tissue];
    % binNet = finalTable(t).wholeNet;
    % tissue = ['TSlinks_' finalTable(t).tissue];
     %%% End of GTEx and TS
    
    binNetT = binNet';
    fullCurrentNet = binNetT + binNet;
    genesIn(:, t) = sum(fullCurrentNet) > 0;
end

file = ['~/data/RDataFiles/general/bl_li_ln_sm_BCM_BCX_BFCX' ...
        '_GTEx_networkGenes.RData']

file = ['~/data/RDataFiles/general/bl_br_li_ln_sm_FDRe-' ...
        '5_Ind0.10_networkGenes.RData']

file = ['~/data/RDataFiles/general/' ...
        'bl_br_li_ln_sm_TSlinks_FDRe-' ...
        '5_Ind0.10_networkGenes.RData']

file = ['~/data/RDataFiles/general/' ...
        'bl_br_li_ln_sm_topSum_FDRe-' ...
        '5_Ind0.10_networkGenes.RData']

saveR(file, 'genesIn')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% >> making the pdist array from top. 

% 1. get the bottom triangle

botTri = full(subTop(logical(tril(ones(size(subTop)), -1))));

botTriDist = ones(1, length(botTri))' - botTri;

% the below doesn't work, cause there is just too many genes. 

% loading the clustering result for a brain dataSet;
book = csvread('~/data/brainClustering.csv');

% which functions or clusters do we find similar? 

% node degree distribution and hubs... do they overlap?

% node degree overlap or sth? do we have common hubs? 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
% loading the binary networks for each tissue 
load('~/networks/tissues/blood/binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat')
bloodNet = binNet;
bloodND = sum(bloodNet) + sum(bloodNet');

load('~/networks/tissues/brain/binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat')
brainNet = binNet;
brainND = sum(brainNet) + sum(brainNet');

load('~/networks/tissues/lung/binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat')
lungNet = binNet;
lungND = sum(lungNet) + sum(lungNet');

load('~/networks/tissues/liver/binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat')
liverNet = binNet;
liverND = sum(liverNet) + sum(liverNet');

load('~/networks/tissues/skeletalMuscle/binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat')
smNet = binNet;
smND = sum(smNet) + sum(smNet');

sumND = (bloodND>0) + (brainND>0) + (smND>0) + (lungND>0) + ...
        (liverND>0);
commonGenes = sumND == 5;

miniA = log10(full(bloodND(commonGenes)));
miniB = log10(full(brainND(commonGenes)));
plot(miniA, miniB, '.')
corr(miniA', miniB')

% set of genes which are highly correlated for the ND

% get a random 1000 genes. select and exchange 100 of them with
% another 100 random genes. If the correlation improved, keep
% % them. If the decreased, do not keep them. 

% % assuming 20% have higher corr

% initialList = rand(1, length(miniA)) > .7;
% firstCorr = corr(miniA(initialList)', miniB(initialList)')
% firstList = initialList;

% for i = 1:10000
%     secondList = firstList;

%     inInd = find(firstList);
%     outInd = find(~firstList);

%     repTo = datasample(inInd, 50);
%     repBy = datasample(outInd, 50);

%     secondList(repTo) = 0;
%     secondList(repBy) = 1;
%     sum(secondList);

%     a = miniA(secondList);
%     b = miniB(secondList);

%     secondCorr = corr(a', b', 'type','Spearman');

%     if(secondCorr > firstCorr)
%         firstList = secondList;
%         firstCorr = secondCorr;
%     end
% end

% plot(miniA(secondList), miniB(secondList), '.')

% myslist = secondList;

% allIn = ((myslist + secondList)  == 2);
% firstCorr = corr(miniA(allIn)', miniB(allIn)')




