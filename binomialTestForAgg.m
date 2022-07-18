
% This is to find the Bernouli probability for observing the links
% in my networks. 

n = 9
bloodProb = log10(binopdf([1:n], n, 0.05));

n = 10
liverProb = log10(binopdf([1:n], n, 0.05));

n = 15
lungProb = log10(binopdf([1:n], n, 0.05));

n = 12
brainProb = log10(binopdf([1:n], n, 0.05));

n = 7
muscleProb = log10(binopdf([1:n], n, 0.05));
%%%

% the possible node count. Load the node count. 
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
tissueDSCount = [9, 12, 10, 15, 7]

% getting the total node count for each tissue to get the FDR

posLinkCount = zeros(1, 5);
for i = 1:length(tissues)
    tissue = tissues{i};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    posLinkCount(i) = sum(expGenesInd) * (sum(expGenesInd) -1) /2;
end

% I will have the count of false discovery in FDCount for each
% tissue

% I have four sets of binary networks: 0.1, 0.05, 0.01, 0.005. For
% each set I can select the FDR.  

%% Step 1: getting the network FDR and densities. 

indProbSet = [0.005, 0.01, 0.05, 0.1]

for p = 1 : 4
    indProb = indProbSet(p);
    for tID = 1: 5
        n = tissueDSCount(tID)
        myProb = 1 - (binocdf([1:n], n, indProb)); % probability of
                                                   % observing links with
                                                   % certain number of
                                                   % repeats in datasets
                                                   % given the initial
                                                   % probability for the
                                                   % individual networks
        FDCount = myProb * posLinkCount(tID) % count of expected links by
                                             % chance in my data

        % now getting the true count for each tissue at each level and FDR
        % for that tissue for each cutoff. 

        t = ceil(abs(log10(indProb)))

        if (t < 2)
            t = 2
        end

        tissue = tissues{tID}
        load(sprintf(['~/networks/tissues/%s/%s' ...
                      'SumNetworkExpThr0.8_%.*f.mat'], tissue, ...
                     tissue, t,indProb))
        FDRatCut = zeros(1, n);
        netDensity = zeros(1, n);
        for i = 1:n
            FDRatCut(i) = sum(FDCount(i:end))/sum(sum(sumNet >= i));
            netDensity(i) = sum(sum(sumNet >= i)) / posLinkCount(tID);
        end
        % write it in a table. 
        netAggInfo(tID).tissue = tissues{tID};
        netAggInfo(tID).binCut = indProb;
        netAggInfo(tID).FDR = FDRatCut;
        netAggInfo(tID).netDensity = netDensity;
    end
    aggInfs(p) = {netAggInfo};
end

% aggInfs(1) = {netAggInfo}; % .10
% aggInfs(2) = {netAggInfo}  % .05
% aggInfs(3) = {netAggInfo}  % .01
% aggInfs(4) = {netAggInfo}  % .005
% aggInfs(5) = {netAggInfo}  % .001

% Step 2: define an FDR cutoff and give me the network density for
% that FDR for each tissue - also the count of 

% just set it to the right netAggInfo:

for n=1:4
    netAggInfo = aggInfs{n}
    binCut = netAggInfo.binCut;
    fdrThrs = [5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3]
    
    l = length(fdrThrs)
    netDensities = zeros(l, 5);
    fdrs = zeros(l, 5);
    fdrInd = zeros(l, 5);

    for f = 1:length(fdrThrs)
        fdrThr = fdrThrs(f);
        for t = 1:5
            tissue = netAggInfo(t)
            sib = netAggInfo(t).FDR;
            smallerFDRInd = (min(find(sib <= fdrThr)))
            % getting the smaller distance for the FDR 
            % biggerFDRInd = smallerFDRInd - 1;
            % disS = abs(fdrThr - sib(smallerFDRInd))
            % disB = abs(fdrThr - sib(biggerFDRInd))
            % if(disS < disB)
            %     fdrInd(t) = smallerFDRInd;
            % else
            %     fdrInd(t) = biggerFDRInd;
            % end
            fdrInd(f, t) = smallerFDRInd;
            fdrs(f, t) = netAggInfo(t).FDR(fdrInd(f, t))
            netDensities(f, t) = netAggInfo(t).netDensity(fdrInd(f, t))
        end
    end
    fdrInd

    addpath('~/codes/MATLAB/myCodes/general')

    h = figure
    heatmap(log10(fdrs(1:end, :)), tissues, fdrThrs, [], 'TickAngle', 45, 'ShowAllTicks', true, 'Colorbar', true, 'colormap', ...
            'pink', 'GridLines', ':')
    title('FDRs for the aggregation of networks at different FDR bounds')

    % each individual thr has one of these tables and results. 
    figFolder = ['~/resultsAndFigures/networkAggregationResults/'];
    fileName = sprintf('FDRheatmap_IndividualNetThr_%.3f_corrected', binCut)
    print(h, '-depsc', [figFolder fileName '.eps']);
    print(h, '-dpdf', [figFolder fileName '.pdf']);

    h = figure
    heatmap(netDensities(1:end, :), tissues, fdrThrs, [], 'TickAngle', 45, 'ShowAllTicks', true, 'Colorbar', true, 'colormap', ...
            'hot', 'GridLines', ':')
    title(['The density of the aggregated networks at different FDR ' ...
           'bounds'])
    figFolder = ['~/resultsAndFigures/networkAggregationResults/' ...
                 '']
    fileName = sprintf('netDensityHeatmap_IndividualNetThr_%.3f', binCut)
    print(h, '-depsc', [figFolder fileName '.eps']);
    print(h, '-dpdf', [figFolder fileName '.pdf']);
end

% print it to the file. 
% just printing. how to print the log power numbers. 

for n = 1:4
    netAggInfo = aggInfs{n}
    binCut = netAggInfo.binCut;
    fdrThrs = [5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3]
    
    l = length(fdrThrs)
    netDensities = zeros(l, 5);
    fdrs = zeros(l, 5);
    fdrInd = zeros(l, 5);

    for f = 1:length(fdrThrs)
        fdrThr = fdrThrs(f);
        for t = 1:5
            tissue = netAggInfo(t)
            sib = netAggInfo(t).FDR;
            smallerFDRInd = (min(find(sib <= fdrThr)))
            % getting the smaller distance for the FDR 
            % biggerFDRInd = smallerFDRInd - 1;
            % disS = abs(fdrThr - sib(smallerFDRInd))
            % disB = abs(fdrThr - sib(biggerFDRInd))
            % if(disS < disB)
            %     fdrInd(t) = smallerFDRInd;
            % else
            %     fdrInd(t) = biggerFDRInd;
            % end
            fdrInd(f, t) = smallerFDRInd;
            fdrs(f, t) = netAggInfo(t).FDR(fdrInd(f, t))
            netDensities(f, t) = netAggInfo(t).netDensity(fdrInd(f, t))
        end
    end
    
    % print them in the file
    figFolder = ['~/resultsAndFigures/networkAggregationResults/'];
    fileName = sprintf('FDRtable_IndividualNetThr_%.3f.txt', binCut)
    file = [figFolder fileName]
    fid = fopen(file, 'w')
    fprintf(fid, ['FDR\t' repmat('%s\t', 1, 5) '\n'], tissues{:})
    for j = 1:length(fdrThrs)
        fprintf(fid, [repmat('%E\t', 1, 6) '\n'], fdrThrs(j), fdrs(j,:))
    end
    fclose(fid)
    
    figFolder = ['~/resultsAndFigures/networkAggregationResults/' ...
                 '']
    fileName = sprintf('netDensityTable_IndividualNetThr_%.3f.txt', ...
                       binCut)
    file = [figFolder fileName]
    fid = fopen(file, 'w')
    fprintf(fid, ['FDR\t' repmat('%s\t', 1, 5) '\n'], tissues{:})
    for j = 1:length(fdrThrs)
        fprintf(fid, ['%E\t' repmat('%.4f\t', 1, 5) '\n'], fdrThrs(j), netDensities(j,:))
    end
    fclose(fid)
    
    figFolder = ['~/resultsAndFigures/networkAggregationResults/' ...
                 '']
    fileName = sprintf('indNetCountTable_IndividualNetThr_%.3f.txt', ...
                       binCut)
    file = [figFolder fileName]
    fid = fopen(file, 'w')
    fprintf(fid, ['FDR\t' repmat('%s\t', 1, 5) '\n'], tissues{:})
    for j = 1:length(fdrThrs)
        fprintf(fid, ['%E\t' repmat('%d\t', 1, 5) '\n'], fdrThrs(j), fdrInd(j,:))
    end
    fclose(fid)
end

% I want to keep some informations in my table:

fileName = 'tissueNetworkAggregatingSettings.txt'
file = ['~/data/general/' fileName]
fid = fopen(file, 'w')
fprintf(fid, [])

% save this as a table 
% one setting: networks 0.05 and the FDR at 5e-4. 
% now getting the networks - comparing them with the old ones. 
for tID = 1: 5
    n = tissueDSCount(tID)
    myProb = 1 - (binocdf([1:n], n, 0.05));
    FDCount = myProb * posLinkCount(tID) % count of expected links by
                                         % chance in my data

    % now getting the true count for each tissue at each level and FDR
    % for that tissue for each cutoff. 

    tissue = tissues{tID}
    load(['~/networks/tissues/' tissue ['/' ...
                        tissue 'SumNetworkExpThr0.8_0.05.mat']])

    sib = sumNet >= fdrInd(tID);

    load(['~/networks/tissues/' tissue ['/' ...
                        'binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat'])

end

% Load different networks based on different thrs for each tissue
% and compare the links? , select an individual threshold. I
% prefere consistency with lower threshold and higher FDR rather
% than individual spikes with lower number of repeats. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

clear sumNet
tissue = 'skeletalMuscle'
t = 2
indProb = 0.1
load(sprintf(['~/networks/tissues/%s/%s' ...
              'SumNetworkExpThr0.8_%.*f.mat'], tissue, ...
             tissue, t,indProb))

sumNet005 = sumNet;
sumNet01 = sumNet;
sumNet05 = sumNet;
sumNet10 = sumNet;

net1 = sumNet005 >= 2;
net2 = sumNet01 >= 1;
net3 = sumNet05 >= 5;
net4 = sumNet10 >= 6;

sum(sum(net4))
net = net2;
sum(sum(net))
book = net4 + net;
sum(sum(book == 2)) 
sum(sum(book == 2)) /sum(sum(net4))

holu = book == 2;

final = final + holu;
final = holu;

% networks are finalized. I can get the overlapping edges. I go for
% the FRD 5e-5 for individual .1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

clear sumNet
tissue = 'skeletalMuscle'
t = 2
indProb = 0.1
load(sprintf(['~/networks/tissues/%s/%s' ...
              'SumNetworkExpThr0.8_%.*f.mat'], tissue, ...
             tissue, t,indProb))

smNet = sumNet; >= 1;
blNet = sumNet; >= 1;
brNet = sumNet; >= 1;
lnNet = sumNet; >= 1;
lvNet = sumNet; >= 1;

sib = smNet + blNet + brNet + lnNet + lvNet;

sib = sib > 1;

% I added this part later, saving another version of the network
% with the two more significant links. 
tissue = 'skeletalMuscle'
binNet = zeros(18494, 18494);
binNet = binNet + (smNet>=7);

save(['~/networks/tissues/' tissue ['/topSumNet_5_0.8Expr_Ind0.10.mat']], ...
     'binNet', '-v7.3')


% I have total of 2,6063 edges in Sib. For those edges, I want to
% get it modeled. Also, I have the TS of the rank test for those
% tissues. 

tissue = 'skeletalMuscle'
load(['~/networks/tissues/' tissue ['/binaryNet_FDR5e-' ...
                    '5_0.8Expr_Ind0.10.mat']])

smBNet = binNet;
blBNet = binNet;
brBNet = binNet;
lnBNet = binNet;
lvBNet = binNet;

holu = smBNet + blBNet + brBNet + lnBNet + lvBNet;

holu = holu > 0;
posLinks = holu;
save('~/networks/linksPresentInOneTissueMin.mat', ...
     'posLinks');

testNet = holu + sib ;

sum(sum(testNet > 0))

temp = testNet > 0;
testNet = temp;
% I am keeping this as list of links which are variable. 
save('~/networks/linksPresenetIn_2Top0.005_OR_Tissues.mat', 'testNet');


