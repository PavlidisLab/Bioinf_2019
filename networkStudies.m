% This is to study the networks only. 
% comparing the hubs for the tissue networks and TS links.
% 0. Count of edges, count of nodes, starting genes, connectivity,
% density
% .0.5 node degree changes 
% 1. node degree and hubs
% 2. study of the hubs for the ts links
% 3. study of the TFs for the links.
% 4. get the common edge info
% 5. HKE studies for the tissue networks
% 6. The GTEx overlap of the networks

% 0. network info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
tissue = tissues{5};
expThr = '0.8'
%ns
for t = 1:5
    tissue = tissues{t}
    ns(t).tissue = tissue;
    
    % Pick one of the sections below: 
    
    load( ['~/networks/tissues/' tissue '/' ...
           'binaryNet_FDR5e-5_' expThr 'Expr_Ind0.10.mat'])
    sum(sum(binNet))
    ns(t).net = binNet;
    
    % load(['~/networks/tissues/' tissue '/' ...
    %      'topSumNet_5_0.8Expr_Ind0.10.mat'])
    % sum(sum(binNet > 0))
    % ns(t).net = binNet > 0;
    
    temp = sum(binNet);
    ns(t).nd = temp + sum(binNet');
    sum(ns(t).nd > 0)
    % figure
    % hist(ns(t).nd(ns(t).nd > 0), 40);
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    ns(t).expGenes = expGenesInd;
end

t = 5
sum(ns(t).nd > 0)
sum(ns(t).expGenes)

FC = '3'
FDR = '0012'
load(sprintf(['~/resultsAndFigures/TSlinks/finalTables/' ...
              'finalTable_CG23_FC%slog_FDR%s.mat'], FC, FDR))

% get the expGenes printed in differentFiles for the diagram. 
for t = 1:5
    tissue = tissues{t}
    m = find(ns(t).nd);
    dlmwrite(['~/resultsAndFigures/tissueGenesIndND_' tissue '.txt'], m')
    m = find(ns(t).expGenes);
    dlmwrite(['~/resultsAndFigures/tissueGenesInd_' tissue '.txt'], ...
             m')
    tissue = finalTable(t).tissue;
    nd = sum(finalTable(t).wholeNet') + sum(finalTable(t).wholeNet);
    m = find(nd > 0);
    dlmwrite(['~/resultsAndFigures/tissueGenesInd_TSNet_' tissue '.txt'], ...
             m')
end

% printing the information on the file
file = ['~/resultsAndFigures/tissueNetworkStudies/' ...
        'generalNetInfo.txt']
fid = fopen(file, 'w')

fprintf(fid, ['attribute\t blood\t brain\t liver\t lung\t skeletalMuscle ' ...
              '\n'])
for t = 1:5
    temp = sum(ns(t).expGenes);
    ns(t).tGC = temp;
    ns(t).nGC = sum(ns(t).nd > 0);
    ns(t).ec = sum(sum(ns(t).net))
    ns(t).netD = sum(sum(ns(t).net)) / (temp*(temp-1)/2);
end

fprintf(fid, ['%s\t' repmat('%d\t', 1, 5) '\n'],...
            'Tissue Genes', [ns(:).tGC])
fprintf(fid, ['%s\t' repmat('%d\t', 1, 5) '\n'],...
            'Network Genes', full([ns(:).nGC])')
fprintf(fid, ['%s\t' repmat('%d\t', 1, 5) '\n'],...
            'Edge Count', full([ns(:).ec])')
fprintf(fid, ['%s\t' repmat('%.3f\t', 1, 5) '\n'],...
            'Network Density', full([ns(:).netD])')

fclose(fid)

% gettign the density distribution and histogram of the links. 
for t = 1:5
    [f, xi] = kds
end

save('~/resultsAndFigures/tissueNetworkStudies/networksInfo_extended.mat', ...
     'ns')

% .0.5 node degree changes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
ND = zeros(18494, 5);

for t = 1:length(tissues)
    tissue = tissues{t};
    load(['~/networks/tissues/' tissue ['/binaryNet_FDR5e-' ...
                        '5_0.8Expr_Ind0.10.mat']])
    nd = sum(binNet) + sum(binNet');
    ND(:, t) = nd;
end

load('~/resultsAndFigures/TSlinks/finalTables/finalTable_CG13_FC3log_FDR0012.mat')
tsND = zeros(18494, 5);
for t = 1:5
    net = finalTable(t).wholeNet;
    sum(sum(net))
    nd = sum(net) + sum(net');
    tsND(:, t) = nd;
end

pND = zeros(18494, 5);
for t = 1:5
    cg = finalTable(t).cgNet;
    not = finalTable(t).noTSNet;
    holu = (cg + not) == 2;
    sum(sum(holu))
    nd = sum(holu) + sum(holu');
    pND(:, t) = nd;
end

allND = [ND tsND pND];
sib = sum(allND'); 
activeGenes = find(sib > 0);

allND = allND(sib > 0, :);

% myGenes 
load(['~/resultsAndFigures/tissueNetworkStudies/' ...
      'the13AllTissueGenes.mat'])
[a, b] = ismember(tempInds, activeGenes);

activeAllExpGenes = tempInds(a);

allND = allND(activeAllExpGenes, :);

rankND = zeros(size(allND));
for t = 1:15
    k = allND(:, t) > 0;
    rankND(k, t) = tiedrank(allND(k), t);
end


t = 1
gin = (tsND(:, t) + ND(:, t)) > 0;
h = figure
plot((tsND(gin,t)), (ND(gin,t)), '.')

cg1 = finalTable(1).cgNet;
not1 = finalTable(1).noTSNet;

holu = cg1 + not1;
holu = holu == 2;
holuND = sum(holu') + sum(holu);
t = 1
plot(holuND(gin), ND(gin, t), '.')

[a, b, c] = find(holu == 2);

[a(1:10), b(1:10)]

% 1. node degree and hubs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1.1 getting the hubs
% ---------------------

clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
tissue = tissues{5};
expThr = '0.8'
%ns
for t = 1:5
    tissue = tissues{t}
    ns(t).tissue = tissue;
    load( ['~/networks/tissues/' tissue '/' ...
           'binaryNet_FDR5e-5_' expThr 'Expr_Ind0.10.mat'])
    sum(sum(binNet))
    ns(t).net = binNet;
    temp = sum(binNet);
    ns(t).nd = temp + sum(binNet');
    sum(ns(t).nd > 0)
    figure
    hist(ns(t).nd(ns(t).nd > 0), 40);
end

% 1.2 studying the overlap of the hubs
% -----------------------------------
for t = 1:5
    ns(t).tissue
    gc = sum(ns(t).nd > 0);
    [a, b] = sort(ns(t).nd, 'descend');
    hubs{t} = b(1:(round(gc/10)));
end

myND = ones(1, 18494);
nums = [2, 3 , 5 , 7, 9];
for t = 1:5
    myND(hubs{t}) = myND(hubs{t}) * nums(t);
end

sum(myND == 9)
superHubs = find(myND == 1890);

% are the hubs connected ?
t = 5
sib = full(ns(t).net(superHubs, superHubs));
sum(sum(sib))
% figure
% view(biograph(sib))
sib = sib + sib';

graphshortestpath(sparse(sib), 1)
% hubs are densely connected. 

% how much do they cover the edges? 
for t = 1:5
    ns(t).tissue
    sib = ns(t).net(superHubs, :) + ns(t).net(:, superHubs)';
    sum1 = sum(sib);
    sum(sum1 > 0)
end

% 2. study of the hubs for TS links
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loading GTEx
load('~/data/general/GPL570GemmaMapNEW.mat')
load('~/networks/GTEx/allNets.mat')
whos GTExNet
gtND = zeros(17, 18494);
for i = 1:17
    gtND(i, :) = GTExNet(i).nd;
end

% checking a particular edge in the GTEx nets 
g1 = 530
g2 = 759
values = zeros(1, 17);
gtnds = zeros(2, 17);
for i = 1:17
    values(1, i) = GTExNet(i).net(g1, g2);
    gtnds(1, i) = GTExNet(i).nd(g1);
    gtnds(2, i) = GTExNet(i).nd(g2);
end
gtnds

% adjusting
book = sum(gtND) > 0;
temp = gtND(:, book);
sib = quantile(temp(:), .95)
adjGTNDs = gtND;
adjGTNDs(adjGTNDs > sib) = sib;

% the TS link network 
load('~/resultsAndFigures/TSlinks/finalTable.mat')
for t = 1:5
    tsn(t).tissue = finalTable(t).tissue;
    tsn(t).net = finalTable(t).wholeNet;
    sum(sum(tsn(t).net))
    tsn(t).nd = sum(finalTable(t).wholeNet) + sum(finalTable(t).wholeNet');
end

% getting the hubs, I guess? 
for t = 1:5
    tsn(t).tissue
    gc = sum(tsn(t).nd > 0);
    [a, b] = sort(tsn(t).nd, 'descend');
    tshubs{t} = b(1:(round(gc/10)));
end

% Q1: what is the correlation of the ND between the tissues. 
tsNdCorr = zeros(5, 5);
ndCorr = zeros(5, 5);
for i =1:5
    tsa = full(tsn(i).nd');
    a = full(ns(i).nd');
    for j = 1:5
        tsb = full(tsn(j).nd');
        b = full(ns(j).nd');
        tsNdCorr(i,j) = corr(tsa, tsb, 'type', 'Spearman');
        ndCorr(i, j) = corr(a, b, 'type', 'Spearman');
    end
end

% Q2: What is the correlation between the two network Nd 
ndCorr = zeros(1, 5);
for i = 1:5
    tsnd = full(tsn(i).nd');
    nsnd = full(ns(i).nd');
    ndCorr(i) = corr(tsnd, nsnd);
end
ndCorr 

% Q3: ok, so how does it look?
addpath('~/codes/MATLAB/myCodes/general/')
% get the nds 
nds = zeros(5, 18494);
tsNDs = zeros(5, 18494);
for i = 1:5
    nds(i,:) = full(ns(i).nd);
    tsNDs(i,:) = full(tsn(i).nd);
end

% getting the hub threshold to help for the heatmap colors
book = sum(nds) > 0;
temp = tsNDs(:, book);
sib = quantile(temp(:), .95)
adjTSNDs = tsNDs;
adjTSNDs(adjTSNDs > sib) = sib;

temp = nds(:, book);
sib = quantile(temp(:), .95)
adjNDs = nds;
adjNDs(adjNDs > sib) = sib;

% tfInds from lines 333 section 4
s = 1
e = 50
h = figure; 
set(h, 'Position', [100, 100, 1500, 900])
subplot(3, 1, 2)
heatmap(adjNDs(:, tfInds(s:e)), tfInds(s:e), {ns(:).tissue},  [], 'TickAngle', 90, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'bone')
%h = figure; 
subplot(3, 1, 1)
heatmap(adjTSNDs(:, tfInds(s:e)), tfInds(s:e), {ns(:).tissue},  [], 'TickAngle', 90, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'bone')
% some observations. I want to look at some of the TFs for now

% looking at GTEx 
subplot(3, 1, 3)
heatmap(adjGTNDs([5, 6, 14:end], tfInds(s:e)), tfInds(s:e), {GTExNet([5, 6, 14:end]).tissue},  [], 'TickAngle', 90, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'bone')
s = s + 50
e = e + 50

% looking at the common genes
load('~/data/general/tissueExpGenes/allExp0.8.mat')

myND = ones(1, 18494);
nums = [2, 3 , 5 , 7, 9];
for t = 1:5
    myND(tshubs{t}) = myND(tshubs{t}) * nums(t);
end

sum(myND == 6)
superTshubs = find(myND == 1890);
% there are no superhubs between the ts hubs, what does that mean.  
myhub = find(myND == 6)

% are the tshubs connected ?
t = 1
sib = full(tsn(t).net(superTshubs, superTshubs));
sum(sum(sib))
% figure
% view(biograph(sib))
sib = sib + sib';

graphshortestpath(sparse(sib), 1)
% tshubs are detsnely connected. 

% how much do they cover the edges? 
for t = 1:5
    tsn(t).tissue
    sib = tsn(t).net(superTshubs, :) + tsn(t).net(:, superTshubs)';
    sum1 = sum(sib);
    sum(sum1 > 0)
end

g = 554
[ns(1).nd(g) ns(2).nd(g) ns(3).nd(g) ns(4).nd(g) ns(5).nd(g)]
[tsn(1).nd(g) tsn(2).nd(g) tsn(3).nd(g) tsn(4).nd(g) tsn(5).nd(g)]

load('~/data/general/GPL570GemmaMapNEW.mat')
gpl570.uniqueSymbols(554)


% 3. TF study for the links 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the TF files, make the newtorks. 
fileName = 'Cell_10.txt' % human TF regulatory networks 
fileName = 'PathwayCommons.txt' % genes with sharing pathways 
fileName = 'Reactome.txt'
fileName = 'Trrust.txt'

file = (['~/data/TFTargetList/' fileName])

fid = fopen(file)

linum = 1;
% I need length of header for the sampleCount
links = textscan(fid, [repmat('%s', 1,2)], ...
                'Headerlines', linum-1, 'Delimiter', ' ');
fclose(fid);
length(links{1})
links{1}(1:10)

% loading my networks and my TS links. 
load('~/data/general/GPL570GemmaMapNEW.mat')
load('~/resultsAndFigures/TSlinks/finalTable.mat')

[ca1, cb1] = ismember(links{1}, gpl570.uniqueSymbols);
[ca2, cb2] = ismember(links{2}, gpl570.uniqueSymbols);

holu = ca1 + ca2;

inLinks = holu == 2;
sum(inLinks)
length(inLinks)

trueInds1 = cb1(inLinks);
trueInds2 = cb2(inLinks);

for i = 1:length(trueInds1)
    if (trueInds1(i) > trueInds2(i))
        temp = trueInds1(i);
        trueInds1(i) = trueInds2(i);
        trueInds2(i) = temp;
    end
end

% making the network from the links 
newNet = sparse(trueInds1, trueInds2, 1, 18494, 18494);

tTSnet = finalTable(5).wholeNet;
sum(sum(tTSnet .* newNet))

% 4. TF studies for the links 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileName = 'nrg2538-s3.txt'

file = (['~/data/TFTargetList/' fileName])

fid = fopen(file)

fline = fgets(fid)

linum = 11;
% I need length of header for the sampleCount
data = textscan(fid, [repmat('%s', 1,7)], ...
                'Headerlines', linum, 'Delimiter', '\t');

[a, b] = ismember(data{6}(1:1364), gpl570.uniqueSymbols);
sum(a)
% switch points [718, 1364, 1391, 1591] - a, b, others, c, x
tfInds = b(a);

load('~/resultsAndFigures/TSlinks/finalTable.mat')
load('~/resultsAndFigures/tissueNetworkStudies/networksInfo.mat')

t = 1
book = finalTable(t).nd(b(a));
sum(book > 0) 
hist(book(book>0), 40)
id = find(book > 100)

myGenesInd = b(a);
id2 = myGenesInd(id)

t = 5
finalTable(t).nd(id2)
gpl570.uniqueSymbols(id2)

% gene ID 6925 is a blood specific transcriptom? 
ID = 6925

for i = 1:5
    finalTable(i).tissue
    finalTable(i).nd(ID)
    ns(i).nd(ID)
end

load('~/resultsAndFigures/TSlinks/TSlinks.mat')
t = 5
kado = sparse(TSlinks{t}.a, TSlinks{t}.b, TSlinks{t}.na, 18494, ...
              18494);
ID1 = 6925
ID2 = 13217
kado(ID1, ID2)
[a, b] = ismember(data{1}, 'b' );

% 5. HKE studies for the tissue networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
load(['~/resultsAndFigures/tissueNetworkStudies/' ...
      'networksInfo_extended.mat'])

load(['~/resultsAndFigures/tissueNetworkStudies/' ...
      'networksInfo.mat'])

load('~/data/general/hkeNet.mat')

sum(sum(hkeNet))

tempNet = zeros(18494, 18494);
for i = 1:5
    tempNet = tempNet + ns(i).net;
end

max((max(tempNet)))

sib = tempNet .* hkeNet';
sib1 = tempNet + (hkeNet' .* 7);

[a, b, c] = find(sib);
kado = [a, b, c];

[ah, bh, ch] = find(hkeNet');

% getting the number of missed links due to the expression level
expMissedLinks = zeros(1, 5);
totalMissedLinks = zeros(1, 5);
totalPresentLinks = zeros(1, 5);
presentLinksToPotLinks = zeros(1, 5);
hkeCount = sum(sum(hkeNet));
for i = 1:5
    [in1, index] = ismember(ah, find(~ns(i).expGenes));
    [in2, index] = ismember(bh, find(~ns(i).expGenes));
    sumin = in1 + in2;
    expMissedLinks(i) = sum(sumin > 0);
    totalMissedLinks(i) = hkeCount - sum(sum(hkeNet' .* ...
                                                     ns(i).net));
    totalPresentLinks(i) = sum(sum(hkeNet' .* ns(i).net));
    presentLinksToPotLinks(i) = totalPresentLinks(i) / (hkeCount - expMissedLinks(i));
end

wholeMat = [expMissedLinks; totalMissedLinks; totalPresentLinks; ...
            presentLinksToPotLinks]
wholeMat1 = [expMissedLinks; totalMissedLinks; totalPresentLinks; presentLinksToPotLinks]

% 5.1  writing this down in a table
file = ['~/resultsAndFigures/tissueNetworkStudies/' ...
        'HKE_rep_extendedNet.txt']
fid = fopen(file, 'w')

fprintf(fid, ['attribute\t blood\t brain\t liver\t lung\t skeletalMuscle ' ...
              '\n'])

fprintf(fid, ['%s\t' repmat('%d\t', 1, 5) '\n'],...
            'Missed links', totalMissedLinks)
fprintf(fid, ['%s\t' repmat('%d\t', 1, 5) '\n'],...
            'Present links', totalPresentLinks)
fprintf(fid, ['%s\t' repmat('%d\t', 1, 5) '\n'],...
            'Missed links due to expression', expMissedLinks)
fprintfn(fid, ['%s\t' repmat('%.3f\t', 1, 5) '\n'],...
            'Percent of potential links present', presentLinksToPotLinks)

fclose(fid)

% 5.2 plot the histogram of percent for before and after networks

myBars = [wholeMat1(end, :); wholeMat(end, :)]

h = figure
c1 = [140, 140, 140];
c2 = [50, 50, 50];
colors = [c1; c2]/256
mybar = bar(myBars')
set(mybar, 'edgecolor', 'none')
set(gca, 'XTick', [1:5])
sudoT = {ns(:).tissue}
sudoT(5) = {'muscle'}
set(gca, 'XTickLabel', sudoT)
ylim([0, 1])
colormap(colors)
legend('Tissue networks with FDR cut-off 5e-5', ['Tissue networks with FDR ' ...
                    'cut-ff 5e-3'])
set(gca, 'FontSize', 14)
ylabel(sprintf('Percent of HKE presenet\n between the expressed genes'))

figFolder = '~/resultsAndFigures/tissueNetworkStudies/figures/HKEpercentPresent'
file = [figFolder '.pdf']
print(h, '-dpdf', file)
file = [figFolder '.eps']
print(h, '-deps', file)

% 5.2 plot the histogram of percent for before and after networks

totLink = sum(sum(hkeNet))
newMat = wholeMat;

myBar(3, :) = (newMat(2, :) - newMat(1, :)) % absent hke
myBar(2, :) = newMat(1, :) % absent due to exp 
myBar(1, :) = newMat(3, :) % total present links

c1 = [140, 140, 140];
c2 = [95, 95, 95];
c3 = [50, 50, 50];

colors = [c1; c2; c3]/256
h = figure
mybar = bar(myBar', 'stacked')

ylim([0, 3200])
set(mybar, 'edgecolor', 'none')
set(gca, 'XTick', [1:5])
sudoT = {ns(:).tissue}
sudoT(5) = {'muscle'}
set(gca, 'XTickLabel', sudoT)
grid on

colormap(colors)
legend('Present HKE', 'HKE absent due to the expression level', 'HKE absent')
set(gca, 'FontSize', 12)
ylabel(sprintf('Count of edges'))

title('Presence of HKE in the extended tissue networks')

figFolder = '~/resultsAndFigures/tissueNetworkStudies/figures/HKEpresence_Extended'
file = [figFolder '.pdf']
print(h, '-dpdf', file)
file = [figFolder '.eps']
print(h, '-deps', file)

% just plotting the  count of datasets
h = figure
mybar = barh([15, 12 ,9, 7 10])
xlim([1, 16])
ylim([0, 6])
set(mybar, 'edgecolor', 'none')
set(gca, 'YTick', [1:5])
set(gca, 'YTickLabel', {'lung', 'brain', 'blood', 'muscle', 'liver'})
grid on

colors = c2/255
colormap(colors)
set(gca, 'FontSize', 12)
ylabel(sprintf('Count of edges'))

title('Count of datasets for eac tissue')

figFolder = '~/resultsAndFigures/tissueNetworkStudies/figures/datasetCount'
file = [figFolder '.pdf']
print(h, '-dpdf', file)
file = [figFolder '.eps']
print(h, '-deps', file)

% 6. The GTEx overlap of the networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

addpath('~/codes/MATLAB/myCodes/general/')
load('~/resultsAndFigures/affyAndGTExComparison/res.mat')

% plot1: Heatmap of the network overlaps
sib = GTExRep.netOverlap
h = figure; 
heatmap(sib([8, 14:17],:), {'brain', 'blood', 'liver', 'lung', 'muscle'}, {'brain', 'blood', 'liver', 'lung', 'muscle'},  [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'gray')
set(gca, 'FontSize', 12)

title('Overlap of the binary tissue networks between GTEx and AffyMetrix datasets')
figFolder = '~/resultsAndFigures/affyAndGTExComparison/GTEx_binaryNetOverlapHM'
file = [figFolder '.pdf']
print(h, '-dpdf', file)
file = [figFolder '.eps']
print(h, '-deps', file)

% plot2: Density of the networks, some sort (edit the title)
h = figure
[xi, fi] = ksdensity(GTExRep.allQs(8, :));
plot(fi, xi)
hold on

for i = 1:5
    [xi, fi] = ksdensity(GTExRep.netQs(ind, :));
    plot(fi, xi)
end


% 7. The Illumina overlap of the networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

addpath('~/codes/MATLAB/myCodes/general/')
load('~/resultsAndFigures/affyAndIllComparison/res.mat')

sib = illRep.netOverlap
h = figure; 
heatmap(sib, {'brain', 'blood', 'liver', 'lung', 'muscle'}, {'blood', 'blood', 'blood', 'brain', 'brain', 'brain'},  [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'gray')
set(gca, 'FontSize', 12)

title(['Overlap of the binary tissue networks between Illumina Bead Chip datasets ' ...
       'and AffyMetrix datasets'])
figFolder = '~/resultsAndFigures/affyAndIllComparison/Ill_binaryNetOverlapHM'
file = [figFolder '.pdf']
print(h, '-dpdf', file)
file = [figFolder '.eps']
print(h, '-deps', file)

% plot2: Density of the networks, some sort (edit the title)
h = figure
set(h, 'Position', [100, 100, 1500, 900])
for i = 1 : 6
    subplot(2, 3, i)
    [xi, fi] = ksdensity(illRep.defaultQs(i, :));
    plot(fi, xi, 'color', 'k')
    hold on
    for j = 1:5
        tColor = 'b'
        ind = (j-1) * 6 + i
        hold on
        if(ismember(i, [1:3]) && (ismember(j, [2, 8, 14, 20, 26])))
            tColor = 'r'
        end
        if(ismember(i, [4:6]) && (ismember(j, [1, 7, 13, 19, 25])))
            tColor = 'r'
        end
        [xi, fi] = ksdensity(illRep.netQs(ind, :));
        plot(fi, xi, 'color' ,tColor)
    end
end

% TODO: plot this plot on the local
