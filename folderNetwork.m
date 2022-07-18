% this file is on making tissue(here I call it network) specific
% coexpression networks. I am doing it using correlation networks. 
% 1. making a corrNetwork from one experiment/ save the sparce matrix
% for all experiments
% 1.1 building the network for a selected set of genes, this I used
% for allExp genes (7835) and also for the tissue ExpGenes
% 2. measure the convergence among the networks build as adding up
% networks. (start with skeletal muscle)
% 3. run the network for at least skeletal muscle. 
% 4. trying to standardize the coexpression data

clear

tissue = 'blood'
%dataFolder = '/space/grp/marjan/data/'
%exprFolder = [dataFolder tissue ['/geneExprMat/gemmaProbeMap/clean/' ...

dataFolder = '~/data/array/'
exprFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/']
fileList = dir([exprFolder '*.mat'])
gCount = 18494;

%gCount = 100;xo
aggMat = zeros(gCount);
aggNet = zeros(gCount);
QSingle = zeros(length(fileList), 1);

%1. the code below saves both the singular and aggregated networks at
%each loop, as the singular networks are added to the aggregated
%one. 
% note: threshold is inside the code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

tissue = 'brain'
dataFolder = '~/data/affyArray/tissues/'
exprFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/']
fileList = dir([exprFolder '*.mat'])
gCount = 18494;
thr = 0.1;

% getting the single binary networks: 
for i = 1: length(fileList)
    tic
    %loading new file
    '1'
    load([exprFolder fileList(i).name]);
    sib = corr(dataSet.mat');
    sib = sib - eye(gCount);
    
    '2'
    %getting the file's network
    upperSingle = sib(logical(triu(ones(size(sib)), 1)));
    
    QSingle(i) = quantile(upperSingle, (1 - thr))
    
    singleNet = sib > QSingle(i);
    sparseSingleNet = sparse(singleNet);
    
    '3'
    %saving the single Network
    [a, b] = regexp(fileList(i).name ,'GSE[0123456789]+');
    Name = fileList(i).name(a:b)
    save(sprintf('~/networks/tissues/%s/singleNet/%s_%.2f.mat', ...
                 tissue, Name, thr), 'sparseSingleNet'); % i is included so I can follow up which experiments were in the i'th aggNet
    
    '4'
    %adding it to the aggregated ranked matrix we had so far and
    %get the results - tied rank takes 10m ?
    % aggMat = aggMat + reshape(tiedrank(sib(:)), gCount, gCount);
    % upperAgg = triu(aggMat, 1);
    % QAgg = quantile(upperAgg(:), (1 - 0.005));
    % aggNet = aggMat > QAgg;
    % sparseAggNet = sparse(aggNet)  
    % '5'
    % %saving the AggNetwork 
    % save(sprintf('%s%s/AggNet/%dagg.mat', dataFolder, tissue, (i)), 'aggNet');
    
    % clear dataSet
    toc
end


save(sprintf('~/networks/tissues/%s/singleNet/QSingle%.2f.mat', tissue, thr), ...
     'QSingle')

% 1.1 filtering based on the exp level amongst all the datasets,
% this is for the common network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'};

tissues = {'liver', 'skeletalMuscle'}
gCount = 18494;
thrArray = [0.10, 0.05, 0.01, 0.005];

for t = 1:length(tissues)
    tissue = tissues{t};
    dataFolder = '~/data/affyArray/tissues/'
    exprFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/']
    fileList = dir([exprFolder '*.mat'])
    % getting the single binary networks: 

    % just loading the set of genes we want to have the network on
    expThr = '0.8'
    sampleCounts = zeros(1, length(fileList))
    QSingle = zeros(1, length(fileList))
    for cthr = 1:length(thrArray)

        thr = thrArray(cthr)
        %thr = .1;

        expTemp = zeros(gCount, gCount);


        %%%%!!!!!! whichever I want!
        % load(['~/data/general/tissueExpGenes/allExp' expThr '.mat'])
        % expTemp(allExp, allExp)  = 1;

        load(['~/data/general/tissueExpGenes/' tissue 'ExpGenes' expThr '.mat'])
        expTemp(expGenesInd, expGenesInd)  = 1;
        
        % expGenesInd = logical(ones(1, 18494));
        % expTemp(expGeneInd, expGeneInd)  = 1;
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        upperTemp = (logical(triu(ones(size(expTemp)), 1)));
        finalTemp = logical(upperTemp .* expTemp);
        for i = 1: length(fileList)
            tic
            %loading new file
            '1'
            load([exprFolder fileList(i).name]);
            sampleCounts(i) = size(dataSet.mat, 2);
            sib = corr(dataSet.mat');
            sib = sib - eye(gCount);
            
            '2'

            %getting the file's network
            upperSingle = sib(finalTemp);
            QSingle(i) = quantile(upperSingle, (1 - thr))        


            sib = sib .* finalTemp;
            singleNet = sib > QSingle(i);
            sparseSingleNet = sparse(singleNet);
            
            '3'
            %saving the single Network
            [a, b] = regexp(fileList(i).name ,'GSE[0123456789]+');
            Name = fileList(i).name(a:b)

            % you have to pick one of them, based on which genes you are
            % looking at
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % save(sprintf('~/networks/tissues/%s/singleNet/%sNetExpThrAll%s_%s_%.2f.mat', ...
            %              tissue ,tissue, expThr, Name, thr), 'sparseSingleNet'); %...


            save(sprintf('~/networks/tissues/%s/singleNet/%sNetExpThr%s_%s_%.3f.mat', ...
                         tissue ,tissue, expThr, Name, thr), ...
                 'sparseSingleNet'); 
            
            
            % save(sprintf('~/networks/tissues/%s/singleNet/%sNetWholeGenes%s_%.3f.mat', ...
            %              tissue ,tissue, Name, thr), ...
            %      'sparseSingleNet'); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            
            '4'
            %adding it to the aggregated ranked matrix we had so far and
            %get the results - tied rank takes 10m ?
            % aggMat = aggMat + reshape(tiedrank(sib(:)), gCount, gCount);
            % upperAgg = triu(aggMat, 1);
            % QAgg = quantile(upperAgg(:), (1 - 0.005));
            % aggNet = aggMat > QAgg;
            % sparseAggNet = sparse(aggNet)  
            % '5'
            % %saving the AggNetwork 
            % save(sprintf('%s%s/AggNet/%dagg.mat', dataFolder, tissue, (i)), 'aggNet');
            
            % clear dataSet
        end
        corrSig(t).tissue = tissue;
        corrSig(t).sc = sampleCounts;
        corrSig(t).corrThr = QSingle;
        toc
    end
end

save('~/resultsAndFigures/tissueNetworkStudies/allGenes_corrSignificanceWholeGenes.mat', ...
     'corrSig')   

save('~/resultsAndFigures/tissueNetworkStudies/corrSignificance.mat', ...
     'corrSig')   

load('~/resultsAndFigures/tissueNetworkStudies/matFiles/allGenes_corrSignificanceWholeGenes.mat')

% get the t statistic significance for the correlation values
clear
load(['~/resultsAndFigures/tissueNetworkStudies/' ...
      'corrSignificance.mat'])
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'};

sib = 0
group1 = 0
group2 = 0
group3 = 0
group4 = 0
for t = 1:5
    tstat = zeros(1, length(corrSig(t).sc));
    pval = zeros(1, length(corrSig(t).sc));
    for j = 1:length(corrSig(t).sc)
        r = corrSig(t).corrThr(j);
        n = corrSig(t).sc(j);
        if ismember(n, [19:29])
            group1 = group1 + 1
        end
        if ismember(n, [30:39])
            group2 = group2 + 1
        end
        if ismember(n, [40:49])
            group3 = group3 + 1
        end
        if(n >= 50)
            group4 = group4 + 1
        end
        tstat(j) = r / sqrt((1 - r^2)/(n - 2));
        pval(j) = 1 - tcdf(tstat(j), (n -2));
        if pval(j) > sib
            sib = pval(j)
            j
            t
        end
    end
    corrSig(t).tstat = tstat;
    corrSig(t).pval = pval;
end

% 1.2 Getting the FDR for edges from the low expressed genes. 
clear
tissues = {'blood', 'lung', 'skeletalMuscle', 'liver', 'brain'}
binNetThr1 = '0.05'
binNetThr2 = '0.050'
expThr = '0.8'
fde = zeros(1, 53);

nlAll = zeros(1,53);
hnlAll = zeros(1,53);
elAll = zeros(1,53);
exhThrLinks = zeros(1,53);

load('~/data/general/linkExprInfo/dataSetExpInf.mat')
c = 0
for t = 1:5  
    tissue = tissues{t};
    fileList1 = dir(['~/networks/tissues/' tissue '/singleNet/' tissue ...
                     '*NetExpThr0.8*' binNetThr1 '.mat'])
    fileList2 = dir(['~/networks/tissues/' tissue '/singleNet/' tissue ...
                     'NetWholeGenes*' binNetThr2 '.mat'])
    
    for f = 1:length(fileList1)
        if(~isempty(strfind(fileList1(f).name, dataSetExpInf(f).name)) ...
           & ~isempty(strfind(fileList2(f).name, ...
                              dataSetExpInf(f).name)))
            qs = quantile(dataSetExpInf(f).meanExp, 2);
            noiseGenes = dataSetExpInf(f).meanExp <= qs(1);
            tempNoiseNet = zeros(18494, 18494);
            tempNoiseNet(noiseGenes, noiseGenes) = 1;
            tempExpNet = zeros(18494, 18494);
            tempExpNet(~noiseGenes, ~noiseGenes) = 1;
            
            tempMixedNet = ones(18494, 18494) - tempNoiseNet;

            load(sprintf('~/networks/tissues/%s/singleNet/%s', ...
                         tissue, fileList1(f).name));
            expThrNet = sparseSingleNet;
            clear sparseSingleNet;
            sum(sum(expThrNet))

            load(sprintf('~/networks/tissues/%s/singleNet/%s', ...
                         tissue, fileList2(f).name));
            wholeGeneNet = sparseSingleNet;
            noiseLinks = tempNoiseNet .* wholeGeneNet;
            expLinks = tempExpNet .* wholeGeneNet;
            sum(sum(noiseLinks))/sum(sum(wholeGeneNet))
            halfNoiseLinks = wholeGeneNet - noiseLinks - expLinks;

            nlAll(f) = sum(sum(noiseLinks));
            hnlAll(f) = sum(sum(halfNoiseLinks == 1));
            elAll(f) = sum(sum(expLinks));
            expThrLinks(f) = sum(sum(expThrNet));
            c = c + 1
        end
    end
end

noiseLinks.nl = nlAll;
noiseLinks.halfNl = hnlAll;
noiseLinks.expNl = elAll;
noiseLinks.expThrLinks = expThrLinks;

save(['~/resultsAndFigures/noiseLinksBetweenLowExp_binNetThr' binNetThr2 ...
     '.mat'], 'noiseLinks')

% 2.  getting the tissue aggregated networks 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

tissue = 'brain'
dataFolder = '~/data/affyArray/tissues/'
exprFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/']
fileList = dir([exprFolder '*.mat'])
saveFolder = ['~/networks/' tissue '/']
gCount = 18494;
aggNet = zeros(gCount);
upSib = zeros(gCount);

% getting the single binary networks: 
for i = 1: length(fileList)
    tic
    %loading new file
    '1'
    load([exprFolder fileList(i).name]);
    sib = corr(dataSet.mat');
    
    %getting the rank values for upper triangle of the net
    upSib = triu(sib, 1);
    [a, b, c] = find(upSib);
    tic
    holu = tiedrank(c);
    toc
    length(holu)
    
    
    %adding the results to the aggNet
    tic
    for j = 1:length(a)
        aggNet(a(j), b(j)) = aggNet(a(j), b(j)) + holu(j);
    end
    toc
end
%getting the 005 of the highest ranks for binary net

allVal = aggNet(logical(triu(ones(size(aggNet)), 1)));
QaggNet = quantile(allVal, (1 - 0.001));

holu = aggNet > QaggNet;
sbinAggNet = sparse(holu);

% save sparcified binary net    
save([saveFolder 'sparseBinaryAggRank001.mat'], 'sbinAggNet')

% 3.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
dataFolder = '/space/grp/marjan/data/';
tissue = 'blood';
fileList = dir([dataFolder tissue '/AggNet/'])
%load([dataFolder tissue '/AggNet/' fileList(9).name])
%fileList(9).name
load(sprintf('%s%s/AggNet/QSingle.mat', dataFolder, tissue))

%find out the number of edges which varies among the networks 

similarity = zeros(1, length(fileList)-4);
for i = 1:(length(fileList)-4)
    %   i = 1
    load(sprintf('%s%s/AggNet/%dagg.mat', dataFolder, tissue, ...
                        i));%aggNet loaded
    tempNeti = full(aggNet);
    clear aggNet
    'net i'
    load(sprintf('%s%s/AggNet/%dagg.mat', dataFolder, tissue, ...
                        i+1));%aggNet loaded
    tempNetiPlus = full(aggNet);
    'net i+'
    total = tempNeti + tempNetiPlus;
    total = total == 2;
    commonEdges = sum(sum(total));
    similarity(i) = commonEdges/sum(sum(tempNeti))
end
clear tempNetiPlus
clear tempNeti
clear total
SMsim = similarity;
lungSim = similarity;
bloodSim = similarity;

% putting on the plot the quantile on which the top 0.5% edges were
% selected - something should be fixed here though - not needed for now
%plotting
xtick = cell(1, length(QSingle)-2)
for i = 1:length(QSingle) - 2
    xtick{i} = sprintf('%0.2f', QSingle(i));
end
xtick


h = figure();
hold all
plot(SMsim , 'Color', 'r')
plot(lungSim , 'Color', 'b')
plot(bloodSim , 'Color', 'g')
grid on
%set(gca, 'Xtick', 1:length(similarity), 'XTickLabel', xtick)
xlabel('number of datasets added')
ylabel('percent of edges in commong between i, i+1 networks')
title('coexpression network convergence for tissues and random network')
fileName = 'corrNetConvergenceTissueRandom02';
figFolder = [dataFolder 'general/']
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);

%% loading the random networks and adding them to the plot
load(['/space/grp/marjan/prioritization Project/results/' ...
      'randomNetworkConv10000.mat'])
hold on
plot(conv(1:20), 'LineWidth', 2)

% 4. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the correlation of Quantile values to the number of samples in
% the experiments
clear
h = figure()
dataFolder = '/space/grp/marjan/data/';
tissue = 'lung';
%load([dataFolder tissue '/AggNet/' fileList(9).name])
%fileList(9).name
load(sprintf('%s%s/AggNet/QSingle.mat', dataFolder, tissue))

exprFolder = [dataFolder tissue '/geneExprMat/clean/']
fileList = dir(exprFolder)
gCount = 21050;
%getting the number of samples

sampleCounts = zeros(1, length(fileList) -2)
for i = 3: length(fileList)
load([exprFolder fileList(i).name]);
sampleCounts(i-2) = size(dataSet.mat, 2);
clear dataSet
end

corr(sampleCounts', QSingle(1: length(sampleCounts)))
hold on
plot(sampleCounts, QSingle(1:length(sampleCounts)), '*', 'Color', 'g')


xlabel('sample count')
ylabel('the 995 permilles')
title(['correlation of the number of samples and the 995 permilles ' ...
       'for different tissues'])
fileName = 'sampleCountPermilles';
figFolder = [dataFolder 'general/']
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);

% 5. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the mean similarity among networks for a tissue

dataFolder = '/space/grp/marjan/data/';
tissue = 'lung';
fileList = dir([dataFolder tissue '/SingleNet/'])

fileList(:).name
similarity = zeros(length(fileList)-2, length(fileList)-2);
for i=3:length(fileList)
    i% = 3
       load(sprintf('%s%s/SingleNet/%s', dataFolder, tissue, ...
                        fileList(i).name));%aggNet load
        tempNeti = full(sparseSingleNet);
        clear sparseSingleNet;
    for j = i+1:length(fileList)
        j %= 4 
        load(sprintf('%s%s/SingleNet/%s', dataFolder, tissue, ...
                        fileList(j).name));%aggNet loaded
        tempNetj = full(sparseSingleNet);
        total = tempNeti + tempNetj;
        total = total == 2;
        commonEdges = sum(sum(total));
        similarity(i-2,j-2) = commonEdges/sum(sum(tempNeti)) ;       
        similarity(j-2,i-2) = commonEdges/sum(sum(tempNeti));
        book = commonEdges/sum(sum(tempNeti))
    end
end

bloodSim = mean(similarity);
lungSim = mean(similarity);
SMSim = mean(similarity);

subres = zeros(1,length(res))
res = [SMSim, bloodSim, lungSim]
wholeLabel = [{'muscle'}, repmat({''}, 1, 8), 'muscle',repmat({''}, 1, 9), ...
              {'lung'},repmat({''}, 1,15 )]

book = [res; 1:length(res)]
h = figure;
subres = zeros(1,length(res))
subres(1:9) = SMSim
bar(subres, 'FaceColor', 'b')
hold on
subres = zeros(1,length(res));
subres(10:19) = bloodSim
bar(subres, 'FaceColor', 'r')
hold on
subres = zeros(1,length(res));
subres(20:end) = lungSim;
bar(subres, 'FaceColor', 'g')
set(gca, 'Xtick', 1:length(res), 'XTickLabel', wholeLabel)
xlabel('tissue')
ylabel(['the mean difference between networks from each dataset and ' ...
        'other datasets in the same tissue'])
fileName = 'networkSimilarityIntraTissue';
figFolder = [dataFolder 'general/figures/']
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);


%heatmap the similarity matrix
boxNames = cell(1, length(fileList) - 2)
for i = 1:length(fileList)-2
    [a,b] = regexp(fileList(i+2).name, 'GSE[0123456789]+');
    boxNames{i} = [fileList(i+2).name(a:b) '-']
    end
h = figure();
%boxplot(similarity,'labels', boxNames)
boxplot(similarity)
% set(gca(), 'XTickLabel', boxNames)
% rotateXLabels(gca(), 45)
ylabel('similarity with other networks-percent of edges in common')
title([tissue ' datasets'])
fileName = 'BPnetworkSimilarity';
figFolder = [dataFolder tissue '/figures/']
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);

%network similarity among the tissues

% 6. heatmap for comparison of edges in my networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
%data folder for the blood 
folder = cell(1, 5)
folder{1} = '~/networks/tissues/blood/singleNet/'
folder{2} = '~/networks/tissues/lung/singleNet/'
folder{3} = '~/networks/tissues/skeletalMuscle/singleNet/'
folder{4} = '~/networks/tissues/liver/singleNet/'
folder{5} = '~/networks/tissues/brain/singleNet/'

ts{1} = 'blood';
ts{2} = 'lung';
ts{3} = 'SM';
ts{4} = 'liver';
ts{5} = 'brain';
 
thr = '01'
 %getting total number of networks
tissueExpCount = zeros(5, 1);
for i = 1:length(folder)
     i
     folder{i}
     fileList = dir([folder{i} 'GSE*_' thr '.mat']);
     tissueExpCount(i) = length(fileList);
     length(fileList)
 end

 % getting the list of network files so I can load them one by one
 networksSim = zeros(sum(tissueExpCount));
 netFiles = cell(sum(tissueExpCount), 1);
 c = 1;
 for i = 1:length(folder)
     fileList = dir([folder{i} 'GSE*_' thr '.mat'])
     for j = 1:length(fileList)
         netFiles{c} = [folder{i} fileList(j).name];
         c = c +1;
     end
 end

 %getting the smiliarity
GSEIDs = cell(1, length(netFiles));
 gCount = 18494
 for i = 1:length(netFiles)
     i
     load(netFiles{i})
     iNet = full(sparseSingleNet);
     iNet = iNet - eye(gCount);
     
     [a, b] = regexp(netFiles{i}, 'GSE[0123456789]+');
     GSEID = netFiles{i}(a:b)

     GSEIDs{i} = GSEID
     
     % check the exp levels here. 
     % get the GSE file ID and read that GSEID netStr file netStr >=
     % 33 OR 27


     for j = i + 1:length(netFiles)
        j
        load(netFiles{j})
        jNet = full(sparseSingleNet);

        % check the exp levels here 
        temp = iNet + jNet;
        similarity = (sum(sum((iNet + jNet) == 2)));
        networksSim(i, j) = similarity;
        networksSim(j, i) = similarity;
    end
end

netSimStruct.networksSim = networksSim;
netSimStruct.GSEIDs = GSEIDs;
netSimStruct.tissueExpCount = tissueExpCount;

file = ['~/networks/networkSim' thr '.mat']
save(file, 'netSimStruct')

book = networksSim ./ sum(sum(jNet));
book = networksSim;

labels = repmat({''}, 1, size(book, 1));
labels{1} = ts{1};
labels{10} = ts{2};
labels{25} = ts{3};
labels{33} = ts{4};
% book = book + eye(length(labels), length(labels));

h = figure; 
heatmap(book, labels, labels, [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'hot')
'heatmap done'
title(['Network similarity among experiments from different ' ...
       'tissues'])
figFolder = '~/'
fileName = ['HMnetworkSim'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);

% 7. getting the same thing but filter for the expression data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
%data folder for the blood 
folder = cell(1, 5)
folder{1} = '~/networks/tissues/blood/singleNet/'
folder{2} = '~/networks/tissues/lung/singleNet/'
folder{3} = '~/networks/tissues/skeletalMuscle/singleNet/'
folder{4} = '~/networks/tissues/liver/singleNet/'
folder{5} = '~/networks/tissues/brain/singleNet/'

ts{1} = 'blood';
ts{2} = 'lung';
ts{3} = 'SM';
ts{4} = 'liver';
ts{5} = 'brain';
 
thr = '01'
%getting total number of networks
tissueExpCount = zeros(5, 1);
for i = 1:length(folder)
    i
    folder{i}
    fileList = dir([folder{i} 'GSE*_' thr '.mat']);
    tissueExpCount(i) = length(fileList);
    length(fileList)
end

% getting the list of network files so I can load them one by one
networksSim = zeros(sum(tissueExpCount));
netFiles = cell(sum(tissueExpCount), 1);
netStrFiles = cell(sum(tissueExpCount), 1);
c = 1;
for i = 1:length(folder)
    fileList = dir([folder{i} 'GSE*_' thr '.mat'])
    for j = 1:length(fileList)
        netFiles{c} = [folder{i} fileList(j).name];

        % here, I want to fill in the netStr file list

        [a, b] = regexp(fileList(j).name, 'GSE[0123456789]+');
        GSEID = fileList(j).name(a:b);

        netStrFiles{c} = [folder{i} 'netStr_' GSEID '.mat'];
        c = c +1;
    end
end

%getting the smiliarity
gCount = 18494
for i = 1:length(netFiles)
    i
    load(netFiles{i})
    iNet = full(sparseSingleNet);
    iNet = iNet - eye(gCount);

     [a, b] = regexp(netFiles{i}, 'GSE[0123456789]+');
     GSEID = netFiles{i}(a:b)

     GSEIDs{i} = GSEID

    
    load(netStrFiles{i})

    % getting links with expr > a certain level
    expNet = abs(netStr.expNet);
    temp1 = expNet >=35;
    temp2 = expNet == 25;
    temp = temp1 + temp2;
    iexpNet = temp;
    clear expNet temp1 temp2 temp;

    iNet = (iexpNet .* iNet) > 0;
    

    for j = i + 1:length(netFiles)
        j
        load(netFiles{j})
        jNet = full(sparseSingleNet);

        load(netStrFiles{j})
        % getting links with expr > a certain level
        expNet = abs(netStr.expNet);
        temp1 = expNet >=35;
        temp2 = expNet == 25;
        temp = temp1 + temp2;
        jexpNet = temp;
        clear expNet temp1 temp2 temp;

        jNet = (jexpNet .* jNet) > 0;

        % check the exp levels here 
        temp = iNet + jNet;
        similarity = (sum(sum((iNet + jNet) == 2)));
        networksSim(i, j) = similarity;
        networksSim(j, i) = similarity;
    end
end


netSimStruct.networksSim = networksSim;
netSimStruct.GSEIDs = GSEIDs;
netSimStruct.tissueExpCount = tissueExpCount;

file = ['~/networks/networkSimExpTHR5_' thr '.mat']
save(file, 'netSimStruct')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8. getting network sim for the aggregated networks
clear 
%data folder for the blood 
folder = cell(1, 4)
file{1} = '/space/grp/marjan/data/blood/AggNet/10agg.mat'
file{2} = '/space/grp/marjan/data/lung/AggNet/17agg.mat'
file{3} = '/space/grp/marjan/data/skeletalMuscle/AggNet/9agg.mat'
file{4} = '/space/grp/marjan/data/liver/AggNet/12agg.mat'

ts{1} = 'blood';
ts{2} = 'lung';
ts{3} = 'SM';
ts{4} = 'liver';

for i = 1:length(file)
    load(file{i})
    iNet = aggNet;
    for j = i + 1:length(file)
        j
        load(file{j})
        jNet = full(aggNet);
        temp = iNet + jNet;
        similarity = (sum(sum((iNet + jNet) == 2)));
        networksSim(i, j) = similarity;
        networksSim(j, i) = similarity;
    end
end

totalEdge = sum(sum(aggNet));

for i = 1:length(file)
    %    networksSim(i,i) = sum(sum(aggNet));
        networksSim(i,i) = nan;
end

h = figure; 
heatmap(networksSim, ts, ts, networksSim, 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'summer')
'heatmap done'
title(sprintf(['Network similarity for tissues - number of edges in common ' ...
               '\n number of edges in each network : %d'], totalEdge))
figFolder = '/space/grp/marjan/data/general/'
fileName = ['HMaggNetSim'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);

% getting the node degree correlation among tissues 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
 %data folder for the blood 
 folder = cell(1, 4)
 folder{1} = '/space/grp/marjan/data/blood/SingleNet/'
 folder{2} = '/space/grp/marjan/data/lung/SingleNet/'
 folder{3} = '/space/grp/marjan/data/skeletalMuscle/SingleNet/'
 folder{4} = '/space/grp/marjan/data/liver/SingleNet/'
 
 ts{1} = 'blood';
 ts{2} = 'lung';
 ts{3} = 'SM';
 ts{4} = 'liver';
 
 %getting total number of networks
tissueExpCount = zeros(4, 1);
for i = 1:length(folder)
     i
     folder{i}
     fileList = dir(folder{i});
     tissueExpCount(i) = length(fileList) - 2;
     length(fileList) - 2
end

% just getting gCount 
fileList = dir(folder{1})
load([folder{1} fileList(3).name])
gCount = size(sparseSingleNet,1);
 
netFilesCount = sum(tissueExpCount);
NDmat = zeros(gCount, netFilesCount);
c = 1;
for i = 1:length(folder)
    fileList = dir(folder{i})
    for j = 3:length(fileList)
        load([folder{i} fileList(j).name])
        NDmat(:, c) = sum(sparseSingleNet);
        c = c +1;
    end
end

% getting the node degree correlation for networks
sib = corr(NDmat);

labels = repmat({''}, 1, size(sib, 1));
labels{1} = ts{1};
labels{11} = ts{2};
labels{28} = ts{3};
labels{37} = ts{4};

h = figure; 
heatmap(sib, labels, labels, [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'hot')
'heatmap done'
title(['Network node degree correlation among datasets from different ' ...
       'tissues'])
figFolder = '/space/grp/marjan/data/general/'
fileName = ['HMnetworkNDCorr'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getting sum of the ranks for the tissue networks 
clear
tissue = 'blood'
gCount = 18377;
dataFolder = '~/data/array/blood/matFiles/geneExprMatGemmaMapBER/outlierRemoved/'

fileList = dir([dataFolder '*.mat'])
rankSum = zeros(gCount, gCount);

for i = 1: length(fileList)
    load([dataFolder fileList(i).name]);
    dataMat = dataSet.mat.^2;
    
    tic
    sib = corr(dataMat);
    tempMat = reshape(tiedrank(sib(:)), gCount, gCount);
    toc
    
    ramkSum = rankSum + tempMat;
    clear tempMat;
    
end

fileName = sprintf('/space/grp/marjan/data/corrVar/rankSum_%s.mat', ...
                   tissue);
save(fileName, 'rankSum', '-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MARCH21st presentation, node degree distribution
clear

h = figure();
hold all


tissue = 'skeletalMuscle'
netFolder = ['~/networks/' tissue '/singleNet/']
netFiles = dir([netFolder '*.mat' ])
for i = 1:length(netFiles)
    load([netFolder netFiles(i).name])     
   Snet = full(sparseSingleNet);
   [ndHist, bins] = networkNDdistPlot(Snet);
   plot(bins, log10(ndHist), 'lineWidth', 1, 'color', [180, 10, 20]/255);
end
xlabel('Node Degree')
ylabel('log10(geneCount)')
title('node degree histograms for all datasets')
print(h, '-depsc', '~/networks/wholeNDdist.eps')
print(h, '-dpdf', '~/networks/wholeNDdist.pdf')

h = figure();
hist(holu, 40)
title(['node degree distribution for random coexpression network - ' ...
       'nd > 10'])
xlabel('node degree')
ylabel('gene counts')
print(h, '-depsc', '~/networks/randomNDdist.eps')
print(h, '-dpdf', '~/networks/randomNDdist.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% p-value specificity for the networks

tissue = 'blood'
%dataFolder = '/space/grp/marjan/data/'
%exprFolder = [dataFolder tissue ['/geneExprMat/gemmaProbeMap/clean/' ...

dataFolder = '~/data/array/'
exprFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/']
fileList = dir([exprFolder '*.mat'])
gCount = 18494;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
samples = zeros(1000000, 53);
dcount = 0;
addpath('~/codes/MATLAB/myCodes/general/')
for t = 1: 5
    tissue = tissues{t};
    dataFolder = '~/data/affyArray/tissues/'
    exprFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/']
    fileList = dir([exprFolder '*.mat'])
    gCount = 18494;

    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])

    for i = 1: length(fileList)
        dcount = dcount + 1;
        %loading new file
        '1'
        dsQNormInf(dcount).tissue = tissue;
        dsQNormInf(dcount).name = ...
            getGSEIDfromStr(fileList(i).name);
        load([exprFolder fileList(i).name]);
        myMat = dataSet.mat(expGenesInd, :);
        sib = corr(myMat');
        sib = sib - eye(sum(expGenesInd));

        '2'
        %getting the file's network
        upperSingle = sib(logical(triu(ones(size(sib)), 1)));
        samples(:, dcount) = datasample(upperSingle, 1000000, ...
                                        'Replace', false);
        dcount
    end 
end

qnorms = quantilenorm(samples);
qqnorms = quantile(qnorms(:, 1), 999);
qsamples = zeros(1000, 53);
for i = 1: 53
    qsamples(:, i) = quantile(samples(:, i), 1000);
end

repValue = zeros(1000, 1);
values = qnorms(qnorms(: ,1) < qqnorms(1), 1);
repValues(1) = mean(values);
for i = 2: 999;
    s1 = qnorms(: ,1) < qqnorms(i);
    s2 = qnorms(: ,1) > qqnorms(i - 1);
    selected = (s1 + s2) > 1;
    values = qnorms(selected, 1);
    repValues(i) = mean(values);
end
values = qnorms(qnorms(: ,1) > qqnorms(i), 1);
repValues(1000) = mean(values);

save('~/data/general/correlationNorm/correlationQuantiles.mat', 'qsamples')
save('~/data/general/correlationNorm/qqnormVector.mat', 'qqnorms')
save('~/data/general/correlationNorm/representedValues.mat', 'repValues')

std1 = std(upperSingle1);
stdzied1 = upperSingle1 /std1;
z1 = zscore(upperSingle1);

hist(upperSingle1, 100)
h1 = figure
hist(stdzied1, 100)

upperSingle2 = sib(logical(triu(ones(size(sib)), 1)));
std2 = std(upperSingle2);
stdzied2 = upperSingle2 /std2;
z2 = zscore(upperSingle2);



h2 = figure;
hist(upperSingle2, 100)
h3 = figure
hist(stdzied2, 100)


QSingle(i) = quantile(upperSingle, (1 - thr))


% getting that plot thing
% >>>>>>>>

quantiles = quantile(upperSingle, 1000);
tic
h1 = figure
hist(upperSingle, 100);
toc


h = figure
tic
zs = zscore(upperSingle);
toc
hist(zs, 100)
% <<<<<<<<
