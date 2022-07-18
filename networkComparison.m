

%%%%% these are the analysis I do after the geneExpression analysis
%%%%% I did in the geneExpressionMain.m file. 
clear 
%data folder for the blood 
addpath '~/codes/MATLAB/myCodes/general/'

folder = cell(1, 4)
folder{1} = '~/data/array/blood/matFiles/geneExprMatGemmaMapBER/'
folder{2} = '~/data/array/lung/matFiles/geneExprMatGemmaMapBER/'
folder{3} = '~/data/array/skeletalMuscle/matFiles/geneExprMatGemmaMapBER/'
folder{4} = '~/data/array/liver/matFiles/geneExprMatGemmaMapBER/'

ts{1} = 'blood';
ts{2} = 'lung';
ts{3} = 'skeletalMuscle';
ts{4} = 'liver';

tissue.dataFolder = folder{1};
for i = 1:4
    tissue(i).dataFolder = folder{i};
    tissue(i).name = ts{i};
    tissue(i).ID = i;
    fileList = dir([folder{i} '*.mat']);
    temp = length(fileList);
    tissue(i).exprCount = temp;
end

gCount = 18494; % or, load a dataset and check it. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I would have 4 list of genes for each tissue : 
% 1. expressed in all datasets, (this one includes the next three -
% union of the three other groups is this one)
% 2. expressed only in all of this tissues' datasets
% 3. expressed in all other datasets 
% 4. expressed here and there

% my interest in the first level is "a particular tissue against
% all other tissues.", so, once working in one tissue, I don't care
% about the others. 

% extracting and saving the justTissueGenes for each dataset and
% tissue. to have a comparison for the datasets, I got the
% quantiles for each dataSet. 

% group 1 genes
muscleGenes 

% group 2 genes
justMuscleGenes 

% group 3 genes 
everyWhereGenes

% group 4 genes
book = muscleGenes - justMuscleGenes - everyWhereGenes;
hereAndThereGenes = book == 1;

% put the group two genes for each tissue here
tissue(1).justGenes = justBloodGenes;
tissue(2).justGenes = justLungGenes;
tissue(3).justGenes = justMuscleGenes;
tissue(4).justGenes = justLiverGenes;


%%%%% 

tissue(1).justTissueNetFolder = ['/space/grp/marjan/tempNetworks/' ...
                    'bloodGenes/'];
tissue(2).justTissueNetFolder = ['/space/grp/marjan/tempNetworks/' ...
                    'lungGenes/']
tissue(3).justTissueNetFolder = ['/space/grp/marjan/tempNetworks/' ...
                    'muscleGenes/'];
tissue(4).justTissueNetFolder = ['/space/grp/marjan/tempNetworks/' ...
                    'liverGenes/']


% open each dataset, build the network, select the muscle genes and
% their connections. 
exprCounter = 0;
quantiles = zeros(totalExprCount, 100);
for i = 1: length(tissue)% tissue Count - this tissue loop is here only
                         % to make me able to read each file
    fileList = dir([folder{i} '*.mat']); 
    tissue(i).name
    for j = 1: length(fileList) 
        exprCounter = exprCounter + 1;
        exprCounter
        % load file
        load([tissue(i).dataFolder fileList(j).name])
        name = getGSEIDfromStr(fileList(j).name)
        
        'data loaded'
        
        % make network
        sib = corr(dataSet.mat');
        sib = sib - eye(gCount);
        
        'corr done'
        
        % rank the edgevalues
        % tic
        % rankNet = reshape(tiedrank(sib(:)), gCount, gCount);
        % toc
        % tempNet = rankNet;
        
        % getting the quantiles
        tic
        book = sib > 0;
        holu = sib(book);
        q100 = quantile(holu, [1:99, 99.5]/100);
        quantiles(exprCounter, :) = q100;
        q100(100)
        toc
                
        
        % now that we have each file, for each file we save the
        % network of its four justTissue genes. 
        for k = 1:length(tissue)
            tissue(k).name
            tissue(k).justTissueNetFolder
            tempGenes = tissue(k).justGenes;
            netSlice = sib(tempGenes, :);
            size(netSlice)
            newNetSlice = zeros(size(netSlice));
            'a'
            for l = 1:99 % to the lenght of quantiles
                newNetSlice(netSlice >= q100(l) & netSlice < q100(l+1)) = l/100;
            end
            'b'
            newNetSlice(netSlice >= q100(100)) = 1;
            saveFolder = [tissue(k).justTissueNetFolder];
            save([saveFolder tissue(i).name '_' name '_netSlice.mat'], 'newNetSlice', ...
                 '-v7.3');
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now loading the networks and comparing them. 
book = sum(exprGenes, 2);
everyWhereGenes = book == 43;
sum(everyWhereGenes);

folder = '/space/grp/marjan/tempNetworks/liverGenes/'
fileList = dir([folder '*.mat'])

currentTissue = 'liver'

load([folder fileList(12).name])
liver01 = newNetSlice;
liver01Bar = hist(liver01(:), 10)
liver01JE = liver01(:, everyWhereGenes);

load([folder fileList(13).name])
liver02 = newNetSlice;
liver02Bar= hist(liver02(:), 10)
liver02JE = liver02(:, everyWhereGenes);

book = currentSlice(:, everyWhereGenes);
size(book)

load([folder fileList(40).name])
muscle01 = newNetSlice;
muscle01Bar = hist(muscle01(:), 10)
muscle021JE = muscle01(:, everyWhereGenes);

bar([muscle01Bar' liver01Bar' liver02Bar'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% making single binary networks 

exprCounter = 0;
for i = 1: length(tissue)% tissue Count - this tissue loop is here only
                         % to make me able to read each file
    fileList = dir([folder{i} '*.mat']); 
    tissue(i).name
    for j = 1: length(fileList) 
        exprCounter = exprCounter + 1
        % load file
        load([tissue(i).dataFolder fileList(j).name])
        
        name = getGSEIDfromStr(fileList(j).name)
        sib = corr(dataSet.mat');
        sib = sib - eye(gCount);
        
        '2'
        %getting the file's network
        upperSingle = triu(sib, 1);
        QSingle(exprCounter) = quantile(upperSingle(:), (1 - 0.005))
        % singleNet = sib > QSingle(exprCounter);
        % sparseSingleNet = sparse(singleNet);
        
        % saveFolder = ['/space/grp/marjan/singleSparceNets/' ...
        %               tissue(i).name '/'];
        % save([saveFolder tissue(i).name '_' name '_sparceNet.mat'], 'sparseSingleNet', ...
        %      '-v7.3');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loading the sparce nets and comparing the node degree of TS genes
% in those networks. 

% find the node degree of muscle genes for each of the datasets. 

exprCounter = 0;
networkFolder = ['/space/grp/marjan/singleSparceNets/']
nd = zeros(gCount, totalExprCount);
currentTissue = 1;
justGeneNets = cell(1, totalExprCount);
for i = 1: length(tissue)                         
    fileList = dir([networkFolder tissue(i).name '/*.mat'])
    tissue(i).name
    
    for j = 1: length(fileList) 
        exprCounter = exprCounter + 1
        % load file
        load([networkFolder tissue(i).name '/' fileList(j).name])
        name = getGSEIDfromStr(fileList(j).name)
        
        % getting the node degree
        % nd(:, exprCounter) = sum(sparseSingleNet);
        % max(nd(:, exprCounter))
        % min(nd(:, exprCounter))
        
        %getting the justgenes networks
        tempNet = full(sparseSingleNet);
        justGeneNets{exprCounter} = ...
            tempNet(tissue(currentTissue).justGenes, tissue(currentTissue).justGenes);
    end
end

msGenesND = nd(tissue(1).justGenes, :);

book = mean(msGenesND(:, 1), 2);
[a, b] = sort(book, 'descend');

sorted = msGenesND(b, :);

    h = figure; 
    heatmap(sorted, [1:43],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
            true, 'Colorbar', true, 'Colormap', 'summer')
    'heatmap done'
    title(sprintf('%s  Number of samples : %d', SName, sampCount))
    fileName = ['cleanGeneSampleCorr' SName];
    print(h, '-depsc', [figFolder fileName '.eps']);
    print(h, '-dpdf', [figFolder fileName '.pdf']);
    
tsNet = zeros(667, 667);
for i = 1:9 
    tsNet = tsNet + justGeneNets{i};
end
tsNet = tsNet ./ 9;

otherNet = zeros(667, 667);
for i = 10:43
    otherNet = otherNet + justGeneNets{i};
end
otherNet = otherNet ./ 34;

    h = figure; 
    heatmap(tsNet, [],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
            true, 'Colorbar', true, 'Colormap', 'summer')
    
figure
bloodC = hist(tsNet(:) , 10)
figure
otherC = hist(otherNet(:), 10)

bar([bloodC' otherC'], 'group')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each dataset, I have the expressed genes (run the code)

% I will open the network, get the node degree of lower expressed
% genes. what is the average? compared to network average? how many of the edges are shared
% between the lower expressed genes? 
expGenes02 = exprGenes;


kado = ones(size(exprGenes));
noExprGenes = kado - exprGenes;

exprCounter = 0;
networkFolder = ['/space/grp/marjan/singleSparceNets/']
nd = zeros(gCount, totalExprCount);
edgeCount = zeros(1, 43);
for i = 1: length(tissue)                         
    fileList = dir([networkFolder tissue(i).name '/*.mat'])
    tissue(i).name
    
    for j = 1: length(fileList) 
        exprCounter = exprCounter + 1
        % load file
        load([networkFolder tissue(i).name '/' fileList(j).name])
        name = getGSEIDfromStr(fileList(j).name)
        
        noGeneNet = sparseSingleNet(logical(noExprGenes(:,exprCounter)), ...
                                    logical(noExprGenes(:,exprCounter)));
        temp = sum(sum(noGeneNet))/2
        edgeCount(exprCounter) = temp;
        % getting the node degree
                 nd(:, exprCounter) = sum(sparseSingleNet);
        % max(nd(:, exprCounter))
        % min(nd(:, exprCounter))
        
        %getting the justgenes networks
        % tempNet = full(sparseSingleNet);
        % justGeneNets{exprCounter} = ...
        %     tempNet(tissue(currentTissue).justGenes,
        %     tissue(currentTissue).justGenes);
    end
end

totalEdge = (gCount * gCount * 0.005) / 2
book = edgeCount/ totalEdge %4.5
hist1 = hist(book, 10);

khiar = edgeCount / totalEdge % 5.25
hist2 = hist(khiar, 10);

holu = edgeCount / totalEdge % 6 
hist3 = hist(holu, 10)
colormap('summer')
h = figure
hist(holu, 10)

h = figure
havij = mean(all, 2)
bar([4, 5, 6]', havij', 'grouped')
xlabel('threshold for the genes')
ylabel('percent of edges in the sparse network for the cut off genes')