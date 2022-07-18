% In this file I have the code for the part "Characterizing the
% links" in the project report. 
% 1. edges shared between every two networks at different
% 2. distribution of number of repeated edges for each tissue. 3*5
% histograms.
% 3. getting sum of all the binary networks for different
% thresholds.  
% 4. saving the summed networks for each tissue
% 5. getting the tissue distance for each edge.

% Single binary networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%
%% 1. edges shared between every two networks at different
%% thresholds
% >>> 1.1 this part is coded in folderNetwork.m, part 6, line 310 ish
% of the file. folderNetwork is where I make the binary networks. 
% >>> 1.2 having the expression involved - same place. 

%%%%
%% 2. distribution of number of repeated edges for each tissue. 3*5
%% histograms. 

clear 
%data folder for the blood 
gCount = 18494;
dataFolder = '~/networks/tissues/'
%tissue = 'blood';

ts{1} = 'blood';
ts{2} = 'lung';
ts{3} = 'skeletalMuscle';
ts{4} = 'liver';
ts{5} = 'brain';

thr = {'01', '001', '005'}

% on tissues
for t = 1:length(ts)
    tissue = ts{t}

    % on thresholds
    for h = 1:length(thr)

        netFolder = [dataFolder tissue '/singleNet/']
        fileList = dir([netFolder 'GSE*_' thr{h} '.mat'])
        fc = length(fileList);
        sumNet = zeros(gCount, gCount);

        GSEIDs = cell(1, length(fileList));
        remainEdgePercent = zeros(1, length(fileList));

        for i = 1:fc
            i
            load([netFolder fileList(i).name]);
            tNet = full(sparseSingleNet);
            tNet = tNet - eye(gCount, gCount);

            % comment the following part if you don't want to
            % filter edges for exp level 
            % >>>>>>>>>>
            [a, b] = regexp(fileList(i).name, 'GSE[0123456789]+');
            GSEID = fileList(i).name(a:b);
            
            GSEIDs{i} = GSEID;

            netStrFile = [netFolder 'netStr_' GSEID '.mat'];
            load(netStrFile);

            expNet = abs(netStr.expNet);
            temp1 = expNet >=35;
            temp2 = expNet == 25;
            temp = temp1 + temp2;
            iexpNet = temp;
            clear expNet temp1 temp2 temp;

            tNet = (iexpNet .* tNet) > 0;
            remainEdgePercent(i) = sum(sum(tNet))/ (gCount * gCount / 2);
            % <<<<<<<<<<<

            sum(sum(tNet))/ (gCount * gCount / 2)
            sumNet = sumNet + tNet;
        end

        % uv = unique(sumNet);
        % edgeSumCount = histc(sumNet(:), uv);

        % saveFolder = [dataFolder tissue '/singleNet/dataFiles/']
        % NOTICE! check the file name, if it is EXP level filter
        % included or not
        % fileName = [saveFolder 'sumNetEdgeCount_ExpThr5_' thr{h} '.mat']
        % save(fileName, 'edgeSumCount');

        % fileName = [saveFolder 'sumNetEdgeCount_' thr{h} '.mat']
        % save(fileName, 'edgeSumCount');


        % >>>> exp filtered info
        % expFilStr.GSEIDs = GSEIDs;
        % expFilStr.edgePercent = remainEdgePercent;
        % fileName = [saveFolder 'expThrRemainEdgePercent_ExpThr5_' thr{h} '.mat']
        % save(fileName, 'expFilStr');
        % <<<< 

        %%%% doing the loop for each dataset to get the dist of
        %%%% edge count for that dataset

        for i = 1:fc
            i
            load([netFolder fileList(i).name]);
            tNet = full(sparseSingleNet);
            tNet = tNet - eye(gCount, gCount);

            % comment the following part if you don't want to
            % filter edges for exp level 
            % >>>>>>>>>>
            [a, b] = regexp(fileList(i).name, 'GSE[0123456789]+');
            GSEID = fileList(i).name(a:b);
            
            GSEIDs{i} = GSEID;

            netStrFile = [netFolder 'netStr_' GSEID '.mat'];
            load(netStrFile);

            expNet = abs(netStr.expNet);
            temp1 = expNet >=35;
            temp2 = expNet == 25;
            temp = temp1 + temp2;
            iexpNet = temp;
            clear expNet temp1 temp2 temp;

            tNet = (iexpNet .* tNet) > 0;
            remainEdgePercent(i) = sum(sum(tNet))/ (gCount * gCount / 2);
            % <<<<<<<<<<<

            tNet = full(tNet);

            tNet = tNet .* sumNet; 

            uv = unique(tNet);
            indEdgeSumCount = histc(tNet(:), uv);

            saveFolder = [dataFolder tissue '/singleNet/dataFiles/']
            % NOTICE! check the file name, if it is EXP level filter
            % included or not
            fileName = [saveFolder 'IndividualSumNetEdgeCount_ExpThr5_' GSEID ...
                        '_' thr{h} '.mat']

            % fileName = [saveFolder 'IndividualSumNetEdgeCount_' thr{h} ...
            %             '_' GSEID '.mat']
            % 
            save(fileName, 'indEdgeSumCount');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.5 exploratory observation for the above results

        saveFolder = [dataFolder tissue '/singleNet/dataFiles/']
        fileName = [saveFolder 'sumNetEdgeCount_ExpThr5_' thr{h} ...
                    '.mat']
        load([saveFolder 'sumNetEdgeCount_01.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. getting sum of all the binary networks for different
% thresholds. 

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
 
thr = '005'
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


GSEIDs = cell(1, length(netStrFiles));
for i = 1:length(netStrFiles)
    [a, b] = regexp(netStrFiles(i), 'GSE[0123456789]+');
    temp = netStrFiles{i};
    GSEIDs{i} = temp(a{1}:b{1});
end

save('~/networks/allGSEIDs.mat', 'GSEIDs')

%getting the added net
gCount = 18494
sumNet = zeros(gCount, gCount);
for i = 1:length(netFiles)
    i
    load(netFiles{i})
    iNet = full(sparseSingleNet);
    iNet = iNet - eye(gCount);

    % [a, b] = regexp(netFiles{i}, 'GSE[0123456789]+');
    % GSEID = netFiles{i}(a:b)

    % GSEIDs{i} = GSEID
    
    load(netStrFiles{i})

    % getting links with expr > a certain level
    expNet = abs(netStr.expNet);
    temp1 = expNet >=35;
    temp2 = expNet == 25;
    temp = temp1 + temp2;
    iexpNet = temp;
    clear expNet temp1 temp2 temp;

    iNet = (iexpNet .* iNet) > 0;
    
    sumNet = sumNet + iNet;
end

thr
sum(sum(sumNet > 2)) / (gCount * gCount /2)

        sumNet = sparse(sumNet);
        dataFolder = '~/networks/tissues/'
        fileName = [dataFolder 'allSumNetworkExpThr5_' thr '.mat'] 
        save(fileName, 'sumNet', '-v7.3');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. saving the summed networks for each tissue

clear 
%data folder for the blood 
gCount = 18494;
dataFolder = '~/networks/tissues/'
%tissue = 'blood';

ts{1} = 'blood';
ts{2} = 'lung';
ts{3} = 'skeletalMuscle';
ts{4} = 'liver';
ts{5} = 'brain';

thr = {'01', '001', '005'}

thr = {'0.05', '0.10'};
thr = {'0.05'}

% on tissues

% 4.5 saving the agg network with already corr filtered links for
% each tissue. EXP FILTER SHOULD BE MENTIONED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
%data folder for the blood 
gCount = 18494;
dataFolder = '~/networks/tissues/'
%tissue = 'blood';

ts{1} = 'blood';
ts{2} = 'lung';
ts{3} = 'skeletalMuscle';
ts{4} = 'liver';
ts{5} = 'brain';

thr = {'0.010', '0.005', '0.050', '0.100'}
expThr = '0.8'

%%%%%%%%%%%%%%
%% IMPORTANT: select if it is the ALL network.
%%%%%%%%%%%%%%

for t = 3:4%length(ts)
    tissue = ts{t}

    % on thresholds
    for h = 1:length(thr)

        netFolder = [dataFolder tissue '/singleNet/']

        %%%%%%%%%%%%%%%%%%%%%%% pick one
        %fileList = dir([netFolder tissue 'NetExpThrAll' expThr '_GSE*_' thr{h} ...
        %               '.mat'])

        fileList = dir([netFolder tissue 'NetExpThr' expThr '_GSE*_' thr{h} ...
                         '.mat'])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fc = length(fileList);
        sumNet = zeros(gCount, gCount);

        GSEIDs = cell(1, length(fileList));

        for i = 1:fc
            i
            load([netFolder fileList(i).name]);
            % sparseSingleNet(9402, 15969)
            % sparseSingleNet(15969, 9402)
            sumNet = sumNet + sparseSingleNet;
        end

        sumNet = sparse(sumNet);

        %%% Pick one
        %        fileName = [dataFolder tissue '/' tissue
        %        'SumNetworkExpThrAll' expThr '_' ...
        %           thr{h} '.mat'] 

        fileName = [dataFolder tissue '/' tissue 'SumNetworkExpThr' expThr '_' ...
                    thr{h} '.mat'] 
        %%%
        
        save(fileName, 'sumNet', '-v7.3');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. getting the tissue distance for each edge. edges are coming
% from part 3 and the distances are coming from part 4. I have 3
% sets of edges based on the coexpression thr from part 3. For each
% of those edges I have values from my tissue sum networks from
% part 4. 

% for each network from part 3, I get the 5 edge values for each
% tissue, which is their distance for that tissue. 

clear 
%data folder for the blood 
gCount = 18494;
dataFolder = '~/networks/tissues/'
%tissue = 'blood';
thr = '_01'
 
ts{1} = 'blood';
ts{2} = 'lung';
ts{3} = 'skeletalMuscle';
ts{4} = 'liver';
ts{5} = 'brain';

% for each sumAllNetwork:

fileName = dir([dataFolder 'allSumNet*' thr '.mat'])

load([dataFolder fileName(1).name]);

kado = sumNet > 0; % all the edges that matter
sum(sum(kado)) / (gCount * gCount /2)

for i = 1 : length(ts)
    disTissue = i;  % tissue is blood

    % getting the distance from the main tissue
    fileName = dir([dataFolder ts{disTissue} '/sumNet*' thr '.mat'])
    load([dataFolder ts{disTissue} '/' fileName.name])
    distSumNet = sumNet;
    mv = max(max(distSumNet))

    distSumNet = 1 - distSumNet / mv;
    distSumNet = distSumNet .* kado;

    for t = 1 : length(ts)
        t
        if (t ~= disTissue)
            fileName = dir([dataFolder ts{t} '/sumNet*' thr '.mat'])
            load([dataFolder ts{t} '/' fileName.name])
            tempNet = sumNet;
            mv = max(max(tempNet));
            tempNet = tempNet/ mv;
            tempNet = tempNet .* kado;
            distSumNet = distSumNet + tempNet;
        end
    end
    sum(sum(distSumNet >0)) / (gCount * gCount / 2)
    % saving the blood distance for kado
    distSumNet = sparse(distSumNet);
    fileName = [dataFolder ts{disTissue} '/tissueDistance_' thr '.mat'] 
    save(fileName, 'distSumNet', '-v7.3');

end

distSumNet(kado) = distSumNet(kado) + 0.001;

book = sparse(distSumNet);

[x, y, v] = find(book);

[s, ind] = sort(v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. load the whole matrix of distances for the selected thr,
% record which tissue has the shortest distance from each link and
% then sort them (for each tissue) based on their distance, and
% then pick that link as that tissues link. THEN , between the
% links of each  tissue, pick the ones which has the shortest
% distance. 

% a matrix of links, their tissue distances, their selected
% tissue. (a matrix of 5 columns)

clear 
%data folder for the blood 
gCount = 18494;
dataFolder = '~/networks/tissues/'
%tissue = 'blood';
thr = '001'
 
ts{1} = 'blood';
ts{2} = 'lung';
ts{3} = 'skeletalMuscle';
ts{4} = 'liver';
ts{5} = 'brain';

% for each sumAllNetwork:

fileName = dir([dataFolder 'allSumNet*' thr '.mat'])

load([dataFolder fileName(1).name]);

kado = sumNet > 0; % all the edges that matter
sum(sum(kado)) / (gCount * gCount /2)
clear sumNet;

edgeCount = sum(sum(kado));

edgeInfo = zeros(edgeCount, 9);

[edgeInfo(:, 1), edgeInfo(:, 2)]  = find(kado);

clear kado
for i = 1:length(ts)

    load([dataFolder ts{i} '/tissueDistance_' thr '.mat'])
    i
    for j = 1: length(edgeInfo)
        a = edgeInfo(j, 1);
        b = edgeInfo(j, 2);
        edgeInfo(j, (i+2)) = distSumNet(a, b);
    end
end

%%%%%%%%%%%
% sub part, checking to see if it is as expected 
ts{1} = 'blood';
ts{2} = 'lung';
ts{3} = 'skeletalMuscle';
ts{4} = 'liver';
ts{5} = 'brain';

k = 12233
edgeInfo(k, 1:2)
edgeInfo(k, 3:end)

book = zeros(5, 1);
for t = 1:5
    fileName = dir([dataFolder ts{t} '/sumNet*' thr '.mat']);
    load([dataFolder ts{t} '/' fileName.name])
    a = edgeInfo(k, 1);
    b = edgeInfo(k, 2);
    book(t) = full(sumNet(a, b));
end

edgeInfo(k, 3:end)
book'

%%%%%%%%

% getting the minDist tissue for each edge, basically, getting ITS
% tissue for each edge and the fold 

[a, b] = min(edgeInfo(:, 3:7)');
edgeInfo(:, 8) = b;

sib = hist(edgeInfo(:,8), [1:5])
ts

%histogram of the distances for each tissue, then I pick which edges ...
%       I like to pick!

t = 5
sib = edgeInfo(:, 8) == t;

temp = edgeInfo(sib, :);
whos temp

[a, b] = sort(temp(:,t + 2));

newblei = temp(b, :);

bloodEdgeInfo = newblei;
lungEdgeInfo = newblei;
msEdgeInfo = newblei;
liverEdgeInfo = newblei;
brainEdgeInfo = newblei;

edgeInfo.blood = bloodEdgeInfo;
edgeInfo.lung = lungEdgeInfo;
edgeInfo.brain = brainEdgeInfo;
edgeInfo.liver = liverEdgeInfo;
edgeInfo.sm = smEdgeInfo;

save(['~/networks/tissues/edgeInfo' thr '.mat'], 'edgeInfo')

% save them. 
% cluster them based on their exp profile
% now, general characteristics of the networks. in what cluster do
% they belong to. what is their characteristics in the network. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%7. a
clear 
%data folder for the blood 
gCount = 18494;
dataFolder = '~/networks/tissues/'
%tissue = 'blood';
thr = '001'
 
ts{1} = 'blood';
ts{2} = 'lung';
ts{3} = 'skeletalMuscle';
ts{4} = 'liver';
ts{5} = 'brain';


load(['~/networks/tissues/edgeInfo' thr '.mat'])

mainTissue = 1;
tempEI = edgeInfo.blood;

%length of edges we are getting
n = 100000;
if (length(tempEI) < n)
    n = length(tempEI)
end

clear edgeInfo
edgeExpInf = zeros(n, 54);
edgeExpInf(:, 1:2) = tempEI(1:n, 1:2);
clear tempEI
fc = 3;
for t=1:length(ts)
    t
    fileList = dir([dataFolder ts{t} '/singleNet/netStr_*.mat'])
    for f = 1:length(fileList)
        f
        load([dataFolder ts{t} '/singleNet/' fileList(f).name])
        explevel = netStr.qExpr;
        clear netStr;
        a = edgeExpInf(1:n, 1);
        b = edgeExpInf(1:n, 2);
        tic
        for j = 1:n
            edgeExpInf(j, fc) = explevel(a(j)) * explevel(b(j));
        end
        toc
        fc = fc + 1;
    end
end

save(['~/networks/tissues/' ts{mainTissue} '/edgeExpInfo' thr '.mat'], 'edgeExpInf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the list of house keeping genes, find their number in my
% gene mapping AND SAVE THEM! 
clear

dataFolder = '~/data/general/'
load([dataFolder 'GPL570GemmaMapNEW.mat'])
gCount = length(gpl570.uniqueSymbols);

file = '~/data/HK_genes.txt'
fid = fopen(file)

data = textscan(fid, repmat('%s', 1, 2))

myMap = containers.Map(gpl570.uniqueSymbols, 1:gCount);

hkInd = zeros(1,length(data{1}));
for i = 1:length(data{1})
    try
        temp = values(myMap, data{1}(i));
        hkInd(i) = temp{1};
    catch
        hkInd(i) = 0;
        continue
    end
end

hkgInd = hkInd(hkInd > 0);
length(hkInd) % 3581 I have 3581 hk genes out of 3804 hk genes
              % provided in the list. 
save('~/data/general/HKGInd.mat', 'hkgInd')


%%%%%%%%%%%%%%%%%%%
% load the list of house keeping edges and find the number of genes
% in my mapping. SAVE THE SUPER SPARSE NETWORK


clear

dataFolder = '~/data/general/'
load([dataFolder 'GPL570GemmaMapNEW.mat'])
gCount = length(gpl570.uniqueSymbols);

file = '~/data/pairs.txt'
fid = fopen(file)

linum = 1
data = textscan(fid, repmat('%s', 1, 3), 'Headerlines', linum-1)

[data{1}(1:10) data{2}(1:10) data{3}(1:10)]


myMap = containers.Map(gpl570.uniqueIDs, 1:gCount);

edgeInd = zeros(length(data{1})-1, 3);
i = 1
for i = 2:length(data{1})
    try
        temp = values(myMap, data{2}(i));
        edgeInd( i-1, 2) = temp{1};
        temp = values(myMap, data{3}(i));
        edgeInd( i-1, 3) = temp{1};
        edgeInd(i-1, 1) = str2num(data{1}{i}(2:(end-1)));
    catch
        edgeInd(i-1, :) = 0;
        continue
    end
end

hkeInd = edgeInd(edgeInd(:, 1) > 0, :);
length(hkeInd) % 2398 of 2670 edges here I have both their genes
              % and therefor can study them. 

% NOTICE all my matrices should be the lower triangle
for i = 1:length(hkeInd)
    if (hkeInd(i, 2) < hkeInd(i, 3))
        temp = hkeInd(i, 2);
        hkeInd(i, 2) = hkeInd(i, 3);
        hkeInd(i, 3) = temp;
    end
end

book = zeros(gCount, gCount);
for i = 1:length(hkeInd)
    book(hkeInd(i,2), hkeInd(i,3)) = 1;
end
hkeNet = sparse(book);
clear book
save('~/data/general/hkeNet.mat' , 'hkeNet')

% checking for the reproducibility of the HKE in my data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the first part is simple repeats 
clear

%set of genes I want 
load('~/data/general/hkeNet.mat')
load('~/data/general/hkgInd.mat')
load('~/data/general/ribosomeIDsGO.mat')

subNet = hkeNet;

thr = '005'
load(['~/networks/tissues/allSumNetworkExpThr5_' thr '.mat'] )

subNet = sumNet(hkgInd, hkgInd);
sib = subNet;

subNet = sumNet(ribIDs, ribIDs);
sib = subNet;

sib = hkeNet .* sumNet;

[r, c, v] = find(sib);


max(v)
hkeHist = zeros(1, max(v));
for i = 1:max(v) 
    hkeHist(i) = sum(v == i);
end
save(['~/data/general/hke_' thr '_repHist.mat'], 'hkeHist')

hkgNetHist = zeros(1, max(v));
for i = 1:max(v) 
    hkgNetHist(i) = sum(v == i);
end
save(['~/data/general/hkgNet_' thr '_repHist.mat'], 'hkgNetHist')

max(v)
ribHist = zeros(1, max(v));
for i = 1:max(v) 
    ribHist(i) = sum(v == i);
end
save(['~/data/general/ribosome_' thr '_repHist.mat'], 'ribHist')


[r, c, v] = find(sumNet);
alleHist = zeros(1, max(v));
for i = 1:max(v) 
    alleHist(i) = sum(v == i);
end
save(['~/networks/tissues/allSumNetExpThr5_' thr '_repHist.mat'], ...
     'alleHist')

% this second part is where I count .6 and .7 for each tissue. 

clear
repThr = 0.4;
load('~/data/general/hkeNet.mat')
thr = '01'
 
ts{1} = 'blood';
ts{2} = 'lung';
ts{3} = 'skeletalMuscle';
ts{4} = 'liver';
ts{5} = 'brain';

gCount = 18494;
totNet = zeros(gCount, gCount);
id = [1 2 4 8 16];

for t = 1:length(ts)
    t
    load(['~/networks/tissues/' ts{t} '/sumNetworkExpThr5_' thr ...
          '.mat'])

    max(max(sumNet));

    sumNet = sumNet ./ max(max(sumNet));

    book = sumNet > repThr;
    %[r, c] = find(book);

    book = book .* id(t);

    totNet = totNet + book;
end

sum(sum(totNet == 31))

hkeNet = hkeNet .* 32;
sib = hkeNet + totNet;

sumSet = [3, 5 , 9, 17, 33, 6, 10, 18, 34, 12, 20, 36, 24, 40, 48, ...
          30, 29, 27, 23, 15, 31, ];

edgeHist = zeros(1, length(sumSet));

for i = 1:length(sumSet)
    i
    edgeHist(i) = sum(sum(sib == sumSet(i))); 
end

%%%%% checking the links for the HKG
clear

load('~/data/general/hkgInd.mat')

repThr = 0.4
thr = '01'
 
ts{1} = 'blood';
ts{2} = 'lung';
ts{3} = 'skeletalMuscle';
ts{4} = 'liver';
ts{5} = 'brain';

gCount = length(hkgInd);
totNet = zeros(gCount, gCount);
id = [1 2 4 8 16];

edgeCount = zeros(1, 5);
for t = 1:length(ts)
    t
    load(['~/networks/tissues/' ts{t} '/sumNetworkExpThr5_' thr ...
          '.mat'])

    sumNet = sumNet(hkgInd, hkgInd);

    sumNet = sumNet ./ max(max(sumNet));

    book = sumNet > repThr;
    %[r, c] = find(book);

    edgeCount(t) = sum(sum(book));

    book = book .* id(t);

    totNet = totNet + book;
end
ts
edgeCount
sum(sum(totNet == 31))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 7.5
% given pairs of IDs, get all the real corrvalues from the data and
% the exp levels, and the HKE +-, and the HKG +-

% load it from the netStr
clear 
%data folder for the blood 
gCount = 18494;
dataFolder = '~/networks/tissues/'
%tissue = 'blood';
thr = '001'
 
ts{1} = 'blood';
ts{2} = 'lung';
ts{3} = 'skeletalMuscle';
ts{4} = 'liver';
ts{5} = 'brain';

load(['~/networks/tissues/edgeInfo' thr '.mat'])

mainTissue = 5;
tempEI = edgeInfo.brain;
n = 3800;

clear edgeInfo
% 2 for index, 52 for corr, 52 for corrBin, 52 for exp, 2 for each gene HK or not,
% 1  for my edge HK or not. 
% also, for each edge all the expressions sorted - based on their cluster
edgeCorrInf = zeros(n + 1, 109+52);
edgeCorrInf(2:end, 1:2) = tempEI(1:n, 1:2);
clear tempEI
fc = 3;
for t=1:length(ts)
    t
    fileList = dir([dataFolder ts{t} '/singleNet/netStr_*.mat'])
    for f = 1:length(fileList)
        f
        load([dataFolder ts{t} '/singleNet/' fileList(f).name])
        edgeCorrInf(1, fc) = t;
        edgeCorrInf(1, fc + 1) = t;
        for e=2:n+1
            a = edgeCorrInf(e, 1);
            b = edgeCorrInf(e, 2);
            edgeCorrInf(e, fc) = netStr.coNt(a, b);
            edgeCorrInf(e, (fc + 1)) = netStr.rankNet(a, b);
            edgeCorrInf(e, (fc + 2)) = netStr.qExpr(a) * ...
                netStr.qExpr(b);
            if (netStr.qExpr(a) < netStr.qExpr(b))
                edgeCorrInf(e, (fc + 2)) = edgeCorrInf(e, (fc + 2)) ...
                    * (-1);
            end
        end
        fc = fc + 3;
    end
end

load('~/data/general/hkeNet.mat')
% filling the HKE info
for i = 2:length(edgeCorrInf)
    if (hkeNet(edgeCorrInf(i, 1), edgeCorrInf(i, 2)) == 1)
        edgeCorrInf(i, 107 + 52) == 1;
    end
end

load('~/data/general/hkgInd.mat')
% filling the HKG info
[a, b] = ismember(edgeCorrInf(2:end, 1), hkgInd);
sum(a)
edgeCorrInf(2:end, 108 + 52) = a;

[a, b] = ismember(edgeCorrInf(2:end, 2), hkgInd);
sum(a)
edgeCorrInf(2:end, 109+52) = a;

% how many HKG do we have?
gset1 = edgeCorrInf(edgeCorrInf(:, 108) == 1, 1);
gset2 = edgeCorrInf(edgeCorrInf(:, 109) == 1, 2);
hkgset = unique([gset1', gset2']);

% the experimentally identified PPIs. REMAINS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster the edges based on their expression levels before saving
%% them. 
% what I do here is I cluster these edges based on their expression
% levels and then sort the clusters in the file. This is to get a
% visualization of the huge table I am going to observe. 
%% the following piece of code is taken from my desktop local code,
%% projectDoctPlots.m. thas is where I plot the heatmap of the
%% expression levels
clear
dataFolder = '~/networks/tissues/'

tissue = 'brain';
files = dir([dataFolder tissue '/edge*.mat'])

% loading the expression levels file for  the edges
load([dataFolder tissue '/' files(1).name])

% loading the edge info for the distances (I already have the edges
% that I want in the edgeCorrInf, but this is how I do it in the projectDoctPlots.m)
files = dir([dataFolder '/edge*.mat'])
load([dataFolder 'edgeInfo001.mat'])

brainEdgeInf = edgeInfo.brain;
clear edgeInfo
edgeCount = 3800
clusterCount = 50;
ts = 5

% clustering
Y = pdist(edgeExpInf(1:edgeCount, 3:end));
Z = linkage(Y, 'complete');

%dendrogram(Z)
T = cluster(Z, 'maxclust', clusterCount);

[a, b] = sort(T);

book = edgeCorrInf(2:end, :);
book = book(b, :); % clustering sort

holu = edgeExpInf(1:edgeCount, 3:end);
clusteredEdgeExpInf = holu(b, :);

% getting the edge distance and adding it to the file as well
edgeDist = brainEdgeInf(1:edgeCount, 7);
edgeColor = edgeDist(b); % clustering sort
size(edgeColor)
edgeColor = [0, edgeColor'];

edgeCorrInf(2:end, :) = book;
edgeCorrInf = [edgeCorrInf, edgeColor'];


% now for each edge I have the gene IDs, the exp level in each data
% set, the correlation in each dataset, (106), the HKE and HKG (109)
% and the distance for that tissue (110) and the correlation bin
csvwrite('~/networks/tissues/brainTFEInfo3800_thr001.csv', ...
         edgeCorrInf)

save('~/networks/tissues/brainTFEINfo3800_thr001.mat', ...
     'edgeCorrInf');

% getting the expression values for genes in my links. 

clear
load('~/networks/tissues/brainTFEINfo3800_thr001.mat');

% or any other gene list
kado =[ edgeCorrInf(2:end, 1)', edgeCorrInf(2:end, 2)']; 
book = sort(unique(kado));

ts{1} = 'blood';
ts{2} = 'lung';
ts{3} = 'skeletalMuscle';
ts{4} = 'liver';
ts{5} = 'brain';

expLevels = zeros(length(book)+2, 5000);
sc = 0; 
for t = 1:length(ts)

    tissue = ts{t};
    files = dir(['~/data/affyArray/tissues/' tissue ['/matFiles/' ...
                        'geneExprMatGemmaMapBER/*.mat']])

    for i = 1:length(files)
        i
        load(['~/data/affyArray/tissues/' tissue '/matFiles/' ...
                            'geneExprMatGemmaMapBER/' ...
                            files(i).name])
        files(i).name
        
        temp = size(dataSet.mat, 2);

        sc + 1
        sc + temp
        expLevels(1, sc+1:sc+temp) = t;
        expLevels(2, sc+1:sc+temp) = i;
        expLevels(3:end, sc+1: sc+temp) = dataSet.mat(book, :);
        sc = sc + temp;
    end
end

finalExpLevels = expLevels(:, 1:sc);

book = [0 , 0, book];
finalExpLevels = [book', finalExpLevels];

csvwrite('~/networks/tissues/brain3800_thr001_expLevels.csv', ...
         finalExpLevels)

save('~/networks/tissues/brain3800_thr001_expLevels.mat', 'finalExpLevels')

% printing the gene symbols in a file
load(['~/data/general/GPL570GemmaMapNEW.mat'])
sib1 = gpl570.uniqueSymbols(edgeCorrInf(2:end, 1));
sib2 = gpl570.uniqueSymbols(edgeCorrInf(2:end, 2));
file = [dataFolder 'brain3800edges_thr001Clustered.txt']
fid = fopen(file, 'w')
for i = 1:edgeCount
   fprintf(fid, [repmat('%s\t', 1, 2) '\n'], sib1{i}, sib2{i}) ;
end

% plotting the heatmap 
addpath('~/codes/MATLAB/myCodes/general/')
figFolder = ['~/networks/tissues/brain/figures/']

h = figure
heatmap(edgeColor(2:end)', [], [], [], 'Colorbar', true, 'Colormap', ...
        summer)

title([tissue 'colorbar for the distance']);
fileName = [tissue '_edgeDistColor001'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);

% add the numbers for edges 
d = length(clusteredEdgeExpInf)/200;
d = floor(d) * 10;
d = 200
labelsInd = [0:d:length(clusteredEdgeExpInf)];
labelsInd(1) = 1;

labels = cell(1, length(clusteredEdgeExpInf));
for i = 1:length(labelsInd)
    labels{labelsInd(i)} = labelsInd(i);
end


h = figure
heatmap(clusteredEdgeExpInf, [], labels, [], 'ShowAllTicks', true, 'Colorbar', true, ...
        'Colormap', bone )

title(sprintf('%s edge expression levels - edges %d', tissue, edgeCount))
fileName = [tissue '_edgeExp001'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% using the .mat files in ~/data/general/linkExprInfo/ for getting
% expression levels - this way I have the probe level expressions
% getting the expression values for genes in my links. 

clear

%%%% getting a gene list
% getting brain specific links genes in book
load('~/networks/tissues/brainTFEINfo3800_thr001.mat');

kado =[ edgeCorrInf(2:end, 1)', edgeCorrInf(2:end, 2)'];
book = sort(unique(kado));


% getting ribosome and histone genes in book 

% ribosome from gene ontology
fid = fopen('~/data/general/ribosomelist.txt')
rlist = textscan(fid, [repmat('%s', 1, 2)])
ribList = rlist{2};

ribList1 = {'RPL3', 'RPL4', 'RPL5', 'RPL6', 'RPL7', 'RPL8', 'RPL9', 'RPL10',  'RPL10A', 'RPL11', 'RPL12', 'RPL13A','RPL14', ...
           'RPL15', 'RPL17','RPL18A', 'RPL19', 'RPL21', 'RPL22', ...
           'RPL23', 'RPL24', 'RPL25', 'RPL26', 'RPL27', 'RPL28', 'RPL29','RPL26L1', 'RPL27', 'RPL30','RL31','RPL32', 'RPL34', ...
           'RPL35', 'RPL35A', 'RPL36', 'RPL36AL', 'RPL37', 'RPL37A', 'RPL38', 'RPL39', 'RPL40', 'RPL41', 'RPLA0', 'RPLA1', 'RPLA2','RPS2', 'RPS3', 'RPS3A', 'RPS4X','RPS4Y1','RPS4Y2','RPS5','RPS6', 'RPS7', 'RPS8', 'RPS9','RPS10','RPS11','RPS12','RPS13', 'RPS14','RPS15','RPS16','RPS17','RPS18','RPS19','RPS20','RPS21','RPS23','RPS24','RPS25','RPS26','RPS27','RPS28','RPS29','RPS30','RPSSA','RPLS6KA3', ...
           'RPS6KB1', 'RPS6KB2', ...
            'RPN1'}

totalRibList = [ribList', ribList1];

histList = {'H1F0', 'H1FNT', 'H1FOO', 'H1FX', 'HIST1H1A', 'HIST1H1B', 'HIST1H1C', ...
            'HIST1H1D', 'HIST1H1E', 'HIST1H1T',  'H2AFB1', 'H2AFB2', 'H2AFB3', ...
            'H2AFJ', 'H2AFV', 'H2AFX', 'H2AFY', 'H2AFY2', 'H2AFZ', 'HIST1H2AA', ...
            'HIST1H2AB', 'HIST1H2AC', 'HIST1H2AD', 'HIST1H2AE', 'HIST1H2AG', ...
           'HIST1H2AI', 'HIST1H2AJ', 'HIST1H2AK', 'HIST1H2AL', 'HIST1H2AM', ...
           'HIST2H2AA3', 'HIST2H2AC', 'H2BFM', 'H2BFS', 'H2BFWT', 'HIST1H2BA', ...
            'HIST1H2BB', 'HIST1H2BC', 'HIST1H2BD', 'HIST1H2BE', 'HIST1H2BF', ...
            'HIST1H2BG', 'HIST1H2BH', 'HIST1H2BI', 'HIST1H2BJ', 'HIST1H2BK', ...
            'HIST1H2BL', 'HIST1H2BM', 'HIST1H2BN', 'HIST1H2BO', 'HIST2H2BE', ...
            'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', ...
           'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST1H4A', 'HIST1H4B', ...
           'HIST1H4C', 'HIST1H4D', 'HIST1H4E', 'HIST1H4F', 'HIST1H4G', 'HIST1H4H', ...
           'HIST1H4I', 'HIST1H4J', 'HIST1H4K', 'HIST1H4L', 'H1F', 'H1H1', 'H2AF', ...
           'H2A1', 'H2A2', 'H2BF', 'H2B1', 'H3A1', 'H3A2', 'H3A3', ...
            'H41', 'H44'}

load('~/data/general/linkExprInfo/wholeExpr.mat')
load('~/data/general/linkExprInfo/probeInd.mat')
load('~/data/general/GPL570GemmaMapNEW.mat')

% getting the gene IDs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a, b] = ismember(totalRibList, gpl570.uniqueSymbols);
sum(a)
ribIDs = b(a);
save('~/data/general/ribosomeIDs.mat', 'ribIDs')
save('~/data/general/ribosomeIDsGO.mat', 'ribIDs')

load('~/data/general/ribosomeIDsGO.mat')
[a, b] = ismember(ribList, gpl570.uniqueSymbols);
sum(a)
ribIDs = b(a);

length(unique(ribIDs))

[a, b] = ismember(histList, gpl570.uniqueSymbols);
sum(a)
histIDs = b(a);
save('~/data/general/histIDs.mat', 'histIDs')
% total number of probe sets for these genes:

book = ribIDs;

book = hkeGenesInd;

book = 1:18494;

sib = probeInd(book, :);
probeCount = sum(sum(sib > 0));
sampleCount = size(wholeExpr, 2);
finalExprMat = zeros(probeCount, sampleCount+2); % geneInd first
                                                 % column, probeInd
                                                 % second column

% for each gene Index in the book, 
% get the probe Indexes, add the expression lines, probeInd and
% gene Ind
probeCounter = 1;
for i = 1:length(book)
    tempProbeInd = probeInd(book(i), :);
    tempProbeInd = tempProbeInd(tempProbeInd > 0);
    geneProbeCount = length(tempProbeInd);
    s = probeCounter;
    e = probeCounter + geneProbeCount - 1;
    finalExprMat(s:e, 1) = book(i);
    finalExprMat(s:e, 2) = tempProbeInd;
    finalExprMat(s:e, 3:end) = wholeExpr(tempProbeInd, :);
    probeCounter = probeCounter + geneProbeCount;
end

fileName = ['~/networks/tissues/' ...
            'brain3800edges_thr001_probeExpLevels.mat'];

fileName = ['~/data/general/wholeExprMat_func.mat']

fileName = '~/networks/tissues/ribosomeExpProbe.mat'

fileName = '~/networks/tissues/hkeGenesExpProbe.mat'

save(fileName, 'finalExprMat')

% gene expression levels 

%% getting the sume networks for a set of genes. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
clear
load('~/data/general/ribosomeIDsGO.mat')

load('~/data/general/histIDs.mat')

geneIDs = ribIDs;
geneIDs = histIDs;

netFolder = '~/networks/tissues/'

load([netFolder 'allSumNetworkExpThr5_01.mat'])

subSumNet = sumNet(geneIDs, geneIDs);

save([netFolder 'GORibosomeSumNet_ExpThr5_01.mat'], 'subSumNet')

% getting tissue networks for ribosome genes. 

netFolder = '~/networks/tissues/'

tissue = 'blood'

thr = '_01'

load([netFolder tissue '/sumNetworkExpThr5' thr '.mat'])

subSumNet = sumNet(geneIDs, geneIDs);

save([netFolder tissue 'GORibosomeSumNet_ExpThr5' thr '.mat'], 'subSumNet')


% sort each cluster based on HKG value, for clusters > 50 edges.

% the whole clustering is over
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the functional terms for these genes and cluster them based
% on their F similarity. I can add edge based on that? or at least
% could check what is happening with funcional similarities! DONE. 

% Draft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('~/data/general/hkgInd.mat')

[a, b] = ismember(ribIDs, hkgInd);

% expression levels of the hkg , the  hke and ribosome genes
% whole geneExpr, wholeprobeExpr, expression bins, different
% networks 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load hkgInd
% load hkeNet
% load ribosomeNet
clear 

load('~/data/general/GPL570GemmaMapNEW.mat')

load('~/data/general/hkgInd.mat')
load('~/data/general/hkeNet.mat')
load('~/data/general/ribosomeIDsGO.mat')

load('~/data/general/linkExprInfo/dataSetProbeInf.mat')
load('~/data/general/linkExprInfo/wholeGeneExpr.mat')

hkgExpr = wholeGeneExpr(hkgInd, :);

length(hkgInd)
sum(sum(hkeNet) > 0)

hkeGenesInd = find((sum(hkeNet) > 0) == 1);

[a, b] = ismember(ribIDs , hkgInd);
sum(a)
[a, b] = ismember(ribIDs , hkeGenesInd);
sum(a)
[a, b] = ismember(hkeGenesInd, hkgInd);
sum(a)

