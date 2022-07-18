% 1.NEW GTEx data 
%%%%%%%%%%%%%%%
% get the sample file and see what samples do we have. each sample
% has a tissue and a sub tissue. Lets see how many of them are
% blood, brain, muscle, liver and lung

% 2.old GTEx data 
%%%%%%%%%%%%%%%%
% this is how I process the GTEx data from seq file to the gene
% expression count
% 1. opening the data
% 2. opening metadata - meta data has the sample tissue information

%%%%%%%%%%%%%%%%
% 6. OBSOLETE Reproducibility of affy netwroks in GTEx
% 6. NEW Reproducibility of affy netwroks in GTEx
% 6.5 saving the GTEx networks.
% 7. Reproducibility of TS and common genes in GTEx
% 8. Reproducibility of the pure links
% 9. GTEx exp genes
% 10. writing the GTEx data into .txt files (Also the new one with
% blood)
% 11. rep box plot 
% 12. GTEx rep of pure links: the actual links with TSS 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% just reading the gene files from GTEx
file = '~/data/GTEx/genes.txt'
fid = fopen(file)

linum = 1
geneList = textscan(fid, [repmat('%s', 1, 1)], ...
                'Headerlines', linum, 'Delimiter', '\t','Bufsize', ...
                100000095);

length(unique(geneList{1}))

[a, b] = ismember(geneList{1}, gpl570.uniqueSymbols);

% there are some genes for which we have multiple values, from
% different transcripts? I am guessing. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the sample annotation file and get the sample IDs for all
% the tissues. 

clear
dataFolder = '~/data/GTEx/'
file = [dataFolder ...
        'GTEx_Data_V6_Annotations_SampleAttributesDS.txt']

fid = fopen(file)
tline = fgets(fid);
strsplit(tline)

linum = 0;

data = textscan(fid, [repmat('%s', 1, 62)], ...
                'Headerlines', linum, 'Delimiter', '\t','Bufsize', ...
                100000095);

col = 7
tissueString = data{col};
c = 0;
book = zeros(1, length(data{col}));

% sample Index for my five tissues
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for i = 1: length(data{col})
    str = data{col}(i);
    str = str{1};
    if(any(strfind(str, 'Blood')))
        c = c + 1;
        book(i) = 1;
    end
end
c

tissueList = data{col}(logical(book));
unique(tissueList)
length(unique(tissueList))

muscleSampleIDs = data{1}(logical(book));
save([dataFolder 'muscleSampleIDs_V6.mat'], 'muscleSampleIDs');

lungSampleIDs = data{1}(logical(book));
save([dataFolder 'lungSampleIDs_V6.mat'], 'lungSampleIDs');

liverSampleIDs = data{1}(logical(book));
save([dataFolder 'liverSampleIDs_V6.mat'], 'liverSampleIDs');

bloodSampleIDs = data{1}(logical(book));
save([dataFolder 'bloodSampleIDs_V6.mat'], 'bloodSampleIDs');

brainSampleIDs.IDs = data{1}(logical(book));
brainSampleIDs.subTS = data{7}(logical(book));
length(unique(data{1}(logical(book))))
save([dataFolder 'brainSampleIDs_V6.mat'], 'brainSampleIDs');

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% getting sampleInds for all the GTEx tissues 
tissues = unique(data{7});
tissues = tissues(2:end);

for i = 1:length(tissues)
    book = ismember(data{7}, tissues{i});
    sampleIDs = data{1}(book);
    samples(i).tissue = tissues{i};
    samples(i).IDs = sampleIDs;
    samples(i).ind = book;
end

save('~/data/GTEx/sampleTissues.mat', 'samples')

% 2. read the file line by line. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 2.1 getting the mapping of transcripts to genes for genes I have
% in affy - I get the mapping from Ensemble. 
clear
dataFolder = '~/data/general/'
load([dataFolder 'GPL570GemmaMapNEW.mat'])
myGeneSym = gpl570.uniqueIDs;

% loading the ensemble mapping file
dataFolder = '~/data/GTEx/'
file = [dataFolder 'mart_export.txt']
fid = fopen(file)
linum = 1

idData = textscan(fid, [repmat('%s', 1, 3)], ...
                'Headerlines', linum, 'Delimiter', ',','Bufsize', ...
                1000000095);

clear fid 

% getting the transcripts with presenet gene in my affy data
length(idData{1})
[a, b] = ismember(idData{3}, gpl570.uniqueSymbols);
ensMapping.trs = idData{1}(a);
ensMapping.genes = idData{3}(a);
clear idData

geneIndMap = containers.Map(gpl570.uniqueSymbols, 1:18494);
ensGeneIDMap = containers.Map(ensMapping.trs, ensMapping.genes);


%%% 2.2 Now getting the right sample counts for each tissue

% the gene level file

dataFolder = '~/data/GTEx/'

%  >>>>>>>>> the gene level file
file = [dataFolder ['GTEx_Analysis_v6p_RNA-seq_RNA-' ...
                    'SeQCv1.1.8_gene_reads.gct']]
fid = fopen(file)
tline = fgets(fid); % sample IDs are at line 3

kado = strsplit(tline);
length(kado) % = 8559, so I have this number of columns, I want to
             % know which columns are the one for each tissue. 
kado(1:10)

% data = textscan(fid, [repmat('%s', 1, 8558)], ...
%                 'Headerlines', linum, 'Delimiter', ',','Bufsize', ...
%                 1000000095);

% >>>>>>>>>> the transcript level file
file = [dataFolder ['GTEx_Analysis_v6_RNA-' ...
                    'seq_Flux1.6_transcript_rpkm.txt']]
fid = fopen(file)
% read the first line and get the sample counts. 
linum = 1

header = textscan(fid, '%s', 1, 'delimiter', '\n', 'headerlines', ...
                  linum - 1, 'Bufsize', 1000000095);

headerStr = header{1}{1};

kado = strsplit(headerStr);
length(kado) % = 8559, so I have this number of columns, I want to
             % know which columns are the one for each tissue. 
kado(1:10)
%%% a little test of ids, cause some are missing! 
% minKado = cell(1, length(kado));
% for i = 5:length(kado)
%     tempID = kado{i};
%     l = length(tempID);
%     minKado{i} = tempID(6:19);
% end
% length(unique(minKado(5:end)))

% tissue = 'blood'
% load([dataFolder tissue 'SampleIDs_V6.mat']);

% tissueIDs = bloodSampleIDs;
% minTSIDs = cell(1, length(tissueIDs));
% for i = 5:length(tissueIDs)
%     tempID = tissueIDs{i};
%     l = length(tempID);
%     minTSIDs{i} = tempID(6:19);
% end

% minTSIDs = minTSIDs(5:end);
% minKado = minKado(5:end);
% [a, samplesInd]= ismember(minTSIDs, minKado);
% sum(a)

% tempc = 0;
% for i = 1:length(tissueIDs)
%     for j = 5:length(kado)
%        sib = strcmp(tissueIDs{i+4}(1:4), kado{j}(1:4));
%        if sib == 1
%           tempc = tempc + 1; 
%        end
%     end
% end

%%% % end of test part

% % this was old
% tissue = 'brain'
% load([dataFolder tissue 'SampleIDs_V6.mat']);
% tissueIDs = brainSampleIDs;
% it is different for BRAIN, cause it has sub tissues too
%[a, samplesInd]= ismember(tissueIDs.IDs, kado);
%[aa, bb] = ismember(kado(samplesInd(a)), tissueIDs.IDs);

% this is the file I got from the annotation file 
load('~/data/GTEx/sampleTissues.mat')

for t = 1:length(samples)
    [a, samplesInd] = ismember(samples(t).IDs, kado);
    sum(a)
    gCount = 18494;
    exprMat = zeros(gCount, sum(a));
    dataSet(t).mat = exprMat;
    dataSet(t).sampleID = kado(samplesInd(a));
    dataSet(t).samplesInd = samplesInd(a);
    dataSet(t).tissue = samples(t).tissue;
end

save('~/data/GTEx/sampleIndInfo_Gene_RPKM_v6p_New.mat', 'dataSet')

load('~/data/GTEx/sampleIndInfo_Gene_RPKM_v6p_New.mat')
%for t = 1:length(samples)

% [a, samplesInd] = ismember(samples(t).IDs, kado);
% sum(a)
% gCount = 18494;
% exprMat = zeros(gCount, sum(a));

% >>>>>>>>>>>>>>>>>>>>>>>>>> this is for reading the transcript file
% reading the next line
addedTrsCount = 0;
lineReadCount = 0;
tline = fgets(fid);

% I want to keep the record of which lines align to which
% transcripts. 
% also recording the sample IDs
GTExTrans = cell(195747, 1);
geneID = zeros(195747, 1);

while ischar(tline)

    lineReadCount = lineReadCount + 1;
    lineSplit = strsplit(tline);

    tline = fgets(fid);

    transID = lineSplit(2);
    
    %    ensID(lineReadCount) = transID;

    % remove the . part of transID
    t = strfind(transID,'.');
    t = t{1};
    if(~isempty(t))
        transID{1} = transID{1}(1:(t-1));
    end

    GTExTrans{lineReadCount} = transID{1};

    % which gene should be trans added to? 
    try
        geneSymbol= values(ensGeneIDMap, transID);
    catch
        continue;
    end
    geneSymbol
    
    addedTrsCount = addedTrsCount + 1;
    geneInd = values(geneIndMap, geneSymbol)
    geneID(lineReadCount) = geneInd{1};

    for t = 1:length(samples)
        % getting the values and adding them
        exprCellStr = lineSplit(dataSet(t).samplesInd);
         exprValues = cellfun(@str2num, exprCellStr);
        dataSet(t).mat(geneInd{1}, :) = dataSet(t).mat(geneInd{1}, :) + exprValues;
    end
    clear geneInd
end
% dataSet(t).mat = exprMat;
% dataSet(t).sampleID = kado(samplesInd(a));
 %dataSet(t).tissue = samples(t).tissue;
%end

save('~/data/GTEx/allGTExData.mat' 'dataSet')

save([dataFolder 'V6_brainDS.mat'], 'dataSet');

load('~/data/GTEx/allGTExData.mat')
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% this is for reading gene level file
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
clear
dataFolder = '~/data/GTEx/'

file = [dataFolder ['GTEx_Analysis_v6p_RNA-seq_RNA-' ...
                    'SeQCv1.1.8_gene_reads.gct']]
fid = fopen(file)
tline = fgets(fid); % sample IDs are at line 3

kado = strsplit(tline);
length(kado) % = 8559, so I have this number of columns, I want to
             % know which columns are the one for each tissue. 
kado(1:10)
load('~/data/GTEx/sampleIndInfo_Gene_RPKM_v6p_New.mat')
samples = dataSet; 
clear dataSet
for t = 1:length(samples)

    [a, samplesInd] = ismember(samples(t).sampleID, kado);
    sum(a)
    gCount = 56238;
    exprMat = zeros(gCount, sum(a));
    dataSet(t).mat = exprMat;
    dataSet(t).sampleID = kado(samplesInd(a));
    dataSet(t).samplesInd = samplesInd(a);
    dataSet(t).tissue = samples(t).tissue;
end

%  >>>>>>>>> the gene level file
% reading the next line
addedTrsCount = 0;
lineReadCount = 0;

% I want to keep the record of which lines align to which
% transcripts. 
% also recording the sample IDs
GTExTrans = cell(56238, 1);
GTExID = cell(56238, 1);

tline = fgets(fid);
t = 1
while ischar(tline)

    lineReadCount
    lineReadCount = lineReadCount + 1;
    lineSplit = strsplit(tline);

    tline = fgets(fid);

    transID = lineSplit(1);
    
    %    ensID(lineReadCount) = transID;

    % remove the . part of transID
    t = strfind(transID,'.');
    t = t{1};
    if(~isempty(t))
        transID{1} = transID{1}(1:(t-1));
    end

    GTExTrans{lineReadCount} = transID{1};
    GTExGenes{lineReadCount} = lineSplit{2} ;

    % which gene should be trans added to? 
    for t = 1:length(samples)
        exprCellStr = lineSplit(dataSet(t).samplesInd);
        exprValues = cellfun(@str2num, exprCellStr);
        dataSet(t).mat(lineReadCount, :) =  exprValues;
    end
end
% dataSet(t).mat = exprMat;
% dataSet(t).sampleID = kado(samplesInd(a));
% dataSet(t).tissue = samples(t).tissue;
%end
GTExDS.dataSets = dataSet;
GTExDS.Genes = GTExGenes;
GTExDS.Trans = GTExTrans;

kado = zeros(1,54301);
for i = 26982:length(GTExDS.Genes)
    [a, b] = ismember(GTExDS.Genes(i), book);
    kado(b) = kado(b) + 1;
end


save('~/data/GTEx/GTExAllDataSetFromGeneLevel_v6p_newWithBlood.mat', 'GTExDS', ...
     '-v7.3')
load(['~/data/GTEx/' ...
      'GTExAllDataSetFromGeneLevel_v6p_newWithBlood.mat'])

% save this as RPM
meanMilFac = zeros(1, 54);
for g = 1:54
    g
    sib = sum(GTExDS.dataSets(g).mat);
    milFac = sib ./ 1e6;
    meanMilFac(g) = mean(milFac);
    tempMat = GTExDS.dataSets(g).mat;
    for i = 1:size(tempMat, 2);
        tempMat(:, i) = tempMat(:, i)./milFac(i);
    end
    newDataSets(g).mat = tempMat;
    newDataSets(g).sampleID = GTExDS.dataSets(g).sampleID;
    newDataSets(g).samplesInd = GTExDS.dataSets(g).samplesInd;
    newDataSets(g).tissue = GTExDS.dataSets(g).tissue;
end

rpmGTExDS.genes = GTExDS.Genes;
rpmGTExDS.Trans = GTExDS.Trans;
rpmGTExDS.dataSets = newDataSets;

save('~/data/GTEx/GTExAllDataSetFromGeneLevel_v6p_newWithBlood_RPM.mat', 'rpmGTExDS', ...
     '-v7.3')

% save('~/data/GTEx/GTExAllDataSetFromGeneLevel_v6p.mat', 'GTExDS', ...
%      '-v7.3')

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% for each tissue, get the sample numbers in the file 

% read the file line by line , for each tissue, and get the
% transcript for that tissue. 

% add the value of that transcript to that gene's expression
% level. 

dataFolder = '~/data/general/'
load([dataFolder 'GPL570GemmaMapNEW.mat'])
myGeneSym = gpl570.uniqueIDs;

% loading the ensemble mapping file
dataFolder = '~/data/GTEx/'
file = [dataFolder 'mart_export.txt']
fid = fopen(file)
linum = 1

idData = textscan(fid, [repmat('%s', 1, 3)], ...
                'Headerlines', linum, 'Delimiter', ',','Bufsize', ...
                1000000095);

clear fid 

% map the tsIDs to the geneSym IDs in the ensemble file
length(unique(idData{1}))

[C, ia, ic] = unique(idData{1});

smallEnsIDs  = idData{1}(ia);
smallGSym = idData{3}(ia);

% my ENSG ids from the dataset have the ensemble version defined by .[n] at the end
% of the ID. I need to cut that part first. 

repENSIDs = dataSet.gIDs;

tic
for i =1:length(repENSIDs)
    tempID = dataSet.gIDs(i);
    tempID = tempID{1};
    t = strfind(tempID,'.');
    if(~isempty(t))
        repENSIDs{i} = tempID(1:(t-1));
    end
end
toc
% repENSIDs IS the IDS from GTEx data.

% get the data lines for my genes 
[a, b] = ismember(repENSIDs, idData{1});

% 183408 transcripts are in my mapping. how many unique genes do I
% have? 

length(unique(idData{3}))
% 49965 unique gene IDs, I am guessing? but these include long
% non-coding RNAs and stuff like that too! 

% let's see how many of them are from my genes! 

[a, b]= ismember(idData{3}, gpl570.uniqueSymbols);
sum(a)

% 132856 of them are from my genes. thise gives the average of
% 7.1837 transcripts. haha, why so many - I should just sum them up

%% 01 Getting all the IDs aligned.

[a, b] = ismember(repENSIDs, idData{1});

dataSetMat01 = dataSet.dataMat(a, :);
GTEXIDs01 = repENSIDs(a);
tsIDS01 = idData{1}(b(a));
geneIDS01 = idData{3}(b(a));

%% 02 Now let's add them to genes

% which gene IDs are in my platform

[a, b] = ismember(geneIDS01, gpl570.uniqueSymbols);

geneIDS01 = geneIDS01(a);
tsIDS01 = tsIDS01(a);
GTEXIDs01 = GTEXIDs01(a);

min(min(dataSet.dataMat))

% 


% 1. opening the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file = [dataFolder fileList(1).name]

 file = ['/home/ptan/Documents/GTEx/GTEx_Analysis_2014-01-17_RNA-seq_Flux1.6_transcript_rpkm_cleaned.txt']

fid = fopen(file)

%    rows: 194845
%columns : 2920
linum = 1

header = textscan(fid, '%s', 1, 'delimiter', '\n', 'headerlines', ...
                  linum - 1, 'Bufsize', 100000095);

%I need length of header for the sampleCount
data = textscan(fid, [repmat('%s', 1, 2) repmat('%f', 1, ...
                                                2918)], ...
                'Headerlines', linum, 'Delimiter', '\t','Bufsize', ...
                100000095);

c = 0;
while ischar(tline)
tline = fgets(fid);
c = c + 1;
end
c 
% 2. opening metadata - meta data has the sample tissue information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% don't know what is happening here

% file = [dataFolder fileList(2).name]

% fid = fopen(file)
% meta = textscan(fid, [repmat('%s', 1, 15) repmat('%f', 1, ...
%                                                 45)], ...
%                 'Headerlines', linum-1, 'Delimiter', '\t','Bufsize', ...
%                 100000095);

% this is the header of the data file
sib = header{1}{1};

[s, e] = regexp(sib, 'GTEX');

sampleIDs = cell( 2916, 1);
geneID = data{2};

% extracting the sample IDs
for i = 1:2915
    ss = strsplit([sib(s(i): s(i+1)-1)]);
    sampleIDs{i} = ss(1);
end

% reading the ts sample id's
% I get the tissue ids from tissue ID files here, then I get those
% columns in the data file with tissue IDs
%%% !!! FIX THIS C ARRAY
c = [3 5 8 9 12]
tissueList = {'blood', 'brain', 'liver', 'lung', 'muscle'}
for i = 1:length(c)
fileList(c(i)).name
    fid = fopen([dataFolder fileList(c(i)).name])
    d = textscan(fid, ['%s'], ...
                 'Headerlines', 0)%, 'Delimiter', '\t','Bufsize', ...
                                  %                    100000095);

    for j = 1:length(d{1})
        dd(j) = d{1}(j);
    end
    
    IDs(i).tissue = tissueList(i)
    IDs(i).ids = dd;
    clear dd
end

% for each tissue getting the sampleIDs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc(k) = c;
for k = 1: 5
    data = IDs(k).ids;
    c = 0;
    k
    IDs(k).sIndex = zeros(cc(k), 1)
    for i = 1:length(data)
        for j = 1:length(sampleIDs)
            if(strcmp(data(i), sampleIDs{j}))
                c = c + 1;
                IDs(k).sIndex(c) = j;
            end
        end
    end
    %    cc(k) = c;
end

save('~/data/GTEx/sampleInfo.mat', 'IDs')

load('~/data/GTEx/sampleInfo.mat')

c = 0
for i = 1:length(sampleIDs)
    c = c + strcmp('GTEX-N7MS-0011-R10A-SM-2HMJK', sampleIDs{i});
    c = c + strcmp('GTEX-N7MS-0011-R8a-SM-2YUMK', sampleIDs{i});
    c = c + strcmp('GTEX-N7MS-0007-SM-26GME', sampleIDs{i});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
dataFolder = '~/data/GTEx/'
file = ['~/data/GTEx/GTEx_Analysis_2014-01-17_RNA-' ...
        'seq_Flux1.6_transcript_rpkm.txt']

file = ['/home/ptan/Documents/GTEx/GTEx_Analysis_2014-01-17_RNA-seq_Flux1.6_transcript_rpkm_cleaned.txt']

fid = fopen(file)

%    rows: 194845
% columns : 2920

linum = 1

data = textscan(fid, [repmat('%s', 1, 2) repmat('%f', 1, ...
                                                2918)], ...
                'Headerlines', linum-1, 'Delimiter', '\t','Bufsize', ...
                1000000095);

data = cell(1,40);
linum = 1
tot = 0
c = 1
while (tot < 194800)
data{c} = textscan(fid, [repmat('%s', 1, 2) repmat('%f', 1, ...
                                                2918)], ...
                'Headerlines', linum, 'Delimiter', '\t','Bufsize', ...
                1000000095);
tot = tot + length(data{c}{1})
c = c + 1
end

% temp = data{1}(5:end);
% tempIDs = data{1}(1:4);
% dataMat = cell2mat(temp);
%load('~/data/GTEx/sampleInfo.mat')

% finding the data file IDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('~/data/GTEx/sampleInfo.mat')
i = 5
dataSet = struct;
tissue = IDs(i).tissue;
ind = IDs(i).sIndex + 4;
ds = data{ind(2)};
length(ds)
dataMat = cell2mat(ds);
tIDs = data{1};
gIDs = data{2};
dataSet.dataMat = dataMat;
dataSet.tissue = tissue;
dataSet.tIDs = tIDs;
dataSet.gIDs = gIDs;

file = ['~/data/GTEx/' tissue{1} '_dataSet.mat']
save(file, 'dataSet', '-v7.3')
dataSet
clear dataSet

clear
load('~/data/GTEx/brain_dataSet.mat')
dataSet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% writing column 2 to the file
file = [dataFolder 'ensembleMapping.txt']
fid = fopen(file, 'w')

for i = 1:length(data{1})
    temp = data{2}(i);
    sib = strtok(temp{1}, '.');
    fprintf(fid, ['%s\n'], sib);
end

% get my genes from the data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load my geneSym file
dataFolder = '~/data/general/'
load([dataFolder 'GPL570GemmaMapNEW.mat'])
myGeneSym = gpl570.uniqueIDs;

% loading the ensemble mapping file
dataFolder = '~/data/GTEx/'
file = [dataFolder 'mart_export.txt']
fid = fopen(file)
linum = 1

idData = textscan(fid, [repmat('%s', 1, 3)], ...
                'Headerlines', linum, 'Delimiter', ',','Bufsize', ...
                1000000095);

% map the tsIDs to the geneSym IDs in the ensemble file
length(unique(idData{1}))

[C, ia, ic] = unique(idData{1});

smallEnsIDs  = idData{1}(ia);
smallGSym = idData{3}(ia);

% my ENSG ids from the dataset have the ensemble version defined by .[n] at the end
% of the ID. I need to cut that part first. 

repENSIDs = dataSet.gIDs;

tic
for i =1:length(repENSIDs)
    tempID = dataSet.gIDs(i);
    tempID = tempID{1};
    t = strfind(tempID,'.');
    if(~isempty(t))
        repENSIDs{i} = tempID(1:(t-1));
    end
end
toc
% repENSIDs IS the IDS from GTEx data.

% get the data lines for my genes 
[a, b] = ismember(repENSIDs, idData{1});

% 183408 transcripts are in my mapping. how many unique genes do I
% have? 

length(unique(idData{3}))
% 49965 unique gene IDs, I am guessing? but these include long
% non-coding RNAs and stuff like that too! 

% let's see how many of them are from my genes! 

[a, b]= ismember(idData{3}, gpl570.uniqueSymbols);
sum(a)

% 132856 of them are from my genes. thise gives the average of
% 7.1837 transcripts. haha, why so many - I should just sum them up

%% 01 Getting all the IDs aligned.

[a, b] = ismember(repENSIDs, idData{1});

dataSetMat01 = dataSet.dataMat(a, :);
GTEXIDs01 = repENSIDs(a);
tsIDS01 = idData{1}(b(a));
geneIDS01 = idData{3}(b(a));

%% 02 Now let's add them to genes

% which gene IDs are in my platform

[a, b] = ismember(geneIDS01, gpl570.uniqueSymbols);

geneIDS01 = geneIDS01(a);
tsIDS01 = tsIDS01(a);
GTEXIDs01 = GTEXIDs01(a);

min(min(dataSet.dataMat))

% get the data for the tsIDs of my genes 
% blood gene marker 

% the new part with checking 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the data for the tsIDs of my genes 
% blood gene marker 

clear

tissues = {'blood, liver, '}

load('~/data/GTEx/V6_blood.mat')
bloodDS = exprMat;
load('~/data/GTEx/V6_brainDS.mat')
brainDS = dataSet.mat;
addpath('~/codes/MATLAB/myCodes/general')

% sample correlation check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sib = corr(brainDS);

h = figure; 
heatmap(sib, [],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'summer')

meanCorr = mean(sib);
toCut = meanCorr >= .70;
sum(toCut)

newExprMat = brainDS(:, toCut);
size(newExprMat)
sib = corr(newExprMat);

h = figure; 
heatmap(sib, [],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'summer')

% cluster and stuff! 

exprMat = brainDS;
Y = pdist(exprMat');
Z = linkage(Y, 'complete');
T = cluster(Z, 'maxclust', 10);
[a, b]= sort(T);

sortedData = exprMat(:,b);
sib = corr(sortedData);

h = figure; 
heatmap(sib, [],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'summer')

h = figure; 
heatmap(log10(sortedData(1:1000, :)), [],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'summer')


% sanity check for the brain data. The sample correlation is very
% low and the gene corr plot might not look right, just need to
% check what is happening, group them by their tissue type maybe?
% brain dataSets are now grouped. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
load('~/data/GTEx/V6_brainDS.mat')
load('~/data/GTEx/brainSampleIDs_V6.mat')

[a, b] = ismember(dataSet.sampleID, brainSampleIDs.IDs);

% getting subIDs for each sample
subIDs = brainSampleIDs.subTS(b(a));

% grouping subIDs together in the dataMatrix and then saving the
% new dataMatrix AFTER sanity check

uniqueSubIDs = unique(subIDs);
sortedBrainDS = zeros(size(brainDS));
sCounter = 1;
figFolder = '~/data/GTEx/figures/'
for i = 1:length(uniqueSubIDs)
    
    % getting the subtissue data
    [a, b] = ismember(subIDs, uniqueSubIDs{i});
    subIDsDS = brainDS(:, a);
    subSIDs = dataSet.sampleID(a);
    
    % getting the correlation matrix before outlier removal
    sib = corr(subIDsDS);
    h = figure; 
    heatmap(sib, [],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
            true, 'Colorbar', true, 'Colormap', 'summer')
    title([uniqueSubIDs{i} 'Before outlier removal'])
    % save figure
    fileName = ['SampleCorrBOR_' uniqueSubIDs{i}];
    print(h, '-depsc', [figFolder fileName '.eps']);
    print(h, '-dpdf', [figFolder fileName '.pdf']);
    
    % adding it to the big dataset before outlier removal
    sortedBrainDS(:, sCounter:sCounter+sum(a)-1) = subIDsDS;
    sCounter = sCounter + sum(a);
    
    % removing the outliers 
    outliers = find(outlierDetectionSBM(subIDsDS, 1, '', '') == 1)
    subIDsDS(:, outliers) = [];
    subSIDs(outliers) = [];
    
    sib = corr(subIDsDS);
    h = figure; 
    heatmap(sib, [],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
            true, 'Colorbar', true, 'Colormap', 'summer')
    title([uniqueSubIDs{i} ' - Outlier removed'])
    size(subIDsDS)
    % save figure
    fileName = ['SampleCorrAOR_' uniqueSubIDs{i}];
    print(h, '-depsc', [figFolder fileName '.eps']);
    print(h, '-dpdf', [figFolder fileName '.pdf']);
    
    % % gene correlation plot
    sib = corr(subIDsDS');
    upCorr = sib(logical(triu(ones(size(sib)), 1)));
    h = figure
    myHist = hist(upCorr, [-1:.05:1])/length(upCorr);
    bar([-1:.05:1],myHist)
    title([uniqueSubIDs{i} ' - Outlier removed'])
    % save correlation plot
    fileName = ['geneCorrCorrAOR_' uniqueSubIDs{i}];
    print(h, '-depsc', [figFolder fileName '.eps']);
    print(h, '-dpdf', [figFolder fileName '.pdf']);
    
    % save the new dataset and sample IDS (not the tissue)
    brainDataSets(i).mat = subIDsDS;
    brainDataSets(i).IDs = subSIDs;
    brainDataSets(i).subTissue = uniqueSubIDs{i};
end

save('~/data/GTEx/brainDataSet_V6.mat', 'brainDataSets')

h = figure; 
heatmap(sib, [],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'summer')


sib = corr(sortedBrainDS);

h = figure; 
heatmap(sib, [],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'summer')

h = figure
heatmap([1:-.01:-1]', [],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', bone)

print(h, '-dpdf', '~/data/general/tempFigures/boneHeatmap.pdf')


data = [340, % blood
        3,
        7
        225,
        37
        90,
        85,
        360,
        75,
        640, % lung 
        335,
        150,
        135,
        412,
        297,
        120,
        292,
        10,
        25,
        191,
        80,
        73,
        1,
        30, 
        550,% muscle
        
        357,
        480,
       ]
length(data)
h = figure
heatmap(, [],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', bone)

print(h, '-dpdf', '~/data/general/tempFigures/edgeHeatmap.pdf')

% checking the TS distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

load('~/data/GTEx/V6_blood.mat')
bloodDS = exprMat;

load('~/data/GTEx/V6_liver.mat')

load('~/data/GTEx/V6_lung.mat')
load('~/data/GTEx/V6_muscle.mat')

% brain datasets are in a structure based on their subtypes, so
% each sub tissue should be regarded as a different dataset
% indices 5-9 are cerebellum, cortex, frontal cortex, hippocampus, hypothalamus
load('~/data/GTEx/brainDataSet_V6.mat')

h = figure
set(h, 'Position', [100, 100, 900, 480])
bar([393, 320, 430, 119, 121, 114, 106])
colormap pink
grid on
set(gca, 'XTickLabel', {'blood', 'lung', 'muscle', 'liver', ...
                    'cerebellum', 'cortex', 'frontal cortex'})
xlim([.2 4.7])
ylim([0, 500])
grid on
set(h, 'PaperOrientation', 'landscape')
set(h, 'PaperPositionMode', 'auto')
print(h, '-dpdf', '~/data/general/tempFigures/GTEXmyData01.pdf')

% just testing
%  load('~/data/general/tissueExpGenes/brainExpGenes0.8.mat')

%  tissue = 'brain'
%     load(['~/networks/tissues/' tissue '/' ...
%           'binaryNet_0.5DSCountCorr_0.8Expr.mat'])
% % holu = kado > 0 + expGenesInd;
% net = tempNet .* brainTS;
% sum(sum(net > 0))
% testTemplate = zeros(18494, 18494);
% testTemplate(expGenesInd, expGenesInd) = 1;
% kado = brainTS .* testTemplate;
% sum(sum(kado > 0))
clear

load('~/data/GTEx/V6_blood.mat')
bloodDS = exprMat;

% brain datasets are in a structure based on their subtypes, so
% each sub tissue should be regarded as a different dataset
% indices 5-9 are cerebellum, cortex, frontal cortex, hippocampus, hypothalamus
load('~/data/GTEx/brainDataSet_V6.mat')

addpath('~/codes/MATLAB/myCodes/general')
load('~/networks/continuousTSLinksResults/brainSpaceMat.mat')
tissue = 'brain'
load(['~/networks/tissues/' tissue '/' ...
      'binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat'])
brainTS = spaceMat >= .58;
brainTS = brainTS .* spaceMat .*binNet;
sum(sum(brainTS > 0))

load('~/networks/continuousTSLinksResults/bloodSpaceMat.mat')
tissue = 'blood'
load(['~/networks/tissues/' tissue '/' ...
      'binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat'])
bloodTS = spaceMat >= .58;
bloodTS = bloodTS .* spaceMat .* binNet;
sum(sum(bloodTS>0))


% using nother network
load(['~/networks/tissues/brain/' ...
      'binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat'])
load('~/networks/continuousTSLinksResults/brainSpaceMat.mat')
temp = spaceMat >= .58;
brainTS = temp .* binNet;

load(['~/networks/tissues/blood/' ...
      'binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat'])
load('~/networks/continuousTSLinksResults/bloodSpaceMat.mat')
temp = spaceMat >= .58;
bloodTS = temp .* binNet;
sum(sum(bloodTS>0))

[a, b, c] = find(brainTS);
% getting brain corr for brainTS links from GTEx brain
for j = 5:7%length(brainDataSets)
    
    brainDS = brainDataSets(j).mat;
    brainSub = brainDataSets(j).subTissue
    size(brainDS)

    % getting brain gene corr plot
    'brain gene corr'
    sib = corr(brainDS');
    upCorr = sib(logical(triu(ones(size(sib)), 1)));
    brainCorr = upCorr;
    % h = figure
    brainAllHist = hist(upCorr, [-1:.05:1])/sum(~isnan(upCorr));
    % plot([-1:.05:1], brainAllHist, 'color', 'k')

    % I need to get the brain specific corr genes here vs blood
    % specific corrs here (from the same data)
    d = zeros(size(a));
    e = zeros(size(a));
    for i = 1:length(a)
        d(i) = corr(brainDS(a(i), :)', brainDS(b(i), :)');
        e(i) = corr(bloodDS(a(i), :)', bloodDS(b(i), :)');
    end

    brainLinks = [a, b, c, d, e];
    
    %%% trying ksdensity
    hh = figure    
    rlist = datasample(brainCorr, 10000000);
    [f, xi]  = ksdensity(rlist);
    plot(xi, f ,'color', [110, 110, 110]/255 ,'LineWidth', 3)
    hold on
    
    tempInd = ~isnan(d);
    [f, xi]  = ksdensity(d(tempInd));
    plot(xi, f,'color', [15 120 220]/256 ,'LineWidth', 3)
    hold on
    
    tempInd = ~isnan(e);    
    [f, xi]  = ksdensity(e(tempInd));
    plot(xi, f, '--','color', [189 0 38]/256,'LineWidth', 3)
    hold on
    
    grid on
    title('correlation of blood and brain TS in brain')
    file = sprintf('~/data/general/tempFigures/ksDen_GTEXBrainDS%d.pdf', ...
                   j)
    print(hh, '-dpdf', file) 
    
    %%%%%%%%%%%%%%%%%%%%%
    % brainTSbrainHist = hist(d, [-1:.05:1])/sum(~isnan(d));
    % brainTSbloodHist = hist(e, [-1:.05:1])/sum(~isnan(e));
    % sum(~isnan(d))
    %     sum(~isnan(e))

    % plot([-1:.05:1], brainTSbrainHist/sum(brainTSbrainHist))
    % hold on
    % plot([-1:.05:1], brainTSbloodHist/sum(brainTSbloodHist), 'color', 'r')

    % getting blood geneCorr plot
    sib = corr(bloodDS');
    upCorr = sib(logical(triu(ones(size(sib)), 1)));
    bloodCorr = upCorr;
    %h = figure
    %    bloodAllHist = hist(upCorr, [-1:.05:1])/sum(~isnan(upCorr));
    %plot([-1:.05:1], bloodAllHist/sum(bloodAllHist), 'color', 'k')

    [ab, bb, cb] = find(bloodTS);
    db = zeros(size(ab));
    eb = zeros(size(ab));
    for i = 1:length(ab)
        db(i) = corr(bloodDS(ab(i), :)', bloodDS(bb(i), :)');
        eb(i) = corr(brainDS(ab(i), :)', brainDS(bb(i), :)');
    end

    bloodLinks = [ab, bb, cb, db, eb];
    
    % with ksdensity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hh = figure    
    rlist = datasample(bloodCorr, 10000000);
    [f, xi]  = ksdensity(rlist);
    plot(xi, f ,'color', [110, 110, 110]/255 ,'LineWidth', 3)
    hold on
    
    tempInd = ~isnan(eb);
    [f, xi]  = ksdensity(eb(tempInd));
    plot(xi, f,'color', [15 120 220]/256 ,'LineWidth', 3)
    hold on
    
    tempInd = ~isnan(db);    
    [f, xi]  = ksdensity(db(tempInd));
    plot(xi, f, '--','color', [189 0 38]/256,'LineWidth', 3)
    hold on
    
    grid on
    title('correlation of blood and brain TS in blood')
    file = sprintf('~/data/general/tempFigures/ksDen_GTEXBloodDS%d.pdf', ...
                   j)
    print(hh, '-dpdf', file) 
    %%%%%%%%%%%%%%%%%%%%%%%
    
    bloodTSbloodHist = hist(db, [-1:.05:1])/sum(~isnan(db))
    bloodTSbrainHist = hist(eb, [-1:.05:1])/sum(~isnan(eb))
    sum(~isnan(db))
    sum(~isnan(eb))
    %plot([-1:.05:1], bloodTSbloodHist/sum(bloodTSbloodHist))
    %hold on
    %plot([-1:.05:1], bloodTSbrainHist/sum(bloodTSbrainHist))

    h = figure
    hold all
    % plot the brain correlation plot
    plot([-1:.05:1], brainAllHist, 'color', 'k', ...
         'LineWidth', 2)
    sum(brainTSbrainHist)
    % plot the hist of brainTS links from brain dataset
    plot([-1:.05:1], brainTSbrainHist, 'color', ...
         'b', 'LineWidth', 2)
    sum(brainTSbrainHist)
    % plot the hist of bloodTS links from brain dataSet
    plot([-1:.05:1], bloodTSbrainHist, 'color', ...
         'r', 'LineWidth', 2)
    sum(bloodTSbrainHist)
    title(['The correlations from brain data - ' brainSub])
    legend(sprintf('brainAll'), ...
           sprintf('brainTS'), ...
           sprintf('bloodTS'))
    
    fileName = sprintf('%s_HistForCorr_GTEx_%d4QDSCounter_0.8Expr_Ind0.05', 'brain', j)
    figFolder = ['~/data/affyArray/figures/TSplots/']
    print(h, '-depsc', [figFolder fileName '.eps']);
    print(h, '-dpdf', [figFolder fileName '.pdf'])

    
    g = figure
    hold all
    plot([-1:.05:1], bloodAllHist, 'color', 'k', ...
         'LineWidth', 2)
    sum(brainTSbrainHist)
    % plot the hist of brainTS links from brain dataset
    plot([-1:.05:1], bloodTSbloodHist, 'color', ...
         'b', 'LineWidth', 2)
    sum(brainTSbrainHist)
    % plot the hist of bloodTS links from brain dataSet
    plot([-1:.05:1], brainTSbloodHist, 'color', ...
         'r', 'LineWidth', 2)
    sum(bloodTSbrainHist)
    title(['The correlations from blood data - ' brainSub])
    legend(sprintf('bloodAll'), ...
           sprintf('bloodTS'), ...
           sprintf('brainTS'))
    fileName = sprintf('%s_HistForCorr_GTEx_%d_4QDSCounter_0.8Expr_Ind0.05', 'blood', j)
    figFolder = ['~/data/affyArray/figures/TSplots/']
    print(g, '-depsc', [figFolder fileName '.eps']);
    print(g, '-dpdf', [figFolder fileName '.pdf'])
end

% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

load('~/data/GTEx/V6_blood.mat')
bloodDS = exprMat;

% exploring the expression and variance of the genes
% get genes which are expressed! 
sib = bloodDS > 0;
testing = sum(sib');
testing = testing / size(bloodDS, 2);
bloodGenes = testing > .8;
sum(bloodGenes)


% brain 5-9 are relevant datasets 
load('~/data/GTEx/brainDataSet_V6.mat')
brainDS6 = brainDataSets(6).mat;
sib = brainDS6 > 0;
testing = sum(sib');
testing = testing / size(brainDS6, 2);
brainGenes = testing > .8;
sum(brainGenes)

addpath('~/codes/MATLAB/myCodes/general')

tissue = 'brain'
load('~/networks/continuousTSLinksResults/brainSpaceMat.mat')
load(['~/networks/tissues/' tissue '/' ...
                    'binaryNet_0.5DSCountCorr_0.8Expr.mat'])
totalNet = spaceMat .* tempNet;
brainTS = totalNet >= .63;
brainTS = brainTS .* totalNet;
sum(sum(brainTS > 0))

tissue = 'blood'
load('~/networks/continuousTSLinksResults/bloodSpaceMat.mat')
load(['~/networks/tissues/' tissue '/' ...
                    'binaryNet_0.5DSCountCorr_0.8Expr.mat'])
totalNet = spaceMat .* tempNet;
bloodTS = totalNet >= .63;
bloodTS = bloodTS .* totalNet;
sum(sum(bloodTS > 0))

% getting the common links
load('~/networks/tissues/commonNetAllTS_0.8_.5.mat')

% OBSOLETE 6. Reproducibility of affy netwroks in GTEx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Get the ovrelap of every network with each other, it is a
% matrix / save the results in a folder like for the other data
% 2. The rankings, like, AUC or what? I donno. 

clear

tissuesAffy = {'blood', 'brain' 'liver' , 'lung', 'skeletalMuscle'}
tissuesGT = {'brainDataSet', 'blood', 'liver', 'lung', 'muscle'}

netOverlap = zeros(17, 5);
edgeCount = zeros(17, 5);
genesInCommon = zeros(17, 5);
GTexGenesExp = zeros(17, 1);
netDensities = zeros(17, 5);
allDens(85, 1) = struct('xi' , zeros(1, 100), 'fi', zeros(1, 100));
netDens(85, 1) = struct('xi' , zeros(1, 100), 'fi', zeros(1, 100));
GTExCorrBins = zeros(17, 1000);

sampleVector = zeros(17, 2, 3500000);
%allQs = zeros(85, 1000);
%netQs = zeros(85, 1000);
qid = 1;
for tg = 1:17 % 13 brain datasets + 4 other tissues 
    tg
    s = 1
    e = 1e6
    if tg <= 13
        load('~/data/GTEx/brainDataSet_V6.mat')
        dsGT = brainDataSets(tg).mat;
        size(dsGT)
    else
        load(['~/data/GTEx/V6_' tissuesGT{tg - 12} '.mat'])
        dsGT = exprMat;
      size(dsGT)
    end
    
    % find the genes which are expressed in the GTEx: 
    holu = dsGT == 0;
    sib = sum(holu');
    GTExExp = sib < (size(dsGT,2)/4);
    
    % get the null network

    finalDS = dsGT(GTExExp, :);
    sib = corr(finalDS');
    sib = sib - eye(sum(GTExExp));

    '2'

    upperSingle = sib(logical(triu(ones(size(sib)), 1)));
    % getting the quantiles. 
    GTExCorrBins(tg, :) = quantile(upperSingle, 1000);
    
    % getting the null network
    sampleVector(tg, 1, s:e) = datasample(upperSingle, 1e6);
    s = e + 1;
    e = e + 5e5;
    for ta = 1:5
        ta
        %bloodDS = exprMat;
        tissue = tissuesAffy{ta};
        load(['~/data/general/tissueExpGenes/' tissue 'ExpGenes0.8.mat'])

        % find the genes which are expressed in the GTEx: 
        holu = dsGT == 0;
        sib = sum(holu');
        GTExExp = sib < (size(dsGT,2)/4);
        %GTExExp = ones(1, 18494);
        
        GTexGenesExp(tg) = sum(GTExExp);

        final = GTExExp + expGenesInd;
        final = final == 2;
        genesInCommon(tg, ta) = sum(final); % these are the genes expressed both in GTEx and
                   % Affy.

        % getting the affy network
        expThr = '0.8'
        load( ['~/networks/tissues/' tissue '/' ...
               'binaryNet_FDR5e-5_' expThr 'Expr_Ind0.10.mat'])

        myNet = binNet(final, final);
        sum(sum(myNet))
        myNetDensity = sum(sum(myNet)) / (sum(final) * (sum(final) ...
                                                        -1)/2)
        
        netDensities(tg, ta) = myNetDensity;
        edgeCount(tg, ta) = sum(sum(myNet));
        finalDS = dsGT(final, :);
        sib = corr(finalDS');
        sib = sib - eye(sum(final));

        '2'
        %getting the network to see the overlap
        upperSingle = sib(logical(triu(ones(size(sib)), 1)));
        thr = myNetDensity;
        QSingle = quantile(upperSingle, (1 - thr))

        singleNet = sib > QSingle;
        sparseSingleNet = sparse(singleNet);
        
        kado = myNet + sparseSingleNet;
        sum(sum(kado == 2))

        netOverlap(tg, ta) = sum(sum(kado == 2))/sum(sum(myNet));
        
        % getting the quantiles for the distribution
        % allQs(qid, :) = quantile(upperSingle, 1000);
        [allDens(qid).fi allDens(qid).xi] = ksdensity(upperSingle, ...
                                                      'Support', [-1,1])
        
        tempNet = sib .* myNet;
        [a, b, c] = find(tempNet);
        
        sampleVector(tg, 1, s:e) = datasample(c, 5e5);
        sampleVector(tg, 2, s:e) = ta;
        s = e + 1;
        e = e + 5e5;
        
        % netQs(qid, :) = quantile(c, 1000);
        [netDens(qid).fi netDens(qid).xi] = ksdensity(c, 'Support', [-1,1])
        qid = qid + 1;
    end
end

save('~/resultsAndFigures/affyAndGTExComparison/netOverlapForViolin.mat', ...
     'sampleVector', '-v7.3')

load('~/resultsAndFigures/affyAndGTExComparison/netOverlapForViolin.mat')

GTExRep.netOverlap = netOverlap;
GTExRep.edgeCount = edgeCount;
GTExRep.genesInCommon = genesInCommon;
GTExRep.GTexGenesExp = GTexGenesExp;
GTExRep.netDensities = netDensities;
%GTExRep.allQs = allQs;
%GTExRep.netQs = netQs;
GTExRep.allDens = allDens;
GTExRep.netDens = netDens;
GTExRep.corrBins = GTExCorrBins;
GTExRep.tissues = tissuesAffy;
GTExRep.brainTissues = {brainDataSets(:).subTissue};
GTExRep.note = 'GTEx was outer loop'
file = '~/resultsAndFigures/affyAndGTExComparison/res.mat'
save(file, 'GTExRep')

load('~/resultsAndFigures/affyAndGTExComparison/res.mat')

h = figure; 
heatmap(GTExRep.netOverlap([17, 6, 14:16], :), [], [], [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'gray')
print(h, '-depsc', ['~/resultsAndFigures/affyAndGTExComparison/ovrelap.eps']);

h = figure; 
heatmap(genesInCommon, [], [], [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'hot')

subplot(1, 2, 1)
heatmap(allQs(61:65, [1:100:1000]), [], [], [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'hot')
subplot(1, 2, 2)
heatmap(netQs(61:65, [1:100:1000]), [], [], [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'hot')

% NEW 6. Reproducibility of affy netwroks in GTEx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Get the ovrelap of every network with each other, it is a
% matrix / save the results in a folder like for the other data
% 2. The rankings, like, AUC or what? I donno. 

clear
tissuesAffy = {'blood', 'brain' 'liver' , 'lung', 'skeletalMuscle'}

% get the GTEx Tissues 
GTExDSInds = [8:20, 36, 37, 39, 54]; 
files = dir(['~/data/GTEx/matFormat_GPL570_rpmFromGeneLevel_V6/*.mat'])

netOverlap = zeros(17, 5);
edgeCount = zeros(17, 5);
genesInCommon = zeros(17, 5);
GTExGenesExp = zeros(17, 1);
netDensities = zeros(17, 5);
allDens(85, 1) = struct('xi' , zeros(1, 100), 'fi', zeros(1, 100));
netDens(85, 1) = struct('xi' , zeros(1, 100), 'fi', zeros(1, 100));
GTExCorrBins = zeros(17, 1000);

sampleVector = zeros(17, 2, 3500000);
%allQs = zeros(85, 1000);
%netQs = zeros(85, 1000);
tissuesGT = cell(1, 17);
qid = 1;
for tg = 1:17 % 13 brain datasets + 4 other tissues 
    tg
    s = 1
    e = 1e6
    % if tg <= 13
    %     load('~/data/GTEx/brainDataSet_V6.mat')
    %     dsGT = brainDataSets(tg).mat;
    %     size(dsGT)
    % else
    %     load(['~/data/GTEx/V6_' tissuesGT{tg - 12} '.mat'])
    %     dsGT = exprMat;
    %   size(dsGT)
    % end
    
    fid = GTExDSInds(tg)
    load(['~/data/GTEx/matFormat_GPL570_rpmFromGeneLevel_V6/' files(fid).name])
    dsGT = dataSet.mat;
    size(dsGT)
    
    fn = files(fid).name
    tissuesGT{tg} = fn(1: end-4)
    
    % find the genes which are expressed in the GTEx: 
    holu = dsGT == 0;
    sib = sum(holu');
    GTExExp = sib < (size(dsGT,2)/4);
    % get the null network

    finalDS = log2(dsGT(GTExExp, :) + 1);
    sib = corr(finalDS');
    sib = sib - eye(sum(GTExExp));

    '2'

    upperSingle = sib(logical(triu(ones(size(sib)), 1)));
    % getting the quantiles. 
    GTExCorrBins(tg, :) = quantile(upperSingle, 1000);
    
    % getting the null network
    sampleVector(tg, 1, s:e) = datasample(upperSingle, 1e6);
    s = e + 1;
    e = e + 5e5;
    for ta = 1:5
        ta
        %bloodDS = exprMat;
        tissue = tissuesAffy{ta};
        load(['~/data/general/tissueExpGenes/' tissue 'ExpGenes0.8.mat'])

        % find the genes which are expressed in the GTEx: 
        % holu = dsGT == 0;
        % sib = sum(holu');
        % GTExExp = sib < (size(dsGT,2)/4);
        % GTExExp = ones(1, 18494);
        
        GTExGenesExp(tg) = sum(GTExExp);

        final = GTExExp + expGenesInd;
        final = final == 2;
        genesInCommon(tg, ta) = sum(final); % these are the genes expressed both in GTEx and
                % Affy.

        % getting the affy network
        expThr = '0.8'
        load( ['~/networks/tissues/' tissue '/' ...
               'binaryNet_FDR5e-5_' expThr 'Expr_Ind0.10.mat'])

        myNet = binNet(final, final);
        sum(sum(myNet))
        myNetDensity = sum(sum(myNet)) / (sum(final) * (sum(final) ...
                                                        -1)/2)
        
        netDensities(tg, ta) = myNetDensity;
        edgeCount(tg, ta) = sum(sum(myNet));
        finalDS = dsGT(final, :);
        sib = corr(finalDS');
        sib = sib - eye(sum(final));

        '2'
        %getting the network to see the overlap
        upperSingle = sib(logical(triu(ones(size(sib)), 1)));
        thr = myNetDensity;
        QSingle = quantile(upperSingle, (1 - thr))

        singleNet = sib > QSingle;
        sparseSingleNet = sparse(singleNet);
        
        kado = myNet + sparseSingleNet;
        sum(sum(kado == 2))

        netOverlap(tg, ta) = sum(sum(kado == 2))/sum(sum(myNet));
        
        % getting the quantiles for the distribution
        % allQs(qid, :) = quantile(upperSingle, 1000);
        [allDens(qid).fi allDens(qid).xi] = ksdensity(upperSingle, ...
                                                      'Support', [-1,1])
        
        tempNet = sib .* myNet;
        [a, b, c] = find(tempNet);
        
        sampleVector(tg, 1, s:e) = datasample(c, 5e5);
        sampleVector(tg, 2, s:e) = ta;
        s = e + 1;
        e = e + 5e5;
        
        % netQs(qid, :) = quantile(c, 1000);
        [netDens(qid).fi netDens(qid).xi] = ksdensity(c, 'Support', [-1,1])
        qid = qid + 1;
    end
end

save('~/resultsAndFigures/affyAndGTExComparison/netOverlapForViolin_geneLevel_rpm.mat', ...
     'sampleVector', '-v7.3')

GTExRep.netOverlap = netOverlap;
GTExRep.edgeCount = edgeCount;
GTExRep.genesInCommon = genesInCommon;
GTExRep.GTExGenesExp = GTExGenesExp;
GTExRep.netDensities = netDensities;
%GTExRep.allQs = allQs;
%GTExRep.netQs = netQs;
GTExRep.allDens = allDens;
GTExRep.netDens = netDens;
GTExRep.corrBins = GTExCorrBins;
GTExRep.AffyTissues = tissuesAffy;
GTExRep.GTExTissues = tissuesGT;
GTExRep.note = 'GTEx was outer loop'
file = '~/resultsAndFigures/affyAndGTExComparison/res_geneLevel_rpm.mat'
save(file, 'GTExRep')

load('~/resultsAndFigures/affyAndGTExComparison/res.mat')

h = figure; 
heatmap(GTExRep.netOverlap([14, 8 , 15:17], :), [], [], [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'gray')
print(h, '-depsc', ['~/resultsAndFigures/affyAndGTExComparison/ovrelap.eps']);

h = figure; 
heatmap(genesInCommon, [], [], [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'hot')

subplot(1, 2, 1)
heatmap(allQs(61:65, [1:100:1000]), [], [], [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'hot')

subplot(1, 2, 2)
heatmap(netQs(61:65, [1:100:1000]), [], [], [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'hot')


% save the tables. 
% get the distribution and save that. 
% do the same for the Illumina selected datasets. 

% 6.5.0 saving the GTEx networks. daf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
tissuesAffy = {'brain', 'blood', 'liver' , 'lung', 'skeletalMuscle'}
tissuesGT = {'brainDataSet', 'blood', 'liver', 'lung', 'muscle'}

load('~/data/GTEx/brainDataSet_V6.mat')
load(['~/data/general/tissueExpGenes/brainExpGenes0.8.mat'])
% getting the affy network
expThr = '0.8'
load( ['~/networks/tissues/brain/' ...
       'binaryNet_FDR5e-5_' expThr 'Expr_Ind0.10.mat'])

tg = 1;
while tg <= 13
    dsGT = brainDataSets(tg).mat;
    size(dsGT)

    % find the genes which are expressed in the GTEx: 
    holu = dsGT == 0;
    sib = sum(holu');
    GTExExp = sib < (size(dsGT,2)/4);
    GTexGenesExp(tg) = sum(GTExExp);

    final = GTExExp + expGenesInd;
    final = final == 2;
    genesInCommon = sum(final); % these are the genes expressed both in GTEx and
                                % Affy.
    
    myNet = binNet(final, final);
    sum(sum(myNet))
    myNetDensity = sum(sum(myNet)) / (sum(final) * (sum(final) ...
                                                    -1)/2)
    finalDS = dsGT(final, :);
    sib = corr(finalDS');
    sib = sib - eye(sum(final));

    '2'
    %getting the network to see the overlap
    upperSingle = sib(logical(triu(ones(size(sib)), 1)));
    thr = myNetDensity;
    QSingle = quantile(upperSingle, (1 - thr))

    singleNet = sib > QSingle;
    
    singleNet = triu(singleNet, 0);
    
    temp1 = zeros(18494, sum(final));
    temp1(final, :) = singleNet;
    temp2 = zeros(18494, 18494);
    temp2(:, final) = temp1;
    
    sparseSingleNet = sparse(temp2);
    
    fileName = sprintf('GTExBinNet_brain_%s.mat', ...
                       brainDataSets(tg).subTissue);
    GTExNet(tg).tissue = brainDataSets(tg).subTissue;
    GTExNet(tg).net = sparseSingleNet;        
    nd = sum(sparseSingleNet) + sum(sparseSingleNet');
    GTExNet(tg).nd = nd;

    tg = tg + 1;
end

tissuesAffy = {'brain', 'blood', 'liver' , 'lung', 'skeletalMuscle'}
tissuesGT = {'brainDataSet', 'blood', 'liver', 'lung', 'muscle'}

for tg = 2:5
    tissue = tissuesAffy{tg};
    load(['~/data/general/tissueExpGenes/' tissue 'ExpGenes0.8.mat'])
    % getting the affy network
    expThr = '0.8'
    load( ['~/networks/tissues/' tissue '/' ...
           'binaryNet_FDR5e-5_' expThr 'Expr_Ind0.10.mat'])
    
    load(['~/data/GTEx/V6_' tissuesGT{tg} '.mat'])
    dsGT = exprMat;
    size(dsGT)
    
        % find the genes which are expressed in the GTEx: 
    holu = dsGT == 0;
    sib = sum(holu');
    GTExExp = sib < (size(dsGT,2)/4);
    GTexGenesExp(tg) = sum(GTExExp);

    final = GTExExp + expGenesInd;
    final = final == 2;
    genesInCommon = sum(final); % these are the genes expressed both in GTEx and
                                % Affy.

    myNet = binNet(final, final);
    sum(sum(myNet))
    myNetDensity = sum(sum(myNet)) / (sum(final) * (sum(final)...
                                                    -1)/2)                                            
    finalDS = dsGT(final, :);
    sib = corr(finalDS');
    sib = sib - eye(sum(final));
    upperSingle = sib(logical(triu(ones(size(sib)), 1)));
    thr = myNetDensity;
    QSingle = quantile(upperSingle, (1 - thr))

    singleNet = sib > QSingle;
    
    singleNet = triu(singleNet, 0);
    % expanding the net to the full net
    temp1 = zeros(18494, sum(final));
    temp1(final, :) = singleNet;
    temp2 = zeros(18494, 18494);
    finalSingleNet(:, final) = temp1;
    
    sparseSingleNet = sparse(finalSingleNet);
    
    fileName = sprintf('GTExBinNet_%s.mat', tissuesGT{tg});
    save(['~/networks/GTEx/' fileName], 'sparseSingleNet')
    
    nd = sum(sparseSingleNet) + sum(sparseSingleNet');
    GTExNet(tg + 12).nd = nd;
    GTExNet(tg + 12).tissue = tissuesGT{tg};
    GTExNet(tg + 12).net = sparseSingleNet;
end
save('~/networks/GTEx/allNets.mat', 'GTExNet')

load('~/networks/GTEx/allNets.mat')

% 6.5.1 saving the GTEx networks: the new GTEx gene files. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
tissuesAffy = {'brain', 'liver' , 'lung', 'skeletalMuscle', 'blood'}

files = dir(['~/data/GTEx/matFormat_GPL570_rpmFromGeneLevel_V6/*.mat'])

GTExDSInds = [8:20, 36, 37, 39, 54]; % I want all the brains and
                                     % the rests 

load(['~/data/general/tissueExpGenes/brainExpGenes0.8.mat'])
% getting the affy network
expThr = '0.8'
load( ['~/networks/tissues/brain/' ...
       'binaryNet_FDR5e-5_' expThr 'Expr_Ind0.10.mat'])

for tg = 1:17
    
    clear expGenesInd
    clear binNet
    if (tg <= 13)
        load(['~/data/general/tissueExpGenes/brainExpGenes0.8.mat'])
        % getting the affy network
        expThr = '0.8'
        load( ['~/networks/tissues/brain/' ...
               'binaryNet_FDR5e-5_' expThr 'Expr_Ind0.10.mat'])
    else
        ta = tissuesAffy{tg - 13};
        load(['~/data/general/tissueExpGenes/' ta 'ExpGenes0.8.mat'])
        % getting the affy network
        expThr = '0.8'
        load( ['~/networks/tissues/' ta '/' ...
               'binaryNet_FDR5e-5_' expThr 'Expr_Ind0.10.mat'])
    end
    
    fName = files(GTExDSInds(tg)).name
    gtTissue = fName(1: end-4)
    load(['~/data/GTEx/matFormat_GPL570_rpmFromGeneLevel_V6/' fName])
    dsGT = dataSet.mat;
    size(dsGT)

    % find the genes which are expressed in the GTEx: 
    holu = dsGT == 0;
    sib = sum(holu');
    GTExExp = sib < (size(dsGT,2)/4);
    GTexGenesExp(tg) = sum(GTExExp);

    final = GTExExp + expGenesInd;
    final = final == 2;
    genesInCommon = sum(final); % these are the genes expressed both in GTEx and
                                % Affy.
    
    myNet = binNet(final, final);
    sum(sum(myNet))
    myNetDensity = sum(sum(myNet)) / (sum(final) * (sum(final) ...
                                                    -1)/2)
    finalDS = log2(dsGT(final, :)+1);
    sib = corr(finalDS');
    sib = sib - eye(sum(final));

    '2'
    %getting the network to see the overlap
    upperSingle = sib(logical(triu(ones(size(sib)), 1)));
    thr = myNetDensity;
    QSingle = quantile(upperSingle, (1 - thr))

    singleNet = sib > QSingle;
    
    singleNet = triu(singleNet, 0);
    
    temp1 = zeros(18494, sum(final));
    temp1(final, :) = singleNet;
    temp2 = zeros(18494, 18494);
    temp2(:, final) = temp1;
    
    sparseSingleNet = sparse(temp2);

    corrQs = quantile(upperSingle, 1000);
    
    % fileName = sprintf('GTExBinNet_brain_%s.mat', ...
    %                    brainDataSets(tg).subTissue);
    % save(['~/networks/GTEx/' fileName], 'sparseSingleNet')
    GTExNet(tg).tissue = gtTissue;
    GTExNet(tg).net = sparseSingleNet;        
    nd = sum(sparseSingleNet) + sum(sparseSingleNet');
    GTExNet(tg).nd = nd;
    GTExNet(tg).corrQs = corrQs;
end

save('~/networks/GTEx/allNets_rpm_geneLevel.mat', 'GTExNet')
load('~/networks/GTEx/allNets_rpm_geneLevel.mat')
load('~/networks/GTEx/allNets.mat')

% 7. Reproducibility of TS and common links in GTEx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% get the GTEx Tissues  and files ...
GTExDSInds = [8:20, 36, 37, 39, 54]; 
GTExDSInds = [13, 36, 37, 39, 54]; 
files = dir(['~/data/GTEx/matFormat_GPL570_rpmFromGeneLevel_V6/*.mat'])

% load my list of genes for affy
load('~/resultsAndFigures/TSlinks/tissues.mat')

tissuesAffy = tissues;
%load('~/resultsAndFigures/TSlinks/TSlinks.mat') % just load one of
                                                % the TSlinks files
                                                %load('~/resultsAndFigures/TSlinks/TSlinks_FDR0012.mat')
load('~/resultsAndFigures/TSlinks/finalTables/finalTable_CG23_FC3log_FDR0012.mat')

%load('~/resultsAndFigures/TSGenes/finalTSGenes_FDR.1_FC1.mat');

% I have the original dist in the GTEx plots. 
% allDist(85, 1) = struct('xi' , zeros(1, 100), 'fi', zeros(1,
% 100));% original GTEX 

% get the GTEx data
%tissuesGT = {'brainDataSet', 'blood', 'liver', 'lung', 'muscle'}
tissuesGT = cell(17,1)

% 1. >>
% count of genes in GTEx 
% (Again, want to see how many genes I am losing) 
GTExGenesExp = zeros(17, 1);

% 2. >>
% total count of genes in common between the two datasets
% (- and therefore the genes available for the TS 
% want to compare it with the genes shared)
genesInCommon = zeros(17, 5);

% 3. >>
% count of genes in affy TS links
affyTSGenes = zeros(5, 1);

% 4. >>
% count of affy TS links (TS links for this tissue)
affyLinkCounts = zeros(5, 1);

% 5. >>
% count of genes in GTEx TS links (how many genes did we lose for
% expression)
GTExTSGenes = zeros(17, 5);

% 6. >>
% count of TS links in GTEx
GTExLinkCount = zeros(17, 5);

% distribution of the links >>  7 & 8
TSDist(85, 1) = struct('xi' , zeros(1, 100), 'fi', zeros(1, 100));
allDist(17, 1) = struct('xi', zeros(1, 100), 'fi', zeros(1, 100));

% NOTE: >> where I am recording the output
qid = 1;
sampleVector = zeros(17, 2, 5000000);
for tg = 1:17 % 13 brain datasets + 4 other tissues 
    tg
    % if tg <= 13
    %     load('~/data/GTEx/brainDataSet_V6.mat')
    %     dsGT = brainDataSets(tg).mat;
    %     size(dsGT)
    % else
    %     load(['~/data/GTEx/V6_' tissuesGT{tg - 12} '.mat'])
    %     dsGT = exprMat;
    %     size(dsGT)
    % end
    
    fid = GTExDSInds(tg)
    load(['~/data/GTEx/matFormat_GPL570_rpmFromGeneLevel_V6/' files(fid).name])
    dsGT = log2(dataSet.mat + 1);
    size(dsGT)
    
    fn = files(fid).name
    tissuesGT{tg} = fn(1: end-4)

    % find the genes which are expressed in the GTEx: 
    holu = dsGT == 0;
    sib = sum(holu');
    GTExExp = sib < (size(dsGT,2)/4);
    %    GTExExp = ones(1, 18494);
    
    % 1. >> 
    GTExGenesExp(tg) = sum(GTExExp);

    % getting the null GTEx network
    nullDS = dsGT(GTExExp, :);
    sib = corr(nullDS');
    sib = sib - eye(sum(GTExExp));
    upperSingle = sib(logical(triu(ones(size(sib)), 1)));
    
    s = 1;
    e = 1e6;
    sampleVector(tg, 1, 1:e) = datasample(upperSingle, e);
    s = s + e
    
    % 8. >>
    % [allDist(tg).fi allDist(tg).xi] = ksdensity(upperSingle, 'Support', ...
    %                                             [-1,1])
    % h = figure;
    % plot(allDist(tg).xi, allDist(tg).fi, 'color', 'b')
    % hold on
    for ta = 1:5
        % I will just go through the affy links.
        ta
        %bloodDS = exprMat;
        tissue = tissuesAffy{ta}
        load(['~/data/general/tissueExpGenes/' tissue 'ExpGenes0.8.mat'])
        final = GTExExp + expGenesInd;
        final = final == 2;
        
        % 2. >>
        genesInCommon(tg, ta) = sum(final); % these are the genes expressed both in GTEx and
                   % Affy.

        % getting the affy network
        % snInd = TSlinks{ta}.na >= TSlinks{ta}.naFDR;
        % binNet = sparse(TSlinks{ta}.a(snInd), TSlinks{ta}.b(snInd), ...
        %                 1, 18494, 18494);
        
        binNet = finalTable(ta).wholeNet;
        
        % 3. >>
        affyTSGenes(ta) = sum((sum(binNet) + sum(binNet'))> 0);
        
        % 4. >>
        affyLinkCounts(ta) = sum(sum(binNet));
        
        % getting the net with Genes in GTEx: GTEx has most of the
        % genes, but the number of edges drop daramatically with
        % missing edges: those genes are TS genes?
        myNet = binNet(final, final);
        sum(sum(myNet))

        % getting the GTEx network
        finalDS = dsGT(final, :);
        sib = corr(finalDS');
        sib = sib - eye(sum(final));
        
        book = myNet .* sib;
        [a, b, c] = find(book);
        
        % 5. >>
        GTExTSGenes(tg, ta) = length(unique([a, b]));
        
        % 6. >>
        GTExLinkCount(tg, ta) = length(c);
        
        % 7. >>
        % [TSDist(qid).fi, TSDist(qid).xi] = ksdensity(c, 'Support', [-1, ...
        %                     1]);
        
        e = s + length(c) -1
        sampleVector(tg, 1, s:e) = c;
        sampleVector(tg, 2, s:e) = ta;
        s = e + 1
        % hold on
        % plot(TSDist(qid).xi, TSDist(qid).fi, 'color', 'k')

        qid = qid + 1;
    end
end

sib = sampleVector(:, 1,:);
kado = sum(sib);
whos kado
e = min(find(kado ==0))

sampleVector = sampleVector(:, :, 1:e);

save('~/resultsAndFigures/TSlinks/GTEx_TS_OverlapForViolin_rpm_geneLevel_V6.mat', ...
     'sampleVector', '-v7.3')
load(['~/resultsAndFigures/TSlinks/' ...
      'GTEx_TS_OverlapForViolin_rpm_geneLevel_V6.mat'])

% save('~/resultsAndFigures/affyAndGTExComparison/netOverlapForViolin.mat', ...
%      'sampleVector', '-v7.3')

res.tissuesAffy = tissuesAffy;
res.tissuesGT = tissuesGT;
res.GTExBrain = {brainDataSets(:).subTissue};
res.GTExGenesExp = GTExGenesExp;
res.genesInCommon = genesInCommon;
res.affyTSGenes = affyTSGenes;
res.affyLinkCounts = affyLinkCounts;
res.GTExTSGenes = GTExTSGenes;
res.GTExLinkCount = GTExLinkCount;
res.TSDist = TSDist;
res.allDist = allDist;
save('~/resultsAndFigures/TSlinks/TSRepGTEx_FDR0012_rpm_geneLevel_V6.mat', 'res')
load('~/resultsAndFigures/TSlinks/TSRepGTEx_FDR0012.mat')

% 8. reproducibility of pure in GTEx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% load my list of genes for affy
load('~/resultsAndFigures/TSlinks/tissues.mat')

tissuesAffy = tissues;

% get the GTEx Tissues  and files ...
GTExDSInds = [8:20, 36, 37, 39, 54]; 
files = dir(['~/data/GTEx/matFormat_GPL570_rpmFromGeneLevel_V6/*.mat'])

%load('~/resultsAndFigures/TSlinks/TSlinks.mat') % just load one of
                                                % the TSlinks files
                                                %load('~/resultsAndFigures/TSlinks/TSlinks_FDR0012.mat')
load('~/resultsAndFigures/TSlinks/finalTables/finalTable_CG13_FC3log_FDR0012.mat')

%load('~/resultsAndFigures/TSGenes/finalTSGenes_FDR.1_FC1.mat');

% I have the original dist in the GTEx plots. 
% allDist(85, 1) = struct('xi' , zeros(1, 100), 'fi', zeros(1,
% 100));% original GTEX 

% get the GTEx data
tissuesGT = {'brainDataSet', 'blood', 'liver', 'lung', 'muscle'}

% 1. >>
% count of genes in GTEx 
% (Again, want to see how many genes I am losing) 
GTExGenesExp = zeros(17, 1);

% 2. >>
% total count of genes in common between the two datasets
% (- and therefore the genes available for the TS nets, cause I
% want to compare it with the genes shared)
genesInCommon = zeros(17, 5);

% 3. >>
% count of genes in affy TS links
affyTSGenes = zeros(5, 1);

% 4. >>
% count of affy TS links (TS links for this tissue)
affyLinkCounts = zeros(5, 1);

% 5. >>
% count of genes in GTEx TS links (how many genes did we lose for
% expression)
GTExTSGenes = zeros(17, 5);

% 6. >>
% count of TS links in GTEx
GTExLinkCount = zeros(17, 5);

% distribution of the links >>  7 & 8
TSDist(85, 1) = struct('xi' , zeros(1, 100), 'fi', zeros(1, 100));
allDist(17, 1) = struct('xi', zeros(1, 100), 'fi', zeros(1, 100));

% NOTE: >> where I am recording the output
qid = 1;
sampleVector = zeros(17, 2, 5000000);
for tg = 1:17 % 13 brain datasets + 4 other tissues 
    % tg
    % if tg <= 13
    %     load('~/data/GTEx/brainDataSet_V6.mat')
    %     dsGT = brainDataSets(tg).mat;
    %     size(dsGT)
    % else
    %     load(['~/data/GTEx/V6_' tissuesGT{tg - 12} '.mat'])
    %     dsGT = exprMat;
    %     size(dsGT)
    % end

    fid = GTExDSInds(tg)
    load(['~/data/GTEx/matFormat_GPL570_rpmFromGeneLevel_V6/' files(fid).name])
    dsGT = log2(dataSet.mat + 1);
    size(dsGT)
    
    fn = files(fid).name
    tissuesGT{tg} = fn(1: end-4)

    
    % find the genes which are expressed in the GTEx: 
    holu = dsGT == 0;
    sib = sum(holu');
    GTExExp = sib < (size(dsGT,2)/4);
    %    GTExExp = ones(1, 18494);
    
    % 1. >> 
    GTExGenesExp(tg) = sum(GTExExp);

    % getting the null GTEx network
    nullDS = dsGT(GTExExp, :);
    sib = corr(nullDS');
    sib = sib - eye(sum(GTExExp));
    upperSingle = sib(logical(triu(ones(size(sib)), 1)));
    
    s = 1;
    e = 1e6;
    sampleVector(tg, 1, 1:e) = datasample(upperSingle, e);
    s = s + e
    
    % 8. >>
    % [allDist(tg).fi allDist(tg).xi] = ksdensity(upperSingle, 'Support', ...
    %                                             [-1,1])
    % h = figure;
    % plot(allDist(tg).xi, allDist(tg).fi, 'color', 'b')
    % hold on
    for ta = 1:5
        % I will just go through the affy links.
        ta
        %bloodDS = exprMat;
        tissue = tissuesAffy{ta}
        load(['~/data/general/tissueExpGenes/' tissue 'ExpGenes0.8.mat'])
        final = GTExExp + expGenesInd;
        final = final == 2;
        
        % 2. >>
        genesInCommon(tg, ta) = sum(final); % these are the genes expressed both in GTEx and
                   % Affy.

        % getting the affy network
        % snInd = TSlinks{ta}.na >= TSlinks{ta}.naFDR;
        % binNet = sparse(TSlinks{ta}.a(snInd), TSlinks{ta}.b(snInd), ...
        %                 1, 18494, 18494);
        
        binNet = (finalTable(ta).cgNet + finalTable(ta).noTSNet) == 2;
        
        % 3. >>
        affyTSGenes(ta) = sum((sum(binNet) + sum(binNet'))> 0);
        
        % 4. >>
        affyLinkCounts(ta) = sum(sum(binNet));
        
        % getting the net with Genes in GTEx: GTEx has most of the
        % genes, but the number of edges drop daramatically with
        % missing edges: those genes are TS genes?
        myNet = binNet(final, final);
        sum(sum(myNet))

        % getting the GTEx network
        finalDS = dsGT(final, :);
        sib = corr(finalDS');
        sib = sib - eye(sum(final));
        
        book = myNet .* sib;
        [a, b, c] = find(book);
        
        % 5. >>
        GTExTSGenes(tg, ta) = length(unique([a, b]));
        
        % 6. >>
        GTExLinkCount(tg, ta) = length(c);
        
        % 7. >>
        % [TSDist(qid).fi, TSDist(qid).xi] = ksdensity(c, 'Support', [-1, ...
        %                     1]);
        
        e = s + length(c) -1
        sampleVector(tg, 1, s:e) = c;
        sampleVector(tg, 2, s:e) = ta;
        s = e + 1
        % hold on
        % plot(TSDist(qid).xi, TSDist(qid).fi, 'color', 'k')

        qid = qid + 1;
    end
end

sib = sampleVector(:, 1,:);
kado = sum(sib);
whos kado
e = min(find(kado ==0))

sampleVector = sampleVector(:, :, 1:e);

save('~/resultsAndFigures/TSlinks/GTEx_TS_pureCG13RSQR_OverlapForViolin_rpmFromGeneLevel_V6.mat', ...
     'sampleVector', '-v7.3')

save('~/resultsAndFigures/affyAndGTExComparison/netOverlapForViolin.mat', ...
     'sampleVector', '-v7.3')

% 9. get the corr bins of given set of links in GTEx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% for the corrBins
%load('~/resultsAndFigures/affyAndGTExComparison/res.mat')

load('~/resultsAndFigures/affyAndGTExComparison/res_geneLevel_rpm.mat')
% load my list of genes for affy
load('~/resultsAndFigures/TSlinks/tissues.mat')

tissuesAffy = tissues;

% get the GTEx Tissues  and files ...
GTExDSInds = [8:20, 36, 37, 39, 54]; 
files = dir(['~/data/GTEx/matFormat_GPL570_rpmFromGeneLevel_V6/*.mat'])

%load('~/resultsAndFigures/TSlinks/TSlinks.mat') % just load one of
                                                % the TSlinks files
                                                %load('~/resultsAndFigures/TSlinks/TSlinks_FDR0012.mat')

% set of links
load('~/resultsAndFigures/TSlinks/pureLinkMat_V0.mat')

%load('~/resultsAndFigures/TSlinks/finalTables/finalTable_CG23_FC3log_FDR0012.mat')
%load('~/resultsAndFigures/TSGenes/finalTSGenes_FDR.1_FC1.mat');

% I have the original dist in the GTEx plots. 
% allDist(85, 1) = struct('xi' , zeros(1, 100), 'fi', zeros(1,
% 100));% original GTEX 

finalMat = zeros(length(pureLinkMat), (3 + 17)); % 3 is the size of
                                                 % pureLinkMat, 17
                                                 % is for the GTEx
                                                 % tissues
finalMat(:, 1:3) = pureLinkMat;
kado = sparse(pureLinkMat(:, 1), pureLinkMat(:, 2), pureLinkMat(:,3),18494, 18494);
for tg = 1:17 % 13 brain datasets + 4 other tissues 

    fid = GTExDSInds(tg)
    load(['~/data/GTEx/matFormat_GPL570/' files(fid).name])
    dsGT = dataSet.mat;
    size(dsGT)
    
    fn = files(fid).name
    tissuesGT{tg} = fn(1: end-4)

    
    % find the genes which are expressed in the GTEx: 
    holu = dsGT == 0;
    sib = sum(holu');
    GTExExp = sib < (size(dsGT,2)/4);
    
    sib = corr(dsGT');
    
    [ga, gb, gc] = find(sib(logical(kado)));
    nans = isnan(gc);
    tempBins = GTExRep.corrBins(tg, :);
    resBins = zeros(length(ga), 1);
    for i = 1:length(ga)
        if(isnan(gc(i)))
            resBins(i) = -1;
        else
            if (gc(i) > tempBins(1000))
                resBins(i) = 1000;
            else
                [ma, mb] = min(find(tempBins > gc(i)));
                resBins(i) = ma;
            end
        end
    end
    finalMat(:, (3+tg)) = resBins;
end

save('~/resultsAndFigures/TSlinks/pureLinkMat_withGTExBins_V0.mat', ...
     'finalMat')

brainFMat = finalMat(pureLinkMat(:, 3) == 2, :);
for t = 4:20
    h = figure
    hist(brainFMat(:, t), 30)
end

% now get the link with lower other tissues and higher brain: brain
% above 900 in either of the tissues : 1 , 2 , 6 , 7 ,8, 10, 
% below 700 in all the other tissues: 14:17

temp = brainFMat(:, 14:17);
templog = temp <= 750; 
lowOthers = sum(templog');
sum(lowOthers >= 3)

temp = brainFMat(:, [1 ,2 ,6, 7, 8, 10] + 3);
templog = temp >= 900;
highBrain = sum(templog');
sum(highBrain > 0)

final = ((highBrain >0 ) + (lowOthers >= 3)) == 2;

selectedLinks = find(final);

% printing the links for Paul: 

% getting the gene symbol to gene name file 
load('~/data/general/GPL570GemmaMapNEW.mat')
fileName = '~/data/general/geneSymbolToGeneName.txt'

fid = fopen(fileName)
myList = textscan(fid, [repmat('%s', 1, 6)], 'headerlines', 1, ...
                  'Delimiter', '\t');

file = ['~/data/linkListForPaul.txt']
fid = fopen(file, 'w')

fprintf(fid, ['Gene1_Symbol\t', 'Gene1_Name\t' 'Gene2_Symbol\t', 'Gene2_Name\n'])

for j = 1:length(selectedLinks)
    g1 = brainFMat(selectedLinks(j), 1);
    g2 = brainFMat(selectedLinks(j), 2);
    sym1 = gpl570.uniqueSymbols{g1};
    [a, b] = ismember(sym1, myList{1});
    if a
        n1 = myList{2}(b);
    else
        n1 = {' '}
    end
    
    sym2 = gpl570.uniqueSymbols{g2};
    [a, b] = ismember(sym2, myList{1});
    if a
        n2 = myList{2}(b);
    else
        n2 = {' '}
    end
    
    fprintf(fid, ['%s\t%s\t%s\t%s\n'], ...
             sym1, n1{1}, sym2, n2{1})
end
fclose(fid)


% 10. writing the GTEx data into .txt files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
load('~/data/GTEx/GTExAllDataSetFromGeneLevel_v6p.mat')
load(['~/data/GTEx/' ...
      'GTExAllDataSetFromGeneLevel_v6p_newWithBlood_RPM.mat'])
GTExDS = rpmGTExDS;
GTExDS.Genes = GTExDS.genes;

dataFolder = '~/data/general/'
load([dataFolder 'GPL570GemmaMapNEW.mat'])

for i = 1:length(GTExDS.dataSets)
    i
    dataSets = GTExDS.dataSets;
    names = strsplit(dataSets(i).tissue, {' ', '-'});
    nameStr = names{1};
    for c = 2:(length(names)-1)
        nameStr = strcat(nameStr, '_', names{c});
    end
    if(length(names) > 1)
        nameStr = strcat(nameStr,'_' ,names{end});
    end
    nameStr

    % printing all the possible genes in the GTEx
        
    % mapping it to the gpl570 genes to have it at the gene level
    [a, b] = ismember(GTExDS.Genes, gpl570.uniqueSymbols);
    mat = dataSets(i).mat;    
    geneMat = zeros(18494, size(mat, 2));
    inds = find(a);
    for j = 1:length(inds)
        localInd = inds(j);
        geneMat(b(localInd), :) = geneMat(b(localInd), :) + mat(localInd, :);
    end 
    
    % % printing the gene level in .mat file. 
    fileName = ['~/data/GTEx/matFormat_GPL570_rpmFromGeneLevel_V6/' nameStr '.mat']
    dataSet.mat = geneMat;
    dataSet.sampleID = dataSets(i).sampleID;
    dataSet.samplesInd = dataSets(i).samplesInd;
    save(fileName, 'dataSet');
    
    % printing the gene level file
    fileName = [nameStr '.txt']
    file = ['~/data/GTEx/txtFormat_GPL570_rpmFromGeneLevel_V6/' fileName]
    fid = fopen(file, 'w')

    sCount = size(dataSets(i).mat, 2)
    header = dataSets(i).sampleID;
    fprintf(fid, ['GeneSymbol\t', repmat('%s\t', 1, sCount) '\n'], header{:});
    for j = 1:18494
        sym = gpl570.uniqueSymbols{j};
        fprintf(fid, ['%s\t' repmat('%f\t', 1, sCount) '\n'], ...
                sym,  geneMat(j, :));
    end
    fclose(fid)
    
    % printing the transcript level file
    % fileName = [nameStr '.txt']
    % file = ['~/data/GTEx/txtFormat_GTExGeneLevel/' fileName]
    % fid = fopen(file, 'w')

    % mat = dataSets(i).mat;
    % sCount = size(dataSets(i).mat, 2)
    % header = dataSets(i).sampleID;
    % fprintf(fid, ['GeneSymbol\t', 'TranscriptID\t', repmat('%s\t', 1, sCount) '\n'], header{:});
    % for j = 1:size(dataSets(i).mat, 1)
    %     sym = GTExDS.Genes{j};
    %     tID = GTExDS.Trans{j};
    %     fprintf(fid, ['%s\t %s\t' repmat('%f\t', 1, sCount) '\n'], ...
    %             sym, tID, mat(j, :));
    % end
    % fclose(fid)
end

% printing the blood
load([dataFolder 'GPL570GemmaMapNEW.mat'])
fileName = ['blood.txt']
file = ['~/data/GTEx/txtFormat_GTExGeneLevel/' fileName]
fid = fopen(file, 'w')

mat = exprMat;
sCount = size(mat, 2)
header = [1:sCount];
fprintf(fid, ['GeneSymbol\t', 'TranscriptID\t', repmat('%d\t', ...
                                                  1, sCount) '\n'], header);
for j = 1:size(mat, 1)
    sym = gpl5670.uniqueSymbols{j};
    fprintf(fid, ['%s\t' repmat('%f\t', 1, sCount) '\n'], ...
            sym, mat(j, :));
end
fclose(fid)

% GTEx vs Affy corr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
tissuesGT = {'brainDataSet', 'blood', 'liver', 'lung', 'muscle'}
tissueAffy = {'brain', 'blood', 'liver', 'lung', 'skeletalMuscle'};

gtex = zeros(18494, 5);
load('~/data/GTEx/brainDataSet_V6.mat')
ds = brainDataSets(4);
ds.subTissue
gtex(:, 1) = median(ds.mat');

for i = 2:5
    tissue = tissuesGT{i};
    load(['~/data/GTEx/V6_' tissue '.mat']);
    gtex(:, i) = median(exprMat');
end

load('~/data/general/tissueExpGenes/allDSExprInfo.mat')

affy = zeros(18494, 5);
for i = 1:5
    tissue = tissueAffy{i};
    tissueID = ismember({allDSExprInfo(:).tissue}, tissue);
    newID = find(tissueID);
    exprMat = zeros(18494, length(newID));
    for j = 1:length(newID)
     expMat(:, j) = allDSExprInfo(newID(j)).exprLevel;
    end
    affy(:, i) = median(expMat');
end

corrs = zeros(1, 18494);

for i = 1:18494
    corrs(i) = corr(affy(i, :)', gtex(i, :)', 'type', 'Spearman');
end

hist(corrs, 20)

% 9. find the overlap of GTEx and affy expressed genes for each
% tissue. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

tissuesAffy = {'blood', 'brain' 'liver' , 'lung', 'skeletalMuscle'}
tissuesGT = {'brainDataSet', 'blood', 'liver', 'lung', 'muscle'}

netOverlap = zeros(17, 5);
edgeCount = zeros(17, 5);
genesInCommon = zeros(17, 5);
GTexGenesExp = zeros(17, 1);
netDensities = zeros(17, 5);
allDens(85, 1) = struct('xi' , zeros(1, 100), 'fi', zeros(1, 100));
netDens(85, 1) = struct('xi' , zeros(1, 100), 'fi', zeros(1, 100));
GTExCorrBins = zeros(17, 1000);

sampleVector = zeros(17, 2, 3500000);
%allQs = zeros(85, 1000);
%netQs = zeros(85, 1000);
qid = 1;
expMat = zeros(18494, (17+ 10));
for tg = 1:17 % 13 brain datasets + 4 other tissues 
    tg
    s = 1
    e = 1e6
    if tg <= 13
        load('~/data/GTEx/brainDataSet_V6.mat')
        dsGT = brainDataSets(tg).mat;
        size(dsGT)
    else
        load(['~/data/GTEx/V6_' tissuesGT{tg - 12} '.mat'])
        dsGT = exprMat;
      size(dsGT)
    end
    
    % find the genes which are expressed in the GTEx: 
    holu = dsGT == 0;
    sib = sum(holu');
    GTExExp = sib < (size(dsGT,2)/4);
    expMat(:, tg) = GTExExp;
end

% get the null network
for ta = 18:22
    ta
    %bloodDS = exprMat;
    tissue = tissuesAffy{ta-17};
    load(['~/data/general/tissueExpGenes/' tissue 'ExpGenes0.8.mat'])
    expMat(:, ta) = expGenesInd;
end

genesOverlap = expMat' * expMat;
save('~/data/general/tissueExpGenes/GtexAffyOverlap.mat', 'genesOverlap');
save('~/data/general/tissueExpGenes/GtexAffyExpGenes.mat', 'expMat');

%  10. get the significance of network overlap for affy and GTEx 

load('~/resultsAndFigures/affyAndGTExComparison/res.mat')

% just get the significance of overlap 

sig = zeros(17, 5);

for i = 1:17
    
    for j = 1:5
    x = GTExRep.netOverlap(i, j) * GTExRep.edgeCount(i, j);
    M = GTExRep.genesInCommon(i,j) * (GTExRep.genesInCommon(i, j) - ...
                                      1)/2;
    k = GTExRep.edgeCount(i,j);
    sig(i,j) = 1 - hygecdf(x, M, k, k);
    end
end

% the overlap of GTEx and other tissues. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% >>>> load the five ATN and TSN networks

clear

FC = '3'
FDR = '0012'
load(sprintf(['~/resultsAndFigures/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))

tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'};
for t = 1 : 5
    tissue = tissues{t};
    load( ['~/networks/tissues/' tissue '/' ...
           'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
    atn{t} = binNet;
    tsnet{t} = finalTable(t).wholeNet;
end

files = dir('~/data/GTEx/matFormat_GPL570/*.mat')

means = zeros(6, 54);
stds =  zeros(6, 54);
qs = zeros(6, 3, 54);

tsMeans = zeros(6, 54);
tsStds =  zeros(6, 54);
tsQs = zeros(6, 3, 54);
for tg = 24:54
    files(tg).name
    load(['~/data/GTEx/matFormat_GPL570/' files(tg).name]);
    dsGT = dataSet.mat;
    size(dsGT)
    
    % find the genes which are expressed in the GTEx: 
    holu = dsGT == 0;
    sib = sum(holu');
    GTExExp = sib < (size(dsGT,2)/4);
    %    GTExExp = ones(1, 18494);
    
    sum(GTExExp)
    % getting the null GTEx network
    nullDS = dsGT(GTExExp, :);
    sib = corr(nullDS');
    sib = sib - eye(sum(GTExExp));
    upperSingle = sib(logical(triu(ones(size(sib)), 1)));
    
    qs(6, :, tg) = quantile(upperSingle, [.25, .5, .75]);
    means(6, tg) = mean(upperSingle);
    stds(6, tg) = std(upperSingle);
    
    tsQs(6, :, tg) = quantile(upperSingle, [.25, .5, .75]);
    tsMeans(6, tg) = mean(upperSingle);
    tsStds(6, tg) = std(upperSingle);
    
    for ta = 1:5
        anet = atn{ta};
        partanet = anet(GTExExp, GTExExp);
        book = sib(partanet);
        qs(ta, :, tg) = quantile(book, [.25, .5, .75]);
        means(ta, tg) = mean(book);
        stds(ta, tg) = std(book);
    end
        
    for ta = 1:5
        anet = tsnet{ta};
        partanet = anet(GTExExp, GTExExp);
        book = sib(logical(partanet));
        tsQs(ta, :, tg) = quantile(book, [.25, .5, .75]);
        tsMeans(ta, tg) = mean(book);
        tsStds(ta, tg) = std(book);
    end
end

disOverlapInfoAll.qs = qs;
disOverlapInfoAll.means = means;
disOverlapInfoAll.stds = stds;
disOverlapInfoAll.tsQs = tsQs;
disOverlapInfoAll.tsMeans = tsMeans;
disOverlapInfoAll.tsStds = tsStds;
disOverlapInfoAll.files = files;

save(['~/resultsAndFigures/affyAndGTExComparison/' ...
      'coexpDist_allGTExTissues.mat'], 'disOverlapInfoAll')

% plotting lines for the dists: 

whos tsMeans
k = 3
book = means(:, k);

plot(book, 'o')
hold all
plot(book - stds(:, k), 'o')
plot(book + stds(:, k), 'o')

h = figure
k = k
h = figure
book = tsMeans(:, k);
plot(book, 'o')
hold on
plot(book - tsStds(:, k), 'o')
hold on
plot(book + tsStds(:, k), 'o')
k = k + 1
title(sprintf('%d', k))

% saving liklihood ratio heatmap for the GTEx. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the network

load('~/networks/GTEx/allNets.mat')
load('~/resultsAndFigures/affyAndGTExComparison/res.mat')


tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'};
for t = 1 : 5
    tissue = tissues{t};
    load( ['~/networks/tissues/' tissue '/' ...
           'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
    atn{t} = binNet;
end

density = GTExRep.netDensities(1, 4)
gLungNet = GTExNet(15).net;
totEdgeCount = sum(sum(gLungNet))/density;
%aLungNet = atn{4};

overlapCount = sum(sum((aLungNet + gLungNet) ==2));
posRatio = overlapCount/sum(sum(gLungNet))

negEdgeCount = totEdgeCount - 476807 + 62587 

negRatio =  (sum(sum(gLungNet)) - overlapCount) /negEdgeCount

ratio = posRatio / negRatio

% 11. rep box plot 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
load(['~/resultsAndFigures/affyAndGTExComparison/' ...
      'netOverlapForViolin.mat'])
tan = sampleVector;

load('~/resultsAndFigures/TSlinks/GTEx_TS_OverlapForViolin.mat')
tsn = sampleVector;

load('~/resultsAndFigures/TSlinks/GTEx_TS_pure13RSQR_OverlapForViolin.mat')
pure = sampleVector;

t = 14
tissue = 'liver';
ttan = reshape(tan(t, : , :), 2, length(tan));
ttsn = reshape(tsn(t, : , :), 2, length(tsn));
tpure = reshape(pure(t, : , :), 2, length(pure));

selected = (ttsn(2, :) > 0);
sum(selected)
ttsnSelected = ttsn(:, selected);
ttsnSelected(2, :) = ttsnSelected(2, :) + 5;

selected = tpure(2, :) > 0;
sum(selected);
tpureSelected = tpure(:, selected);
tpureSelected(2, :) = tpureSelected(2, :) + 10;

whos ttan
final = [ttan, ttsnSelected, tpureSelected];
h = figure
boxplot(final(1, :), final(2, :))
set(gca, 'yTick', [-.5, 0, .5, 1])
set(gca, 'xTickLabel', [-.5, 0, .5, 1])
figFolder = '~/resultsAndFigures/affyAndGTExComparison/figures/'
file = [tissue '_GTExBoxplot']
print(h, '-depsc', [figFolder file '.eps']);
print(h, '-dpdf', [figFolder file '.pdf']);


% 12. GTEx rep of pure links: the actual links with TSS : for the
% TS links ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
% load my list of genes for affy
load('~/resultsAndFigures/TSlinks/tissues.mat')
tissuesAffy = tissues;

% get the pure nets with their tissue: 
load(['~/resultsAndFigures/TSlinks/finalTables/' ...
      'finalTable_CG13_FC3log_FDR0012.mat'])

sumNet = zeros(18494, 18494);
for t = 1:5
    sumNet = sumNet + ((finalTable(t).cgNet + finalTable(t).noTSNet) ...
             == 2) * t;
end
[a, b, c] = find(sumNet);
pureNets = [a, b, c];


sumNet = zeros(18494, 18494);
for t = 1:5
    sumNet = sumNet + ((finalTable(t).wholeNet)) * t;
end
[a, b, c] = find(sumNet);
pureNets = [a, b, c];

% get the GTEx Tissues  and files: brain cortex, liver, lung, sm, blood
GTExDSInds = [54, 13, 36, 37, 39]; 
files = dir(['~/data/GTEx/matFormat_GPL570_rpmFromGeneLevel_V6/*.mat'])

% get the GTEx data
tissuesGT = {'brainDataSet', 'blood', 'liver', 'lung', 'muscle'}

% NOTE: >> where I am recording the output
results = zeros(length(a), 8);
results(:, 1:3) = pureNets;
GTExExpMat = zeros(18494, 5);
GTExMeanExpMat = zeros(18494, 5);
for tg = 1:5 
    tg
    fid = GTExDSInds(tg)
    load(['~/data/GTEx/matFormat_GPL570_rpmFromGeneLevel_V6/' files(fid).name])
    dsGT = log2(dataSet.mat + 1);
    GTExMeanExpMat(:, tg) = mean(dsGT');
    size(dsGT)

    fn = files(fid).name
    tissuesGT{tg} = fn(1: end-4)
    
    % find the genes which are expressed in the GTEx: 
    holu = dsGT == 0;
    sib = sum(holu');
    GTExExp = sib < (size(dsGT,2)/4);
    
    GTExExpMat(:, tg) = GTExExp;

    % getting the null GTEx network
    nullDS = log2(dsGT(GTExExp, :)+1);
    sib = corr(nullDS');
    sib = sib - eye(sum(GTExExp));
    upperSingle = sib(logical(triu(ones(size(sib)), 1)));

    % >> get the coexpression bins 
    bins = quantile(upperSingle, 1000);
    
   'done'
    temp = zeros(18494, sum(GTExExp));
    temp(GTExExp, :) = sib;
    GTExCorr = ones(18494, 18494) * -2;
    GTExCorr(:, GTExExp) = temp;

    % >> get the coexp rank for all the pure links for this network
    for i = 1:length(pureNets)
        corrV = GTExCorr(results(i, 1), results(i, 2));
        if corrV >  -1.5
            if corrV < bins(1000)
                results(i, tg+3) = min(find(corrV < bins));
            else
                results(i, tg+3) = 1000;
            end
        end
    end
end

[a, b] = sort(results(:, 3));
sResults = results(b, :);
% now get the TSS for the pure links: the average from all the
% other tissues. 
GTExTSS = zeros(length(results), 1);
for i = 1:length(results)
    book = results(i, 4:end);
    if (sum(book>0) == 5)
        main = book(results(i, 3));
        book(results(i, 3)) = [];
        book(book < 500) = 500;
        kado = main - book;
        kado(kado < 0) = 0;
        GTExTSS(i) = (sum(kado)/4);
    end
end

GTExTSS = GTExTSS ./ 500;

% get the average expression level for the links, quantile
% normalized: 
% I do the quantile normalization for those present in all
% tissues:
sib = sum(GTExExpMat') == 5;
whos GTExMeanExpMat
unionExpGenes = (sum(GTExExpMat')> 0);
sib = quantilenorm(GTExMeanExpMat(unionExpGenes, :));

qnGTExMeanExpMat = zeros(18494, 5);
qnGTExMeanExpMat(unionExpGenes, :) = sib;
qnExp = qnGTExMeanExpMat;

% save each of the files for each of the tissues: tables
sum(GTExTSS> .58)
tissuesGT = {'blood', 'brain_cortex', 'liver', 'lung', 'muscle'}
load('~/data/general/GPL570GemmaMapNEW.mat')
for t = 1:5
    f = find(results(:, 3) == t);
    crs = zeros(length(f), 1);
    printRes = [results(f, :), GTExTSS(f)];
    
    file = sprintf(['~/resultsAndFigures/affyAndGTExComparison/' ...
                    'TSlinksTSS/%s_pureLinks_exp_corrs_rpm_geneLevel_v6.csv'], tissuesGT{t})
    fid = fopen(file, 'w')
    fprintf(fid, ['gene_A,', 'gene_B,', 'blood_corr,', 'brain_cortex_corr,', ...
                 'liver_corr,', 'lung_corr,', 'skeletalMuscle_corr,', ...
                  'GTExTSS, gA_bloodExp, gB_bloodExp, gA_brainExp, ' ...
                   'gB_brainExp, gA_liverExp, gB_liverExp, gA_lungExp, ' ...
                   'gB_lungExp, gA_muscleExp, gB_muscleExp\n'])    
    for i = 1:length(printRes)
        g1ID = printRes(i, 1);
        g2ID = printRes(i, 2);
        
        % getting the correlation for min exp and thing. 
        corrs = printRes(i, 4:8);
        mins = zeros(1, 5);
        for j = 1:5
            mins(j) = min(qnExp(g1ID, j), qnExp(g2ID, j));
        end
        sib = corr(mins', corrs');
        g1Sym = gpl570.uniqueSymbols(g1ID);
        g2Sym = gpl570.uniqueSymbols(g2ID);
        fprintf(fid, ['%s,%s, %d, %d, %d, %d, %d, %.3f, %.2f, %.2f, %.2f, ' ...
                      '%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n'], g1Sym{1}, g2Sym{1}, ...
                printRes(i, 4), printRes(i, 5), printRes(i, 6), ...
                printRes(i, 7), printRes(i, 8), printRes(i, 9), ...
                qnExp(g1ID, 1), qnExp(g2ID, 1), qnExp(g1ID, 2), ...
                qnExp(g2ID, 2), qnExp(g1ID, 3), qnExp(g2ID, 3), ...
                qnExp(g1ID, 4), qnExp(g2ID, 4), qnExp(g1ID, 5));
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('~/resultsAndFigures/affyAndGTExComparison/res.mat')

