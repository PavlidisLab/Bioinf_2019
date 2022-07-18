% This file is for getting the gene expression matrix from the cel
% files from GEO, sanity check and heatmap/boxplots of the
% samples. This code is Gemma free (unlike arrayMain02.m and arrayMain01.m)
% THE MUST RUN PART
% Step 0. open the dataset info file and read the names
% Step 1. getting the .CEL files from GEO
% Step 2. making the probe matrix from .cel files
% Step 3. loading the platform file
% Step 4. making the gene matrix
% Step 5. sanity check and correlation heatmap plots 
% Step 6. removing the outliers 
% Step 7. correlation of all samples together for one tissue. 
% Step 8. getting the correlation between list of tissues

% THE MUST RUN PART
%%%%%%%%%%%%%%%%%%%
addpath('/space/grp/marjan/prioritization Project/MATLABcodes/strings')

% folders
celFolder = '/space/grp/marjan/data/skeletalMuscle/celFiles/';
dataFolder = '/space/grp/marjan/data/';
tissue = 'lung';

% Step 0. open the dataset info file and read the names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = [dataFolder tissue '/DataSetInfo_38lung_30-500.txt']
fid = fopen(fileName)

dataSetInf = textscan(fid, ['%s', repmat('%d', 1, 4), '%s', '%d'], ...
                      'headerlines', 1, 'Delimiter', '\t');

% Step 1. getting the .CEL files from GEO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% skipping the following command, I am gettin files from my
% system since downloading from server takes time. 

% for i = 1:length(dataSetInf{1})
%     book = dataSetInf{1}(i); 
%     templ = length(book{1});
%     system(['wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/' book{1}(1:(templ - 3)), ...
%             'nnn/' book{1} '/suppl/' book{1} '_RAW.tar'])
%     i 
%     'Done'
% end

% unzipping the .tar files. 

% getting the list of all files in a folder, for each file : 
celFolder = [dataFolder tissue '/celFiles/']
tarList = dir(celFolder);

for i = 3:length(tarList)
    i
    if(tarList(i).isdir == 0)
        tarList(i).name
        system(['mkdir ' celFolder tarList(i).name(1:(end-4))]);
        system(['tar -xvf ' celFolder tarList(i).name ' -C ' ...
                celFolder tarList(i).name(1:(end-4))], '-echo')
    end

end


% Step 2. making the probe matrix from .cel files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this step is done for each experiment/folder 
celFolder = [dataFolder tissue '/celFiles/'];
celFolder = '~/data/affyArray/tissues/blood/celFiles/'
fileList = dir(celFolder);
libFile = '/space/grp/marjan/data/HG-U133_Plus_2.CDF';

for i = 3:length(fileList)
    i
    tic
    if(fileList(i).isdir == 1) 
        experimentFolder = [celFolder fileList(i).name]
        
        % experimentFolder = [celFolder 'GSE27489_RAW']
        
        %unzip all files here 
        system(['gunzip ' experimentFolder '/*.gz'], '-echo');
        
        %getting list of the cel files here
        tempCelFileList = dir(fullfile(experimentFolder ,'*.CEL'));
        celFileList = cell(1, length(tempCelFileList));
        for j = 1:length(tempCelFileList)
            celFileList{j} = tempCelFileList(j).name;
        end
        
        expr = affyrma(celFileList, libFile, 'CELPath', ...
                       experimentFolder);
        
        tempName = fileList(i).name;
        name = tempName(1:(end - 4))
        dmwrite(expr, [dataFolder tissue '/probeExprMat/' name '.txt' ...
                      ]);
        clear expr
     end
    toc
end

% Step 3. loading the platform file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen('/space/grp/marjan/data/GPL570-13270.txt')
gpl = textscan(fid, [repmat('%s', 1, 16)], 'Delimiter', '\t', ...
               'Headerlines',17, ...
               'Bufsize', 100000095);
%this buffersize just seem to be enough!

%GPL570 probe ID's
gpl570.probeID = gpl{1};

%GPL570 gene ID's
gpl570.geneID = gpl{11};

%GPL570 gene symbols
gpl570.symbols = unique(gpl{11});
gCount = length(gpl570.symbols);
 
geneSymbols = gpl570.symbols;
save(sprintf('%sGPL570geneSymbols.mat', dataFolder), ...
     'geneSymbols');

%loading the geneSymbols
load([dataFolder 'GPL570geneSymbols.mat'])

% Step 4. making the gene matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reading the probe text file
probeFolder = [dataFolder tissue '/probeExprMat/'];
fileList = dir(probeFolder);

%TODO: what is wrong with file 18 20
for i = 3:length(fileList)
    i 
    Name = fileList(i).name
    linum = 1;
    fileName = [probeFolder Name]
    fid = fopen(fileName);
    
    header = textscan(fid, '%s', 1, 'delimiter', '\n', 'headerlines', ...
                      linum - 1, 'Bufsize', 100000095);
    
    %getting the header:
    tempHeader = header{1}{1};
    headers = strsplit(tempHeader, ' ');
    sampleCount = length(headers);
    
    %getting the GOEID from a string. 
    GEOID = cell(1, sampleCount);
    expression = 'GSM[0123456789]+';
    for j = 1:sampleCount
        GEOID(j) = regexp(headers{j}, expression, 'match');
    end
        
    %I need length of header for the sampleCount
    data = textscan(fid, [repmat('%s', 1,1) repmat('%f', 1, ...
                                                   sampleCount)], ...
                    'Headerlines', linum, 'Delimiter', '\t','Bufsize', ...
                    100000095);
    fclose(fid);
    
    pCount = length(data{1});
    dataExpr = cell2mat(data(2:end));
    
    [a, b] = ismember(data{1}, gpl{1});
    geneID = gpl{11}(b(a));
    
    geneMap = containers.Map(gpl570.symbols, [1:gCount]);
    probeMap = containers.Map(data{1}, geneID);
    finalMat = zeros(gCount, sampleCount);
    divMat = zeros(gCount, sampleCount);
    
    tic
    for j = 1:pCount
        sib = values(geneMap, values(probeMap, data{1}(j)));
        finalMat(sib{1}, :) = finalMat(sib{1}, :) + dataExpr(j, :);
        divMat(sib{1}, :) = divMat(sib{1}, :) + 1;
    end
    toc
    
    finalMat = finalMat./divMat;
    
    dataSet.mat = finalMat;
    dataSet.GEOID = GEOID;
    dataSet.Gemma = 'null';
     
    save([dataFolder tissue '/geneExprMat/' Name(1:(end-4)) '.mat'], 'dataSet');
    
end


% Step 5. sanity check and correlation heatmap plots 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%heatmap for each file and saving them
exprFolder = [dataFolder tissue '/geneExprMat/clean/'];
exprFileList = dir(exprFolder);
figFolder = [dataFolder tissue '/figures/'];

%TODO for lung tissue : what is wrong with experiment 19
for i = 3:length(exprFileList)
    i
    Name = exprFileList(i).name
 
    load([exprFolder Name]);
    
    'fileloaded'

    exprData = dataSet.mat;
    sampCount = size(exprData, 2);
    sib = corr(exprData);
    
    % book = cell(1, sampCount);
    % for j = 1: sampCount
    %     book{j} = dataSet.GEOID{j}{1};
    % end
    
    [a,b] = regexp(Name, 'GSE[0123456789]+')
    SName = Name(a:b)
     
    'corr'

    % make and save heatmap
    h = figure; 
    heatmap(sib, dataSet.GEOID,  dataSet.GEOID,  [], 'TickAngle', 45, 'ShowAllTicks', ...
            true, 'Colorbar', true, 'Colormap', 'summer')
    'heatmap done'
    title(sprintf('%s  Number of samples : %d', SName, sampCount))
    fileName = ['cleanGeneSampleCorr' SName];
    print(h, '-depsc', [figFolder fileName '.eps']);
    print(h, '-dpdf', [figFolder fileName '.pdf']);

end


%heatmap for allfiles together and saving them

% Step 6. removing the outliers 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for each experiment: 
%get the result of "Basic algorithm"
%get the result of "Sort-by-median"

% get the intersect and remove them using the code in arrayMain01.m
% for removing outliers. Make and save the heatmap. 

%use the function removeSamples02.m for removing samples and making
%the heatmaps
tissue = 'skeletalMuscle'
dataFolder = '/space/grp/marjan/data/'
exprFolder = [dataFolder tissue '/geneExprMat/']
fileList = dir(exprFolder)

% phase 1 : exploration.
%%%%%%% this loop is for exploring the two outlier detection algorithms
for i = 3: length(fileList)
    i = 3

    Name = fileList(i).name

    load([exprFolder Name]);

    samplesBasic = find(outlierDetectionBasic(dataSet.mat) == 1)
    size(dataSet.mat, 2)
    
    BPviewOnly = 1; %to save or not to save the box plot
    BPfigFolder = [dataFolder tissue '/figures/sampleCorrBP/']
    samplesSBM = find(outlierDetectionSBM(dataSet.mat, BPviewOnly, BPfigFolder,Name) == 1)

    samples = intersect(samplesBasic, samplesSBM)

    %get the intersect of two samples here
    viewOnly = 1;
    titleExtra = 'Basic'
    figFolder = [dataFolder tissue '/figures/basicClean/']
    removeSamples03(exprFolder, Name, figFolder, samplesBasic, titleExtra, ...
                    viewOnly);
    viewOnly = 1;
    titleExtra = 'SBM'
    figFolder = [dataFolder tissue '/figures/SBMClean/']
    removeSamples03(exprFolder, Name, figFolder, samplesSBM, titleExtra, ...
                    viewOnly);

    titleExtra = 'intersect'
    figFolder = [dataFolder tissue '/figures/intersectClean/']
    removeSamples03(exprFolder, Name, figFolder, samples, titleExtra, 0);

end

% phase 2. removing samples. 
%%%%%%%% this loop is for actually removing samples

% after exploring I am just removing samples
for i = 3:length(fileList)
    viewOnly = 0;
    Name = fileList(i).name
    load([exprFolder Name]);
    samplesSBM = find(outlierDetectionSBM(dataSet.mat, BPviewOnly, BPfigFolder,Name) == 1)
    removeSamples02(exprFolder, Name, samplesSBM, viewOnly);

end

    exprData = dataSet.mat;
    sampCount = size(exprData, 2);
    

b = 1
e = sampCount

    sib = corr(exprData(:, b:e));
    
    'Corr'

    % make and save heatmap
    h = figure; 
    heatmap(sib, dataSet.GEOID(b:e), dataSet.GEOID(b:e),  [], 'TickAngle', 45, 'ShowAllTicks', ...
            true, 'Colorbar', true, 'Colormap', 'summer')

% Step 7. correlation of all samples together for one tissue. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get all the samples in the directory
tissue = 'lung'
geneExprFolder = [dataFolder tissue '/geneExprMat/clean/']
fileList = dir(geneExprFolder);

%finding the total number of samples
totalSampleCount = 0;
for i = 3:length(fileList)
    Name = fileList(i).name
    load([geneExprFolder Name]);
    'loaded'
    exprData = dataSet.mat;
    sampCount = size(exprData, 2);
    gCount = size(exprData, 1);
    totalSampleCount = totalSampleCount + sampCount;
end

wholeExpr = zeros(gCount, totalSampleCount);

c = 1;
label = repmat({''}, 1, totalSampleCount);

for i = 3:length(fileList)
    Name = fileList(i).name
    
    load([geneExprFolder Name]);
    'loaded'
    exprData = dataSet.mat;
    
    sampCount = size(exprData, 2);
    wholeExpr(:, c:(c + sampCount -1)) = exprData;
    
    label(c) = {Name};
    c = c + sampCount
    
end

sib = corr(wholeExpr);


h = figure; 
heatmap(sib, label, label,  [], 'TickAngle', 45, 'ShowAllTicks', ...
            true, 'Colorbar', true, 'Colormap', 'hot')
    'heatmap done'
    
temp = length(fileList) - 2
temptitle = sprintf(['%s samples from %d experiments-total of %d ' ...
                        'samples'], tissue, temp, c)
title(temptitle)
    
fileName = [tissue 'AllSamplesCorr']
figFolder = [dataFolder tissue '/figures/']
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);
    
%Step 8. getting the correlation between list of tissues
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 
 %data folder for the blood 
 folder = cell(1, 3)
 folder{1} = '/space/grp/marjan/data/blood/geneExprMat/clean/'
 folder{2} = '/space/grp/marjan/data/lung/geneExprMat/clean/'
 folder{3} = '/space/grp/marjan/data/skeletalMuscle/geneExprMat/clean/'
 
 tissueExpCount = zeros(3, 1);
 %getting the total number of samples
totalSampleCount = 0;
 for i = 1:length(folder)
     i
     fileList = dir(folder{i});
     tissueExpCount(i) = length(fileList) - 2;
     for j = 3:length(fileList)
         j
         Name = fileList(j).name
         load([folder{i} Name]);
         'loaded'
         exprData = dataSet.mat;
         sampCount = size(exprData, 2);
         gCount = size(exprData, 1);
         totalSampleCount = totalSampleCount + sampCount;
     end
 end
 

wholeExpr = zeros(gCount, totalSampleCount);

c = 1;
label = repmat({''}, 1, totalSampleCount);

for i = 1:length(folder)
    fileList = dir(folder{i})
    for j = 3:length(fileList)
        Name = fileList(j).name
        
        load([folder{i} Name]);
        'loaded'
        exprData = dataSet.mat;
        
        sampCount = size(exprData, 2);
        wholeExpr(:, c:(c + sampCount -1)) = exprData;
        
    [a,b] = regexp(Name, 'GSE[0123456789]+');
    label{c} = [Name(a:b)];
    %        label(c) = {Name};
        c = c + sampCount
    end
end

sib = corr(wholeExpr);

h = figure; 
heatmap(sib, label, label,  [], 'TickAngle', 45, 'ShowAllTicks', ...
            true, 'Colorbar', true, 'Colormap', 'hot')
    'heatmap done'
    
    
tissue = 'blood, lung, skeletalMuscle'    
temp = sprintf('%d, %d, %d', tissueExpCount(1), tissueExpCount(2), tissueExpCount(3))
temptitle = sprintf(['%s samples from %s experiments-total of %d ' ...
                        'samples'], tissue, temp, c)
title(temptitle)
    
fileName = ['bloodLungSMuscleAllSamplesCorr']
figFolder = [dataFolder 'general/']
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this is just some draft Code
Name = 'something!'
load([exprFolder Name])
dataSet.mat = dataSet.mat(:, 45:end);
dataSet.GEOID(1:44) = [];

save([exprFolder Name], 'dataSet')
 
% for the next dataset: 
% I am using the file extra to detect the non muscle samples. 

% 1. are the sample orders in my dataSet file as are in the GEO?
% yes. 
Name = 'GSE13070.mat'
%open the file, read line by line, for each line, record the GSM id
%and if it has muscle or not. 
isMuscle = zeros(1, length(dataSet.GEOID));
%walking on the file: 
fid = fopen([dataFolder tissue '/extra.txt'])
myDataSetWalker = 1;

while myDataSetWalker <= length(dataSet.GEOID)
    myDataSetWalker
    tline = fgetl(fid);
    [a, b] = regexp(tline, 'GSM[0123456789]+');
    if(strcmp(dataSet.GEOID{myDataSetWalker}, tline(a:b)))
        if(logical(strfind(tline, 'muscle')))
           isMuscle(myDataSetWalker) = 1;
        end
        myDataSetWalker = myDataSetWalker + 1;
    end
end

toKillIndex = find(isMuscle == 0);
dataSet.GEOID(toKillIndex) = [];
dataSet.mat(:, toKillIndex) = [];
save([exprFolder Name], 'dataSet')


 

 
