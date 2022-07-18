% same as arrayMain03, but using files :
% /space/grp/marjan/data/GPL570GemmaMap.txt and GPL570AffyMap.txt
% probes with multiple genes and probes with no genes are deleted
% from probe list in GSEthese two files, also, instead of gene symbol, I
% am using NCBI gene ID in each file. These two mappings share
% 18000 NCBI gene IDs
% This file is for getting the gene expression matrix from the cel
% files from GEO, sanity check and heatmap/boxplots of the
% samples. This code is Gemma free (unlike arrayMain02.m and arrayMain01.m)
% THE MUST RUN PART
% Step 0. open the dataset info file and read the names - obsolete
% Step 1. getting the .CEL files from GEO - obsolete 
% Step 2. making the probe matrix from .cel files 
% Step 3. loading the platform file
% Step 4. making the gene matrix 
% Step 5. sanity check and correlation heatmap plots 
% Step 6. removing the outliers 
% Step 7. correlation of all samples together for one tissue. 
% Step 8. getting the correlation between list of tissues

% THE MUST RUN PART
%%%%%%%%%%%%%%%%%%%
addpath('~/codes/MATLAB/myCodes/strings')
addpath('~/codes/MATLAB/myCodes/general')

% folders
% celFolder = '/space/grp/marjan/data/skeletalMuscle/celFiles/';
% dataFolder = '/space/grp/marjan/data/';
dataFolder = '~/data/array/'
tissue = 'brain';

% Step 0. open the dataset info file and read the names
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fileName = [dataFolder tissue '/DataSetInfo_38lung_30-500.txt']
% fid = fopen(fileame)

% dataSetInf = textscan(fid, ['%s', repmat('%d', 1, 4), '%s', '%d'], ...
%                       'headerlines', 1, 'Delimiter', '\t');

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
tarList = dir([celFolder '*.tar*']);

for i = 4:length(tarList)
    i
    if(tarList(i).isdir == 0)
        tarList(i).name
        %making the new directory for the .CEL files
        system(['mkdir ' celFolder tarList(i).name(1:(end-4))]);
        %unzipping them 
        system(['tar -xvf ' celFolder tarList(i).name ' -C ' ...
                celFolder tarList(i).name(1:(end-4))], '-echo')
        %unzipping them second level 
        system(['gunzip ' celFolder tarList(i).name(1:(end-4)) '/*.gz'])
    end
end


% Step 1.5. reading the batch info data from the files and saving
% it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for each folder, 
% for each .CEL file in that folder , read the particular line, 


% line which has the batch info - maybe get it automatically later
l = 14;
folderList = dir([celFolder 'GSE*RAW']);

for i = 11:length(folderList)
    files = dir([celFolder folderList(i).name '/*.CEL'])
    celStruct = struct;
    for j = 1:length(files)
        j
        celStruct(j).name = files(j).name(1:end-4);
        
        % get the line number
        [book, book1] = system(['awk ''/DatHeader/{print NR - 1}'' ' celFolder ...
                               folderList(i).name '/' celStruct(j).name ...
                            '.CEL']);
        
        l = str2num(book1(1:2));
        l = l + 1;
        l = num2str(l)
        
        
        sib = ['''' l 'q;d' '''' ];
        [book, book2] = system(['sed ' sib ' ' celFolder folderList(i).name '/' files(j).name]);
        % extract the date
        [a, b] = regexp(book2, '\d\d/\d\d/\d\d');
        date = book2(a:b);
        celStruct(j).date = date; 
        % collecting the batchInfo
    end
    % extracting the batches
    batchdate = unique({celStruct(:).date});
    batchMap = containers.Map(batchdate, [1:length(batchdate)]);
    
    fileName = [folderList(i).name(1:end-4) '_sampleInfo.txt']
    file = [dataFolder tissue '/batchInfo/' fileName]
    file = [celFolder 'batchInfo/' fileName]
    
    fid = fopen(file, 'w')
    fprintf(fid, ['sampleID\t date\t batch \n']);
    for j = 1:length(files)
        celStruct(j).batchNumber = values(batchMap, ...
                                          {celStruct(j).date})
        fprintf(fid, ['%s\t %s\t %d\n'], celStruct(j).name, ...
                celStruct(j).date, celStruct(j).batchNumber{1})
    end
    fclose(fid)
    % write the batchInfo file
end

% Step 2. making the probe matrix from .cel files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this step is done for each experiment/folder 
clear

dataFolder = '~/data/array/'
tissue = 'lung'
celFolder = [dataFolder tissue '/celFiles/'];
fileList = dir([celFolder '*_RAW'])
libFile = '/space/grp/marjan/data/HG-U133_Plus_2.cdf';

for i = 1:length(fileList)
    tic
    if(fileList(i).isdir == 1) 
        experimentFolder = [celFolder fileList(i).name]
        
        % experimentFolder = [celFolder 'GSE27489_RAW']
        
        %unzip all files here 
        %   system(['gunzip ' experimentFolder '/*.gz'], '-echo');
        
        %getting list of the cel files here
        tempCelFileList = dir([experimentFolder '/*.CEL']);
        celFileList = cell(1, length(tempCelFileList));
        for j = 1:length(tempCelFileList)
                celFileList{j} = tempCelFileList(j).name;
        end
        
        % % next line doesn't worke because it used *~/* as home in
        % experimentFolder path 
        % expr = affyrma(celFileList, libFile, 'CELPath', ...
        %                    experimentFolder);
        
        
        expr = affyrma(celFileList, libFile, 'CELPath', ...
                           ['/home/mfarahbod/data/array/' tissue ...
                            '/celFiles/' fileList(i).name]);
        
        
        tempName = fileList(i).name;
        name = tempName(1:(end - 4))
        dmwrite(expr, [dataFolder tissue '/textFiles/probeTxtFromCEL/' name '.txt' ...
                      ]);
        clear expr
     end
    toc
end

% Step 3. loading the platform file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen('/space/grp/marjan/data/GPL570GemmaMap.txt')
gpl = textscan(fid, [repmat('%s', 1, 3)], 'Delimiter', '\t', ...
               'Headerlines',1, ...
               'Bufsize', 100000095);
%this buffersize just seem to be enough!

%GPL570 probe ID's
gpl570.probeID = gpl{1};

%GPL570 gene symbolss
gpl570.geneSymbols = gpl{2};

%GPL570 ncbi gene ID
gpl570.ncbiID = gpl{3};

%getting uniqueSymbolsg
gpl570.uniqueSymbols = unique(gpl{2});

%getting the corresponding uniqueIDs
holu = containers.Map(gpl570.geneSymbols, gpl570.ncbiID);
gpl570.uniqueIDs = values(holu, gpl570.uniqueSymbols);

gCount = length(gpl570.uniqueSymbols);
 
geneSymbols = gpl570.uniqueSymbols;
save(sprintf('%sGPL570GemmaMap.mat', dataFolder), ...
     'gpl570');

%loading the geneSymbols
dataFolder = '/space/grp/marjan/data/'
load([dataFolder 'GPL570GemmaMap.mat'])

% Step 4. making the gene matrix - batche effect removed in R using
% batch data ~/codes/R/batchEffectRemoval/BERmove02.R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reading the probe text file
clear
dataFolder = '~/data/general/'
 load([dataFolder 'GPL570GemmaMapNEW.mat'])

dataFolder = '~/data/affyArray/tissues/'
addpath('~/codes/MATLAB/myCodes/strings/')
addpath('~/codes/MATLAB/myCodes/general/')
tissue = 'brain'

% temp addresses
%probeFolder = [dataFolder tissue '/lungTumorSamples/textFiles/probeComBatBatchEffectRemoved/'];
probeFolder = [dataFolder tissue '/'] 

% permanent address
probeFolder = [dataFolder tissue ['/textFiles/' ...
                    'probeComBatBatchEffectRemoved/']];
               
fileList = dir([probeFolder '*.txt'])
gCount = length(gpl570.uniqueSymbols)
thr = zeros(1, length(fileList));

for i = 1:length(fileList)
    Name = fileList(i).name
    linum = 1;
    fileName = [probeFolder Name]
    fid = fopen(fileName)
    
    header = textscan(fid, '%s', 1, 'delimiter', '\n', 'headerlines', ...
                      linum - 1, 'Bufsize', 100000095);
    
    % getting the header:
    tempHeader = header{1}{1};
    headers = strsplit(tempHeader, '\t');
    %        headers = strsplit(tempHeader, '\t');
    sampleCount = length(headers)

    
    % getting the GOEID from a string. 
    GEOID = cell(1, sampleCount);
    expression = 'GSM[0123456789]+';
    for j = 1:sampleCount
         GEOID(j) = regexp(headers{j}, expression, 'match');
    end
    
    % I need length of header for the sampleCount
    data = textscan(fid, [repmat('%s', 1,1) repmat('%f', 1, ...
                                                   sampleCount)], ...
                    'Headerlines', linum-1, 'Delimiter', '\t','Bufsize', ...
                    100000095);
    fclose(fid);
    
    pCount = length(data{1});
    dataExpr = cell2mat(data(2:end));
    %        dataExpr = cell2mat(data(7:end));
    
    %     % %%    removing outliers - yes, I am doing it from the probe matrix
%     BPviewOnly = 1;
%     BPfigFolder = ''
%     %outlier removal here
%     samplesSBM = find(outlierDetectionSBM(dataExpr, BPviewOnly, ...
%                                           BPfigFolder,Name) == 1)
%     dataExpr(:, samplesSBM) = [];
%     GEOID(samplesSBM) = [];
%     size(dataExpr)
%     size(GEOID)
    
%     % ORProbeDataSet = removeSamples04(probeDataSet ,samplesSBM);
%     % ORSampleCount = size(ORProbeDataSet.mat, 2);
    
        %removing the extra ""
    for j = 1:length(data{1})
        if(strcmp(data{1}{j}(1), '"'))
            data{1}{j} = data{1}{j}(2:end);
        end
        if(strcmp(data{1}{j}(end), '"'))
            data{1}{j} = data{1}{j}(1:(end-1));
        end
    end
        
%     % save the outlier removed probeMat into .txt files 
%     saveFolder = ['~/data/array/' tissue ['/textFiles/' ...
%                         'probeBatchRemOutRem/']];
%     fileName = Name;
%     file = [saveFolder Name]
%     fid = fopen(file, 'w')
%     fprintf(fid, ['probeID\t' repmat('%s\t', 1, size(dataExpr, 2)) '\n'], GEOID{:});
%     for j = 1:length(data{1})
%         book = data{1}(j);
%         fprintf(fid, ['%s\t' repmat('%f\t', 1, size(dataExpr, 2)) '\n'], ...
%                 book{1}, dataExpr(j,:));
%     end
%     fclose(fid)
% end
        
    [a, b] = ismember(data{1}, gpl570.probeID);
    geneID = gpl570.ncbiID(b(a)); %gene ID 33k
    geneSymbols = gpl570.geneSymbols(b(a)); %gene NCBI ID 33k
    myProbes = data{1}(a); %probe names - 33k
    myDataExpr = dataExpr(a ,:); %probe values
    pCount = length(myProbes); 
    
    geneMap = containers.Map(gpl570.uniqueSymbols, [1:gCount]);
    probeMap = containers.Map(myProbes, geneSymbols);%  NCBI geneID
    %                                                   map to
    %                                                   probes
        % the new probe gene mapping with filter - these codes are
    % taken from the probeFiltering.m
    % >>>>>>>>>>>>>>>>>>>>>
    
    % making the probeDataSet
    probeDataSet.mat = myDataExpr;
    probeDataSet.GEOID = GEOID;
    %        probeDataSet.GEOID = headers(7:end);
    probeDataSet.Gemma = 'null';
    
    % %%    removing outliers - yes, I am doing it from the probe matrix
    BPviewOnly = 1;
    BPfigFolder = ''
    %outlier removal here
    samplesSBM = find(outlierDetectionSBM(probeDataSet.mat, BPviewOnly, BPfigFolder,Name) == 1)
    ORProbeDataSet = removeSamples04(probeDataSet ,samplesSBM);
    ORSampleCount = size(ORProbeDataSet.mat, 2);
    
    %%%% moved to arrayMain04sub.m
    % %%%% here, get the probe correlation dist - save all the dist
    % %%%% functions for a tissue. 
    % tic
    % holu = corr(probeDataSet.mat');
    % toc
    
    % upHolu = holu(logical(triu(ones(size(holu)), 1)));
    % X = randi(length(upHolu), [1, floor(length(upHolu)/10)]);
    % tic
    % [f, xi] = ksdensity(upHolu(X));
    % toc


    % % getting the index of probes for each gene - this is not based
    % % on data, but on the platform
    divMat = zeros(gCount, ORSampleCount);

    %    adding the maps
    tic
    for j = 1:pCount
        sib = values(geneMap, values(probeMap, myProbes(j)));
        divMat(sib{1}, :) = divMat(sib{1}, :) + 1;
    end
    toc

    %getting the porbe list for each gene
    %max probe count for each gene
    probeInd = zeros(gCount, max(max(divMat))); 
    for j = 1:pCount
        sib = values(geneMap, values(probeMap, myProbes(j)));
        temp = find(probeInd(sib{1}, :) == 0);
        pInd = temp(1);
        probeInd(sib{1}, pInd) = j;
    end

    % % probe filtering, this would be the next step - choosing the
    % filter thr with the sex related genes
    %%%%%%%%%%%%%%%%%%%%%%% >>>>>>>>>>>
    % sryID = values(geneMap, {'SRY'}) % this gene has only one
    %                                      % probe
    % sryProbeList = probeInd(sryID{1}, :)
    % % do I have XIST-femaleMarker and KDM5D-maleMarker ? yes
    % kdm5dID = values(geneMap, {'KDM5D'}) % this gene has only one
    %                                      % probe
    % kdm5dProbeList = probeInd(kdm5dID{1}, :)
    % xistID = values(geneMap, {'XIST'}) % has many probes
    % xistProbeList = probeInd(xistID{1}, :)
    % %       filtering 
    % thr = 6; % threshold is gotten by comparing two sex specific
    %          % genes (probeFiltering.m - lines 150 - 474)
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    probeMean = mean(ORProbeDataSet.mat, 2);
    a = mean(probeMean)
    newGeneExpr = zeros(gCount, ORSampleCount);
    affectedGeneCount = 0;
    
    % for each dataset, set thr as the 2/3 of the expression level
    
    % get the 2/3 quantile, 
    q = quantile(ORProbeDataSet.mat(:), 2)
    thr(i) = q(1);
    
    % put the thr for that 
    
    % for each gene
    for j = 1: length(probeInd)
        pCount = sum(probeInd(j,:)>0);
        %if there are any with exp > thr, get the mean of bigger ones
        largeProbes = find(probeMean(probeInd(j,1:pCount))> thr(i));

        if(length(largeProbes) > 0)
            largeProbesInd = probeInd(j, largeProbes);
            newGeneExpr(j,:) = mean(ORProbeDataSet.mat(largeProbesInd, ...
                                                       :), 1);
    
            if(pCount > length(largeProbes))
                affectedGeneCount = affectedGeneCount + 1;
            end
        else
            newGeneExpr(j, :) = mean(ORProbeDataSet.mat(probeInd(j, ...
                                                              1:pCount), ...
                                                        :), 1);
        end
    end
    
     finalMat = newGeneExpr;
    
    % the following old lines were commendted, I did not filter the
    % probes in them
    %%%%%%%%%%% >>>>>>>>>>>>>>>>>>>>>>>
    % finalMat = zeros(gCount, sampleCount);
    % divMat = zeros(gCount, sampleCount);
    
    % tic
    % for j = 1:pCount
    %     sib = values(geneMap, values(probeMap, myProbes(j)));
    %     finalMat(sib{1}, :) = finalMat(sib{1}, :) + myDataExpr(j, :);
    %     divMat(sib{1}, :) = divMat(sib{1}, :) + 1;
    % end
    % toc
    
    %finalMat = finalMat./divMat;
    
    
    % probeInd = zeros(gCount, max(max(divMat)));
    % for j = 1:pCount
    %     sib = values(geneMap, values(probeMap, myProbes(j)));
    %     temp = find(probeInd(sib{1}, :) == 0);
    %     pInd = temp(1);
    %     probeInd(sib{1}, pInd) = j;
    % end
    
    % book = find(divMat(:,1) == 3);
    % twoProbeGenesPIDs = probeInd(book, 1:3);
    % twoProbeGenesPexp = zeros(length(book), 3);
    % for j = 1: length(book)
    %     for k = 1:3
    %         twoProbeGenesPexp(j,k) = mean(myDataExpr(twoProbesGenesPIDs(j,k), :));
    %     end
    % end
    
    
    % plot(twoProbeGenesPexp(:,1), twoProbeGenesPexp(:,3), '.')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% <<<<<<<<
    % removing the new outliers with the z-score
    % tic
    % book = clustergram(finalMat, 'Standardize', 'Row');
    % toc
    
    % % this should be fixed, cause that is not always the case
    % toRemove = clusterGroup(book, (ORSampleCount-2), 1, 'InfoOnly', ...
    %                         'true');
    
    % toRemove = toRemove.ColumnNodeNames;
    % removeInd = zeros(1, length(toRemove));
    % for(j = 1:length(toRemove))
    %     removeInd(j) = str2num(toRemove{j});
    % end
    
    
    % finalMat(:, removeInd) = [];
    % ORProbeDataSet.GEOID(removeInd) = [];
    
    dataSet.mat = finalMat;
    dataSet.GEOID = ORProbeDataSet.GEOID;
    dataSet.Gemma = 'null';
    
    %temp address
        % save([dataFolder tissue '/lungTumorSamples/matFiles/geneExprMatGemmaMapBER/' Name(1:(end-4)) '.mat'], 'dataSet');
    
        %permanent address
            save([dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/' Name(1:(end-4)) '.mat'], 'dataSet');
    
end
 

% Step 5. sanity check and correlation heatmap plots 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%heatmap for each file and saving them
clear
addpath('~/codes/MATLAB/myCodes/general/')
dataFolder = '~/data/array/'
tissue = 'blood';
exprFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/']
exprFileList = dir([exprFolder '*.mat'])
figFolder = [dataFolder tissue '/figures/'];

%TODO for lung tissue : what is wrong with experiment 19
for i = 1:length(exprFileList)

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
    
    % holu = find(sib(:,1) <0.93);
    % t = ones(1, size(dataSet.mat, 2));
    % t(holu) = 0;
    % dataSet.mat = dataSet.mat(:, logical(t));
    % dataSet.GEOID = dataSet.GEOID(logical(t));
    % book  = dataSet.GEOID(logical(t));
    % sib = corr(book);
    
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
% get the intersect and remove them using the code in arrayMain01.m
% for removing outliers. Make and save the heatmap. 

%use the function removeSamples02.m for removing samples and making
%the heatmap
tissue = 'liver'
dataFolder = '/space/grp/marjan/data/'
exprFolder = [dataFolder tissue '/geneExprMat/gemmaProbeMap/BER/']
fileList = dir(exprFolder)

tempfile = [exprFolder Name]
save(tempfile, 'dataSet');

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
clear
addpath('~/codes/MATLAB/myCodes/general/')
dataFolder = '/space/grp/marjan/data/'
newDataFolder = '~/data/array/'
tissue = 'skeletalMuscle'
BPviewOnly = 0;

exprFolder = [newDataFolder tissue '/matFiles/geneExprMatGemmaMapBER/']
BPfigFolder = [newDataFolder tissue '/figures/outlierDetectionBoxplots/']
fileList = dir([exprFolder '*.mat'])
% after exploring I am just removing samples
for i = 1:length(fileList)
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
            true, 'Colorbar', true, 'Colormap', 'map')

% Step 7. correlation of all samples together for one tissue. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get all the samples in the directory

%dataFolder = '/space/grp/marjan/data/'
clear
dataFolder = '~/data/affyArray/tissues/'
tissue = 'blood'

geneExprFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/']

%geneExprFolder = [dataFolder tissue ['/geneExprMat/gemmaProbeMap/' ...
%                    'clean/']
                  
fileList = dir([geneExprFolder '*.mat'])

%finding the total number of samples
totalSampleCount = 0;
for i = 1:length(fileList)
    Name = fileList(i).name
    load([geneExprFolder Name]);
    'loaded'
    exprData = dataSet.mat;
    
    sib = mean(exprData');
    quantile(sib, [.5, .6, .7, .8, .9])
end

    sampCount = size(exprData, 2);
    gCount = size(exprData, 1);
    totalSampleCount = totalSampleCount + sampCount;
end

wholeExpr = zeros(gCount, totalSampleCount);

c = 1;
label = repmat({''}, 1, totalSampleCount);

%NOTE : I am doing the inverse log transfer in the following code,
%just checking how the plot looks - I like it better with log
%transferred 
for i = 1:length(fileList)
    Name = fileList(i).name
    
    load([geneExprFolder Name]);
    'loaded'
    exprData = 2.^(dataSet.mat);
    
    sampCount = size(exprData, 2);
    wholeExpr(:, c:(c + sampCount -1)) = exprData;
    
    
    [a,b] = regexp(Name, 'GSE[0123456789]+');
    label{c} = [Name(a:b)];
    %   
    
    %    label(c) = {Name};
    c = c + sampCount
    
end

sib = corr(wholeExpr, 'type', 'Pearson');
cmap = makingRedBlueColorMap();

h = figure;
colormap('cmap');
heatmap(sib)

h = figure; 
heatmap(sib, label, label,  [], 'TickAngle', 45, 'ShowAllTicks', ...
            true, 'Colorbar', true, 'Colormap', 'bone')
    'heatmap done'
    
temp = length(fileList)
temptitle = sprintf(['%s samples from %d experiments-total of %d ' ...
                        'samples - after removing batch effect from experiments'], tissue, temp, c)
title(temptitle)
    
fileName = [tissue 'AllSamplesCorrABERPearsonBone']
figFolder = ['~/data/array/' tissue '/figures/sampleCorrHeatmapsGenelevelFinal/']
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf'])

%8. getting the correlation between list of tissues
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
clear    
%Step(h, '-dpdf', [figFolder fileName '.pdf'])clear 
 %data folder for the blood 
 folder = cell(1, 4)
 folder{1} = '/space/grp/marjan/data/blood/geneExprMat/gemmaProbeMap/clean/'
 folder{2} = '/space/grp/marjan/data/lung/geneExprMat/gemmaProbeMap/clean/'
 folder{3} = '/space/grp/marjan/data/skeletalMuscle/geneExprMat/gemmaProbeMap/clean/'
 folder{4} = '/space/grp/marjan/data/liver/geneExprMat/gemmaProbeMap/clean/'

 % data with batch effect removed!
 folder{1} = '~/data/affyArray/tissues/blood/matFiles/geneExprMatGemmaMapBER/'
 folder{2} = '~/data/affyArray/tissues/lung/matFiles/geneExprMatGemmaMapBER/'
 folder{3} = '~/data/affyArray/tissues/skeletalMuscle/matFiles/geneExprMatGemmaMapBER/'
 folder{4} = '~/data/affyArray/tissues/liver/matFiles/geneExprMatGemmaMapBER/'
 folder{5} = '~/data/affyArray/tissues/brain/matFiles/geneExprMatGemmaMapBER/'

 
 ts{1} = 'blood';
 ts{2} = 'lung';
 ts{3} = 'SM';
 ts{4} = 'liver';
 ts{5} = 'brain';
 
 tissueExpCount = zeros(5, 1);
 tissueSampleCount = zeros(5, 1);
 tissueSample = zeros(5, 15);
 %getting the total number of samples
totalSampleCount = 0;
 for i = 1:length(folder)
     i
     'tissue'
     fileList = dir([folder{i} '*.mat']);
     tissueExpCount(i) = length(fileList);
     for j = 1:length(fileList)
         j
         Name = fileList(j).name
         load([folder{i} Name]);
         'loaded'
         exprData = dataSet.mat;
         sampCount = size(exprData, 2);
         gCount = size(exprData, 1);
         totalSampleCount = totalSampleCount + sampCount;
         tissueSampleCount(i) = tissueSampleCount(i) + sampCount;
         tissueSample(i, j) = sampCount;
     end
 end
 

wholeExpr = zeros(gCount, totalSampleCount);

% this is an extra piece of code here just to record the dataset
% names 


for i = 1:length(folder)
    fileList = dir([folder{i} '*.mat'])
    tissue(i).fileList = fileList(:).name;
end

c = 1;
label = repmat({''}, 1, totalSampleCount);

for i = 1:length(folder)
    fileList = dir([folder{i} '*.mat']);
    for j = 1:length(fileList)
        Name = fileList(j).name
        
        load([folder{i} Name]);
        'loaded'
       exprData = (dataSet.mat);
        max(max(dataSet.mat))
        max(max(exprData))
        
        sampCount = size(exprData, 2);
        wholeExpr(:, c:(c + sampCount -1)) = exprData;
        
    [a,b] = regexp(Name, 'GSE[0123456789]+');
    %    label{c} = [ts{i} '-' Name(a:b)];
    label{c} = [Name(a:b)];
    %        label(c) = {Name};
        c = c + sampCount;
    end
end

sib = corr(wholeExpr);

%tissue labels:
label{1} = ts{1};
label{tissueSampleCount(1) + 1} = ts{2};
label{sum(tissueSampleCount(1:2)) + 1} = 'skeletal Muscle';
label{sum(tissueSampleCount(1:3)) + 1} = ts{4};
label{sum(tissueSampleCount(1:4)) + 1} = ts{5};

h = figure; 
heatmap(sib, label, label,  [], 'TickAngle', 45, 'ShowAllTicks', ...
            true, 'Colorbar', true, 'Colormap', 'bone')
    'heatmap done'
    
    
tissue = 'blood, lung, skeletalMuscle, liver'    
temp = sprintf('%d, %d, \n %d, %d', tissueExpCount(1), tissueExpCount(2), ...
               tissueExpCount(3), tissueExpCount(4))
temptitle = sprintf(['%s samples from %s experiments-total of %d ' ...
                        'samples from %d datasets'], tissue, temp, ...
                    length(label), sum(tissueExpCount))
title(temptitle)

dataFolder = '~/data/array/'
fileName = ['bloodLungSMuscleLiverAllSamplesCorrPearsonBERbone']
figFolder = [dataFolder]
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

% plotting the sample counts for each tissue from datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getting the sample counts for each tissue. 

%get all the samples in the directory
dataFolder = '~/data/array/'
tissue = 'skeletalMuscle'
geneExprFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/']
fileList = dir([geneExprFolder, '*.mat'])
sampCount = zeros(length(fileList),1)

%finding the total number of samples
totalSampleCount = 0;
for i = 1:length(fileList)
    Name = fileList(i).name
    load([geneExprFolder Name]);
    'loaded'
    exprData = dataSet.mat;
    sampCount(i) = size(exprData, 2);
    gCount = size(exprData, 1);
    totalSampleCount = totalSampleCount + sampCount;
end

SMS = sampCount;
lungS = sampCount;
bloodS = sampCount;
liverS = sampCount;

% stack plot for sample counts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;

bloodHSMat = [0, 109; 
              61, 59;
              45, 20;
              81, 64;
              11, 33;
              52, 0;
              31, 131;
              20, 10;
              1, 30];

[a, b] = sort(sum(bloodHSMat, 2), 'descend');
bloodHSMat(:,1) = bloodHSMat(b, 1) 
bloodHSMat(:, 2) = bloodHSMat(b, 2)

t = bloodHSMat(:,1);
bloodHSMat(:, 1) = bloodHSMat(:, 2);
bloodHSMat(:, 2) = t;

bar(bloodHSMat, 'stacked')
colormap(summer)

lungHSMat = [58, 0;
             72, 0;
             60, 0;
             75, 0;
             0, 38;
             46, 45;
             23, 6;
             25, 25;
             15, 15;
             40, 0;
             100, 0;
             293, 14;
             246, 0;
             0, 12;
             30, 0];

[a, b] = sort(sum(lungHSMat, 2), 'descend');
lungHSMat(:,1) = lungHSMat(b, 1) 
lungHSMat(:, 2) = lungHSMat(b, 2)

t = lungHSMat(:,1);
lungHSMat(:, 1) = lungHSMat(:, 2);
lungHSMat(:, 2) = t;


hold all
bar([12: 26], lungHSMat, 'stacked')

SMHSMat = [52, 0;
           30, 6;
           0, 68;
           44, 18;
           0, 42;
           50, 0;
           0, 24;
           0, 89];


[a, b] = sort(sum(SMHSMat, 2), 'descend');
SMHSMat(:,1) = SMHSMat(b, 1) 
SMHSMat(:, 2) = SMHSMat(b, 2)

t = SMHSMat(:,1);
SMHSMat(:, 1) = SMHSMat(:, 2);
SMHSMat(:, 2) = t;


bar([29:36], SMHSMat, 'stacked')

liverHSMat = [0, 49;
              20, 8;
              10, 10;
              52, 0;
              12, 12;
              47, 0;
              12, 8;
              63, 0;
              65, 10;
              108, 0;
              15, 7];

[a, b] = sort(sum(liverHSMat, 2), 'descend');
liverHSMat(:,1) = liverHSMat(b, 1) 
liverHSMat(:, 2) = liverHSMat(b, 2)

t = liverHSMat(:,1);
liverHSMat(:, 1) = liverHSMat(:, 2);
liverHSMat(:, 2) = t;

bar([39: 49], liverHSMat, 'stacked')

labels{1} = 'blood';
labels{2} = 'lung';
labels{3} = 'skeletal muscle';
labels{4} = 'liver';


set(gca, 'Xtick', [5, 20, 33, 44], 'XTickLabel', labels)
legend('healthy', 'condition')
ylabel('sample count')

figFolder = [dataFolder 'general/'];
fileName = 'healthyCondition'
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);


% data for healthy and conditioned tissue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bloodHSCount = [109, 0,  ]

% the general sample count
%%%%%%%%%%%%%%%%%%%%%%%%%

lungS = lungS';
lungS = sort(lungS, 'descend');
SMS = [SMS', 0 0 0 0 0 0 0 ]
SMS = sort(SMS, 'descend');
bloodS = [bloodS', 0 0 0 0 0]
bloodS = sort(bloodS, 'descend');
liverS = [liverS', 0 0 0 0]
liverS = sort(liverS, 'descend');
Smat = [liverS; SMS; bloodS; lungS]

labels{1} = 'liver';
labels{2} = 'skeletal muscle';
labels{3} = 'blood';
labels{4} = 'lung';

% defining the colormap to be used in the plot
myColorMap = zeros(15, 3);
oi = [1:2:15]
oe = [2:2:14]
for i = 1:8
myColorMap(oi(i), :) = [100, 20, 20]/255;
myColorMap(oe(i), :) = [230, 150, 150]/255;
end

h = figure;
bar(Smat, 'stacked','EdgeColor', 'none')%, 'LineWidth', 1)
set(gca, 'Xtick', 1:4, 'XTickLabel', labels)
colormap(myColorMap)
ylabel(['Stacked sample counts from all experiments for each ' ...
        'tissue'])
title('samples for each tissue')



figFolder = [dataFolder 'general/'];
fileName = 'sampleCountsTotal'
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);

%the not stack plot for samplecounts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

totSM = sum(SMS)
totLung = sum(lungS)
totBlood = sum(bloodS)
totLiver = sum(liverS)


h = figure;
hold all
plot([0 1 2 3 4 5], [0 0 0 0 0 0 ], '.', 'Color', 'k')
book = ones(size(SMS));
plot(book, SMS, '*')

book = ones(size(lungS));
book = book*2;
plot(book, lungS, '*')

book = ones(size(bloodS));
book = book*3;
plot(book, bloodS, '*')

book = ones(size(liverS));
book = book*4;
plot(book, liverS, '*')

legened({'total number of samples', totSM, totLung, totBlood, totLiver})

label{1} = 'Muscle'
label{2} = 'Lung'
label{3} = 'Blood'
label{4} = 'Liver'
set(gca, 'Xtick', 1:4, 'XTickLabel', label)
xlabel('tissue')
ylabel('number of samples per experiment')

fileName = 'sampleCountsPerTissue';
figFolder = [dataFolder 'general/']
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);

hold all