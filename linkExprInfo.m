% this is a file which has the information for studying the link
% expression 
% 1. savign the probe expression data and it's related
% information. 
% 2. savinge the gene expression data. the information saved for
% probe expression can be shared 

%reading the probe text file
clear

addpath('~/codes/MATLAB/myCodes/general/')
dataFolder = '~/data/general/'
load([dataFolder 'GPL570GemmaMapNEW.mat'])

dataFolder = '~/data/affyArray/tissues/'
addpath('~/codes/MATLAB/myCodes/strings/')
addpath('~/codes/MATLAB/myCodes/general/')

ts{1} = 'blood';
ts{2} = 'lung';
ts{3} = 'skeletalMuscle';
ts{4} = 'liver';
ts{5} = 'brain';

% total sample count: I am getting it from arrayMain04.m, it is
% 3790 for the batch effect removed files. I am getting 4k as
% sample count

totalSampleCount = 4000;
probeCount = length(gpl570.probeID);

wholeExpr = zeros(probeCount, totalSampleCount);
sampleCounter = 1;
dataSetCounter = 1;

for t = 1:length(ts)
    tissue = ts{t}
    probeFolder = [dataFolder tissue ['/textFiles/' ...
                        'probeComBatBatchEffectRemoved/']]

    fileList = dir([probeFolder '*.txt'])
    length(fileList)
    gCount = length(gpl570.uniqueSymbols)

    for i = 1:length(fileList)
    Name = fileList(i).name
    linum = 1;
    fileName = [probeFolder Name]
    fid = fopen(fileName)
    
    header = textscan(fid, '%s', 1, 'delimiter', '\n', 'headerlines', ...
                      linum - 1, 'Bufsize', 100000095);
    
    %    getting the header:
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
                    'Headerlines', linum-1, 'Delimiter', '\t','Bufsize', ...
                    100000095);
    fclose(fid);
        
    pCount = length(data{1});
    dataExpr = cell2mat(data(2:end));

    for j = 1:length(data{1})
        if(strcmp(data{1}{j}(1), '"'))
            data{1}{j} = data{1}{j}(2:end);
        end
        if(strcmp(data{1}{j}(end), '"'))
            data{1}{j} = data{1}{j}(1:(end-1));
        end
    end

    [a, b] = ismember(data{1}, gpl570.probeID);
    geneID = gpl570.ncbiID(b(a)); %gene ID 33k
    geneSymbols = gpl570.geneSymbols(b(a)); %gene NCBI ID 33k
    myProbes = data{1}(a); %probe names - 33k
    myDataExpr = dataExpr(a ,:); %probe values
    pCount = length(myProbes); 
    
    geneMap = containers.Map(gpl570.uniqueSymbols, [1:gCount]);
    probeMap = containers.Map(myProbes, geneSymbols);%  NCBI geneID

    probeDataSet.mat = myDataExpr;
    probeDataSet.GEOID = GEOID;
    probeDataSet.Gemma = 'null';
    
    % %%    removing outliers - yes, I am doing it from the probe matrix
    BPviewOnly = 1;
    BPfigFolder = ''
    %outlier removal here
    samplesSBM = find(outlierDetectionSBM(probeDataSet.mat, BPviewOnly, ...
                                          BPfigFolder,Name) == 1)

    if(sampleCount - length(samplesSBM) > 3)
        ORProbeDataSet = removeSamples04(probeDataSet ,samplesSBM);
        ORSampleCount = size(ORProbeDataSet.mat, 2);
    else
        ORProbeDataSet = probeDataSet;
        ORSampleCount = size(probeDataSet.mat, 2)
    end


    % so far, I have the probe dataset for this dataset, I add it
    % to the expr folder. 
    wholeExpr(:, sampleCounter:(sampleCounter + ORSampleCount - 1)) = ...
        ORProbeDataSet.mat;
    sampleCounter = sampleCounter + ORSampleCount;


    divMat = zeros(gCount, ORSampleCount);

    %    adding the maps
    tic
    for j = 1:pCount
        sib = values(geneMap, values(probeMap, myProbes(j)));
        divMat(sib{1}, :) = divMat(sib{1}, :) + 1;
    end
    toc

    % divmat has the same columns. each column of divmat has the
    % value that tells which probes are being used for each gene
    % originally. 

    % max probe count for each gene - I am not using divMat cause I
    % am going to filter the probes based on their expression levels
    probeInd = zeros(gCount, max(max(divMat))); 
    for j = 1:pCount
        sib = values(geneMap, values(probeMap, myProbes(j)));
        temp = find(probeInd(sib{1}, :) == 0);
        pInd = temp(1);
        probeInd(sib{1}, pInd) = j;
    end

    
    probeMean = mean(ORProbeDataSet.mat, 2);
    a = mean(probeMean)
    newGeneExpr = zeros(gCount, ORSampleCount);
    affectedGeneCount = 0;
    
    % for each dataset, set thr as the 2/3 of the expression level
    % get the 2/3 quantile, 
    q = quantile(ORProbeDataSet.mat(:), 2)
    thr = q(1);

    % TODO HERE: in the dataSet specific vector which is size of
    % the probeCount, mark which probes are being used in this
    % dataSet. 

    probePresenceVector = zeros(1, probeCount);

    for j = 1: length(probeInd) % this is the same as geneCount
        pCount = sum(probeInd(j,:)>0); % how many probes does that
                                       % gene have

        %if there are any with exp > thr, get the mean of bigger ones
        largeProbes = find(probeMean(probeInd(j,1:pCount))> thr);

        if(length(largeProbes) > 0)
            largeProbesInd = probeInd(j, largeProbes);
            newGeneExpr(j,:) = mean(ORProbeDataSet.mat(largeProbesInd, ...
                                                       :), 1);

            
            probePresenceVector(largeProbesInd) = 1;
    
            if(pCount > length(largeProbes))
                affectedGeneCount = affectedGeneCount + 1;
            end
        else
            newGeneExpr(j, :) = mean(ORProbeDataSet.mat(probeInd(j, ...
                                                              1:pCount), ...
                                                        :), 1);

            probePresenceVector(probeInd(j, 1:pCount)) = 1;
        end
    end
    
    sum(probePresenceVector)
    %     finalMat = newGeneExpr;

    % for each dataSet, save the name, the sampleInd, the index and
    % the vector
    dataSetProbeInf(dataSetCounter).name = getGSEIDfromStr(fileName);
    dataSetProbeInf(dataSetCounter).sampleInd = [sampleCounter - ORSampleCount ,sampleCounter - 1];
    dataSetProbeInf(dataSetCounter).probePresenceVector = ...
        probePresenceVector;
    dataSetProbeInf(dataSetCounter).tissue = tissue;
   
    dataSetCounter = dataSetCounter + 1
    
    end
end

wholeExpr = wholeExpr(:, 1:sampleCounter);

% save probeInd, dataSetProbeInf and wholeExpr

save('~/data/general/linkExprInfo/probeInd.mat', 'probeInd')
save('~/data/general/linkExprInfo/dataSetProbeInf.mat', ...
     'dataSetProbeInf')
save('~/data/general/linkExprInfo/wholeExpr.mat', 'wholeExpr')

% 2. collecting the expression data at gene level 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

addpath('~/codes/MATLAB/myCodes/general/')
dataFolder = '~/data/general/'
load([dataFolder 'GPL570GemmaMapNEW.mat'])

dataFolder = '~/data/affyArray/tissues/'
addpath('~/codes/MATLAB/myCodes/strings/')
addpath('~/codes/MATLAB/myCodes/general/')

ts{1} = 'blood';
ts{2} = 'lung';
ts{3} = 'skeletalMuscle';
ts{4} = 'liver';
ts{5} = 'brain';

gCount = 18494;

wholeGeneExpr = zeros(gCount, 3600);
sCounter = 1;
for t=1:length(ts)
    tissue = ts{t};
    fileList = dir([dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/*.mat'])
    for i = 1:length(fileList)
        load([dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/' ...
              fileList(i).name]);
        fileList(i).name
        thisSCount = size(dataSet.mat, 2);
        thisSCount
        wholeGeneExpr(:, sCounter:sCounter + thisSCount -1) = ...
            dataSet.mat;
        sCounter = sCounter + thisSCount;
    end
end

wholeGeneExpr = wholeGeneExpr(:, 1:3563);

save('~/data/general/linkExprInfo/wholeGeneExpr.mat', ...
     'wholeGeneExpr')

% making and saving the gene Expr and STD for the 10th time!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('~/data/general/linkExprInfo/wholeGeneExpr.mat')
load('~/data/general/linkExprInfo/dataSetProbeInf.mat')

for i = 1:length(dataSetProbeInf)
    s = dataSetProbeInf(i).sampleInd(1);
    e = dataSetProbeInf(i).sampleInd(2);
    book = wholeGeneExpr(:, s:e);
    meanExp = mean(book');
    dataSetExpInf(i).meanExp = meanExp;
    dataSetExpInf(i).name = dataSetProbeInf(i).name;
    dataSetExpInf(i).sampleInd = dataSetProbeInf(i).sampleInd;
    dataSetExpInf(i).tissue = dataSetProbeInf(i).tissue;
end

    
    