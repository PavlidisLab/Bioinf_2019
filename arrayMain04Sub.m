% I chopped this part of arrayMain04.m code for making the probe
% correlation density plots. I save the results for each tissue in
% a file named [tissue]_probeCorrDistribution.mat 

clear
dataFolder = '/space/grp/marjan/data/'
load([dataFolder 'GPL570GemmaMapNEW.mat'])

dataFolder = '~/data/array/'
addpath('~/codes/MATLAB/myCodes/strings/')
addpath('~/codes/MATLAB/myCodes/general/')
tissue = 'lung'
probeCorr = struct
gCount = length(gpl570.uniqueSymbols)

%temp folder
%probeFolder = [dataFolder tissue ['/lungTumorSamples/textFiles/' ...
                    'probeComBatBatchEffectRemoved/']];

%permanent folder 
probeFolder = [dataFolder tissue ['/textFiles/' ...
                    'probeComBatBatchEffectRemoved/']];
               
fileList = dir([probeFolder '*.txt'])

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
    
    
    [a,b] = regexp(Name, 'GSE[0123456789]+')
    SName = Name(a:b)
    probeCorr(i).dataSet = SName;
    
    %I need length of header for the sampleCount
    data = textscan(fid, [repmat('%s', 1,1) repmat('%f', 1, ...
                                                   sampleCount)], ...
                    'Headerlines', linum-1, 'Delimiter', '\t','Bufsize', ...
                    100000095);
    fclose(fid);
    
    pCount = length(data{1});
    dataExpr = cell2mat(data(2:end));
    
    %removing the extra ""
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
    %                                                   map to
    %                                                   probes
        % the new probe gene mapping with filter - these codes are
    % taken from the probeFiltering.m
    % >>>>>>>>>>>>>>>>>>>>>
    
    % making the probeDataSet
    probeDataSet.mat = myDataExpr;
    probeDataSet.GEOID = GEOID;
    probeDataSet.Gemma = 'null';
    
    % %%    removing outliers - yes, I am doing it from the probe matrix
    BPviewOnly = 1;
    BPfigFolder = ''
    %outlier removal here
    samplesSBM = find(outlierDetectionSBM(probeDataSet.mat, BPviewOnly, BPfigFolder,Name) == 1)
    ORProbeDataSet = removeSamples04(probeDataSet ,samplesSBM);
    ORSampleCount = size(ORProbeDataSet.mat, 2);
    
    %%%% here, get the probe correlation dist - save all the dist
    %%%% functions for a tissue. 
    tic
    holu = corr(probeDataSet.mat');
    toc
    
    upHolu = holu(logical(triu(ones(size(holu)), 1)));
    X = randi(length(upHolu), [1, floor(length(upHolu)/10)]);
    tic
    [f, xi] = ksdensity(upHolu(X));
    toc
    i
    
    probeCorr(i).dist = [xi;f];

end

tissue

% temp folder
%save([dataFolder tissue '/lungTumorSamples/matFiles/probeCorrelationDistribution.mat'], ...
     'probeCorr')

% permanent folder
save([dataFolder tissue '/matFiles/probeCorrelationDistribution.mat'], ...
     'probeCorr')