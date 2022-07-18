% getting the 3 plots of correlation distribution. 

clear
dataFolder = '~/data/general/'
load([dataFolder 'GPL570GemmaMapNEW.mat'])

dataFolder = '~/data/array/'
addpath('~/codes/MATLAB/myCodes/strings/')
addpath('~/codes/MATLAB/myCodes/general/')

gCount = length(gpl570.uniqueSymbols)

folderList = dir(dataFolder);
tissues = {folderList([4 5 7 8 9]).name}
tissueCorrStruct = struct;

for j = 1: length(tissues)

    tissue = tissues{j};
    probeFolder = [dataFolder tissue ['/textFiles/' ...
                        'probeComBatBatchEffectRemoved/']];
    
    fileList = dir([probeFolder '*.txt'])

    % reading the files for each tissue
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
        for k = 1:sampleCount
            GEOID(k) = regexp(headers{k}, expression, 'match');
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
        for k = 1:length(data{1})
            if(strcmp(data{1}{k}(1), '"'))
                data{1}{k} = data{1}{k}(2:end);
            end
            if(strcmp(data{1}{k}(end), '"'))
                data{1}{k} = data{1}{k}(1:(end-1));
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
        
        tic
        holu = corr(probeDataSet.mat');
        toc
        
        upHolu = holu(logical(triu(ones(size(holu)), 1)));
        y = randsample(length(upHolu), floor(length(upHolu)/1000));
        q92per(i) = quantile(upHolu(y),0.92);
        q94per(i) = quantile(upHolu(y),0.94);
        q95per(i) = quantile(upHolu(y),0.95);
        q96per(i) = quantile(upHolu(y),0.96);
        q97per(i) = quantile(upHolu(y),0.97);
        q98per(i) = quantile(upHolu(y),0.98);
        q99per(i) = quantile(upHolu(y),0.99);
    end

    tissueCorrStruct(j).tissue = tissue;
    tissueCorrStruct(j).q92per = q92per;
    tissueCorrStruct(j).q94per = q94per;
    tissueCorrStruct(j).q95per = q95per;
    tissueCorrStruct(j).q96per = q96per;
    tissueCorrStruct(j).q97per = q97per;
    tissueCorrStruct(j).q98per = q98per;
    tissueCorrStruct(j).q99per = q99per;
    tissueCorrStruct(j).fileList = fileList; 

end


save([dataFolder 'corrQAllTissues.mat'], 'tissueCorrStruct')

% printing it in a .txt file 

file = [dataFolder 'corrQ92AllTissues.txt'];
fid = fopen(file, 'w')
fprintf(fid, ['tissue\t experiment\t Q92\t Q94\t Q95\t Q96\t ' ...
              'Q97\t Q99\n']);
for i = 1:5
    for j = 1: length(tissueCorrStruct(i).fileList)
        fileName = tissueCorrStruct(i).fileList(j).name;

        [a,b] = regexp(fileName, 'GSE[0123456789]+')
        GSEID = fileName(a:b)
                
        fprintf(fid, ['%s\t %s\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n'], ...
                tissueCorrStruct(i).tissue, ...
                GSEID, ...
                tissueCorrStruct(i).q92per(j), ...
                tissueCorrStruct(i).q94per(j), ...
                tissueCorrStruct(i).q95per(j), ...
                tissueCorrStruct(i).q96per(j), ...
                tissueCorrStruct(i).q97per(j), ...
                tissueCorrStruct(i).q98per(j), ...
                tissueCorrStruct(i).q99per(j));
        i
        j
    end
end
fclose(fid)


