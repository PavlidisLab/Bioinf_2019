% this file is to get the batch info for experiments from their
% .CEL files

clear

addpath('~/codes/MATLAB/myCodes/general/')
dataFolder = '~/data/array/'

tissue = 'skeletalMuscle'
wholeExpFolder = dir([dataFolder tissue '/celFiles/'])

%13 remians
for i = 3: length(wholeExpFolder)
    i = 3
    fileList = dir([dataFolder tissue '/celFiles/' ...
                    wholeExpFolder(i).name] )

    date = cell(length(fileList)-2, 1);
    sampleName = cell(length(fileList) -2, 1);
    expName = getGSEIDfromStr(wholeExpFolder(i).name)
    
    expFolder = [dataFolder tissue '/celFiles/' wholeExpFolder(i).name ...
                '/']
    
    % to get the dates from a folder of .CEL files
    
    %    [status, cmdout] = system(sprintf('head -n 18 %s | awk ''BEGIN{FS=" "}{if(NF>4)print $8}''', expFolder ))
  
    %    dates = textscan();
    %    this was the old one, I got the new command from Sanja
    for j = 3: length(fileList)
         [status, cmdout] = system(sprintf(['grep -n -a ''^DatHeader'' %s ' ...
                            '|cut -f1 -d:'], [expFolder ...
                            fileList(j).name]));
        
         if(isempty(cmdout))
            'file is empty'
            date(j-2) = {'null'};
            nameExp = 'GSM[0-9]+'
            sampleName(j-2) = regexp(fileList(j).name, nameExp, ...
                                     'match');
            
        else
            dateLine = str2num(cmdout)
            

            fid = fopen([dataFolder tissue '/celFiles/' wholeExpFolder(i).name ...
                         '/' fileList(j).name])
            
            
            header = textscan(fid, '%s', 1, 'delimiter', '\n', 'headerlines', dateLine - 1, ...
                              'Bufsize', 1000000095)
            
            headStr = header{1}{1};

            expression = '\d*/\d*/\d*';
            book = regexp(header{1}, expression, 'match')

            date(j-2) = book{1};

            nameExp = 'GSM[0-9]+'
            sampleName(j-2) = regexp(fileList(j).name, nameExp, ...
                                     'match');
        end
    end
    
    sample.date = date;
    sample.name = sampleName;
    
    % get the batch vector 
    batchMap = containers.Map(unique(date), [1:length(unique(date))]);
    batchArray = cell2mat(values(batchMap, date));
    sample.batchArray = batchArray;
    
    batchInfoFile = [dataFolder tissue '/' expName 'batchInfo.mat']
    save(batchInfoFile, 'sample');
    
    % saving it into a .text file
    fileName = [dataFolder tissue '/' expName '_sampleInfo.txt']
    fid = fopen(fileName, 'w')
    fprintf(fid, 'sampleID\tdate\tbatch\n')
    for k = 1:length(fileList)-2
        fprintf(fid, '%s\t%s\t%d\n', sampleName{k}, date{k}, batchArray(k))
    end
    fclose(fid)
end

%lung 13 remains
% liver 9
%sm 8
