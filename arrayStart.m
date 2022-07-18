% http://www.mathworks.com/products/demos/bioinfo/affydemo/affydemo.html
% This file has the steps from array experiment in GEO, to the matrix of raw probe/sample in a
% .txt file
% % step 1. found the desired experiment in GEO, load all the .CEL
% files in a folder. Also, load the platform .cdf file in the
% parent directory. Sample experiment : GSE33223

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Exploration
  
celStruct = affyread('GSE7753_RAW/GSM187584.CEL');
celStruct.ProbeColumnNames
celStruct.Probes(1, :)

cdfStruct = affyread('/space/grp/marjan/data/blood/HG-U133_Plus_2.cdf');
cdfStruct
cdfStruct.ProbeSets
cdfStruct.ProbeSets(100)

%to test the normalization method on files : experiment : GSE7753
  
expr = affyrma('GSE7753_RAW/GSM187584.CEL',['/space/grp/marjan/' ...
                    'data/blood/HG-U133_Plus_2.cdf']);

celFileList = dir(['/space/grp/marjan/prioritization Project/' ...
                   'MATLABcodes/GSE7753_RAW/']);
%getting the file names into cell array 
book = cell(1, length(celFileList));
for i = 1:length(celFileList)
    book{i} = ['GSE7753_RAW/' celFileList(i).name];
end

book = book(3:end);

expr = affyrma(book, ['/space/grp/marjan/data/blood/HG-' ...
                    'U133_Plus_2.cdf']);

book = double(expr);
sib = corr(book);
heatmap(sib);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for the blood data getting the probe expression files and sanity
%check the matrix

%open the expr file and read the names. 


fileName = '/space/grp/marjan/data/blood/theBloodDataInfoClean.txt';
fid = fopen(fileName);

expInf = textscan(fid, ['%s', '%d', '%s', '%d', '%d'], 'headerlines', 1, ...
                  'Delimiter', '\t')
fclose(fid)
experimentList = expInf{1};
% experimentCount = length(experimentList);
% infCol = expInf{2};
% infStr = expInf{3};
% infColCount = expInf{4};
% hlineCount = expInf{5};

ofor i = 1:length(experimentList)
    celFolder = [experimentList{i} '_RAW'];
    celPath = ['/space/grp/marjan/data/blood/celData' celFolder];
    libFile = '/space/grp/marjan/data/blood/HG-' ...
              'U133_Plus_2.cdf';

    tempCelFileList = dir([celPath celFolder]);
    celFileList = cell(1, length(tempCelFileList));
    for i = 1:length(tempCelFileList)
        celFileList{i} = [tempCelFileList(i).name];
    end
    celFileList = celFileList(3:end);

    expr = affyrma(celFileList, libFile, 'CELPath', celPath);
    
    %saving the .txt file with the header and row names. NOTE : I
    %need to change the header linecounts in the main array
end