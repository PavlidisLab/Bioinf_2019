% you can just use this as is, if you give it the right inputs

clear

% the general folder for my data, you can just change it to the
% path where your data is
dataFolder = '~/data/array/'
tissue = 'lung'

% the general folder for where all the cel files from different
% experiments are
celFolder = [dataFolder tissue '/celFiles/'];

% Folder name where the CEL files for one
% experiment are. The format I had on file was:
% GSE10041_RAW 
fileList = dir([celFolder '*_RAW'])

% still on the server, no need ot change that
libFile = '/space/grp/marjan/data/HG-U133_Plus_2.cdf';

% for each dataset folder
for i = 1:length(fileList)
    tic
    if(fileList(i).isdir == 1) 
        experimentFolder = [celFolder fileList(i).name]
        
        %getting list of the cel files here
        tempCelFileList = dir([experimentFolder '/*.CEL']);
        celFileList = cell(1, length(tempCelFileList));
        for j = 1:length(tempCelFileList)
                celFileList{j} = tempCelFileList(j).name;
        end
        
        expr = affyrma(celFileList, libFile, 'CELPath', ...
                           ['/home/mfarahbod/data/array/' tissue ...
                            '/celFiles/' fileList(i).name]);
        
        % just making the file name for each dataset
        tempName = fileList(i).name;
        name = tempName(1:(end - 4))

        % writing it into .txt file (instead of .mat file), so I
        % can load it in R
        dmwrite(expr, [dataFolder tissue '/textFiles/probeTxtFromCEL/' name '.txt' ...
                      ]);
        clear expr
     end
    toc
end
