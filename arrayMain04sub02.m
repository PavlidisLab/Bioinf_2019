% get the sample count for each dataset for each tissue, save them
% in a file. 

clear
dataFolder = '~/data/array/';
tissue = 'blood'
exprFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/'];
files = dir([exprFolder '*.mat'])

exprInfo = struct
for i = 1:length(files)
    Name = files(i).name;
    [a,b] = regexp(Name, 'GSE[0123456789]+')
    SName = Name(a:b);
    exprInfo(i).name = SName;
    
    load([exprFolder Name])
    exprInfo(i).sampleCount = size(dataExpr, 2)
end

saveFolder = [dataFolder tissue '/matFiles/']
save([saveFolder tissue '_exprInfo.mat'], 'exprInfo')

