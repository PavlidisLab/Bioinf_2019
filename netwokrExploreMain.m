% here I work with the files from the infNetwork.m. This is a
% main file cause I am going to call functions from here, which
% work on the whole networks. 
% 1. first part: getting the mean and std for correlation values at
% each expression level of the edges. files are saved as: netExpr01_[GSEID].mat
% 2. second part: getting the average expression level of the genes
% for each tissue, with their STD. 

% the common for loop for loading the networks
clear
mainFolder = '~/networks/tissues/'
tissueList = dir(mainFolder)

% 1. first part: getting the mean and std for correlation values at
% each expression level of the edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getting the number of edges at each expression level and the
% average coexpression value for each file. 
% for every tissue 
for t = 3:length(tissueList)
    tissue = tissueList(t).name;
    fileList = dir([mainFolder tissue '/singleNet/netStr*.mat'])

    for i = 1:length(fileList)
        load([mainFolder tissue '/singleNet/' ...
              fileList(i).name]);
        
        kado = netStr.expNet;
        sib = [4 6 10 14 18 22 9 15 21 27 33 25 35 45 55 49 63 77 ...
               81 99 121]
        counts = zeros(1, length(sib));
        meanCorr = zeros(1, length(sib));
        for j = 1: length(sib)
            
            % counts of the edges with each value
            counts(j) = sum(sum(kado == sib(j)));
            
            % getting all the edges which have the current value
            tempNet = netStr.coNt(kado == sib(j));
            
            % getting the correlation values > 0
            tt = tempNet(tempNet>0);
            meanCorr(j) = mean(tt(:));
            stdCorr(j) = std(tt(:));
        end
        
        % saving it for the file
        netInf.counts = counts;
        netInf.meanCorr = meanCorr;
        netInf.stdCorr = stdCorr;
        
        [a, b] = regexp(fileList(i).name ,'GSE[0123456789]+');
        Name = fileList(i).name(a:b)
        save(sprintf('~/networks/tissues/%s/netInfo/netExpr01_%s.mat', ...
                     tissue, Name), 'netInf')
    end
end

% 2. second part: getting the average expression level of the genes
% for each tissue, with their STD. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getting the number of edges at each expression level and the
% average coexpression value for each file. 
% for every tissue 
clear
mainFolder = '~/networks/tissues/'
tissueList = dir(mainFolder)
gCount = 18498;
expLevel = zeros(gCount, 5);
expStd = zeros(gCount, 5);

for t = 3:length(tissueList)
    tissue = tissueList(t).name;
    fileList = dir([mainFolder tissue '/singleNet/netStr*.mat'])
    
    expLevles = zeros(gCount, length(fileList));
    for i = 1:length(fileList)
        load([mainFolder tissue '/singleNet/' ...
              fileList(i).name]);
        
        qCount = 6;
        kado = [0:1/6:1]
        Qs = quantile(netStr.avgExpr, kado)
        qSize= floor(gCount/qCount)
        qExpr = ones(1, gCount);
        for j = 2:qCount
            tempInd = netStr.avgExpr >= Qs(j);
            qExpr(tempInd) = qExpr(tempInd) + 1;
        end
        
        expLevels(:, i) = qExpr;
    end
    
    tempExpLevel = mean(expLevels');
    tempExpStd = std(expLevels');
    tissues{t-2} = tissue;
    expLevel(:, t-2) = tempExpLevel;
    expStd(:, t-2) = tempExpStd;
end

tisStr.expLevel = expLevel;
tisStr.expStd = expStd;
tisStr.tissues = tissues;

% data structure is different in the two files:
save('~/data/general/tissueExprLevelInfo01.mat','tisStr')
save('~/data/general/tissueExprLevelInfo.mat','tisStr')