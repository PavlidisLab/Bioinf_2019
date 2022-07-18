%this code is for exploring the array data, starting with the
%expression matrix which is already normalized. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%getting the 0.05
%open the expression file

%AIM
% 1. convergence for the highest expression level in each dataset
% 2. convergence for the highest expression level in each tissue
% 3. investigating the highest expression values 

%Further exploration: 
% do we have a group of genes which would cluster in a tissue? >>
% clustering genes in a tissue, and recognizing different
% expression groups within each tissue. >> for this part, I need to
% find a way to normalize the number of samples we are taking from
% each dataset. maybe cut the hierarchial clustering somewhere? 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the code below gets the top 1000 expressed genes among all
% samples within datasets, and also all datasets within a tissue. 
% saved are a cell containg the top ranking genes for each experiment(experimentsTotHigh),
% the convergence in a experiment(wholeIncomm), the convergence between
% experiments within a tissue(tissueIncomm)
% saved in file named topExpConv , topExpressedGeneConvergence_[tissue].mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the ranked top values for all samples 
'------------------------------------------'
clear
tissue = 'skeletalMuscle'
dataFolder = '/space/grp/marjan/data/'
exprFolder = [dataFolder tissue '/geneExprMat/clean/']
fileList = dir(exprFolder)
gCount = 21050;


wholeIncomm = cell(1, length(fileList)-2);
Q = 0.99; % the quantile
tissueTotal = zeros(gCount, 1);
tissueIncomm = zeros(length(fileList)-3, 1);
topGeneCount = round((1-Q)*gCount) + 1;
experimentsTotHigh = zeros(length(fileList)-2, (round((1-Q) * ...
                                                  gCount)+1)); % keeping the top ranking genes from each experiment
                                                
size(experimentsTotHigh)
% for each dataset:
for i = 3: length(fileList)  
    i
    load([exprFolder fileList(i).name]);
    sampleCount = size(dataSet.mat, 2)
    expIncomm = zeros(sampleCount-1, 1); %keeping percent of the
                                         %genes incommon among the
                                         %samples within an experiment
    totalExp = tiedrank(dataSet.mat(:,1)); % keeping  the total
                                           % gene ranking for each experiment
    q = quantile(totalExp, Q);
    totHighExp = find(totalExp > q); %keeping the top ranking genes
                                     %index for each experiment
    
    'first'

    % for each sample: 
    for j = 2: sampleCount
        %        j
        oldTotHighExp = totHighExp; % the old set of high ranking genes
        
        %get the current sample's high expressed genes 
        sampRank = tiedrank(dataSet.mat(:, j));
        q = quantile(sampRank, Q);
        sampHighExp = find(sampRank > q);
        
        % sum+ the ranked values with the dataset total
        totalExp = totalExp + sampRank;
        q = quantile(totalExp, Q);
        totHighExp = find(totalExp > q);
                
        expIncomm(j-1) = length(intersect(oldTotHighExp, totHighExp))/((1 ...
                                                          - Q) * ...
                                                          gCount);

    end

            'DONE experiment'
    if (i == 3) % this is the first loop
        totHighTissue = totHighExp; 
        experimentsTotHigh(i-2, (1:length(totHighExp))) = totHighExp;
        totalTissue = totalExp/sampleCount; % of the current
                                            % experiment
        wholeIncomm{i - 2} = expIncomm;
        'second'
    else
        
        oldTotHighTissue = totHighTissue;
        'third'

        %saving the experiments high ranked genes
        experimentsTotHigh(i-2, (1:length(totHighExp))) = totHighExp; 
        
        'fourth'
        
        % saving the whole experiment incomon
        wholeIncomm{i - 2} = expIncomm;
        
        %adding totalTissue to the mean samplerank of the experiment
        totalTissue = totalTissue + (totalExp/sampleCount); 
        
        q = quantile(totalTissue, Q);
        totHighTissue = find(totalTissue > q);
        
        tissueIncomm(i - 3) = length(intersect(oldTotHighTissue, ...
                                               totHighTissue))/((1-Q)*gCount);
    end
end
% end
topExpConv.WholeConv = wholeIncomm;
topExpConv.tissueConv = tissueIncomm;
topExpConv.expTotHigh = experimentsTotHigh;
topExpConv.totHigh = totHighTissue;
fileName = sprintf('top%dExpressedGeneConvergence_%s.mat', topGeneCount, ...
                   tissue)
folder = [dataFolder tissue '/generalResults/']
save([folder fileName], 'topExpConv')
% saved in file named topExpConv ,
% topExpressedGeneConvergence_[tissue].mat

%plotting the results from the above code : 
% data:
%1. the wholeConv - vector of convergence values among samples
%within experiment
%2. the experimentsTotHigh(this is list of genes, I want to compare
%it within different tissues)
%3. the tissueConv - vector of convergence values among experiments

clear
tissue = 'lung'
topGeneCount = 212;
dataFolder = '/space/grp/marjan/data/';
%load datafile for the tissue,
load(sprintf('%s%s/generalResults/top%dExpressedGeneConvergence_%s.mat', ...
             dataFolder, tissue, topGeneCount, tissue));
top1 = topExpConv;
clear topExpConv

topGeneCount = 1054;
load(sprintf('%s%s/generalResults/top%dExpressedGeneConvergence_%s.mat', ...
             dataFolder, tissue, topGeneCount, tissue));
top2 = topExpConv;
clear topExpConv

% Plotting
%%%%%%%%%%%%%%%%%% 

%%%% 1.
%plotting the convergence within samples for each experiment within
%tissue : I need three of these plots, since I have three tissues
h = figure();

%hold all
for i = 1: 5%length(topExpConv.WholeConv)
    subplot(5, 1, i)
    x = top1.WholeConv{i};
    y = top2.WholeConv{i};
    hold all
    plot(x, 'LineWidth', 2)
    plot(y, 'LineWidth', 2)
end

%get the experiment name, 
% put the title
% put the general xlabel and general y label
% save them

%%%% 2.
% plot the convergence among tissues, 3 of these plots too. 
hold all
plot(top1.tissueConv, 'LineWidth', 2)
plot(top2.tissueConv, 'LineWidt', 2)

%title, labels, save, legend

%%% 3. 
%incommong among them three
clear
tissue = 'lung'
topGeneCount = 1054;
dataFolder = '/space/grp/marjan/data/';
%load datafile for the tissue,
load(sprintf('%s%s/generalResults/top%dExpressedGeneConvergence_%s.mat', ...
             dataFolder, tissue, topGeneCount, tissue));
toplung = topExpConv.totHigh;
clear topExpConv

tissue = 'blood'
%topGeneCount = 212;
dataFolder = '/space/grp/marjan/data/';
%load datafile for the tissue,
load(sprintf('%s%s/generalResults/top%dExpressedGeneConvergence_%s.mat', ...
             dataFolder, tissue, topGeneCount, tissue));
topblood = topExpConv.totHigh;
clear topExpConv

tissue = 'skeletalMuscle'
%topGeneCount = 212;
dataFolder = '/space/grp/marjan/data/';
%load datafile for the tissue,
load(sprintf('%s%s/generalResults/top%dExpressedGeneConvergence_%s.mat', ...
             dataFolder, tissue, topGeneCount, tissue));
topSM = topExpConv.totHigh;
clear topExpConv

length(intersect(topblood, toplung))
length(intersect(toplung, topSM))

%TODO plots

% Comparison between these and high degree nodes
% then the clustering among tissues and their similarities. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%getting the mean expression level among experiments and tissues
clear
tissue = 'lung'
dataFolder = '/space/grp/marjan/data/'
exprFolder = [dataFolder tissue '/geneExprMat/clean/']
fileList = dir(exprFolder)
gCount = 21050;

total = zeros(1, gCount);
% for each dataset:
for i = 3: length(fileList)  
    i
    load([exprFolder fileList(i).name]);
    sampleCount = size(dataSet.mat, 2)
    expMean = mean(dataSet.mat, 2);
    expMean = tiedrank(expMean);
    total = total + expMean';
    
end
bloodresult = total./(length(fileList) - 2);
lungresult  = total./(length(fileList) - 2);
SMresult = total./(length(fileList) - 2);

%% loading the three networks
% it is loaded running the networkExploration.m code

corr(bloodNdegree', bloodresult') % 0.190
corr(lungNdegree', lungresult') % 0.028 
corr(SMNdegree', SMresult') % 0.173

%%% clustering and consistency among expression levels 
%GET THE PRINTER
%
