% this is to study regression : 
% 3. Getting regression for the list of TS links. HUH. 
% 4. finding the Thr for the rsqrd

% for a given link, get the expression level and correlation rank. 

ID1 = 234 
ID2 = 427

clear
load('~/data/general/linkExprInfo/wholeGeneExpr.mat')
load('~/data/general/linkExprInfo/dataSetProbeInf.mat')

% for each dataset

dsCount = length(dataSetProbeInf)
n = 53
featureMat = zeros(n, 2);

featureMat = zeros(n, 7);
YMat = zeros(n, 1);
c = zeros(n, 1);

ID1 = 4 
ID2 = 74
for i = 1:53

    % ID1 = randi(18494, 1);
    % ID2 = randi(18494, 1);
    ds = i
    s = dataSetProbeInf(ds).sampleInd(1);
    e = dataSetProbeInf(ds).sampleInd(2);
    exp1 = wholeGeneExpr(ID1, s:e);
    exp2 = wholeGeneExpr(ID2, s:e);
    
    c(i) = corr(exp1', exp2');

    featureMat(i, 2) = e - s + 1;

    YMat(i) = corrQuantileReturn(dataSetProbeInf(ds).name, ...
       c(i))/1000;


    expMat = wholeGeneExpr(:, s:e);
    avgExp = mean(expMat, 2);
    avgExpC = avgExp - mean(avgExp);

    % featureMat(i, 1) = avgExpC(ID1);
    % featureMat(i, 2) = avgExpC(ID2);

    % lexp = min(avgExpC(ID1), avgExpC(ID2));
    
     tissue = dataSetProbeInf(ds).tissue;
     kado = unique({dataSetProbeInf(:).tissue});
     [a, b] = ismember(tissue, kado);
     
     b + 2
     featureMat(i, b + 2) = 1;
        featureMat(i, 1) = min(avgExpC(ID1), avgExpC(ID2));
        %    featureMat(i) = min(avgExpC(ID1), avgExpC(ID2));
end

tic
book = LinearModel.fit(featureMat, c)
toc

err = zeros(1, n);
for i = 1:n
    err(i) = c(i) - (-4.14e-05 * featureMat(i, 2) +...
                     -0.016 * featureMat(i, 3) +...
                     0.357 * featureMat(i, 4)+...
                     -0.388 * featureMat(i, 5)+...
                     0.144 * featureMat(i, 6));
end
h = figure
hist(err, [-1:.1:1])
hold on

sib = hist(c, [-1:.1:1]);

plot([-1:.1:1], sib, '*', 'Color', 'r')

he = (featureMat > -1) && (featureMat < 0);
corr(featureMat(he), YMat(he))
corr(featureMat(he), c(he))

h = figure
plot(featureMat, YMat, 'o')
corr(featureMat, YMat)

h = figure
plot(featureMat(:, 1), c, 'o')
corr(featureMat, c)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
load('~/data/general/linkExprInfo/wholeGeneExpr.mat')
load('~/data/general/linkExprInfo/dataSetProbeInf.mat')

% for each dataset

dsCount = length(dataSetProbeInf)
n = 53
featureMat = zeros(n, 2);

featureMat = zeros(n, 7);
YMat = zeros(n, 1);
c = zeros(n, 1);


ID1 = 5914
ID2 = 15659
tic
for i = 1:53

    % ID1 = randi(18494, 1);
    % ID2 = randi(18494, 1);
    %    ds = i
    
    
    s = dataSetProbeInf(ds).sampleInd(1);
    e = dataSetProbeInf(ds).sampleInd(2);
    exp1 = wholeGeneExpr(ID1, s:e);
    exp2 = wholeGeneExpr(ID2, s:e);
    
    c(i) = corr(exp1', exp2');

    YMat(i) = corrQuantileReturn(dataSetProbeInf(ds).name, ...
       c(i))/1000;

    expMat = wholeGeneExpr(:, s:e);
    avgExp = mean(expMat, 2);
    avgExpC = avgExp - mean(avgExp);

    avgExpC = avgExpC + abs(min(avgExpC));

    featureMat(i, 1) = avgExpC(ID1) + avgExpC(ID2);
    featureMat(i, 2) = abs(avgExpC(ID1) - avgExpC(ID2));

    % featureMat(i, 1) = avgExpC(ID1);
    % featureMat(i, 2) = avgExpC(ID2);

    % lexp = min(avgExpC(ID1), avgExpC(ID2));
    
     tissue = dataSetProbeInf(ds).tissue;
     kado = unique({dataSetProbeInf(:).tissue});
     [a, b] = ismember(tissue, kado);
     
     %     b + 2;
     featureMat(i, b + 2) = 1;
        featureMat(i, 1) = min(avgExpC(ID1), avgExpC(ID2));
        %    featureMat(i) = min(avgExpC(ID1), avgExpC(ID2));
end

% getting the regression model
toc
tic
lm1 = LinearModel.fit(featureMat, YMat, 'Intercept', false)

lm2 = LinearModel.fit(featureMat(:, 1:2), YMat, 'Intercept', false)
toc

% do the code for 1000 possible brain links? how does it look?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make the matrix. 
for i = 1:53

    % ID1 = randi(18494, 1);
    % ID2 = randi(18494, 1);
    %    ds = i
    
    s = dataSetProbeInf(ds).sampleInd(1);
    e = dataSetProbeInf(ds).sampleInd(2);
    expMat =  wholeGeneExpr(:, s:e);
    
    corrMat = corr(expMat');
    upOnes = corrMat(logical(triu(ones(size(corrMat)), 1)));
    clear corrMat
    
    % get a chunk of matrix that you want 

    YMat(i, :) = corrQuantileReturn(dataSetProbeInf(ds).name, ...
       upOnes(1:1000))/1000;

    expMat = wholeGeneExpr(:, s:e);
    avgExp = mean(expMat, 2);
    avgExpC = avgExp - mean(avgExp);

    avgExpC = avgExpC + abs(min(avgExpC));

    featureMat(i, 1) = avgExpC(ID1) + avgExpC(ID2);
    featureMat(i, 2) = abs(avgExpC(ID1) - avgExpC(ID2));

    % featureMat(i, 1) = avgExpC(ID1);
    % featureMat(i, 2) = avgExpC(ID2);

    % lexp = min(avgExpC(ID1), avgExpC(ID2));
    
     tissue = dataSetProbeInf(ds).tissue;
     kado = unique({dataSetProbeInf(:).tissue});
     [a, b] = ismember(tissue, kado);
     
     %     b + 2;
     featureMat(i, b + 2) = 1;
        featureMat(i, 1) = min(avgExpC(ID1), avgExpC(ID2));
        %    featureMat(i) = min(avgExpC(ID1), avgExpC(ID2));
end


% The regression with the finalized tissue links. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

clear sumNet
tissue = 'brain'
t = 2
indProb = 0.1
load(sprintf(['~/networks/tissues/%s/%s' ...
              'SumNetworkExpThr0.8_%.*f.mat'], tissue, ...
             tissue, t,indProb))


smNet = sumNet >= 6;
blNet = sumNet >= 7;
brNet = sumNet >= 8;
lnNet = sumNet >= 9;
lvNet = sumNet >= 7;

sib = smNet + blNet + brNet + lnNet + lvNet; 

[ID1s, ID2s, corrs] = find(sib); 
temp = sib > 0;

load('~/data/general/tissueExpGenes/allExp0.8.mat')

% make the matrix.

% Get the genes which are expressed at least in one tissue: get
% some links

myGenes = zeros(1, 18494);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    myGenes = myGenes + expGenesInd;
end
expGenes = myGenes > 0;


inds = find(allExp);

inds = find(expGenes);
ID1s = datasample(inds, 10000);
ID2s = datasample(inds, 10000);
holu = ID1s - ID2s;
sum(holu == 0)

sib = zeros(18494, 18494);

ID1s = kado(:, 1);
ID2s = kado(:, 2);
for i = 1:length(ID1s)
    sib(ID1s(i), ID2s(i)) = 1;
end

% number of links

% 3. Getting regression for all the links
% >>>>>>>>>>>>>>>> Getting the links with some variety for corr
clear
load('~/data/general/linkExprInfo/wholeGeneExpr.mat')
load('~/data/general/linkExprInfo/dataSetProbeInf.mat')

load('~/networks/linksPresenetIn_2Top0.005_OR_Tissues.mat');
load('~/resultsAndFigures/TSlinks/TSlinksExtended.mat') % 'TSlinks'

load(['~/resultsAndFigures/TSlinks/' ...
      'TSlinks_FDR0012Union33Inter66.mat'])

fdr = '0012'
load(['~/resultsAndFigures/TSlinks/TSlinks_FDR' fdr '.mat'])

fdr = '001'
load(['~/resultsAndFigures/TSlinks/TSlinks_FDR' fdr '_NA_EA_withNextTissueValues.mat'])

load('~/resultsAndFigures/TSlinks/tissues.mat') % 'TSlinks'
tID = 4
tissues
sib = sparse(TSlinks{tID}.a', TSlinks{tID}.b',  ones(1, ...
                                                  length(TSlinks{tID}.a) ...
                                                  ), 18494, 18494);
TSlinks{tID}

%[a, b, c] = find(testNet);
n = length(TSlinks{tID}.a)
% the model matrix defining
featureMat = zeros(53, 2, n);
YMat = zeros(53, 1, n);

% sib = sparse([a(10001:12000)', 18494], [b(10001:12000)', 18494], c(1:n));
%sib = sparse([kado(:,1)', gCount], [kado(:,2)', gCount], c(1:n));

for i = 1:53
    i
    % ID1 = randi(18494, 1);
    % ID2 = randi(18494, 1);
    ds = i
    
    % getting the correlation for the dataset
    s = dataSetProbeInf(ds).sampleInd(1);
    e = dataSetProbeInf(ds).sampleInd(2);
    expMat =  wholeGeneExpr(:, s:e);
    
    tic
    corrMat = corr(expMat');
    corrMat = corrMat + .00001; 
    sibCorr = sib .* corrMat;
    [ID1s, ID2s, corrs] = find(sibCorr);
    toc
    
    % clear corrMat
    % % get correlation quantiles for the chunk of matrix that you want 
    % YMat(i, 1, :) = corrQuantileReturn(dataSetProbeInf(ds).name, ...
    %    corrs(1:n))/1000;
    
    YMat(i, 1, :) = corrs;

    expMat = wholeGeneExpr(:, s:e);
    
    % % >>>>>>>> I am going to use zscore instead of these
    % avgExp = mean(expMat, 2);
    % avgExpC = avgExp - mean(avgExp);
    % avgExpC = avgExpC + abs(min(avgExpC));
    % <<<<<<<<
    
    avgExp = mean(expMat, 2);
    avgExpC = zscore(avgExp);
    featureMat(i, 1, :) = avgExpC(ID1s(1:n)) + avgExpC(ID2s(1:n));
    featureMat(i, 2, :) = abs(avgExpC(ID1s(1:n)) - avgExpC(ID2s(1:n)));
    
    % featureMat(i, 1) = avgExpC(ID1);
    % featureMat(i, 2) = avgExpC(ID2);

    % lexp = min(avgExpC(ID1), avgExpC(ID2));
    
    % tissue = dataSetProbeInf(ds).tissue;
    % kado = unique({dataSetProbeInf(:).tissue});
    % [a, b] = ismember(tissue, kado);
    
    % %     b + 2;
    % featureMat(i, b + 2, c) = 1;
    %    featureMat(i, 1) = min(avgExpC(ID1), avgExpC(ID2));
    %    %    featureMat(i) = min(avgExpC(ID1), avgExpC(ID2));
end

% featureMat(1:9, 3, :) = 1;
% featureMat(10:24, 4, :) = 1;
% featureMat(25:31, 5, :) = 1;
% featureMat(32:41, 6, :) = 1;
% featureMat(42:53, 7, :) = 1;

% I made the feature matrix for 100 edges. 

% i = 6
% kado = featureMat(:, : ,i);

% r = datasample([1:53], 53, 'Replace', false);
% temp = kado(r, 3:end);
% kado(:, 3:end) = temp;

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> just testing
i = 3
tic
kado = featureMat(:, : ,i);
lm2 = LinearModel.fit(kado(:, 1:2), YMat(:, 1, i));
toc
Lm2.Rsquared

% x1int1 = lm2.Coefficients(1,1);
% x1int2 = dataset2struct(x1int1);
% inter = x1int2.Estimate;
% x1int1 = lm2.Coefficients(2,1);
% x1int2 = dataset2struct(x1int1);
% coef1 = x1int2.Estimate;
% x2int1 = lm2.Coefficients(3,1);
% x2int2 = dataset2struct(x2int1);
% coef2 = x2int2.Estimate;

tres = dataset2struct(lm2.Residuals);
Res = [tres(:).Raw];
plot(lm2.Fitted, Res, '*')
sib = find(Res > .4)
figure
qqplot([tres(:).Standardized])


(YMat(1,1,1) - (kado(1, 1) * coef1) - kado(1,2) * coef2 - inter)

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< just testing

% THIS
shuffle = true
n = 2000

% OR THIS
shuffle = false
RSEd = zeros(1, n);
SSE = zeros(1, n);
SSR  = zeros(1, n);
SST  = zeros(1, n);
coef1 = zeros(1, n);
coef2 = zeros(1, n);

% NOTE: this loop is modified fro the shuffled data, fix it. 
%ris = randi([1, n], 1, 10000);
for i = 1:n
    ri = i;
    kado = featureMat(:, : ,ri);
    Y = (YMat(:,1, ri));
    if(shuffle)
        shuffInd = datasample([1:53], 53, 'Replace', false);
        YS = Y(shuffInd);
        Y = YS;
        'sh'
    end
    lm2 = LinearModel.fit(kado(:, 1:2), Y);
    RSEd(i) = lm2.Rsquared.Ordinary;
    SSE(i) = lm2.SSE;
    SSR(i) = lm2.SSR;
    SST(i) = lm2.SST;
    x1int1 = lm2.Coefficients(2,1);
    x1int2 = dataset2struct(x1int1);
    coef1(i) = x1int2.Estimate;
    x2int1 = lm2.Coefficients(3,1);
    x2int2 = dataset2struct(x2int1);
    coef2(i) = x2int2.Estimate;
end

tissues(tID)
TSReg.RSEd = RSEd;
TSReg.SSE = SSE;
TSReg.SSR = SSR;
TSReg.SST = SST;
TSReg.coef1 = coef1;
TSReg.coef2 = coef2;

tissues{tID}
save(sprintf(['~/resultsAndFigures/expAndCorrRegression/' ...
      '%sTSReg_UnionFDR_' fdr '_NA_EA.mat'], tissues{tID}), ...
     'TSReg')

save(sprintf(['~/resultsAndFigures/expAndCorrRegression/' ...
      '%sTSReg_UnionFDR_' fdr '_shuffled_NA_EA.mat'], tissues{tID}), 'TSReg')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('~/resultsAndFigures/TSlinks/tissues.mat') % 'TSlinks'
load(sprintf(['~/resultsAndFigures/expAndCorrRegression/' ...
      '%sTSReg_UnionFDR001.mat'], tissues{tID}))



save('~/resultsAndFigures/expAndCorrRegression/brain_10000_SH_TSReg_extended.mat', ...
     'TSReg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
load('~/resultsAndFigures/TSlinks/TSlinks.mat')

load('~/resultsAndFigures/expAndCorrRegression/brainTSReg.mat')
load('~/resultsAndFigures/expAndCorrRegression/bloodTSReg.mat')
load('~/resultsAndFigures/expAndCorrRegression/liverTSReg.mat')
load('~/resultsAndFigures/expAndCorrRegression/lungTSReg.mat')
load('~/resultsAndFigures/expAndCorrRegression/muscleTSReg.mat')

holu = brainTSReg.RSEd > .24;
h = scatter(brainTSReg.coef2(holu), brainTSReg.RSEd(holu), 'k')
grid on
xlabel('sum')
ylabel('diff')

figure
h = scatter(brainTSReg.coef2, brainTSReg.RSEd, 'k')
heatscatter(brainTSReg.coef2(1:1000), brainTSReg.RSEd(1:1000), '~/', 'heatS.pdf')

h = figure
plot(brainTSReg.RSEd(sam), brainTSReg.coef1(sam), 'o', 'color', [200, ...
                   200, 200]/256)

holu = brainTSReg.coef2< -.8;
as = TSlinks{2}.a(holu);
bs = TSlinks{2}.b(holu);
i = 1
[as(i), bs(i)]

h = lsline
set(h(1), 'color', 'r')

sam = datasample([1:length(brainTSReg.RSEd)], 10000, 'Replace', false);
f = fittype('poly22')
myFit = fit(brainTSReg.coef2(sam)', brainTSReg.coef1(sam)', 'poly2')
plot(myFit, brainTSReg.coef2(sam), brainTSReg.coef1(sam), 'o')

a = hist(bloodTSReg.RSEd, 40)/length(bloodTSReg.RSEd);
bar(a)

colors = [228,26,28;
          55,126,184;
          77,175,74;
          152,78,163;
          10, 10,10]/256

h = figure
[f, xi] = ksdensity(bloodTSReg.RSEd);
plot(xi, f, 'color', colors(1,:), 'linewidth', 2)
hold all
[f, xi] = ksdensity(liverTSReg.RSEd);
plot(xi, f, 'color', colors(2,:), 'linewidth', 2)
[f, xi] = ksdensity(lungTSReg.RSEd);
plot(xi, f, 'color', colors(3,:), 'linewidth', 2)
[f, xi] = ksdensity(muscleTSReg.RSEd);
plot(xi, f ,'color', colors(4,:), 'linewidth', 2)
[f, xi] = ksdensity(brainTSReg.RSEd);
plot(xi, f, 'color', colors(5,:), 'linewidth', 2)

legend('blood - 24,058 links', 'liver - 14,890 links', ['lung - 19,' ...
                    '094 links'], 'muscle - 2,727 links', ['brain - ' ...
       '209,240 links'])
grid on

title(['Density of the R-squared values from regression model for the ' ...
       'TS links'])

figFolder = '~/resultsAndFigures/expAndCorrRegression/'
fileName = 'RsquaredDensity'
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf'])

% >>>>>>>>>>>>>>>>>>>>> A plot I made last minute
h = figure
[f, xi] = ksdensity(mytest);
plot(xi, f ,'color', colors(4,:), 'linewidth', 2)
hold on
[f, xi] = ksdensity(RSEd(1:45000));
plot(xi, f, 'color', colors(5,:), 'linewidth', 2)
legend('true labels', 'shuffled labels')
ylim([0 22])
grid on

title('Brain RSquared vs shuffled labels')
figFolder = '~/resultsAndFigures/expAndCorrRegression/'
fileName = 'BrainShuff'
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf'])
% <<<<<<<<<<<<<<<<<<<<<<<<<< 

length(muscleTSReg.RSEd)


load('~/resultsAndFigures/TSlinks/TSlinks.mat') % 'TSlinks'

lowRSQ = find(brainTSReg.RS



Ed < 0.0001);
length(lowRSQ)
l = 14
[TSlinks{2}.a(lowRSQ(1:l)) TSlinks{2}.b(lowRSQ(1:l))]

h = figure
hist(liverTSReg.RSEd, 40)


holu = find(bloodTSReg.RSEd < 0.1);
i = holu(2)
[TSlinks{tID}.a(i) TSlinks{tID}.b(i)]

% checking the SSR for the low RSEd
plot(SSE, RSEd, '.')
all.RSEd = RSEd;
all.SSE = SSE;
all.SSR = SSR;
all.SST = SST;

% Sometimes, I tend to have lower RSEd in my data. Small RSEd means
% that the SSE/SST is big, means that SSE is actually very similar
% to SST but that similarity might just be my expected noise, which
% I don't care about; I am not trying to predict
% that anyway. Now, if for small RSEd I also have a
% relatively small SSR (let's call it the noise SSR), it shows that RSEd is not a good
% represenation of my performance, since I am making a good
% prediction (considering the expected noise), but my data is
% actually very much mean-centered anyway. 
% how do I model the noise then? 

h = figure
plot(all.SSE, all.RSEd, '.')

sum(SSE < 2)

% test for the error:

inter1 = lm2.Coefficients(1,1);
inter2 = dataset2struct(inter1);
intercept = inter2.Estimate

x1int1 = lm2.Coefficients(2,1);
x1int2 = dataset2struct(x1int1);
x1 = x1int2.Estimate

x2int1 = lm2.Coefficients(3,1);
x2int2 = dataset2struct(x2int1);
x2 = x2int2.Estimate

p = x1 .* kado(:, 1) + x2 .* kado(:, 2) + intercept;
predictedErr = p - YMat(:, 1, 1);
SSres = sum(predictedErr .^ 2);
SStot = sum((YMat(:, 1, 1) - mean(YMat(:, 1 ,1))) .^ 2)
1 - (SSres/SStot)

h = figure; 
hist(RSor, 40)

smallError = RSor < 0.02;
sum(smallError)
smE = find(smallError);

RSnull = zeros(1, 10000);
for i = 1:10000
    nullError = 0.0194 - YMat(:, 1, i);
    RSnull(i) = sqrt(sum(nullError .^ 2))/53;
end
h = figure
hist(RSnull, 40)

fileName = sprintf('~/data/general/tempFigures/regressionError')
print(h, '-depsc', [fileName '.eps']);
print(h, '-dpdf', [fileName '.pdf']);

sib(ID1s(i), ID2s(i))

% I will need to get the variables for the two models, it is not
% really much. 

symbols = gpl570.uniqueSymbols(hkgInd);
% About 10million edges that we studied randomly... 

file = '~/data/HKG.txt'
fid = fopen(file, 'w')
for j = 1:length(symbols)
    fprintf(fid, ['%s\n'], symbols{j});
end
fclose(fid)

load('~/resultsAndFigures/TSlinks/TSlinks.mat') % 'TSlinks'
i = 10935
kado = [TSlinks{2}.a(10001:12000), TSlinks{2}.b(10001:12000)];

RSorOfAll = RSor;

h = figure;
hist(RSorOfAll, 40);

save('~/resultsAndFigures/expAndCorrRegression/Rsqred_10000OfAll.mat', 'RSorOfAll')
save('~/resultsAndFigures/expAndCorrRegression/Rsqred_2000BrainTS_Starts10001.mat', 'RSor')
save(['~/resultsAndFigures/expAndCorrRegression/' ...
      'the4SSsAll_10000_12000.mat'], 'all')
save(['~/resultsAndFigures/expAndCorrRegression/' ...
      'the4SSsBrain_10000_12000.mat'], 'brain')

(RSor < .1)

h = figure
hist(brain.SSR, 40)

h = figure
plot(all.SST, all.RSEd, '.')
title('all')
hist(all.RSEd(highVarAllInd), 10)

highVarAllInd = all.SST>4;

h = figure
plot(brain.SST, brain.RSEd, '.')
title('brain')
pureBrainInd = find(brain.RSEd<.1);
sum(pureBrainInd)

t = 887
kado = [TSlinks{2}.a(10000 + t), TSlinks{2}.b(10000 + t)]



