% Here my goal is just to get the shared neighbour thing and the
% shared JCC thing. 

clear

tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
% get the genes expressed in more than one tissue. 
expMat = zeros(18494, 5);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    expMat(:, t) = expGenesInd;
end
 
teCount = sum(expMat');
% for those genes, get this whole exp profile. 

% making the matrix exp template for links which could be present
% in more than one tissue
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
% get the genes expressed in more than one tissue. 
expNetTemp = zeros(18494, 18494);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    expNetTemp(expGenesInd, expGenesInd) = expNetTemp(expGenesInd, expGenesInd) + 1;
end

expNetTempMTO = expNetTemp>1; % the link is expressed in more than
                              % one tissue


% get the function of these genes in each of the tissues. 
myGenes = find(teCount > 0);

load(['~/resultsAndFigures/commonAndPure/' ...
     'multiTissueExpGenesFunctions_GO07_nonComp_expInOneT.mat'])

FDR = '0012'
FC = '3'

load(sprintf(['~/resultsAndFigures/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))

% >>>> making an alternative template for the cgnet, where I have
% the in more than one tissue
expTemplate = zeros(18494, 18494);
moreThanOneExp = teCount > 1;
expTemplate(moreThanOneExp, moreThanOneExp) = 1;
% get all the genes in the pure net.
sumNet = zeros(18494, 18494);
for t = 1: 5
    
    % links with genes expressed in more than one tissue
    %tempPure = ((expNetTempMTO + finalTable(t).noTSNet) == 2) ...
    %        .* t;
     
    %links with genes expressed in all the tissuesn
     tempPure = ((finalTable(t).cgNet + finalTable(t).noTSNet) == 2) ...
          .* t;

    sumNet = sumNet + tempPure;
end
max(max(sumNet))
sum(sum(sumNet > 0))
book = sum(sumNet) + sum(sumNet');

% getting the ATND and TSND
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'};
atnds = zeros(5, 18494);
tsnds = zeros(5, 18494);
for t = 1 : 5
    tissue = tissues{t};
    load( ['~/networks/tissues/' tissue '/' ...
           'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
    atn{t} = binNet;
    atnds(t, :) = sum(binNet) + sum(binNet');
    tsnet = finalTable(t).wholeNet;
    tsnds(t, :) = sum(tsnet) + sum(tsnet');
end

% load pure links and get that neighbour count 
load('~/resultsAndFigures/tempResults/pureData_allThings.mat')

% for each pure link, get the overlap of the partners in TAN:
overlapRec = zeros(80637, 8);
[a, b, c] = find(sumNet);
overlapRec(:, 1) = a;
overlapRec(:, 2) = b;
overlapRec(:, 3) = c;
for t = 1:5
    thisTAN = atn{t};
    fullThisTAN = thisTAN + thisTAN';
    for i = 1:length(a)
        v1 = fullThisTAN(a(i), :);
        v2 = fullThisTAN(b(i), :);
        insct = v1 * v2';
        unn = sum(v1) + sum(v2) - insct;
        overlapRec(i, (t+3)) = insct / (unn + .01);
    end
end

save('~/resultsAndFigures/tempResults/pureData_neighborJCC.mat', ...
     'overlapRec')

pureTOP = zeros(1, 80637);
maxOTOP = zeros(1, 80637);
minOTOP = zeros(1, 80637);
for i = 1:80637
   pureT = overlapRec(i, 3);
   book = overlapRec(i, 4:end);
   pureTOP(i) = book(pureT);
   book(pureT) = [];
   maxOTOP(i) = max(book);
   minOTOP(i) = min(book);
end

top.pure = pureTOP;
top.max = maxOTOP;
top.min = minOTOP;
save('~/resultsAndFigures/tempResults/pureData_TOP_pureMaxMin.mat', ...
     'top')

h = figure
hist(pureTOP)

h = figure
hist(minOTOP)

hh = figure
hist(minOTOP)
plot([1:10])
% for each link, get the overlap of their neighbourhoods in each of
% the tissues for their TSNs and TANs. This is a loop, just code
% it. 

% The max is the max, find the minimum. box plot them. 

% load semi pure and get that neighbour count 
load('~/resultsAndFigures/tempResults/semiPureData_allThings.mat')

% same for here. 

% load common and get that neighbour count 
load('~/resultsAndFigures/tempResults/sparseNet_null.mat')

overlapRec = zeros(size(sparseNets));
overlapRec(:, 1) = sparseNets(:, 1);
overlapRec(:, 2) = sparseNets(:, 2);
for t = 1:5
    t
    thisTAN = atn{t};
    fullThisTAN = thisTAN + thisTAN';
    overlap = fullThisTAN * fullThisTAN;
    nds = sum(fullThisTAN);
    for i = 1:length(sparseNets)
        i
        g1 = sparseNets(i, 1);
        g2 = sparseNets(i, 2);
        insct = overlap(g1, g2);
        unn = nds(g1) + nds(g2) - insct;
        overlapRec(i, (t+2)) = insct / (unn + .01);
    end
end

save('~/resultsAndFigures/tempResults/nullData_neighborJCC.mat', ...
     'overlapRec')

mainTOP = zeros(1, 368411);
maxOTOP = zeros(1, 368411);
minOTOP = zeros(1, 368411);
for i = 1:368411
   book = overlapRec(i, 3:end);
   book = book(logical(sparseNets(i, 3:end)));
   [mainTOP(i), b] = max(book);
   book(b) = [];
   maxOTOP(i) = max(book);
   minOTOP(i) = min(book);
end

top.main = mainTOP;
top.max = maxOTOP;
top.min = minOTOP;
save('~/resultsAndFigures/tempResults/nullData_TOP_mainMaxMin_logical.mat', ...
     'top')

thisTan = atn{1};
fullThisTan = thisTan + thisTan';

overlap = fullThisTan * fullThisTan;

whos overlap

