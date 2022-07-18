% This file is to load the data from Gemma, get the control samples
% and put them together in one expression matrix. Starting is a
% list of experiments for a tissue with the tissue name. for each experiment,
% file is loaded in the datafolder, the columns needed are
% extracted
% Having the experiment file:
% 1. getting the files from Gemma
% 2. getting the platform file 
% 3. getting normal samples
% 4. sanity check and sample correlation heatmap plots
% 5. removing the bad samples
% 6. making the network 
%NOTE: if you wanted to add the name of biassays in the
%finalMat.mat file, save finalMat in a text file and add the
%variable normalName as the heading. 
%NOTE : there is a bug making the network, I should get the top .5%
%from the triangle ranking matrix, not the whole matrix. - fixed in arrayMain02.m

%global variables
dataFolder = '/space/grp/marjan/data';
global rawDataFolder;
global filDataFolder;
global designFolder;
global tissue;

rawDataFolder = '/space/grp/marjan/data/blood/rawData/';
filDataFolder = '/space/grp/marjan/data/blood/filteredData/';
designFolder = '/space/grp/marjan/data/blood/design/';
exprFolder = '/space/grp/marjan/data/blood/normalGeneExpression/';

try
createClassFromWsdl('http://www.chibi.ubc.ca/Gemma/ws/gemma.wsdl')
catch 
end 

% Step 0. loading the experiment file and info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loading the text file created for the tissue. 
tissue = 'blood'
%fileName = '/space/grp/marjan/data/blood/theBloodDataInfo.txt';
fileName = '/space/grp/marjan/data/blood/theBloodDataInfoClean.txt';
fid = fopen(fileName);

expInf = textscan(fid, ['%s', '%d', '%s', '%d', '%d'], 'headerlines', 1, ...
                  'Delimiter', '\t')
fclose(fid)
experimentList = expInf{1};
experimentCount = length(experimentList);
infCol = expInf{2};
infStr = expInf{3};
infColCount = expInf{4};
hlineCount = expInf{5};
%these words are manually curated for each data set from the design
%matrix, so that I can pick the control/normal samples. 

%Step 1. getting files from GEMMA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:experimentCount
    Name = experimentList{i};

       'raw'
    %getting the rawfile
    system(['gemmaCli.sh ubic.gemma.apps.ExpressionDataMatrixWriterCLI ' ...
            '-e ' Name '  -o ' rawDataFolder Name '.txt']);
    
    % 'filtered' 
    % %getting the filteredData
    % system(['gemmaCli.sh ubic.gemma.apps.ExpressionDataMatrixWriterCLI ' ...
    %         '-e ' Name '  -o ' filDataFolder Name 'fil.txt -filter']);
    
    'design'
    %getting the design file
    system(['gemmaCli.sh ' ...
            'ubic.gemma.datastructure.matrix.ExperimentalDesignWriterCLI ' ...
            '-e ' Name ' -o ' designFolder Name 'design.txt']);
    
end

% Step 2. reading the gpl file
%%%%%%%%%% reading gpl file

fid = fopen('/space/grp/marjan/data/GPL570-13270.txt')
gpl = textscan(fid, [repmat('%s', 1, 16)], 'Delimiter', '\t', ...
               'Headerlines',17, ...
               'Bufsize', 100000095);
%this buffersize just seem to be enough!

%GPL570 probe ID's
gpl570.probeID = gpl{1};

%GPL570 gene ID's
gpl570.geneID = gpl{11};

%GPL570 gene symbols
gpl570.symbols = unique(gpl{11});
gCount = length(gpl570.symbols);
 
geneSymbols = gpl570.symbols;
save(sprintf('%sGPL570geneSymbols.mat', exprFolder), ...
     'geneSymbols');

%loading the geneSymbols
load([exprFolder 'GPL570geneSymbols.mat'])

% Step 3.reading the files 
% having the list of files in the data location, FOR EACH FILE: 
% 0. read the file. 
% .5 get the normal samples. 
% 1. convert the probe matrix to the gene matrix
% 2. save the gene matrix
%%%%%%%%%%%%%%%%% getting the gene expression matrix 
% get the raw gene expressions and then the network from there. 

%the section below is for getting the expression matrix. 

for i = 1:experimentCount
     
    Name = experimentList{i}
    
    %3.1 getting normal samples
    %getting the number of samples for the file
    tempExpNum = experimentId(gemmaService, Name);
    expNum = tempExpNum.ee_id;
    tempSampCount = experimentNumSamples(gemmaService, expNum);
    sampCount = tempSampCount.eeNumSample_id;
    sampCount = str2num(sampCount);    
    
    linum = 7;
    fileName = [rawDataFolder Name '.txt'];
    fid = fopen(fileName);
    
    header = textscan(fid, '%s', 1, 'delimiter', '\n', 'headerlines', ...
                         linum - 1, 'Bufsize', 100000095);
    
    data = textscan(fid, [repmat('%s', 1,2) repmat('%f', 1, sampCount)], ...
                    'Headerlines', 7, 'Delimiter', '\t', 'Bufsize', ...
                    100000095);
    fclose(fid);
    
    [normalName, normalInd] = getNormalSamplesIndex(Name, infCol(i), ...
                                                          infStr(i), ...
                                                          infColCount(i), ...
                                                          sampCount, ...
                                                          header, hlineCount(i), ...
                                                          'raw');
    
     normalCount = length(normalInd);
     
     %getting the probe expression matrix
     pCount = length(data{1});        
     dataExpr = cell2mat(data(normalInd + 2));

     %%%%%%%%% mapping the probes to the genes. 
     %%% NOTE : in each expression matrix, genes are as they are in
     %%% the gpl570.symbols.

     [a, b] = ismember(data{1}, gpl{1});
     geneID = gpl{11}(b(a));

     %3.2 converting the probe matrix to the gene matrix
     geneMap = containers.Map(gpl570.symbols, [1:gCount]);
     probeMap = containers.Map(data{1}, geneID);
     finalMat = zeros(gCount, normalCount);
     divMat = zeros(gCount, normalCount);

     tic
     for i = 1:pCount
         sib = (values(geneMap, values(probeMap, data{1}(i))));
         finalMat(sib{1}, :) = finalMat(sib{1}, :) + dataExpr(i, :);
         divMat(sib{1}, :) = divMat(sib{1}, :) + 1;
     end
     toc
     finalMat = finalMat./divMat;

     exprIDs = zeros(1, normalCount);

     
     %getting gemma experiment IDs for the samples
    for j = 1:normalCount
        holu = normalName{j}; 
        tempb = strfind(holu, 'Id=');
        tempe = strfind(holu, 'Name');
        exprIDs(j) = str2num(holu(tempb+3:tempe - 1));
    end
    
    dataSet.mat = finalMat;
    dataSet.IDs = exprIDs;
    dataSet.names = normalName;

     save(sprintf('%sgeneExprDataSet%s.mat', ...
                  exprFolder, Name), 'dataSet');
     
end

% 4. getting the sanity check heatmap plots. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% making the heatmap sample correlation plots
for i = 1:experimentCount
    
    Name = experimentList{i}
 
    load(sprintf('%scleangeneExprDataSet%s.mat', exprFolder, Name));
    
    'fileloaded'

    exprData = dataSet.mat;
    sampCount = size(exprData, 2);
    sib = corr(exprData);
    
    'corr'

    % make and save heatmap
    h = figure; 
    heatmap(sib, dataSet.IDs, dataSet.IDs,  [], 'TickAngle', 45, 'ShowAllTicks', ...
            true, 'Colorbar', true)
    'heatmap done'
    title(sprintf('%s  Number of samples : %d', Name, sampCount))
     fileName = ['cleanRawSampleCorr' Name];
     figFolder = '/space/grp/marjan/data/blood/figures/';
    print(h, '-depsc', [figFolder fileName '.eps']);
    print(h, '-dpdf', [figFolder fileName '.pdf']);

end

% 5. removing the bad samples, based on the heatmap plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%based on the heatmap plots, removing the bad samples. 
% I am picking the samples manually here, maybe later I automate 
i = 2;
viewOnly = 0;
Name = experimentList{i}
samples = [175250];
fileName = ['geneExprDataSet' Name '.mat'];
removeSamples(exprFolder, fileName, samples, viewOnly);


b = 50
e = 80

    load(sprintf('%sgeneExprDataSet%s.mat', exprFolder, Name));
    
    'fileloaded'

    exprData = dataset.mat;
    sampCount = size(exprData, 2);
    sib = corr(exprData(:, b:e));
    
    'corr'

    % make and save heatmap
    h = figure; 
    heatmap(sib, dataset.IDs(b:e), dataset.IDs(b:e),  [], 'TickAngle', 45, 'ShowAllTicks', ...
            true, 'Colorbar', true)
    
    %5.5 getting the over all sample rank correlation 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% putting togehter all the sample files, I want to get the
% correlation for all the samples together. 
 wholeBlood = zeros(gCount, 329);
 c = 1;
 label = repmat({''}, 1, 329);
  
for i = 1:experimentCount
    
    Name = experimentList{i};

    load(sprintf('%scleangeneExprDataSet%s.mat', exprFolder, Name));
    
    sampCount = size(dataSet.mat, 2);
    if(max(max(dataSet.mat)) > 1000)
     wholeBlood(:, c:(c+sampCount -1)) = log2(dataSet.mat);
    else
     wholeBlood(:, c:(c+sampCount -1)) = dataSet.mat;
   end
   
    label(c) = {Name};
    max(max(dataSet.mat))
    c = c + sampCount;
    
end

sib = corr(wholeBlood);

 h = figure; 
    heatmap(sib, label, label,  [], 'TickAngle', 45, 'ShowAllTicks', ...
            true, 'Colorbar', true, 'Colormap', 'hot')
    'heatmap done'
    title('blood01- Blood datasets, total of 329 samples')
     fileName = ['blood01AllSamplesCorr' ];
     figFolder = '/space/grp/marjan/data/blood/figures/';
    print(h, '-depsc', [figFolder fileName '.eps']);
    print(h, '-dpdf', [figFolder fileName '.pdf']);
    
  
    % 6. getting the correlation matrix and making the network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gCount = 21050;

%the one network memory allocation
finalCorrMat = zeros(gCount);

%the two network memory allocation
finalPosCorrMat = zeros(gCount);  %only for neg and pos network
finalNegCorrMat = zeros(gCount); % only for neg and pos network


networkFolder = '/space/grp/marjan/data/blood/network/' % to save
                                                        % the network in
exprFolder = '/space/grp/marjan/data/blood/normalGeneExpression/' % to read the expression from
selectionName = 'blood01'


%define gCount
for i =1: experimentCount
    i
    tic
    
    Name = experimentList{i};
    load(sprintf('%scleangeneExprDataSet%s.mat', exprFolder, Name));
    
    [ ' 1.expression file loaded for ' Name]

    exprCount = size(dataSet.mat, 2);

    corrMat = corr(dataSet.mat', 'type', 'Spearman');
    tempMat = -diag(ones(1, gCount), 0); % putting it zero so it is
                                         % not picked in the ranking
    corrMat = corrMat + tempMat;
    size(corrMat)
    
    clear tempMat 
    clear corrExpr

       [ ' 2. coexpression made for ' Name]
    % I am not saving the coexpression matrix. NO USE. 

    %ranking the whole matrix. 
    tr = reshape(corrMat, 1, gCount* gCount) ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %single network
    rankedTR = tiedrank(tr);
    clear tr
    
    tempRankCorrMat = reshape(rankedTR, gCount, gCount);
    finalCorrMat = finalCorrMat + tempRankCorrMat;
    clear tempRankCorrMat
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neg and pos networks
    
    % if you decided to keep pos and negative separately, uncomment
    % the folloing 
    % Negtr = tr;
    % Negtr = corrMat;
    % Negtr(Negtr > 0) = 0; %just keeping the negative correlation
    
    % Postr = tr;
    % Postr = corrMat;
    % Postr(Postr < 0) = 0; %just keeping the positive correlation
    
    % %each tied rank takes 10m. 
    % [ ' 3. starting the ranking for ' Name]
    % ptr = tiedrank(Postr);  
    % clear Postr
    
    % [ ' 4. postivie rank finished']
    
    % ntr = tiedrank(Negtr);
    % clear Negtr
    
    % [ '5. negative rank finished']
    
    % tempPosCorrMat = reshape(ptr, gCount, gCount);
    % tempNegCorrMat = reshape(ntr, gCount, gCount);
    
    % clear ntr
    % clear ptr

    % [ '6. adding to the final mat']
    
    % finalPosCorrMat = finalPosCorrMat + tempPosCorrMat;
    % finalNegCorrMat = finalNegCorrMat + tempNegCorrMat;
    
    % % finalPosCorrMat = finalPosCorrMat + Postr;
    % % finalNegCorrMat = finalNegCorrMat + Negtr;

    % clear tempPosCorrMat
    % clear tempNegCorrMat
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    whos
    toc
    %End of things I want todo for each file
end

%one network 
%*just code here!* % save both the rawRankNetwork and the sparsePosNetwork

%saving raw rank
save([networkFolder selectionName 'RawRank.mat'], 'finalCorrMat', '-v7.3')

%getting 0.5% top ranks for the network 
book = reshape(finalCorrMat, 1, gCount * gCount);
Qlist = quantile(book, 200);
Q = Qlist(200);
posNet = finalCorrMat;
posNet(posNet < Q) = 0;
posNet(posNet >= Q) = 1;
posNet = sparse(posNet);
save([networkFolder selectionName 'sparsePosNetwork.mat'], ...
     'posNet');

Q = Qlist(1);
negNet = finalCorrMat;
negNet(negNet <= Q) = 1;
negNet(negNet > Q) = 0;
negNet = sparse(negNet);
save([networkFolder selectionName 'sparseNegNetwork.mat'], 'negNet');


    % save([networkFolder selectionName 'RawRankNeg.mat'], ...
    %      'finalNegCorrMat', '-v7.3');
    % save([networkFolder selectionName 'RawRankPos.mat'], ...
    %      'finalPosCorrMat', '-v7.3');
 
    % %getting the 0.05% top ranks. 
    % ptr = reshape(finalPosCorrMat, 1, gCount*gCount);
    % Qlist = quantile(ptr, 2000);
    % Q = Qlist(2000);

    % finalPosCorrMat(finalPosCorrMat <= Q) = 0;
    % finalPosCorrMat(finalPosCorrMat > Q) = 1;

    % ntr = reshape(finalNegCorrMat, 1, gCount*gCount);
    % Qlist = quantile(ntr, 2000);
    % Q = Qlist(1);


    % finalNegCorrMat(finalNegCorrMat < Q) = 1;
    % finalNegCorrMat(finalNegCorrMat >= Q) = 0;
    
    % %saving the network
    % sparseNegCoexprMat = sparse(finalNegCorrMat);
    % sparsePosCoexprMat = sparse(finalPosCorrMat);
    
    % save([networkFolder selectionName 'FinalNegNetwork.mat'], 'sparseNegCoexprMat')
    % save([networkFolder selectionName 'FinalPosNetwork.mat'],
    % 'sparsePosCoexprMat')
    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %loading the network and sanity checks
    
   NegCorrNet = full(sparseNegCoexprMat);
   PosCorrNet = full(sparsePosCoexprMat);
   
   posNodeDegree = sum(PosCorrNet);
   hist(sum(PosCorrNet), 40) 
   over250NodeDegreeInd = find(posNodeDegree > 250);
   
   negNodeDegree = sum(NegCorrNet);
   hist(sum(NegCorrNet), 40)
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %analysis of negative and positive correlations among the gens,
   %for finalPosCorrMat and final NegCorrMat
   
   negCorrVector = reshape(finalNegCorrMat, 1, gCount*gCount);
   hist(negCorrVector, 40)% 0.019% are smaller than -3 
                          %  smaller than -4      vs random 0.6%

   negQlist = quantile(negCorrVector, 2000); %Qlist(1) = -3.31  
   negCoExpr = finalNegCorrMat;
   negCoExpr(negCoExpr > -3.31) = 0;
   negCoExpr(negCoExpr < -3.31) = 1;
   
   
   posCorrVector = reshape(finalPosCorrMat, 1, gCount*gCount);
   hist(posCorrVector , 40)% 0.1% are bigger than 3
                           % 8.63e-04 greater than 4 vs random 0.6%
                           % what is the 2000 (0.0005 of data) quantile?
   posQlist = quantile(posCorrVector, 2000); % Qlist(2000) = 4.23
   posCoExpr = finalPosCorrMat;
   posCoExpr(posCoExpr < 4.23) = 0;
   posCoExpr(posCoExpr > 4.23) = 1;
   
   
   mixVec = posCorrVector + negCorrVector;
   mixQlist = quantile(mixVec, 2000); %Qlist(2000) = 4.2030,
                                      %Qlist(1) = -3.2743;
   mixNegCoExpr = finalPosCorrMat + finalNegCorrMat;   
   mixNegCoExpr(mixNegCoExpr > -3.2743) = 0;
   mixNegCoExpr(mixNegCoExpr < -3.2743) = 1;

   
gCount = 21050;
finalPosCorrMat = zeros(gCount);  
finalNegCorrMat = zeros(gCount);
networkFolder = '/space/grp/marjan/data/blood/network/' % to save
                                                        % the network in
exprFolder = '/space/grp/marjan/data/blood/normalGeneExpression/' % to read the expression from
selectionName = 'blood01'

   %trying a random sample:
   for i = 1 : experimentCount
   
       %1. create a sample of random values
       Name = experimentList(i)
    load(sprintf('%scleangeneExprDataSet%s.mat', exprFolder, Name));
    
    [ ' 1.expression file loaded for ' Name]

    exprCount = size(dataSet.mat, 2);

    sample = rand(gCount, exprCount);
    
    
    corrMat = corr(sample', 'type', 'Spearman');
    tempMat = -diag(ones(1, gCount), 0);
    corrMat = corrMat + tempMat;
    size(corrMat)
    
    clear tempMat 
    clear corrExpr

    
    [ ' 2. coexpression made for ' Name]
    % I am not saving the coexpression matrix. NO USE. 

    %ranking the whole matrix. 
    % tr = reshape(corrMat, 1, gCount* gCount) ;
    % Negtr = tr;
    Negtr = corrMat;
    Negtr(Negtr > 0) = 0; %just keeping the negative correlation
    
    % Postr = tr;
    Postr = corrMat;
    Postr(Postr < 0) = 0; %just keeping the positive correlation
    
    
    % finalPosCorrMat = finalPosCorrMat + tempPosCorrMat;
    % finalNegCorrMat = finalNegCorrMat + tempNegCorrMat;
    
    finalPosCorrMat = finalPosCorrMat + Postr;
    finalNegCorrMat = finalNegCorrMat + Negtr;
       
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % making the data ready for the finalizedMain.m (pp)

   %kernel
   load([networkFolder selectionName 'sparsePosNetwork.mat']);
   data.kernel = full(sparsePosNetwork);
   
   load([exprFolder 'GPL570geneSymbols.mat'])
   data.geneSymbols =  geneSymbols;
   
   load([dataFolder 'GOdata_GPL570.mat'])
   data.matF = GOdata.matF;
   data.matP = GOdata.matP;
   data.matC = GOdata.matC;
   data.GOTerms = GOdata.GOTerms;
   data.GOID = GOdata.GOID;
   
   % save([networkFolder selectionName 'RawRankNeg.mat'], ...
    %      'finalNegCorrMat', '-v7.3');
   
   %labelMatrix
   %network
   %function terms
   %function id's
   %gene symbols