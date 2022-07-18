# Bioinf_2019
Farahbod &amp; Pavlidis 2019 Bioinformatics 


These three functions are used for measuring the TS value for each of the links:

    continuousTSLinksFunction_NA_qnorms.m - meaasures the TS value for a given tissue, given list of links

    continuousTSLinksFunction_RandomRep_NA_qnorms.m - calculates the TS value for pseudo tissues

    corrQuantileReturn.m - returns the measure quantiles for given sets of correlation values from a dataaset


>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
LIST OF ALL FILES AND THEIR HEADERS IN THE REPOSITORY
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

arrayStart.m

% http://www.mathworks.com/products/demos/bioinfo/affydemo/affydemo.html
% This file has the steps from array experiment in GEO, to the matrix of raw probe/sample in a
% .txt file
% % step 1. found the desired experiment in GEO, load all the .CEL
% files in a folder. Also, load the platform .cdf file in the
% parent directory. Sample experiment : GSE33223

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

arrayExplorationMain.m

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

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

arrayMain01.m

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

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

arrayMain03.m

% This file is for getting the gene expression matrix from the cel
% files from GEO, sanity check and heatmap/boxplots of the
% samples. This code is Gemma free (unlike arrayMain02.m and arrayMain01.m)
% THE MUST RUN PART
% Step 0. open the dataset info file and read the names
% Step 1. getting the .CEL files from GEO
% Step 2. making the probe matrix from .cel files
% Step 3. loading the platform file
% Step 4. making the gene matrix
% Step 5. sanity check and correlation heatmap plots 
% Step 6. removing the outliers 
% Step 7. correlation of all samples together for one tissue. 
% Step 8. getting the correlation between list of tissues

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

arrayMain04.m

% same as arrayMain03, but using files :
% /space/grp/marjan/data/GPL570GemmaMap.txt and GPL570AffyMap.txt
% probes with multiple genes and probes with no genes are deleted
% from probe list in GSEthese two files, also, instead of gene symbol, I
% am using NCBI gene ID in each file. These two mappings share
% 18000 NCBI gene IDs
% This file is for getting the gene expression matrix from the cel
% files from GEO, sanity check and heatmap/boxplots of the
% samples. This code is Gemma free (unlike arrayMain02.m and arrayMain01.m)
% THE MUST RUN PART
% Step 0. open the dataset info file and read the names - obsolete
% Step 1. getting the .CEL files from GEO - obsolete 
% Step 2. making the probe matrix from .cel files 
% Step 3. loading the platform file
% Step 4. making the gene matrix 
% Step 5. sanity check and correlation heatmap plots 
% Step 6. removing the outliers 
% Step 7. correlation of all samples together for one tissue. 
% Step 8. getting the correlation between list of tissues

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

arrayMain04Sub.m

% I chopped this part of arrayMain04.m code for making the probe
% correlation density plots. I save the results for each tissue in
% a file named [tissue]_probeCorrDistribution.mat 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


arrayMain04Sub02.m

% get the sample count for each dataset for each tissue, save them
% in a file. 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

arrayMain04Sub04.m

% getting the 3 plots of correlation distribution. 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

binomialTestForAgg.m

% This is to find the Bernouli probability for observing the links
% in my networks. 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

CELtest.m

% this file is to get the batch info for experiments from their
% .CEL files

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

CELtoProbeTxt.m

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


clusterStudies.m

% in this file I study clusterings of the genes. How dense each
% cluster is and how connected the clusterings are with each
% other. My inputs are the clusterings and the networks (tissue
% networks). Questions I want to find the answer to are: 
% 1. How dense is the cluster? 
% 2. How is it "clustered"? ND distribution and cluster coefficient
% maybe? Does it have like, sub clusters or genes are smoothly
% clustered?
% 3. how is each cluster connected to the other clusters? 


>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


continuousTSLinksFunction.m

% this is the main function for running the code in parallel and
% saving the data. 
% where tissue is the tissue and s is the starting point for the
% correlation matrix and n is the number of edges to be checked
% (cause I am running in parallel, so I choose how long each run takes)
%function[] = continuousTSLinksFunction(tissue, s, n)
% s is the starting point (cause I get in limited number of edges
% each time) - n is the count of edges , I have ~57e+6 edges, for
% brain, for example. 1e+6 takes about half an hour to run. I give
% it 10mil at a time, I guess. 10e+6 takes 5 hours, and probably
% 10G of RAM. 


>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

continuousTSLinksFunction_NA_qnorm.m

% this is the main function for running the code in parallel and
% saving the data.
% where tissue is the tissue and s is the starting point for the
% correlation matrix and n is the number of edges to be checked
% (cause I am running in parallel, so I choose how long each run takes)
%function[] = continuousTSLinksFunction(tissue, s, n)
% s is the starting point (cause I get in limited number of edges
% each time) - n is the count of edges , I have ~57e+6 edges, for
% brain, for example. 1e+6 takes about half an hour to run. I give
% it 10mil at a time, I guess. 10e+6 takes 5 hours, and probably
% 10G of RAM. ~


>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

contiunousTSLinksFunction_negAdjusted_ExpAdjusted.m

% this is the main function for running the code in parallel and
% saving the data. 
% where tissue is the tissue and s is the starting point for the
% correlation matrix and n is the number of edges to be checked
% (cause I am running in parallel, so I choose how long each run takes)
%function[] = continuousTSLinksFunction(tissue, s, n)
% s is the starting point (cause I get in limited number of edges
% each time) - n is the count of edges , I have ~57e+6 edges, for
% brain, for example. 1e+6 takes about half an hour to run. I give
% it 10mil at a time, I guess. 10e+6 takes 5 hours, and probably
% 10G of RAM. 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

continuousTSLinksFunction_negativeAdjusted.m

% this is the main function for running the code in parallel and
% saving the data. 
% where tissue is the tissue and s is the starting point for the
% correlation matrix and n is the number of edges to be checked
% (cause I am running in parallel, so I choose how long each run takes)
%function[] = continuousTSLinksFunction(tissue, s, n)
% s is the starting point (cause I get in limited number of edges
% each time) - n is the count of edges , I have ~57e+6 edges, for
% brain, for example. 1e+6 takes about half an hour to run. I give
% it 10mil at a time, I guess. 10e+6 takes 5 hours, and probably
% 10G of RAM. 


>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

continuousTSLinksFunction_random.m

% this function is a modification of the function:
% continuousTSLinksFunction() 
% this is for random testing, so it takes a list of datasets as the
% "tissue" (which are supposed to be a mixture of all other
% tissues) and the rest of datasets are random. The dataSetList
% should be a list of dataSet IDs. 
% s is the starting point (cause I get in limited number of edges
% each time) - n is the count of edges , I have ~57e+6 edges, for
% brain, for example. 1e+6 takes about half an hour to run. I give
% it 10mil at a time, I guess. 10e+6 takes 5 hours, and probably
% 10G of RAM.

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

continuousTSLinksFunction_randomRep.m

% this function is a modification of the function:
% continuousTSLinksFunction() 
% this is for random testing, so it takes a list of datasets as the
% "tissue" (which are supposed to be a mixture of all other
% tissues) and the rest of datasets are random. The dataSetList
% should be a list of dataSet IDs. 
% s is the starting point (cause I get in limited number of edges
% each time) - n is the count of edges , I have ~57e+6 edges, for
% brain, for example. 1e+6 takes about half an hour to run. I give
% it 10mil at a time, I guess. 10e+6 takes 5 hours, and probably
% 10G of RAM. 


>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


continuousTSLinksFunction_randomRep_WRV.m

% this function is a modification of the function:
% continuousTSLinksFunction() 
% this is for random testing, so it takes a list of datasets as the
% "tissue" (which are supposed to be a mixture of all other
% tissues) and the rest of datasets are random. The dataSetList
% should be a list of dataSet IDs. 
% s is the starting point (cause I get in limited number of edges
% each time) - n is the count of edges , I have ~57e+6 edges, for
% brain, for example. 1e+6 takes about half an hour to run. I give
% it 10mil at a time, I guess. 10e+6 takes 5 hours, and probably
% 10G of RAM. 


>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

continuousTSLinksFunction_randomRep_NA_qnorms.m

% this is the main function for running the code in parallel and
% saving the data. 
% where tissue is the tissue and s is the starting point for the
% correlation matrix and n is the number of edges to be checked
% (cause I am running in parallel, so I choose how long each run takes)
%function[] = continuousTSLinksFunction(tissue, s, n)
% s is the starting point (cause I get in limited number of edges
% each time) - n is the count of edges , I have ~57e+6 edges, for
% brain, for example. 1e+6 takes about half an hour to run. I give
% it 10mil at a time, I guess. 10e+6 takes 5 hours, and probably
% 10G of RAM. 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

continuousTSLinksFunction_randomRep_WRV_oneSided.m

% this function is a modification of the function:
% continuousTSLinksFunction() 
% this is for random testing, so it takes a list of datasets as the
% "tissue" (which are supposed to be a mixture of all other
% tissues) and the rest of datasets are random. The dataSetList
% should be a list of dataSet IDs. 
% s is the starting point (cause I get in limited number of edges
% each time) - n is the count of edges , I have ~57e+6 edges, for
% brain, for example. 1e+6 takes about half an hour to run. I give
% it 10mil at a time, I guess. 10e+6 takes 5 hours, and probably
% 10G of RAM.

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

continuousTSLinksFunction_WilcoxonRankValue.m

% this is where I find the tissue specificity for each link with the wilcoxon
% Rank test. 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

correlationHomoTest.m

% this is basically similar to correlationVarFunc02 for the input
% variables and datastructure, meant to run on different chunks of
% the correlation matrix 

% this function gets the mean and variance of the correlation
% values between the given index of genes from the correlation
% matrix of the datasets for one tissue. I had to write a function
% for it, limiting the input matrix since loading all the datasets
% on the ram at a time was not possible and also I can do it in
% parallel using the function
%
% function[] = correlationVarFunc02(ib, ie, jb, je, tissue)
%
% where ib, ie, jb and je define the sub-matrix of the correlation
% matrix and tissue is a string, name of the tissue(which is name
% of the folder for tissue)
% Notes:
% 1. this function is faster than the original
% correlationVarFunc.m, since it uses MATLAB 3dim matrix 
% 2. this function uses the rank normalization. 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

correlationInvestigation.m

% the code is to investigate the correlation among genes and it's
% correlation to expression level
% Paul's idea:
% I'm interested in the relationship between mean expression level and correlation. If you take the correlation matrix for one data set, and make distributions of correlations for the following cases:
% - between genes that are in the bottom 1/3 of expression
% - between genes that are in the top 2/3 of expression levels
% - between those two groups
% My hypothesis is that the correlations among the bottom 1/3 should be closer to zero, and the correlations between the high and low groups should be intermediate (on average).
% Also, when you get your thresholded network from a single data
% set, what iXs the distribution of links among those three
% categories (compared to expectation).
% 5. Find the correlation significance for each of the 1000
% correlation ranks, using the shuffled data. 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

correlationInvestigationSub01.m

% this file is to get the absolute mean correlation among the
% genes. I am just interested to see how much of the network we are
% building with the current method involves the genes which are not
% expressed. (bottom 1/3)

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

continuousTSLinksRandomFunction.m

% this is the main function for running the code in parallel and
% saving the data. 
% where tissue is the tissue and s is the starting point for the
% correlation matrix and n is the number of edges to be checked
% (cause I am running in parallel, so I choose how long each run takes)
%function[] = continuousTSLinksFunction(tissue, s, n)
% s is the starting point (cause I get in limited number of edges
% each time) - n is the count of edges , I have ~57e+6 edges, for
% brain, for example. 1e+6 takes about half an hour to run. I giveo
% it 10mil at a time, I guess. 10e+6 takes 5 hours, and probably
% 10G of RAM. 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

correlationVar.m

% this is the main file for getting the correlation variation among
% different datasets for a tissue. The functions for correlation
% variation are correlationVarFunc.m and a newer version
% correlationVarFunc01.m which is much faster. 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

correlationVarFunc01.m

% this function gets the mean and variance of the correlation
% values between the given index of genes from the correlation
% matrix of the datasets for one tissue. I had to write a function
% for it, limiting the input matrix since loading all the datasets
% on the ram at a time was not possible and also I can do it in
% parallel using the function
%
% function[] = correlationVarFunc01(ib, ie, jb, je, tissue)
%
% where ib, ie, jb and je define the sub-matrix of the correlation
% matrix and tissue is a string, name of the tissue(which is name
% of the folder for tissue)
% Notes:
% 1. this function is faster than the original
% correlationVarFunc.m, since it uses MATLAB 3dim matrix 
% 2. for normalization, this function doesn't use the ranking but
% transfers the correlation values between [0,1], this is faster
% than rank normalization but might not be the best normalization

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

correlationVarFunc02.m

% this function gets the mean and variance of the correlation
% values between the given index of genes from the correlation
% matrix of the datasets for one tissue. I had to write a function
% for it, limiting the input matrix since loading all the datasets
% on the ram at a time was not possible and also I can do it in
% parallel using the function
%
% function[] = correlationVarFunc02(ib, ie, jb, je, tissue)
%
% where ib, ie, jb and je define the sub-matrix of the correlation
% matrix and tissue is a string, name of the tissue(which is name
% of the folder for tissue)
% Notes:
% 1. this function is faster than the original
% correlationVarFunc.m, since it uses MATLAB 3dim matrix 
% 2. this function uses the rank normalization. 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

corrHomoStatsAndPlots.m

% in this file I work on getting some statistics of my aggregated
% tissue networks obtained from the homotest and the average rank
% comparison.

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

corrQNormReturn.m

% dsName: name of the dataset
% corrArray: the list of correlation values
% for the given dataset, this function returns the normalized
% ranked values of the correlation values in the corrArray

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

corrQuantileReturn.m

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

dataSetCorrQReturn.m

% this is a function to return the quantile vector of a dataset

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

edgeInfo.m

% In this file I have the code for the part "Characterizing the
% links" in the project report. 
% 1. edges shared between every two networks at different
% 2. distribution of number of repeated edges for each tissue. 3*5
% histograms.
% 3. getting sum of all the binary networks for different
% thresholds.  
% 4. saving the summed networks for each tissue
% 5. getting the tissue distance for each edge.

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

folderNetwork.m

% this file is on making tissue(here I call it network) specific
% coexpression networks. I am doing it using correlation networks. 
% 1. making a corrNetwork from one experiment/ save the sparce matrix
% for all experiments
% 1.1 building the network for a selected set of genes, this I used
% for allExp genes (7835) and also for the tissue ExpGenes
% 2. measure the convergence among the networks build as adding up
% networks. (start with skeletal muscle)
% 3. run the network for at least skeletal muscle. 
% 4. trying to standardize the coexpression data


>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


getExpRank.m

% getting the ranked expression levels for the genes

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

getNetwork01.m

%given the address of a dataset, this file gets the correlation
%network, having the correlation method

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

gettingTissueExpressedGenes.m

% In this file I am getting the genes which are expressed in 80% of
% the datasets in a tissue. 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

gettingTissueNetworks.m

% this file is to get networks for each tissue. I have the binary
% networks the tissues, here I want to:
% 
% 1. get the links between the expressed genes. I want to study the
% mean and median of the top %.5 edges of the genes expressed in a
% tissue. I want to know where my links stand. The result of this
% part, is for, 5,000,000 links, I keep the corr ranks, median and
% mean? I select from them. I think I can prove if the median is
% dropping, I am not missing any links based on other thrs.
%
% For each tissue:
% between the links which are expressed everywhere, 
% I want those which are highly coexpressed in the tissue. 
% Therefore:
% For each dataset, I rank the corrMatrix, 
% I add the exp gene links to each other. 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

GTExdata.m

% 1.NEW GTEx data 
%%%%%%%%%%%%%%%
% get the sample file and see what samples do we have. each sample
% has a tissue and a sub tissue. Lets see how many of them are
% blood, brain, muscle, liver and lung

% 2.old GTEx data 
%%%%%%%%%%%%%%%%
% this is how I process the GTEx data from seq file to the gene
% expression count
% 1. opening the data
% 2. opening metadata - meta data has the sample tissue information

%%%%%%%%%%%%%%%%
% 6. OBSOLETE Reproducibility of affy netwroks in GTEx
% 6. NEW Reproducibility of affy netwroks in GTEx
% 6.5 saving the GTEx networks.
% 7. Reproducibility of TS and common genes in GTEx
% 8. Reproducibility of the pure links
% 9. GTEx exp genes
% 10. writing the GTEx data into .txt files (Also the new one with
% blood)
% 11. rep box plot 
% 12. GTEx rep of pure links: the actual links with TSS 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

infNetworks.m

% for each dataset, build the sparse coexpression networ with the
% -2% of edges and 2% of edges for correlation. Keep the average
% gene expression level for each gene as well as its STD among the
% samples, as well as its expression rank which is a value between
% 1:6
%
% The goal after this process is to explore the characteristics of
% the TS links VS the common links. 
% 1. getting the netStr file for each of the data sets
% 2. checking the percent of edges in every coexpression category for
% some of the networks
% 3. recording the correlation quantiles
% 4. saving the qExpr (the expression level separately )

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

linkExprInfo.m

% this is a file which has the information for studying the link
% expression 
% 1. savign the probe expression data and it's related
% information. 
% 2. savinge the gene expression data. the information saved for
% probe expression can be shared 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


netOverlap_function.m

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

networkComparison.m

%%%%% these are the analysis I do after the geneExpression analysis
%%%%% I did in the geneExpressionMain.m file. 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

networkExploreMain.m

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


networkStudies.m

% This is to study the networks only. 
% comparing the hubs for the tissue networks and TS links.
% 0. Count of edges, count of nodes, starting genes, connectivity,
% density
% .0.5 node degree changes 
% 1. node degree and hubs
% 2. study of the hubs for the ts links
% 3. study of the TFs for the links.
% 4. get the common edge info
% 5. HKE studies for the tissue networks
% 6. The GTEx overlap of the networks

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

nodeDegreeDist.m

% here I make node degree distributions for all my 20 networks:
% TSN, TAN, GTEx, pure

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

outlierDetectionBasic.m

%getting the sample matrix, it returns array including the index-or
%names- of the samples to be removed

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

outlierDetectionSBM.m

% detecting outliers based on correlation and sort-by-median
% algorithm 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

pureLinksAndCommonGenes.m

% load the functional terms for each of them 
% load the pure edges: common functionality of the two genes in one
% tissue, versus their functionality in other tissues: How many of
% them have list of functions in other tissues, and they should not
% be overlapped. 
% This code is for the last part of the results and studies three
% things:
% 1. shift of functioanlity for the confidently highly expressed
% genes in ALL tissues : subset of the 3.
% 1.5 shift of functionality for confidently highly expressed in
% more than one tissue. subset of the 3. 
% 2. shift of functionality for the partners of the pure links. 
% 3. shift of functionality for any of the genes expressed in more
% than one tisssue. 
% 3.2 using JCC distance: here is to get genes whith little SOF
% 4. get the histograms. 
% 6. I check the pureLinkMat (with shift of functionality) in
% GTEx. The code is in GTEx, part 9. 
% 7. To test the clusters
% 8. the plots for the functional identity of the genes >> those
% volcano like plots are in this part, 1, 2, 3.
part, 1, 2, 3.
% 8.35 studying the selected instances:
% 9. examining the moderate JCC instances for the cases
% I can use almost all the above code
% 10. getting the tsn for all vs things : also, this is for JCC <
% 0.01, I should do it for JCC < .1 really and then it is done.
% 11. that tiny bar plot of count of genes expressed in multiple
% tissues
% 12. The violin plot of the JCC and FI# for the pure, semi-pure and
% null
% 13. getting the heatmap of XYWZ
% 14. Individual examples in the last two result parts (SOF and
% pure)
% 15. links that represent functional similarity in ATN
% 16. links that represent functional similarity in TSN
% 17. pure links that represenet functional similarity
% 18. functional implication of different groups (networks) of
%links: I find how well functions are represented in the networks
% 19. Which functions are represented in WHICh networks
% 20. The network numbers: how the TS link portion decrease
% 21. The node degree distribution for genes which have pure link
%     vs those who don't, between the TS links
% 22. The functional enrichment for genes which have pure links,
%      versus those who do not
% 23. TANs and TSNs node degree distributions
% 24. STD of the expression for partners of pure links between
% tissues: examples of low STD with examples of high STD.
% 25.
% 26.
% 27.
% 28. the functions enriched in tissues
% 29.
% 30. response to Bioinf Reviewer: random links for TSN and pure
% (copy of 19)~

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

pureNeighboursShared.m

% Here my goal is just to get the shared neighbour thing and the
% shared JCC thing. 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

puretTSlinkStudy.m

% This is to look at the expression level of the coexpressed links,
% make some plots and check the reproducibility of the links. 
% 
% 1. get the count of pure links as those expressed in the
% datasets, compare their TS distribution with others and make dist
% plots
% 
% 10. testing the expression level for the TS links and finding
% examples
% 11. comparing the pValues vs the distances
% 11.5 getting the TS links - obsolete
% 11.75. getting the TS links - new
% 12. new pure link study: getting the everywhere expressed genes
% using ANOVA
% 12.5 doing the regression for different tissue expression
% 13. Getting the TS genes using ANOVA
% 15. get the final table
%15.1 get the genral info for ts
% 15.2 get the count of pure and other ts links
% 16: get the expression level for the regression identified pure
% links. get the ATN for those genes. 
% 17: get the ATN and TSN for the top 1/3 links.
% Draft

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

rankSum.m

% this function is to get sum of the ranked correlation values for
% all genes in a tissue, given the tissue

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

regressionLink.m

% this is to study regression : 
% 3. Getting regression for the list of TS links. HUH. 
% 4. finding the Thr for the rsqrd

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

removeSamples02.m

% function [] = removeSamples(folder, fileName, samples)
% reads the dataset from the folder/file, remove the samples from
% the 'dataset' structure read from the file and save it back
% again, adding 'clean' to the beginning of the fileName
% samples : index of the desired removing samples

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

removeSamples03.m

% NOTE : removeSample03 does not actually save the changed
% matrix. removeSample02 would. 
% function [] = removeSamples03(folder, fileName, samples)
% reads the dataset from the folder/file, remove the samples from
% the file for the heatmap plot. 
% samples : index of the desired removing samples

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

removeSamples04.m

% function [] = removeSamples(folder, fileName, samples)
% reads the dataset from the folder/file, remove the samples from
% the 'dataset' structure read from the file and save it back
% again, adding 'clean' to the beginning of the fileName
% samples : index of the desired removing samples

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

runCorrelationHomoTest.m

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


sortedGeneIndex.m

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

testSpecificCases.m

% 1. Testing individual cases in GTEx data, in GO data and also
% closer looks into networks. 

% 1. Testing individual cases in GTEx data, in GO data and also
% closer looks into networks. 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

theCommonNetwork.m

% basically, the common network and other things. 
% here I find the common network amongst the genes that are
% expressed in 80% of the datasets for each tissue. I got these
% genes from the file theContiunousTSlinkSelection. 
% 1. get the common network from individual thr
% 2. get the common network from the overal thr for each tissue 
% ?. get the common network amongs all the datasets and with one
% thr only. 
% 3. the distribution of TS values amongst all the tissues
% 3.5 distance for the next TS of the links
% 4. the stat test for the common network and also based on the
% node degree. 
% 5. Comparison of links in different thresholds for individual
% networks
% 6. Genes shared between the networks

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


theContinuousTSlinkSelection.m

% In this code I do all the side tests and explorations for
% continuous ts link selection. 
%
%% Part 1: pilot phases for the measurement and measuring the TS
%% for the links
% 
% 1. First step: trying out my algorithm for some random TS
% links. Just wanted to check if it behaves as I expected and code
% the core of the algorithm. 
%
% 2. Trying the code for some real links between HKG - also a
% shoortcut for making random finalMat - it is different than
% the
% 1st part's random finalMat, that you don't get to choose values
% for each tissue. I used it for checking the "everywhere
% coexpressed" or "everywhere absent" space value for links. 
%
% 3. Trying the same code differently - I get the correlations by
% for randomly selected pairs of genes. (not the first part!)
%
% 4. assembling the results and some other checkings - this is
% where I study the everywhere potential links
%
% 5. finding final mat and space for a given pair of genes and
% tissue
%
% 6. getting links expressed in all the tissues - between them,
% getting the TS links.
%
% 6.5. getting the TS links for each tissue, the links under the
% three categoreis
% 
% 7. finding final mat for a given set of g1 and g2
%
% 77. Doing a function for a whole tissue
%
%% Part 2: Determining the THR value for identification of the
%% links and putting the cut off for the tissues. 
%
% 8. TS for a random list of networks: generating random list of datasets for testing and testing
% them. (code for calling the function is here)
%
% 9. Cut off the Threshold: For each link passed the thr, what is
% the TS value for the next best tissue? Get this for all the tissues.
% 
% 10. load the extended TS links and the previous ts links and
% check the correlation of the two values. 
% % 11. put together the results from psuedo tissues and do sth. 
%
% 12. Get the hist for the various neg adjusted
% 
% 13.% 13. Get the tissue negAdjusted and pvalues. Also, modified pvalues
%
% 14. get the thresholds for the union. 
% 15. get the average values

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

tissueNetworkComparisons.m

% here is the file to compare the tissue networks, this is where I
% do the clustering. For clustering, I use WGCNA in R, here I
% prepare the data for WGCNA. 
% 0. Just loading the networks
% 1. getting the expressed gens 
% 2. the ribosome link check
% 3. clustering
% 4. getting the genes for the R clustering




