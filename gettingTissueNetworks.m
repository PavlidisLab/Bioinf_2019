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


clear


dataFolder = '~/data/affyArray/tissues/'
folderName = dir(dataFolder)

counter = 0
gCount = 18494;

for t = 3:length(folderName)

    tissue = folderName(t).name;
    exprFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/']
    fileList = dir([exprFolder '*.mat'])

    % getting the expressed genes for that tissue
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    tempMat = zeros(gCount, gCount);
    tempMat(expGenesInd, expGenesInd) = 1;
    
    % getting the single binary networks: 
    % count of gene quantiles
    for i = 1: length(fileList)
        counter = counter + 1
        tic
        %loading new file
        '1'
        load([exprFolder fileList(i).name]);
        
        sib = corr(dataSet.mat');

        % getting the genes which are labeled as expressed in the
        % tissue
        sib = sib - eye(gCount);
        sib = sib .* tempMat;
        
        '2'
        %getting the sparse raw values network
        upperSingle = sib(logical(triu(ones(size(sib)), 1)));
        Q = quantile(upperSingle, (1 - .1))
        sNet1 = sib > Qplus;
        sNet1(1:gCount+1:gCount*gCount) = 0; % putting the upper
                                             % triangle zero
        net = sib .* sNet1;
        sNet2 = sib < Qminus;
        net = net + sib .* sNet2;

        [r, c] = find(net);

        % getting the sparse corr p-values? or sth similar? 
        % for the .1 upper quantile and .6 lower quantile, we
        % divide it into 100. 
        q = 100;
        bookp = upperSingle(upperSingle > Qplus);
        qbookp = quantile(bookp, q);
        bookm = upperSingle(upperSingle < Qminus);
        qbookm = quantile(bookm, q);

        % find the quantile for each of the net values in the final
        % net.
        
        rankNet = zeros(gCount, gCount);
        qs = [qbookm , qbookp];
        tic
        for ii = 1:length(r)
            temp = find(qs < net(r(ii), c(ii)));
            % if (net(r(i), c(i))) < 0
            %     '<0'
            % end
            if length(temp) < 100 % the corr is negative
                rankNet(r(ii), c(ii)) = -(100 - length(temp));
            else % the corr is positive
                rankNet(r(ii), c(ii)) = temp(end) - 100;
            end
        end
        toc

        % getting the sparse expression values
        sparseExpNet = zeros(gCount, gCount);
        for ii=1:length(r)
            g1 = r(ii);
            g2 = c(ii);
            if(qExpr(g1) < qExpr(g2))
                sparseExpNet(g1,g2) = -(qExpr(g1) * qExpr(g2));            
            else
                sparseExpNet(g1,g2) = qExpr(g1) * qExpr(g2);                        
            end
        end
        
        netStr.expNet = sparse(sparseExpNet .* tril(ones(gCount,gCount),-1));
        netStr.coNt = sparse(net .* tril(ones(gCount,gCount),-1));
        netStr.rankNet = rankNet; 
        netStr.qExpr = qExpr;
        netStr.avgExpr = avgExpr;
        netStr.stdExpr = stdExpr;

        
        '3'
        %saving the single Network
        [a, b] = regexp(fileList(i).name ,'GSE[0123456789]+');
        Name = fileList(i).name(a:b)
        save(sprintf('~/networks/tissues/%s/singleNet/netStr_%s.mat', ...
                     tissue, Name), 'netStr', '-v7.3'); % i is included so I can follow up which experiments were in the i'th aggNet
        
        '4'
        %adding it to the aggregated ranked matrix we had so far and
        %get the results - tied rank takes 10m ?
        % aggMat = aggMat + reshape(tiedrank(sib(:)), gCount, gCount);
        % upperAgg = triu(aggMat, 1);
        % QAgg = quantile(upperAgg(:), (1 - 0.005));
        % aggNet = aggMat > QAgg;df
        % sparseAggNet = sparse(aggNet);
        
        % '5'
        % %saving the AggNetwork 
        % save(sprintf('%s%s/AggNet/%dagg.mat', dataFolder, tissue, (i)), 'aggNet');
        
        % clear dataSet
        toc
    end
end



