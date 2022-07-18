
function[a, b] = netOverlap_function(ID1, ID2)
FC = '3'
FDR = '0012'
load(sprintf(['~/resultsAndFigures/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))


tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
geneID = ID1
% 70. get the partners in the two networks
atpartners = zeros(5, 18494);
tspartners = zeros(5, 18494);
for t = 1 : 5
    tissue = tissues{t};
    load( ['~/networks/tissues/' tissue '/' ...
           'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
    atn{t} = binNet;
    atpartners(t, :) = (binNet(geneID, :)) + (binNet(:, geneID)');
    tsnet = finalTable(t).wholeNet;
    tspartners(t, :) = (tsnet(geneID, :)) + (tsnet(:, geneID)');
end
first.atp = atpartners;
first.tsp = tspartners;


geneID = ID2
% 70. get the partners in the two networks
atpartners = zeros(5, 18494);
tspartners = zeros(5, 18494);
for t = 1 : 5
    tissue = tissues{t};
    load( ['~/networks/tissues/' tissue '/' ...
           'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
    atn{t} = binNet;
    atpartners(t, :) = (binNet(geneID, :)) + (binNet(:, geneID)');
    tsnet = finalTable(t).wholeNet;
    tspartners(t, :) = (tsnet(geneID, :)) + (tsnet(:, geneID)');
end
second.atp = atpartners;
second.tsp = tspartners;

sib = first.atp + second.atp;
a = sum(sib' == 2);

sib = first.tsp + second.tsp;
b = sum(sib' == 2);
end

