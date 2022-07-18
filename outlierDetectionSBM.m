% detecting outliers based on correlation and sort-by-median
% algorithm 

function[outlierArray] = outlierDetectionSBM(geneMat, BPviewOnly, ...
                                             figFolder, Name)

%get the correlation among samples
sib = corr(geneMat);
sampleCount = length(sib);
sampleCount

%for each sample, get the 1st, 2nd and 3rd quantile:
Q1st = zeros(sampleCount, 1);
Q2nd = zeros(sampleCount, 1);
Q3rd = zeros(sampleCount, 1);
for i = 1: sampleCount
    quartiles = quantile(sib(i,:), 4);
    Q1st(i) = quartiles(1);
    Q2nd(i) = quartiles(2);
    Q3rd(i) = quartiles(3);
end

%checking the boxplot
outlierArray = zeros(sampleCount, 1);
[a, b] = sort(Q2nd);

% heatmap(sib(:, b), [], [], [], 'Colorbar', true, 'Colormap', ...
%         'summer')

for i = 1:(sampleCount-1)
    Q3rd(b(i)) - Q1st(b(i+1));
    if(Q2nd(b(i)) < Q1st(b(i+1)))
        outlierArray(b(1:i)) = 1;
    end
end

% h = figure;
% boxplot(sib(:, b))
if(BPviewOnly)
    return
else
    %saving the BP
    
    title(sprintf('%s  Number of samples : %d', Name(1:(end-4)), sampleCount))
    fileName = ['sampleCorrBoxplot' Name(1:(end-4))];
    print(h, '-depsc', [figFolder fileName '.eps']);
    print(h, '-dpdf', [figFolder fileName '.pdf']);

end

end