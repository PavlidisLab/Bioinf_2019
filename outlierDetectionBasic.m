%getting the sample matrix, it returns array including the index-or
%names- of the samples to be removed

function[outlierArray] = outlierDetectionBasic(geneMat)

%get the correlation among samples
sib = corr(geneMat);
sampleCount = length(sib);

lmat = tril(sib, -1);
larray = lmat(lmat>0);

%checking
%length(larray)

%corrTresh: get the 15% quntile among the correlations
corrQ = quantile(larray, 100);
treshold = corrQ(15);

%for each sample, if it has more than 90% of it's correlations less
%than the threshold 
outlierArray = zeros(1, sampleCount);
for i = 1:length(sib)
    %    i
    temp = sum(sib(i,:) < treshold);
    if((temp / sampleCount) > 0.9)
        outlierArray(i) = 1;
    end
end
end