% function [] = removeSamples(folder, fileName, samples)
% reads the dataset from the folder/file, remove the samples from
% the 'dataset' structure read from the file and save it back
% again, adding 'clean' to the beginning of the fileName
% samples : index of the desired removing samples

function[] = removeSamples02(folder, fileName, samples, viewOnly)

load([folder fileName]);


    toKillIndex = samples; 
    dataSet.GEOID(toKillIndex) = [];
    dataSet.mat(:, toKillIndex) = [];

    exprData = dataSet.mat;
    sampCount = size(exprData, 2);
    sib = corr(exprData);
    
    'corr'
    
    % book = cell(1, sampCount);
    % for j = 1: sampCount
    %     book{j} = dataSet.GEOID{j}{1};
    % end
    
    book = dataSet.GEOID;
    
    [a, b] = regexp(fileName, 'GSE[0123456789]*');
    Name = fileName(a:b);
    
    

    % make and save heatmap
    h = figure; 
    heatmap(sib, book, book,  [], 'TickAngle', 45, 'ShowAllTicks', ...
            true, 'Colorbar', true, 'Colormap', 'summer')
    'heatmap done'
      title(sprintf('%s  Number of samples : %d', Name, sampCount))
    if(viewOnly)
        return
else
    save([folder 'outlierRemoved/' fileName], 'dataSet');
    'saved'
end

end