% NOTE : removeSample03 does not actually save the changed
% matrix. removeSample02 would. 
% function [] = removeSamples03(folder, fileName, samples)
% reads the dataset from the folder/file, remove the samples from
% the file for the heatmap plot. 
% samples : index of the desired removing samples


function[] = removeSamples03(folder, fileName, figFolder, samples, ...
                             titleExtra, viewOnly)

load([folder fileName]);

    toKillIndex = samples; 
    dataSet.GEOID(toKillIndex) = [];
    dataSet.mat(:, toKillIndex) = [];
    
    outliersCount = length(samples);

    exprData = dataSet.mat;
    sampCount = size(exprData, 2);
    sib = corr(exprData);
    
    'corr'
    
    
    book = cell(1, sampCount);
    for j = 1: sampCount
        book{j} = dataSet.GEOID{j}{1};
    end
    

    % make and save heatmap
    h = figure; 
    heatmap(sib, dataSet.GEOID, dataSet.GEOID,  [], 'TickAngle', 45, 'ShowAllTicks', ...
            true, 'Colorbar', true, 'Colormap', 'summer')
    'heatmap done'
    titleExtra
      title(sprintf('%s %s  Number of samples : %d, outliers: %d', titleExtra, ...
                    fileName, sampCount, outliersCount))
      if(viewOnly)
          return
      else
          fileName = [titleExtra 'geneSampleCorr' fileName(1:(end-4))];
          print(h, '-depsc', [figFolder fileName '.eps']);
          print(h, '-dpdf', [figFolder fileName '.pdf']);
          
      end

end
