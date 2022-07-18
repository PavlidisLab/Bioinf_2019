% function [] = removeSamples(folder, fileName, samples)
% reads the dataset from the folder/file, remove the samples from
% the 'dataset' structure read from the file and save it back
% again, adding 'clean' to the beginning of the fileName
% samples : index of the desired removing samples

function[dataSet] = removeSamples04(dataSet, samples)

    toKillIndex = samples; 
    dataSet.GEOID(toKillIndex) = [];
    dataSet.mat(:, toKillIndex) = [];

end