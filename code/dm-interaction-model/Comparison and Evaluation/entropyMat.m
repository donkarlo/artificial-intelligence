function  entropyOut = entropyMat(freqMat)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
transMat = freqMat./repmat(sum(freqMat,2),1,size(freqMat,2));
bitPerData = log2(size(transMat,2));                                        %   Maximum entropy that could be reached
logMat = log2(2*(transMat == 0) + transMat);
logMat = logMat.*(logMat ~= 1);
listEntropy = -sum(transMat.*logMat,2)/bitPerData;
meanEntropy = mean(listEntropy);
stdEntropy = std(listEntropy);

entropyOut.List = listEntropy;
entropyOut.mean =  meanEntropy;
entropyOut.std = stdEntropy;
end

