%   Fourth algorithm: Generation of dictionary
set(0,'defaultfigurecolor',[1 1 1])
clc
clear
close all
curDir = pwd;
MB = true;
scenarioTrain = 1; % 1(PM): training data
Case = 1; %1:  Less numbers, 2: Best case, 3: High numbers
%% Directories
if MB == true
    dirMat = ['C:\Users\mohamad.baydoun\Desktop\Interaction journal (2019)',...
        '\Generated Data\' num2str(scenarioTrain)];
    dirResultsTrain = ['C:\Users\mohamad.baydoun\Desktop',...
        '\Interaction journal (2019)\Separately Technique',...
        '\ResultsTraining\' num2str(scenarioTrain) '\C' num2str(Case)];
    
else
    dirMat = ['C:\Users\damian.campo\Desktop\Interaction journal (2019)',...
        '\Generated Data\' num2str(scenarioTrain)];
    dirResultsTrain = ['C:\Users\damian.campo\Desktop',...
        '\Interaction journal (2019)\Separately Technique',...
        '\ResultsTraining\' num2str(scenarioTrain) '\C' num2str(Case)];
end
cd(dirMat)
load('genSyncdata.mat')
cd(dirResultsTrain)
load('VocabulariesFinalControl.mat')
load('VocabulariesFinalOdometric.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Dictionary Learning
vocOdometricSize = netOdometric.N;
vocControlSize = netcontrol.N;
dataCode1 = trainingTracksOdometric.tracksDiscrete;
dataCode2 = trainingTrackscontrol.tracksDiscrete;

eTot = vocOdometricSize*vocControlSize;                                     %   Total number of events
label = zeros(eTot,4);                                                      %   Combination of all possible events
for i = 1:vocOdometricSize
    label(vocControlSize*(i-1)+1:vocControlSize*(i),1) = i;                 %   Initialize the first column
    label(vocControlSize*(i-1)+1:vocControlSize*(i),2) = 1:vocControlSize;  %   Initialize the second column
end
label(:,3) = 1:eTot;
dataCode = [dataCode1 dataCode2];
label2 = unique(dataCode, 'rows');                                           % remove repetitions (and set in sorted order)
newlabel = 1:size(label2);
newlabel = newlabel';
label2 = horzcat(label2, newlabel);
for j = 1:size(label2,1)
    label(((label2(j,1)-1)*vocControlSize + label2(j,2)),4) = label2(j,3);
end
transMatDictionary = zeros(size(label2,1));
for i = 1:size(dataCode,1)-1
    coderow = (dataCode(i,1)-1)*vocControlSize + dataCode(i,2);
    codecol = (dataCode(i+1,1)-1)*vocControlSize + dataCode(i+1,2);
    rowMat = label(coderow,4);                                             %   Row of transition matrix is first node
    colMat = label(codecol,4);
    transMatDictionary(rowMat,colMat) = transMatDictionary(rowMat,colMat) + 1;                 %   Increase occurrence of transition
end
%% Normalize relative frequency
normVect = sum(transMatDictionary,2) + (sum(transMatDictionary,2) == 0);
transMatDictionary = transMatDictionary./repmat(normVect,1,size(transMatDictionary,1));
%% visualize the transition matrix
figure
imagesc(transMatDictionary);
colorbar;
colormap jet
title('Transtion matrix for superstates')
cd(dirResultsTrain)
save('Dictionary.mat','transMatDictionary','label')
cd(curDir)