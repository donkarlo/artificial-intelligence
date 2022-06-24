%   Sixth algorithm: Comparison and Evaluation in training process (Var,
% connections, entropy)

set(0,'defaultfigurecolor',[1 1 1])
clc
clear
close all
curDir = pwd;
MB = true;
scenarioTrain = 1; % 1 (PM): training data
Case = 1; % 1: Less numbers, 2: Best case, 3: High numbers
%% Directories of training results (Separate)
if MB == true
    dirResultsTrain = ['C:\Users\mohamad.baydoun\Desktop',...
        '\Interaction journal (2019)\Separately Technique',...
        '\ResultsTraining\' num2str(scenarioTrain) '\C' num2str(Case)];
    dirComparisonResultsTrain = ['C:\Users\mohamad.baydoun\Desktop',...
        '\Interaction journal (2019)\Comparison and Evaluation',...
        '\TrainingComparisonResuls\' num2str(scenarioTrain) '\C' num2str(Case)];
else
    dirResultsTrain = ['C:\Users\damian.campo\Desktop',...
        '\Interaction journal (2019)\Separately Technique',...
        '\ResultsTraining\' num2str(scenarioTrain) '\C' num2str(Case)];
    dirComparisonResultsTrain = ['C:\Users\damian.campo\Desktop',...
        '\Interaction journal (2019)\Comparison and Evaluation',...
        '\TrainingComparisonResuls\' num2str(scenarioTrain) '\C' num2str(Case)];
end
cd(dirResultsTrain)
load('Dictionary.mat')
load('VocabulariesFinalOdometric.mat')
load('VocabulariesFinalControl.mat')

validLabSep = label(find(label(:,4)~=0),:);
posNodes = netOdometric.datanodes;
velNodes = netOdometric.CompDataNodes;
SVNodes = netcontrol.datanodes;
velSVNodes = netcontrol.CompDataNodes;

std1 = [];
std2 = [];
std3 = [];
std4 = [];

for i = 1:size(validLabSep,1)
    indOdom = validLabSep(i,1);
    indCont = validLabSep(i,2);
    std1 = [std1; std(posNodes{1,indOdom})];
    std2 = [std2; std(velNodes{1,indOdom})];
    std3 = [std3; std(SVNodes{1,indCont})];
    std4 = [std4; std(velSVNodes{1,indCont})];
end

stdSeparately= [std1 std2 std3 std4];
ComparisonTraining.stdSeparatelyMean = mean(stdSeparately); % Var of Separate

TranMatSeparately = transMatDictionary;
ComparisonTraining.sizeMatSepararely = size(TranMatSeparately,1);
ComparisonTraining.numbConnectSeparartely = sum(sum(TranMatSeparately> 0,1));

%% Directories of training results (Joint)
if MB == true
    dirResultsTrain = ['C:\Users\mohamad.baydoun\Desktop',...
        '\Interaction journal (2019)\Jointly Technique',...
        '\ResultsTraining\' num2str(scenarioTrain) '\C' num2str(Case)];
else
    dirResultsTrain = ['C:\Users\damian.campo\Desktop',...
        '\Interaction journal (2019)\Jointly Technique',...
        '\ResultsTraining\' num2str(scenarioTrain) '\C' num2str(Case)];
end
cd(dirResultsTrain)
load('Dictionary.mat')
load('VocabulariesFinalOdometric.mat')

validLabJoint = label(find(label(:,4)~=0),:);
stdJoint = [];

for i = 1:size(validLabJoint,1)
    stdJoint = [stdJoint; std(netOdometric.datanodes{1,validLabJoint(i,1)})];
end
ComparisonTraining.stdJointMean = mean(stdJoint);

TranMatJointOdm = transitionInfoOdometric.transitionMat;
ComparisonTraining.sizeMatJoint = size(TranMatJointOdm,1);
ComparisonTraining.numbConnectJoint= sum(sum(TranMatJointOdm> 0,1));
cd(curDir)

ComparisonTraining.mismatchesNnodes = abs(ComparisonTraining.sizeMatSepararely - ComparisonTraining.sizeMatJoint);
ComparisonTraining.mismatchesNconnect = abs(ComparisonTraining.numbConnectSeparartely - ComparisonTraining.numbConnectJoint);

ComparisonTraining.entropySeparately = entropyMat(TranMatSeparately);
ComparisonTraining.entropyJointOdm = entropyMat(TranMatJointOdm);
cd(dirComparisonResultsTrain)
save(['ComparisonTrainingT' num2str(scenarioTrain)],'ComparisonTraining')
cd(curDir)


