%   Fifth algorithm: State estimation and abnormality detection
set(0,'defaultfigurecolor',[1 1 1])
clc
clear
close all
curDir = pwd;
MB = true;
scenarioTrain = 1;
scenarioTesting = 2;% 1(PM): training data, 2(PA), 3(ES), 4(Uturn): testing data
plotBool = true;
Case = 1; % 1: Less numbers, 2: Best case, 3: High numbers
N = 40;%10,100,80,60, 40,20,5 number of particles

%% Directories
if MB == true
    dirMat = ['C:\Users\mohamad.baydoun\Desktop\Interaction journal (2019)',...
        '\Generated Data\' num2str(scenarioTrain)];
    dirResultsTrain = ['C:\Users\mohamad.baydoun\Desktop',...
        '\Interaction journal (2019)\Separately Technique',...
        '\ResultsTraining\' num2str(scenarioTrain) '\C' num2str(Case)];
    dirResultsTesting = ['C:\Users\mohamad.baydoun\Desktop',...
        '\Interaction journal (2019)\Separately Technique',...
        '\ResultsTesting\' num2str(scenarioTesting) '\C' num2str(Case)];
else
    dirMat = ['C:\Users\damian.campo\Desktop\Interaction journal (2019)',...
        '\Generated Data\' num2str(scenarioTrain)];
    dirResultsTrain = ['C:\Users\damian.campo\Desktop',...
        '\Interaction journal (2019)\Separately Technique',...
        '\ResultsTraining\' num2str(scenarioTrain) '\C' num2str(Case)];
    dirResultsTesting = ['C:\Users\damian.campo\Desktop',...
        '\Interaction journal (2019)\Separately Technique',...
        '\ResultsTesting\' num2str(scenarioTesting) '\C' num2str(Case)];
end
cd(dirMat)
load('genSyncdata.mat')
cd(dirResultsTrain)
load('VocabulariesFinalControl.mat')
load('VocabulariesFinalOdometric.mat')
load('Dictionary.mat')
cd(curDir)

%% Mean and Covariance (training data)
% odometric
odometric.averageStateOdm = netOdometric.nodesMean;                                              %   Mean neurons of position data
odometric.averageDivOdm = netOdometric.nodesMeanComp;
odometric.radiusStateOdm = netOdometric.nodesRadAccept*2;                                        %   Acceptance neuron radius of position data
odometric.radiusDivOdm = netOdometric.nodesRadAcceptComp*2;
covarianceOdm = ones(4,4,length(netOdometric.datanodes))*1e+100;
for jj = 1:length(netOdometric.datanodes)
    dataOdmnodes {1,jj} =  [netOdometric.datanodes{1, jj},netOdometric.CompDataNodes{1, jj}];
    covarianceOdm(:,:,jj) = cov(dataOdmnodes{1,jj});
end
odometric.CovarianceOdm = covarianceOdm;
% control
control.averageStateCon = netcontrol.nodesMean;                                              %   Mean neurons of position data
control.averageDivCon = netcontrol.nodesMeanComp;
control.radiusStateCon = netcontrol.nodesRadAccept*2;                                        %   Acceptance neuron radius of position data
control.radiusDivCon = netcontrol.nodesRadAcceptComp*2;
covarianceCont = ones(4,4,length(netOdometric.datanodes))*1e+100;

for kk = 1:length(netcontrol.datanodes)
    dataContnodes {1,kk} =  [netcontrol.datanodes{1,kk},netcontrol.CompDataNodes{1,kk}];
    covarianceCont(:,:,kk) = cov(dataContnodes{1,kk});
end
control.CovarianceCont = covarianceCont;
%%% For Normalization
odometric.minMAxNormOdm = [structSyncData.NormMinMax.xPos,...
    structSyncData.NormMinMax.yPos,structSyncData.NormMinMax.xPosdiv,...
    structSyncData.NormMinMax.yPosdiv]';
control.minMAxNormControl = [structSyncData.NormMinMax.S,structSyncData.NormMinMax.V,...
    structSyncData.NormMinMax.Sdiv,structSyncData.NormMinMax.Vdiv]';
%%%%%%%%%%%%%%%%%%%%%%%
codebook = label;
transitionMat = transMatDictionary;
estimationAbn = MJPF_Separately(odometric, control,...
    transitionMat,codebook,curDir, scenarioTesting,MB,plotBool,dirResultsTesting,scenarioTrain,Case,N);
EstimationAbn.estimationAbn = estimationAbn;
cd(dirResultsTesting)
save(['EstimationAbnInteractionT' num2str(scenarioTrain) 'Par' num2str(N)],'EstimationAbn'); % T: Training
cd(curDir)