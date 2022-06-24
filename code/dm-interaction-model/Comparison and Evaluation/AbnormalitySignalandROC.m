% Eighth algorithm: Plotting abnormality signals and Roc Curves
set(0,'defaultfigurecolor',[1 1 1])
clc
clear
close all
curDir = pwd;
MB = true;
numberParticle = 10;
scenarioTesting = 2;% 1(PM): training data, 2(PA), 3(ES), 4(Uturn): testing data
scenarioTrain = 1; % 1 (PM): training data
Case = 2; % 1: Less numbers, 2: Best case, 3: High numbers

%% Directorie of training results (Joint) 
if MB == true
    dirMatTrain = ['C:\Users\mohamad.baydoun\Desktop\Interaction journal (2019)',...
        '\Generated Data\' num2str(scenarioTrain)];
    dirResultsTesting = ['C:\Users\mohamad.baydoun\Desktop',...
        '\Interaction journal (2019)\Jointly Technique',...
        '\ResultsTesting\' num2str(scenarioTesting) '\C' num2str(Case)];
    
else
    dirMatTrain = ['C:\Users\damian.campo\Desktop\Interaction journal (2019)',...
        '\Generated Data\' num2str(scenarioTrain)];
    dirResultsTesting = ['C:\Users\damian.campo\Desktop',...
        '\Interaction journal (2019)\Jointly Technique',...
        '\ResultsTesting\' num2str(scenarioTesting) '\C' num2str(Case)];
end
cd(dirMatTrain)
load('genSyncdata.mat')
inputTrainData = [structSyncData.Filtered.xPos,structSyncData.Filtered.yPos];
inputTrainData = inputTrainData(100:end-100,:);
cd(dirResultsTesting)
load(['EstimationAbnInteractionT' num2str(scenarioTrain) 'Par' num2str(numberParticle) '.mat'])
AbnormalSignalJointly = EstimationAbn.estimationAbn;

%% Directorie of training results (Separate)
if MB == true
    dirResultsTesting = ['C:\Users\mohamad.baydoun\Desktop',...
        '\Interaction journal (2019)\Separately Technique',...
        '\ResultsTesting\' num2str(scenarioTesting) '\C' num2str(Case)];
    dirComparisonResultsTest = ['C:\Users\mohamad.baydoun\Desktop',...
        '\Interaction journal (2019)\Comparison and Evaluation',...
        '\TestingAbnormalityRocResults\' num2str(scenarioTesting) '\C' num2str(Case)];
else
    dirResultsTesting = ['C:\Users\damian.campo\Desktop',...
        '\Interaction journal (2019)\Separately Technique',...
        '\ResultsTesting\' num2str(scenarioTesting) '\C' num2str(Case)];
    dirComparisonResultsTest = ['C:\Users\damian.campo\Desktop',...
        '\Interaction journal (2019)\Comparison and Evaluation',...
        '\TestingAbnormalityRocResults\' num2str(scenarioTesting) '\C' num2str(Case)];
end
cd(dirResultsTesting)
load(['EstimationAbnInteractionT' num2str(scenarioTrain) 'Par' num2str(numberParticle) '.mat'])
AbnormalSignalSeparately = EstimationAbn.estimationAbn;
cd(curDir)

%% scatter scenario
h = figure;
h.Position = [538   227   717   593];%[2023 168 810 661];
hold on
plot1 = scatter(inputTrainData(1:718,1),inputTrainData(1:718,2),7,'k','fill');

x = AbnormalSignalSeparately.dataTestOdometric(1,:);
y = AbnormalSignalSeparately.dataTestOdometric(2,:);
plot2 = scatter(x,y,16,'b','fill');

plot3 = scatter(x(66:112),y(66:112),13,'r','fill');
scatter(x(395:457),y(395:457),13,'r','fill')

plot4 = scatter(x(1),y(1),150,'g','*');
plot5 = scatter(x(end),y(end),150,'c','+');

xlab = xlabel('$x$','interpreter','latex');
xlab.FontSize = 16;
ylab = ylabel('$y$','interpreter','latex');
ylab.FontSize = 16;
axis([-22 22 0 40])
grid on
lgg = legend([plot1, plot2, plot3,plot4, plot5], {'Training Data',...
    'Observed data','Abnormal area','Start point','End point'},'interpreter','latex');
lgg.FontSize = 11;
lgg.Position = [0.1004    0.9432    0.8605    0.0354];  
lgg.Orientation = 'horizontal';
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

%%  Plot abnormality Signals for both Approaches
%% Case 1:
% Abnormality measurement of S-A at discrite level: Odometry
h1 = figure;
tit1 = title('Abnormality measurement of S-A at discrite level: Odometry','interpreter','latex');
tit1.FontSize = 15;
plotAbnormalitySignal(AbnormalSignalSeparately.abnormdb1Odm_Probabilistic,...
    AbnormalSignalSeparately.abnormdb1Odm_ProbabilisticSmooth,h1);
% Abnormality measurement of S-A at discrite level: Control

% Abnormality measurement of S-J at continous level: Odometry

% Abnormality measurement of S-J at discrite level: Control

%% Case 2:
% Abnormality measurement of S-A at discrite level: Odometry
% Abnormality measurement of S-A at discrite level: Control
% Abnormality measurement of S-J at continous level: Odometry
% Abnormality measurement of S-J at discrite level: Control

%% Case 3:
% Abnormality measurement of S-A at discrite level: Odometry
% Abnormality measurement of S-A at discrite level: Control
% Abnormality measurement of S-J at continous level: Odometry
% Abnormality measurement of S-J at discrite level: Control
%% ROC performances (AUC) And ACC
AbnSepOdm = AbnormalSignalSeparately.abnormdb1Odm_Probabilistic(2:end);
AbnSepOdmSmooth = AbnormalSignalSeparately.abnormdb1Odm_ProbabilisticSmooth(2:end);
AbnSepCont = AbnormalSignalSeparately.abnormdb1Cont_Probabilistic(2:end);
AbnSepContSmooth = AbnormalSignalSeparately.abnormdb1Cont_ProbabilisticSmooth(2:end);

AbnJoinOdm = AbnormalSignalJointly.abnormdb1Odm_Probabilistic(2:end);
AbnJoinOdmSmooth = AbnormalSignalJointly.abnormdb1Odm_ProbabilisticSmooth(2:end);
AbnJoinCont = AbnormalSignalJointly.abnormdb1Cont_Probabilistic(2:end);
AbnJoinContSmooth = AbnormalSignalJointly.abnormdb1Cont_ProbabilisticSmooth(2:end); 

AbnSepOdmSmooth = normalization(AbnSepOdmSmooth);
AbnJoinOdmSmooth = normalization (AbnJoinOdmSmooth);

AbnSepContSmooth = normalization(AbnSepContSmooth);
AbnJoinContSmooth = normalization (AbnJoinContSmooth);

part1 = zeros(1,65); part2 = ones(1,47); part3 = zeros(1,283);
part4 = ones(1,63);part5 = zeros(1,268);
Ground_truth = [part1,part2,part3,part4,part5];

[Roc_OdmSep ,AUC_OdmSep,ACC_OdmSep] = Roc_calculation(AbnSepOdmSmooth,Ground_truth);
[Roc_OdmJoin ,AUC_OdmJoin,ACC_OdmSJoin] = Roc_calculation(AbnJoinOdmSmooth,Ground_truth);

[Roc_ContSep ,AUC_ContSep,ACC_ContSep] = Roc_calculation(AbnSepContSmooth,Ground_truth);
[Roc_ContJoin ,AUC_ContJoin,ACC_ContSJoin] = Roc_calculation(AbnJoinContSmooth,Ground_truth);


%% plotting
figure;
R1 = plot(Roc_OdmSep(1,:),Roc_OdmSep(2,:),'b','LineWidth',2);
hold on
R2 = plot(Roc_OdmJoin(1,:),Roc_OdmJoin(2,:),'r','LineWidth',2);
lgg1 = legend([R1 R2],{'ROC "Odometry case of Separately approache"','ROC "Odometry case of Jointly approache"'});
lgg1.Position = [0.3542    0.1925    0.5446    0.0869];
xlabel('False Positive');
ylabel('True Positive');
title('ROC Features')
grid on
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

figure;
R1 = plot(Roc_ContSep(1,:),Roc_ContSep(2,:),'b','LineWidth',2);
hold on
R2 = plot(Roc_ContJoin(1,:),Roc_ContJoin(2,:),'r','LineWidth',2);
lgg1 = legend([R1 R2],{'ROC "Control case of Separately approache"','ROC "Control case of Jointly approache"'});
lgg1.Position = [0.3542    0.1925    0.5446    0.0869];
xlabel('False Positive');
ylabel('True Positive');
title('ROC Features')
grid on

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%export_fig rocfinale2.pdf -pdf -transparent


