% seventh algorithm: Comparisons: Estimation error and space coverage
set(0,'defaultfigurecolor',[1 1 1])
clc
clear
close all
curDir = pwd;
MB = true;
numberParticle = 40;
scenarioTesting = 1;% 1(PM): training data, 2(PA), 3(ES), 4(Uturn): testing data
scenarioTrain = 1; % 1 (PM): training data
Case = 2; % 1: Less numbers, 2: Best case, 3: High numbers
plotting = true;

%% Directorie of training results (Joint) 
if MB == true
    dirResultsTesting = ['C:\Users\isip40\Desktop',...
        '\Interaction journal (2019)\Jointly Technique',...
        '\ResultsTesting\' num2str(scenarioTesting) '\C' num2str(Case)];
    
else
    dirResultsTesting = ['C:\Users\damian.campo\Desktop',...
        '\Interaction journal (2019)\Jointly Technique',...
        '\ResultsTesting\' num2str(scenarioTesting) '\C' num2str(Case)];
end
cd(dirResultsTesting)

load(['EstimationAbnInteractionT' num2str(scenarioTrain) 'Par' num2str(numberParticle) '.mat'])
CountOutsideJoint = EstimationAbn.estimationAbn.meanCountOutside; %State space coverage
errorJoint = EstimationAbn.estimationAbn.errorLocal; %Estimation error

%% Directorie of training results (Separate)
if MB == true
    dirResultsTesting = ['C:\Users\isip40\Desktop',...
        '\Interaction journal (2019)\Separately Technique',...
        '\ResultsTesting\' num2str(scenarioTesting) '\C' num2str(Case)];
    dirComparisonResultsTest = ['C:\Users\isip40\Desktop',...
        '\Interaction journal (2019)\Comparison and Evaluation',...
        '\TestingComparisonResuls\' num2str(scenarioTesting) '\C' num2str(Case)];
else
    dirResultsTesting = ['C:\Users\damian.campo\Desktop',...
        '\Interaction journal (2019)\Separately Technique',...
        '\ResultsTesting\' num2str(scenarioTesting) '\C' num2str(Case)];
    dirComparisonResultsTest = ['C:\Users\damian.campo\Desktop',...
        '\Interaction journal (2019)\Comparison and Evaluation',...
        '\TestingComparisonResuls\' num2str(scenarioTesting) '\C' num2str(Case)];
end
cd(dirResultsTesting)
load(['EstimationAbnInteractionT' num2str(scenarioTrain) 'Par' num2str(numberParticle) '.mat'] )
CountOutsideSeparat = EstimationAbn.estimationAbn.meanCountOutside; %State space coverage
errorSeparat = EstimationAbn.estimationAbn.errorLocal; %Estimation error
cd(curDir)

%% Global evaluation to select the best techniques                                                          
for i = 1:size(CountOutsideJoint,1)
    if  CountOutsideJoint(i,1) == CountOutsideSeparat(i,1)
        countFinal(i,1) = 0;
    elseif CountOutsideJoint(i,1) > CountOutsideSeparat(i,1)
        countFinal(i,1) = 1;                                                %   voting separate is coded as 1
    else
        countFinal(i,1) = 2;                                                %   voting joint is coded as 2
    end
end
countFinalIni = countFinal;
countFinal = countFinal(find(countFinal~=0));
votesSeparate = sum(countFinal==1)/size(countFinal,1);
votesJoint = sum(countFinal==2)/size(countFinal,1);
display(['Separate-technique Best coverage: ', num2str(votesSeparate*100)])
display(['Joint-technique Best coverage: ', num2str(votesJoint*100)])
display('//////////////////')
ComparisonTesting.SeparateTechniqueCoverage= votesSeparate*100; 
ComparisonTesting.JointTechniqueCoverage = votesJoint*100;

%% Global performances
errorGlobal = (errorSeparat-errorJoint)./(errorSeparat+errorJoint);
errorTotal = sum(errorGlobal,2);
errorTotal = errorTotal(2:end);

beta = 0.9;
beta2 = 0.4;
x_t = 0;
x_t2 = 0;

errorTotalnew = [];
errorTotalnew2 = [];
for i=1:size(errorTotal,1)
x_t = beta *x_t + (1- beta)*errorTotal(i);
x_t2 = beta2 *x_t2 + (1- beta2)*errorTotal(i);
errorTotalnew = [errorTotalnew x_t];
errorTotalnew2 = [errorTotalnew2 x_t2];
end

threshold = mean(errorTotal)-0.3;                                              % Threshold for selecting best technique
globalPerformancejoint = sum(errorTotal>threshold)/size(errorTotal,1) ;
globalPerformanceSeparate = sum(errorTotal<-threshold)/size(errorTotal,1);
globalEqualPerformance = 1 - (globalPerformancejoint + globalPerformanceSeparate);
globalPerformance = mean(errorTotal); 

display(['Separate-technique best performance: ', num2str(globalPerformanceSeparate*100)])
display(['Joint-technique best performance: ', num2str(globalPerformancejoint*100)])
display(['Equal performance for both techniques: ', num2str(globalEqualPerformance*100)])
display(['Average performance comparison (Threshold): ', num2str(globalPerformance)])
ComparisonTesting.SeparateTechniquePerformance = globalPerformanceSeparate*100;
ComparisonTesting.JointTechniquePerformance = globalPerformancejoint*100;
ComparisonTesting.EqualPerformanceForBothTechniques = globalEqualPerformance*100;

cd(dirComparisonResultsTest)
if plotting == true
%% Plotting: space coverage
% % % h1 = figure;
% % % hold on
% % % xlab = xlabel('Time instants ($k$)','interpreter','latex');
% % % xlab.FontSize = 16;
% % % ylab = ylabel('$\% N_{Out}$','interpreter','latex');
% % % ylab.FontSize = 16;
% % % plot1 = plot(CountOutsideJoint,'r','Linewidth',1.3);
% % % plot2 = plot(CountOutsideSeparat,'b','Linewidth',1.3);
% % % h1.Position = [1693         571        1663         365];
% % % axis([0  size(CountOutsideSeparat,1) 0 1])
% % % grid on
% % % lg = legend([plot1 plot2],{'Joint case','Separate case'});
% % % lg.FontSize = 12;
% % % lg.Position = [0.7888    0.7694    0.1341    0.1274];
% % % title('Comparison between state space coverage','interpreter','latex','FontSize',18)
% % % ax = gca;
% % % outerpos = ax.OuterPosition;
% % % ti = ax.TightInset; 
% % % left = outerpos(1) + ti(1);
% % % bottom = outerpos(2) + ti(2);
% % % ax_width = outerpos(3) - ti(1) - ti(3);
% % % ax_height = outerpos(4) - ti(2) - ti(4);
% % % ax.Position = [left bottom ax_width ax_height];
% % % set(h1,'Units','Inches');
% % % pos = get(h1,'Position');
% % % set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% In average, in each iteration, CountOutsideSeparatMean*N particles are 
%outside of the model for the seprate experiment,
%whereas CountOutsideSeparatMean*N particles go outside in the Joint experiment
%% Plotting: Estimation error
h2 = figure;
hold on
xlab = xlabel('Time instants ($k$)','interpreter','latex');
xlab.FontSize = 16;
ylab = ylabel('Error comparison ($E$)','interpreter','latex');
ylab.FontSize = 16;

rec1 = rectangle('Position',[0 -3 size(errorTotal,1)  -threshold+3],'Curvature',0);
rec1.EdgeColor = 'none';
rec1.FaceColor = [0 1 0 0.1];

rec2 = rectangle('Position',[-3 -threshold  size(errorTotal,1) 2*threshold],'Curvature',0);
rec2.EdgeColor = 'none';
rec2.FaceColor = [1 0 0 0.1];

rec3 = rectangle('Position',[-3 threshold  size(errorTotal,1) 4-threshold],'Curvature',0);
rec3.EdgeColor = 'none';
rec3.FaceColor = [1 1 0 0.1];
title('Comparison between error performances','interpreter','latex','FontSize',18)

plot2 = plot([1:size(errorTotal,1)],ones(size(errorTotal,1),1)*threshold,'--r','LineWidth',2);

plot3 = plot([1:size(errorTotal,1)],ones(size(errorTotal,1),1)*threshold*-1,'--r','LineWidth',2);

h2a = plot(errorTotalnew2,'b', 'LineWidth',12); h2a.Color(4)=0.1;  % 70% transparent
hold on
plot1 = plot(errorTotalnew,'b');

% % plot1 = plot(errorTotal,'b');

h2.Position = [1693         571        1663         365];
axis([0  size(errorTotal,1) -3 4])
grid on
lg = legend([plot1 plot2],{'Error comparison','Equal performance area'});
lg.FontSize = 12;
lg.Position = [0.7888    0.7694    0.1341    0.1274];

tx1 = text( 5,-2.4,'High performance for separate clusters');
tx1.FontSize = 13;
tx1.Interpreter = 'latex';

tx2 = text( 5,3.2,'High performance for joint clusters');
tx2.FontSize = 13;
tx2.Interpreter = 'latex';
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
set(h2,'Units','Inches');
pos = get(h2,'Position');
set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

%% Saving plots
% % % print(h1,'-djpeg','-r300',['coverageT' num2str(scenarioTrain) 'Par' num2str(numberParticle) 'C' num2str(Case)])
% % % %print(h1,['coverageT' num2str(scenarioTrain) 'Par' num2str(numberParticle) 'C' num2str(Case)],'-dpdf','-r0')
% % % export_fig(h1,['coverageT' num2str(scenarioTrain) 'Par' num2str(numberParticle) 'C' num2str(Case)],'-pdf','-transparent')
print(h2,'-djpeg','-r300',['PerformanceTMB3' num2str(scenarioTrain) 'Par' num2str(numberParticle) 'C' num2str(Case)])
print(h2,['PerformanceTMB3' num2str(scenarioTrain) 'Par' num2str(numberParticle) 'C' num2str(Case)],'-dpdf','-r0')
% % export_fig(h2,['PerformanceTMB' num2str(scenarioTrain) 'Par' num2str(numberParticle) 'C' num2str(Case)],'-pdf','-transparent')
end

%% Saving results
save(['ComparisonTestingT' num2str(scenarioTrain) 'Par' num2str(numberParticle)],'ComparisonTesting'); % T: Training
cd(curDir)