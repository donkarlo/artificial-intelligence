%   Third algorithm: Generation of vocabularies and its properties
set(0,'defaultfigurecolor',[1 1 1])
clc
clear
close all
curDir = pwd;

%%
MB = true;
scenarioTrain = 1; % 1 (PM): training data
plotflagGNG = true;
plotMap = true;
Case = 1; % 1: Less numbers, 2: Best case, 3: High numbers
%% Modalities
odometric = true;
if   odometric == true
    control = false;
else
    control = true;
end
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
for counter = 1:1
    %     counter
    if odometric == true
        %% GNG parameters :Best case/Less numbers/High numbers
        if Case == 1 %Less numbers
            params.L_growing = 1200;
            params.L_decay =  40;
            params.k = 0.6;
        elseif Case == 2 %Best case
            params.L_growing = 750;
            params.L_decay =  90;
            params.k = 0.6;
        elseif Case == 3%High numbers
            params.L_growing = 700;
            params.L_decay =  200;
            params.k = 100;
        end
        params.N = 300;% 300
        params.MaxIt = 5; %5 %5%
        %params.L_growing = 700; %750 %1200 % 700
        params.epsilon_b = 0.05*1.5; %0.05*1.5
        params.epsilon_n = 0.0006*1.5; %0.0006*1.5
        params.alpha = 0.5; %0.5
        params.delta = 0.995; %0.995
        params.T = 100; %100                                                            %    It could be a function of params.L_growing, e.g., params.LDecay = 2*params.L_growing
        %         params.L_decay =  200; %90 %40 %200
        %         params.k = 100; %0.6%0.6%100
        params.seedvector = 1;% 1
        
        %% GNG processing
        inputData = [structSyncData.Filtered.xPos,structSyncData.Filtered.yPos];
        inputData = inputData(100:end-100,:);
        inputDataDeriv = [structSyncData.Filtered.divxPos,structSyncData.Filtered.divyPos];
        inputDataDeriv = inputDataDeriv(100:end-100,:);
        
        compData = inputDataDeriv;
        cd(curDir)
        netOdometric = GrowingNeuralGasNetwork(inputData, compData, params, plotflagGNG);
        %         timeLocal = netOdometric.timeLocal;
        timeGlobal(counter,1) = netOdometric.timeGlobal;
        mycolors = colorcube;
        mycolors = [mycolors;mycolors];
        if plotMap == true
            h = figure;
            hold on
            scatter(inputData(:,1),inputData(:,2),60,mycolors(netOdometric.dataColorNode,:),'.')
            scatter(netOdometric.w(:,1),netOdometric.w(:,2),250,'+','k','linewidth',2)
            h1 = figure;
            hold on
            scatter(inputData(:,1),inputData(:,2),60,mycolors(netOdometric.dataColorNode,:),'.')
            scatter(netOdometric.w(:,1),netOdometric.w(:,2),200'.','k','linewidth',2)
            quiver(netOdometric.w(:,1), netOdometric.w(:,2), netOdometric.nodesMeanComp(:,1),netOdometric.nodesMeanComp(:,2),'k','Autoscale','off','LineWidth',1)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%  Caculation of transition matrices
        nodesInTime = netOdometric.dataColorNode;
        %%  Transition matrix
        transitionMatChanges = zeros(netOdometric.N,netOdometric.N);
        transitionMat = zeros(netOdometric.N,netOdometric.N);
        timesSpent = [];
        currLength = length(inputData);
        trackDicrete = nodesInTime(1:currLength);
        ind = find(diff(trackDicrete) ~= 0);
        for k = 1:size(ind,1)
            transitionMatChanges(trackDicrete(ind(k),1),trackDicrete(ind(k)+1,1)) =...
                transitionMatChanges(trackDicrete(ind(k),1),trackDicrete(ind(k)+1,1)) + 1;
        end
        %%%%%
        for k = 1:currLength-1
            transitionMat(trackDicrete(k,1),trackDicrete(k+1,1)) =...
                transitionMat(trackDicrete(k,1),trackDicrete(k+1,1)) + 1;
        end
        %%%%%
        codeInd = [0; ind];
        tspentTran = diff(codeInd);
        %%%
        for k = 1:size(tspentTran,1)
            if size(unique([timesSpent;tspentTran(k,1)]),1) ~= size(unique(timesSpent),1)
                timeMats{1,tspentTran(k)} = zeros(netOdometric.N,netOdometric.N);
                timesSpent = [timesSpent; tspentTran(k)];
            end
            timeMats{1,tspentTran(k)}(trackDicrete(ind(k),1),trackDicrete(ind(k)+1,1)) =...
                timeMats{1,tspentTran(k)}(trackDicrete(ind(k),1),trackDicrete(ind(k)+1,1)) + 1;
        end
        %%%%%
        ind2 = find(diff(trackDicrete) == 0);
        tspentSame = 1;
        for k = 1:size(ind2,1)
            if k > 1
                if ind2(k) == ind2(k-1) + 1
                    tspentSame = tspentSame + 1;
                else
                    tspentSame = 1;
                end
            end
            
            if size(unique([timesSpent;tspentSame]),1) ~= size(unique(timesSpent),1)
                timeMats{1,tspentSame} = zeros(netOdometric.N,netOdometric.N);
                timesSpent = [timesSpent; tspentSame];
            end
            
            timeMats{1,tspentSame}(trackDicrete(ind2(k),1),trackDicrete(ind2(k)+1,1)) =...
                timeMats{1,tspentSame}(trackDicrete(ind2(k),1),trackDicrete(ind2(k)+1,1)) + 1;
        end
        
        transitionMatChanges = transitionMatChanges./repmat(sum(transitionMatChanges,2) + (sum(transitionMatChanges,2)==0),1,netOdometric.N);
        transitionMat = transitionMat./repmat(sum(transitionMat,2) + (sum(transitionMat,2)==0),1,netOdometric.N);
        
        transitionInfoOdometric.transitionMat = transitionMat;
        transitionInfoOdometric.TimeMats = timeMats;
        transitionInfoOdometric.transitionMatChanges = transitionMatChanges;                 %   We do not use this in MPJF
        if plotflagGNG
            figure
            imagesc(transitionMat);
            colorbar;
            colormap jet
            title('Transtion matrix for superstates')
        end
        trainingTracksOdometric.tracksDiscrete = trackDicrete;
        cd(dirResultsTrain)
        save('VocabulariesFinalOdometric.mat','netOdometric','transitionInfoOdometric','trainingTracksOdometric')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif control == true
        %% GNG processing
        %% GNG parameters   Best case/Less numbers/High numbers
        if Case == 1 %Less numbers
            params.L_growing = 750;
            params.L_decay =  100;
            params.k = 1.1;
        elseif Case == 2 %Best case
            params.L_growing = 750;
            params.L_decay =  100;
            params.k = 3;
        elseif Case == 3%High numbers
            params.L_growing = 700;
            params.L_decay =  100;
            params.k = 100;
        end
        params.N = 300;% 300
        params.MaxIt = 5; %5
        %params.L_growing = 700; %750%750% 700
        params.epsilon_b = 0.05*1.5;%0.05*1.5
        params.epsilon_n = 0.0006*1.5;% 0.0006*1.5
        params.alpha = 0.5;%0.5
        params.delta = 0.995;%0.995
        params.T = 100;%100                                                             %    It could be a function of params.L_growing, e.g., params.LDecay = 2*params.L_growing
        %params.L_decay =  100; %100
        %params.k = 100; %3 % 1.1%100
        params.seedvector = 1; %1
        
        %% GNG processing
        inputData = [structSyncData.Filtered.S,structSyncData.Filtered.V];
        inputData = inputData(100:end-100,:);
        inputDataDeriv = [structSyncData.Filtered.divS,structSyncData.Filtered.divV];
        inputDataDeriv = inputDataDeriv(100:end-100,:);
        compData = inputDataDeriv;
        cd(curDir)
        netcontrol = GrowingNeuralGasNetwork(inputData, compData, params, plotflagGNG);
        timeGlobal(counter,1) = netcontrol.timeGlobal;
        % %         timeLocal = netOdometric.timeLocal;
        mycolors = colorcube;
        mycolors = [mycolors;mycolors];
        if plotMap == true
            h = figure;
            hold on
            scatter(inputData(:,1),inputData(:,2),60,mycolors(netcontrol.dataColorNode,:),'.')
            scatter(netcontrol.w(:,1),netcontrol.w(:,2),250,'+','k','linewidth',2)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%  Caculation of transition matrices
        nodesInTime = netcontrol.dataColorNode;
        %%  Transition matrix
        transitionMatChanges = zeros(netcontrol.N,netcontrol.N);
        transitionMat = zeros(netcontrol.N,netcontrol.N);
        timesSpent = [];
        currLength = length(inputData);
        trackDicrete = nodesInTime(1:currLength);
        ind = find(diff(trackDicrete) ~= 0);
        for k = 1:size(ind,1)
            transitionMatChanges(trackDicrete(ind(k),1),trackDicrete(ind(k)+1,1)) =...
                transitionMatChanges(trackDicrete(ind(k),1),trackDicrete(ind(k)+1,1)) + 1;
        end
        %%%%%
        for k = 1:currLength-1
            transitionMat(trackDicrete(k,1),trackDicrete(k+1,1)) =...
                transitionMat(trackDicrete(k,1),trackDicrete(k+1,1)) + 1;
        end
        %%%%%
        codeInd = [0; ind];
        tspentTran = diff(codeInd);
        %%%
        for k = 1:size(tspentTran,1)
            if size(unique([timesSpent;tspentTran(k,1)]),1) ~= size(unique(timesSpent),1)
                timeMats{1,tspentTran(k)} = zeros(netcontrol.N,netcontrol.N);
                timesSpent = [timesSpent; tspentTran(k)];
            end
            timeMats{1,tspentTran(k)}(trackDicrete(ind(k),1),trackDicrete(ind(k)+1,1)) =...
                timeMats{1,tspentTran(k)}(trackDicrete(ind(k),1),trackDicrete(ind(k)+1,1)) + 1;
        end
        %%%%%
        ind2 = find(diff(trackDicrete) == 0);
        tspentSame = 1;
        for k = 1:size(ind2,1)
            if k > 1
                if ind2(k) == ind2(k-1) + 1
                    tspentSame = tspentSame + 1;
                else
                    tspentSame = 1;
                end
            end
            
            if size(unique([timesSpent;tspentSame]),1) ~= size(unique(timesSpent),1)
                timeMats{1,tspentSame} = zeros(netcontrol.N,netcontrol.N);
                timesSpent = [timesSpent; tspentSame];
            end
            
            timeMats{1,tspentSame}(trackDicrete(ind2(k),1),trackDicrete(ind2(k)+1,1)) =...
                timeMats{1,tspentSame}(trackDicrete(ind2(k),1),trackDicrete(ind2(k)+1,1)) + 1;
        end
        
        transitionMatChanges = transitionMatChanges./repmat(sum(transitionMatChanges,2) + (sum(transitionMatChanges,2)==0),1,netcontrol.N);
        transitionMat = transitionMat./repmat(sum(transitionMat,2) + (sum(transitionMat,2)==0),1,netcontrol.N);
        
        transitionInfocontrol.transitionMat = transitionMat;
        transitionInfocontrol.TimeMats = timeMats;
        transitionInfocontrol.transitionMatChanges = transitionMatChanges;                 %   We do not use this in MPJF
        if plotflagGNG
            figure
            imagesc(transitionMat);
            colorbar;
            colormap jet
            title('Transtion matrix for superstates')
        end
        trainingTrackscontrol.tracksDiscrete = trackDicrete;
        cd(dirResultsTrain)
        save('VocabulariesFinalControl.mat','netcontrol','transitionInfocontrol','trainingTrackscontrol')
    end
    cd(curDir)
end
%export_fig rocfinale2.pdf -pdf -transparent
%print(figure(3),'-djpeg','-r300','BestOdm')
