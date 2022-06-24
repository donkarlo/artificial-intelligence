function estimationAbn = MJPF_NewJoint(odometric, control,...
    transitionMat,codebook,curDir, scenarioTesting,MB,plotBool,dirResultsTesting,scenarioTrain,Case,N)
nomalizationBool = true;                                                   %   False when testing with training data
if MB == true
    dirMatTesting = ['C:\Users\mohamad.baydoun\Desktop\Interaction journal (2019)',...
        '\Generated Data\' num2str(scenarioTesting)];
else
    dirMatTesting = ['C:\Users\damian.campo\Desktop\Interaction journal (2019)',...
        '\Generated Data\' num2str(scenarioTesting)];
end
cd(dirMatTesting)
load('genSyncdata.mat')
dataTestOdometric = [structSyncData.Filtered.xPos,...
    structSyncData.Filtered.yPos,structSyncData.Filtered.divxPos,...
    structSyncData.Filtered.divyPos]';
dataTestOdometric = dataTestOdometric(:,200:925);
dataTestcontrol = [structSyncData.Filtered.S,structSyncData.Filtered.V,...
    structSyncData.Filtered.divS,structSyncData.Filtered.divV]';
dataTestcontrol = dataTestcontrol(:,200:925);

%% Definition of the parameters (training)
% Transition matrix for the continous-time system.
numState = 2;
A = [eye(numState),zeros(numState);zeros(numState),zeros(numState)];
%   Measurement model.
H = [eye(numState),zeros(numState)];
%   Control input
B = [eye(numState);eye(numState)];
% Variance in the measurements.
r1 = 1e-2;
r2 = 1e-8;
R = eye(numState)*r1*2;
R1 = eye(2*numState)*r1*2;

DatasetSize = length(dataTestOdometric);
vocSizeOdm = length(odometric.averageStateOdm);
vocSizeCont = length(control.averageStateCon);
dictionarySize = length(transitionMat);
% Radius of superstates
radiusOdm = odometric.radiusStateOdm;
radiusCont = control.radiusStateCon;
% Mean and covariance
UOdm = odometric.averageStateOdm(:,1:4);
UCont = control.averageStateCon(:,5:8);
covarianceOdm = odometric.CovarianceOdm;
covarianceCont = control.CovarianceCont;
% Normalization values
if nomalizationBool == true
    if MB == true
        dirResultsTestingNorm = ['C:\Users\mohamad.baydoun\Desktop',...
            '\Interaction journal (2019)\Jointly Technique',...
            '\ResultsTesting\' num2str(scenarioTrain) '\C' num2str(Case)];
    else
        dirResultsTestingNorm = ['C:\Users\damian.campo\Desktop',...
            '\Interaction journal (2019)\Jointly Technique',...
            '\ResultsTesting\' num2str(scenarioTesting) '\C' num2str(Case)];
    end
    cd(dirResultsTestingNorm)
    load(['EstimationAbnInteractionT' num2str(scenarioTrain) 'Par40.mat'])%num2str(N)
    meanOdmAbn1 = mean(EstimationAbn.estimationAbn.abnormdb1Odm_Probabilistic);
    minOdmAbn1 = min(EstimationAbn.estimationAbn.abnormdb1Odm_Probabilistic);
    stdOdmAbn1 = std(EstimationAbn.estimationAbn.abnormdb1Odm_Probabilistic);
    normOdmAbn1 = meanOdmAbn1 + 2*stdOdmAbn1;
    normOdmAbn1 = normOdmAbn1 - minOdmAbn1;
    
    meanContAbn1 = mean(EstimationAbn.estimationAbn.abnormdb1Cont_Probabilistic);
    minContAbn1 = min(EstimationAbn.estimationAbn.abnormdb1Cont_Probabilistic);
    stdContAbn1 = std(EstimationAbn.estimationAbn.abnormdb1Cont_Probabilistic);
    normContAbn1 = meanContAbn1 + 2* stdContAbn1;
    normContAbn1 = normContAbn1 - minContAbn1;
    
    minOdmAbn2 = min(EstimationAbn.estimationAbn.abnormdb2Odm_Probabilistic);
    meanOdmAbn2 = mean(EstimationAbn.estimationAbn.abnormdb2Odm_Probabilistic);
    stdOdmAbn2 = std(EstimationAbn.estimationAbn.abnormdb2Odm_Probabilistic);
    normOdmAbn2 = meanOdmAbn2 + 2*stdOdmAbn2;
    normOdmAbn2 = normOdmAbn2 - minOdmAbn2;
    
    meanContAbn2 = mean(EstimationAbn.estimationAbn.abnormdb2Cont_Probabilistic);
    minContAbn2 = min(EstimationAbn.estimationAbn.abnormdb2Cont_Probabilistic);
    stdContAbn2 = std(EstimationAbn.estimationAbn.abnormdb2Cont_Probabilistic);
    normContAbn2 = meanContAbn2 + 2* stdContAbn2;
    normContAbn2 = normContAbn2 - minContAbn2;
    
else
    normOdmAbn1 = 1;
    normContAbn1 = 1;
    minOdmAbn1 = 0;
    minContAbn1 = 0;

    normOdmAbn2 = 1;
    normContAbn2 = 1;
    minOdmAbn2 = 0;
    minContAbn2 = 0;
    
    stdOdmAbn1 = 0.3;
    stdContAbn1 = 0.3;    
    stdOdmAbn2 = 0.3;
    stdContAbn2 = 0.3;    
end
cd(curDir)
addpath('./ekfukf')

%% Initialization Filtering steps.
% Initial guesses for the state mean and covariance.
Pinit = eye(numState*2)*r2;
Q = eye(numState*2)*r1*100;                                                % Prediction noise equal for each superstate
weightsPar = zeros(1,N);                                                   % Weight of particles (1x100 matrix that contains all zeros)

%% Definition of parameters (MJPF)
statepredOdm = zeros(4,DatasetSize,N);                                        % predict state 4D Odm
PpredOdm = zeros(4,4,DatasetSize,N);                                          % predict covariance matrix Odm
statepredCont = zeros(4,DatasetSize,N);                                       % predict state 4D Cont
PpredCont = zeros(4,4,DatasetSize,N);                                         % predict covariance matrix cont
stateUpdatedOdm = zeros(4,DatasetSize,N);
updatedPOdm = zeros(4,4,DatasetSize,N);
stateUpdatedCont = zeros(4,DatasetSize,N);                                   % stato aggiornato
updatedPCont = zeros(4,4,DatasetSize,N);                                     % P aggiornato
weightscoeff = zeros(N,DatasetSize);                                         % pesi per ogni particle
sOdm = zeros(N,DatasetSize);                                                 % superstato del Odm
sCont = zeros(N,DatasetSize);                                                % superstato del Cont
d1Odm_Probabilistic = zeros(N,DatasetSize-1);
d1Cont_Probabilistic = zeros(N,DatasetSize-1);
% d1Overall_Probabilistic = zeros(N,DatasetSize);
abnormdb1Odm_Probabilistic = zeros(1,DatasetSize-1);
abnormdb1Cont_Probabilistic = zeros(1,DatasetSize-1);
% abnormdb1Overall_Probabilistic = zeros(1,DatasetSize);

% % d1Odm_Deterministic = zeros(N,DatasetSize);
% % d1Cont_Deterministic = zeros(N,DatasetSize);
% % d1Overall_Deterministic = zeros(N,DatasetSize);
% % abnormdb1Odm_Deterministic = zeros(1,DatasetSize);
% % abnormdb1Cont_Deterministic = zeros(1,DatasetSize);
% % abnormdb1Overall_Deterministic = zeros(1,DatasetSize);
indexabnormdb1Odm = zeros(N,DatasetSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%Abnormality measurement at continuous level
d2Odm_Probabilistic = zeros(N,DatasetSize-1);
d2Cont_Probabilistic = zeros(N,DatasetSize-1);
abnormdb2Odm_Probabilistic = zeros(1,DatasetSize-1);
abnormdb2Cont_Probabilistic = zeros(1,DatasetSize-1);
errorLocal = zeros(DatasetSize,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
activeLabels = zeros(N,DatasetSize);
firstIt = 1;
cd(dirResultsTesting)
v1 = VideoWriter(['Abnormality detectionT' num2str(scenarioTrain) 'Par' num2str(N) '.avi']); %   Initialize a video data
open (v1)
cd(curDir)
h = figure;
h.Position = [561    34   560   931];

%% main loop
%Initialization of states/update KFs' info
for i = 1:DatasetSize
    counOutside = 0;
    for n = 1:N
        if i == firstIt
            currStatePos = awgn(dataTestOdometric(1:2,i),50);
            currStateVel = awgn(dataTestOdometric(3:4,i),100);
            currStateOdm = [currStatePos;currStateVel];
            currStateSV = awgn(dataTestcontrol(1:2,i),50);
            currStateDiSV = awgn(dataTestcontrol(3:4,i),100);
            currStateCont = [currStateSV;currStateDiSV];
            currPOdm = Pinit;
            currPCont = Pinit;
            stateUpdatedOdm(:,i,n) = currStateOdm;
            stateUpdatedCont(:,i,n) = currStateCont;
            updatedPOdm(:,:,i,n) = currPOdm;
            updatedPCont(:,:,i,n) = currPCont;
        else
            %   NEXT MEASUREMENT APPEARS
            % UPDATE
            [stateUpdatedOdm(:,i,n),updatedPOdm(:,:,i,n),~] =...
                kf_update(statepredOdm(:,i-1,n),PpredOdm(:,:,i-1,n),...
                dataTestOdometric(1:2,i),H,R);%update Odm
            [stateUpdatedCont(:,i,n),updatedPCont(:,:,i,n),~] =...
                kf_update(statepredCont(:,i-1,n),PpredCont(:,:,i-1,n),...
                dataTestcontrol(1:2,i),H,R);%update Cont
            %   Association of updated states to variables
            currStateOdm = stateUpdatedOdm(:,i,n);
            currStateCont = stateUpdatedCont(:,i,n);
            currPOdm = updatedPOdm(:,:,i,n);
            currPCont = updatedPCont(:,:,i,n);
        end
        
        %% Calculation of current superstates
        [probDistOdm,~] =...
            nodesProb_New(currStateOdm(1:2,1)',UOdm(:,1:2),radiusOdm);  %sum(radiusOdm(1:2,:))/2      %   Find probability of being in each superstate of voc1
        [probDistCont,~] =...
            nodesProb_New(currStateCont(1:2,1)',UCont(:,1:2),radiusCont);%sum(radiusCont(1:2,:))/2
        probDistSortOdm = [];
        propIndSortOdm = [];
        probDistDynamic = probDistOdm;
        for jj =1:size(find(probDistOdm),2)
            probDistT = makedist('Multinomial','Probabilities',probDistDynamic);             %   Multinomial distribution to pick multiple likely particles
            indC = probDistT.random(1,1);
            probDistSortOdm = [probDistSortOdm probDistOdm(1,indC)];
            propIndSortOdm = [propIndSortOdm indC];
            probDistDynamic(1,indC) = 0;
            probDistDynamic = probDistDynamic/sum(probDistDynamic);
        end
        probDistSortOdm = [probDistSortOdm zeros(1,size(probDistOdm,2)- size(find(probDistOdm),2))];
        propIndSortOdm = [propIndSortOdm, setdiff(1:size(probDistOdm,2),propIndSortOdm) ];
        probDistSortCont = [];
        propIndSortCont = [];
        probDistDynamic = probDistCont;
        for jj =1:size(find(probDistCont),2)
            probDistT = makedist('Multinomial','Probabilities',probDistDynamic);             %   Multinomial distribution to pick multiple likely particles
            indC = probDistT.random(1,1);
            probDistSortCont = [probDistSortCont probDistCont(1,indC)];
            propIndSortCont = [propIndSortCont indC];
            probDistDynamic(1,indC) = 0;
            probDistDynamic = probDistDynamic/sum(probDistDynamic);
        end
        probDistSortCont = [probDistSortCont zeros(1,size(probDistCont,2)- size(find(probDistCont),2))];
        propIndSortCont = [propIndSortCont, setdiff(1:size(probDistCont,2),propIndSortCont) ];
        %   Normalization of nodes and current states (Odometry)
        currStateNormOdm = currStateOdm(1:2,1)' - [odometric.minMAxNormOdm(1,1) odometric.minMAxNormOdm(3,1)];
        currStateNormOdm = currStateNormOdm./([odometric.minMAxNormOdm(2,1) odometric.minMAxNormOdm(4,1)]);
        UOdmNorm = UOdm(:,1:2) - repmat([odometric.minMAxNormOdm(1,1) odometric.minMAxNormOdm(3,1)],size(UOdm,1),1);
        UOdmNorm = UOdmNorm./repmat([odometric.minMAxNormOdm(2,1) odometric.minMAxNormOdm(4,1)],size(UOdm,1),1);
        distancesNormOdm = pdist2(currStateNormOdm,UOdmNorm);                                  % Calculation of norm 1 distance from current positions to each neuron
        %   Normalization of nodes and current states (Control)
        currStateNormCont = currStateCont(1:2,1)' - [control.minMAxNormControl(1,1) control.minMAxNormControl(3,1)];
        currStateNormCont = currStateNormCont./([control.minMAxNormControl(2,1) control.minMAxNormControl(4,1)]);
        UControlNorm = UCont(:,1:2) - repmat([control.minMAxNormControl(1,1) control.minMAxNormControl(3,1)],size(UCont,1),1);
        UControlNorm = UControlNorm./repmat([control.minMAxNormControl(2,1) control.minMAxNormControl(4,1)],size(UCont,1),1);
        distancesNormCont = pdist2(currStateNormCont,UControlNorm);                                  % Calculation of norm 1 distance from current positions to each neuron
        %   TRANSFORMATION OF STATES INTO SUPERSTATES
        zeroProbBool = false;
        outsideBool = false;
        i1 = 1;
        i2 = 1;
        %   TRANSFORMATION OF STATES INTO SUPERSTATES
        while(true)
            emptyNBool(i,n) = false;
            if i1 > vocSizeOdm || i2 > vocSizeCont
                %   All combinations are tried and no label is identified
                emptyNBool(i,n) = true;
                zeroProbBool = false;
                outsideBool = true;
                break
            end
            ind1Curr = propIndSortOdm(i1);
            ind2Curr = propIndSortCont(i2);
            if ind1Curr > vocSizeOdm || ind2Curr > vocSizeCont
                %   The closest label to the current state is the empty neuron
                emptyNBool(i,n) = true;
                zeroProbBool = false;
                outsideBool = true;
                break
            end
            codeEval = codebook((ind1Curr-1)*vocSizeCont+ind2Curr,4) ~= 0;
            s1Curr = probDistSortOdm(i1);
            s2Curr = probDistSortCont(i2);
            if codeEval == true
                sOdm(n,i) = ind1Curr;
                sCont(n,i) = ind2Curr;
                minDistLabel = codebook((sOdm(n,i)-1)*...
                    vocSizeCont + sCont(n,i),4);                        % look if this couple has been observed
                break
            else
                sTot  = [s1Curr + probDistSortCont(i2+1),...
                    s2Curr + probDistSortOdm(i1+1)];
                [~,sIndMax] = max(sTot);
                if sIndMax == 1
                    i2 = i2+1;
                else
                    i1 = i1+1;
                end
            end
            if s1Curr == 0 || s2Curr == 0
                %   The probability of at least one vocabulary is null
                emptyNBool(i,n) = true;
                zeroProbBool = true;
                outsideBool = false;
                break
            end
        end
        %   find minimum sum
        if zeroProbBool == 1 && outsideBool == 0
            validCodebook = codebook(find(codebook(:,4)~=0),:);
            disOdmCodebook = distancesNormOdm(validCodebook(:,1));
            disContCodebook = distancesNormCont(validCodebook(:,2));
            [~, minDistLabel] = min(disOdmCodebook + disContCodebook);
        elseif zeroProbBool == 0 && outsideBool == 1
            validCodebook = codebook(find(codebook(:,4)~=0),:);
            disOdmCodebook = distancesNormOdm(validCodebook(:,1));
            disContCodebook = distancesNormCont(validCodebook(:,2));
            disProbCodebook = disOdmCodebook + disContCodebook;
            disProbCodebook = disProbCodebook/sum(disProbCodebook);
            probDistT = makedist('Multinomial','Probabilities',disProbCodebook);             %   Multinomial distribution to pick multiple likely particles
            minDistLabel = probDistT.random(1,1);
            counOutside = counOutside + 1;
        end
        
        activeLabels(n,i) = minDistLabel;
        weightscoeff(n,i) = 1/N;                                            %   Weight of particles
        indActiveLab = find(codebook(:,4) == activeLabels(n,i));
        nodeIndOdm = codebook(indActiveLab,1);
        nodeIndCont = codebook(indActiveLab,2);
        indexabnormdb1Odm(n,i) = nodeIndOdm;
        if i > firstIt
            %   CALCULATION OF ABNORMALITY MEASUREMENTS
            %% Bhattacharyya distance (Discrete level)
            d1Odm_Probabilistic(n,i-1) = max([(bhattacharyyadistance(statepredOdm(:,i-1,n)',...
                UOdm(nodeIndOdm,:),PpredOdm(:,:,i-1,n),...                   %   measure bhattacharrya distance between p(xk/zk-1) and p(xk/sk) object 1
                positivedefinite(covarianceOdm(1:4,1:4,nodeIndOdm))) - minOdmAbn1)/normOdmAbn1,0]);
            d1Cont_Probabilistic(n,i-1)  = max([(bhattacharyyadistance(statepredCont(:,i-1,n)',...
                UCont(nodeIndCont,:),PpredCont(:,:,i-1,n),...                % measure bhattacharrya distance between p(xk/zk-1) and p(xk/sk) object 2
                positivedefinite(covarianceCont(5:8,5:8,nodeIndCont)))- minContAbn1)/normContAbn1,0]);          
            %% Bhattacharyya distance (Continuous level)
            d2Odm_Probabilistic(n,i-1) = max([(bhattacharyyadistance(statepredOdm(:,i-1,n)',...
                dataTestOdometric(:,i)',PpredOdm(:,:,i-1,n),R1) - minOdmAbn2)/normOdmAbn2,0]);
            d2Cont_Probabilistic(n,i-1) = max([(bhattacharyyadistance(statepredCont(:,i-1,n)',...
                dataTestcontrol(:,i)',PpredCont(:,:,i-1,n),R1) - minContAbn2)/normContAbn2,0]);
            
            error(n,1) = pdist2(statepredOdm(1,i-1,n),dataTestOdometric(1,i))/abs(statepredOdm(1,i-1,n) + dataTestOdometric(1,i));
            error(n,2) = pdist2(statepredOdm(2,i-1,n),dataTestOdometric(2,i))/abs(statepredOdm(2,i-1,n) + dataTestOdometric(2,i));
            error(n,3) = pdist2(statepredCont(1,i-1,n),dataTestcontrol(1,i))/abs(statepredCont(1,i-1,n) + dataTestcontrol(1,i));
            error(n,4) = pdist2(statepredCont(2,i-1,n),dataTestcontrol(2,i))/abs(statepredCont(2,i-1,n) + dataTestcontrol(2,i));
            weightsPar(n) = weightscoeff(n,i-1)/sum(error(n,:));
            
            % RESAMPLING OF PARTICLES
            if n == N
                [~,minParticle] = min(sum(error,2));
                errorLocal(i,:) = error(minParticle,:);
                
                indexAbnormdb1Odm = indexabnormdb1Odm(minParticle,i);
                abnormdb1Odm_Probabilistic(i-1) = d1Odm_Probabilistic(minParticle,i-1);
                abnormdb1Cont_Probabilistic(i-1) = d1Cont_Probabilistic(minParticle,i-1);               
                %%%%%%%%%%%%%%%%%%%%%%Abnormality measurement at continuous level
                abnormdb2Odm_Probabilistic(i-1) = d2Odm_Probabilistic(minParticle,i-1);
                abnormdb2Cont_Probabilistic(i-1) = d2Cont_Probabilistic(minParticle,i-1);
                %%%%%%%%%%%%%%%%%%%%%%                        
                weightsPar = weightsPar/sum(weightsPar);                                               %   Normalize weights in such a way that they all sum 1
                weightscoeff(:,i)= weightsPar';
                pd = makedist('Multinomial','Probabilities',weightsPar);             %   Multinomial distribution to pick multiple likely particles
                wRes = pd.random(1,N);
                for ij = 1:N
                    % REPLACEMENT OF CORRECTED DATA DEPENDING ON
                    % SURVIVING NEURONS
                    stateUpdatedOdm(:,i,ij) = stateUpdatedOdm(:,i,wRes(ij));
                    updatedPOdm(:,:,i,ij) = updatedPOdm(:,:,i,wRes(ij));
                    stateUpdatedCont(:,i,ij) = stateUpdatedCont(:,i,wRes(ij));
                    updatedPCont(:,:,i,ij) = updatedPCont(:,:,i,wRes(ij));
                    %   Association of updated states to variables
                    currStateOdm = stateUpdatedOdm(:,i,ij);
                    currStateCont = stateUpdatedCont(:,i,ij);
                    currPOdm = updatedPOdm(:,:,i,ij);
                    currPCont = updatedPCont(:,:,i,ij);
                    emptyNBool(i,ij) = emptyNBool(i,wRes(ij));
                    activeLabels(ij,i) = activeLabels(wRes(ij),i);
                    weightscoeff(ij,i) = 1/N;
                    transitionCurr = zeros(1,dictionarySize);
                    transitionCurr(1,1:dictionarySize) = transitionMat(activeLabels(ij,i),:);
                    labPredictProbs = makedist('Multinomial','Probabilities',...    %   probability of label of two neurons
                        transitionCurr);
                    labPredict(ij,i) = labPredictProbs.random(1,1);
                    indPred = find(codebook(:,4)==labPredict(ij,i));
                    node1Pred(ij,i) = codebook(indPred,1);
                    node2Pred(ij,i) = codebook(indPred,2);
                    U1 = UOdm(node1Pred(ij,i),3:4)';                                   %   From som1
                    U2 = UCont(node2Pred(ij,i),3:4)';
                    [statepredOdm(:,i,ij),PpredOdm(:,:,i,ij)] =...
                        kf_predict(currStateOdm,currPOdm, A, Q, B, U1);                 %predicition Odmetric
                    [statepredCont(:,i,ij),PpredCont(:,:,i,ij)] =...
                        kf_predict(currStateCont,currPCont, A, Q, B, U2);                 %predicition Control
                   %%%%    Caluclation of outside samples (mean)
                    meanCountOutside(i,1) = counOutside/N;
                    
                end
            end
        end
        if i == firstIt
            transitionCurr = zeros(1,dictionarySize);
            transitionCurr(1,1:dictionarySize) = transitionMat(activeLabels(n,i),:);
            labPredictProbs = makedist('Multinomial','Probabilities',...    %   probability of label of two neurons
                transitionCurr);
            labPredict(n,i) = labPredictProbs.random(1,1);
            indPred = find(codebook(:,4)==labPredict(n,i));
            node1Pred(n,i) = codebook(indPred,1);
            node2Pred(n,i) = codebook(indPred,2);
            U1 = UOdm(node1Pred(n,i),3:4)';                                   %   From som1
            U2 = UCont(node2Pred(n,i),3:4)';
            [statepredOdm(:,i,n),PpredOdm(:,:,i,n)] =...
                kf_predict(currStateOdm,currPOdm, A, Q, B, U1);                 %predicition Odmetric
            PpredOdm(:,:,i,n) = [1.0196         0         0         0
         0    1.0196         0         0
         0         0    1.0000         0
         0         0         0    1.0000];
            [statepredCont(:,i,n),PpredCont(:,:,i,n)] =...
                kf_predict(currStateCont,currPCont, A, Q, B, U2);                 %predicition Control
        end
    end
    if plotBool == true
        if i > firstIt
            subplot(3,1,1);
            scatter(UOdm(indexAbnormdb1Odm,1),UOdm(indexAbnormdb1Odm,2),'k','o')
            hold on
            scatter(dataTestOdometric(1,1:i),dataTestOdometric(2,1:i),'r','.')
            axis([-25 25 -5 40])
            xlabel('$x$','interpreter','latex')
            ylabel('$y$','interpreter','latex')
            grid on
            if i == 2
                %   Prediction
                %db1
                [abnormdb1OdmPred,PpredOdmPred1] =...
                    kf_predict(abnormdb1Odm_Probabilistic(i-1),3*stdOdmAbn1);
                [abnormdb1ContPred,PpredContPred1] =...
                    kf_predict(abnormdb1Cont_Probabilistic(i-1),3*stdContAbn1);
                abnormdb1Odm_ProbabilisticSmooth(i-1) = abnormdb1Odm_Probabilistic(i-1);
                abnormdb1Cont_ProbabilisticSmooth(i-1) = abnormdb1Cont_Probabilistic(i-1);
                
                %db2
                [abnormdb2OdmPred,PpredOdmPred2] =...
                    kf_predict(abnormdb2Odm_Probabilistic(i-1),3*stdOdmAbn2);
                [abnormdb2ContPred,PpredContPred2] =...
                    kf_predict(abnormdb2Cont_Probabilistic(i-1),3*stdContAbn2);
                abnormdb2Odm_ProbabilisticSmooth(i-1) = abnormdb2Odm_Probabilistic(i-1);
                abnormdb2Cont_ProbabilisticSmooth(i-1) = abnormdb2Cont_Probabilistic(i-1);
            else
                %   Update
                %db1
                [abnormdb1Odm_ProbabilisticSmooth(i-1),~,~] =...
                    kf_update(abnormdb1OdmPred,PpredOdmPred1, abnormdb1Odm_Probabilistic(i-1),1,10*stdOdmAbn1);
                [abnormdb1Cont_ProbabilisticSmooth(i-1),~,~] =...
                    kf_update(abnormdb1ContPred,PpredContPred1, abnormdb1Cont_Probabilistic(i-1),1,10*stdContAbn1);
                % db2
                [abnormdb2Odm_ProbabilisticSmooth(i-1),~,~] =...
                    kf_update(abnormdb2OdmPred,PpredOdmPred2, abnormdb2Odm_Probabilistic(i-1),1,10*stdOdmAbn2);
                [abnormdb2Cont_ProbabilisticSmooth(i-1),~,~] =...
                    kf_update(abnormdb2ContPred,PpredContPred2, abnormdb2Cont_Probabilistic(i-1),1,10*stdContAbn2);
                
                %   Prediction
                %db1
                [abnormdb1OdmPred,PpredOdmPred1] =...
                    kf_predict(abnormdb1Odm_ProbabilisticSmooth(i-1),3*stdOdmAbn1);
                [abnormdb1ContPred,PpredContPred1] =...
                    kf_predict(abnormdb1Cont_ProbabilisticSmooth(i-1),3*stdContAbn1);
                %db2
                [abnormdb2OdmPred,PpredOdmPred2] =...
                    kf_predict(abnormdb2Odm_ProbabilisticSmooth(i-1),3*stdOdmAbn2);
                [abnormdb2ContPred,PpredContPred2] =...
                    kf_predict(abnormdb2Cont_ProbabilisticSmooth(i-1),3*stdContAbn2);                
            end
            subplot(3,1,2);
            hold on
            plot(abnormdb2Odm_Probabilistic(1:i-1) ,'--b','LineWidth',1)
            plot(abnormdb2Odm_ProbabilisticSmooth(1:i-1) ,'-k','LineWidth',2)
%             plot(0:i,ones(1,i+1) ,'--g','LineWidth',2)
            xlabel('Time seconds (S)','interpreter','latex')
            ylabel('Abnormality measurments','interpreter','latex')
            title('Abnormality measurement (Odometry)','interpreter','latex')
            grid minor
            
            subplot(3,1,3);
            hold on
            plot(abnormdb2Cont_Probabilistic(1:i-1) ,'--r','LineWidth',1)
            plot(abnormdb2Cont_ProbabilisticSmooth(1:i-1) ,'-k','LineWidth',2)
%             plot(0:i,ones(1,i+1) ,'--b','LineWidth',2)
            xlabel('Time seconds (S)','interpreter','latex')
            ylabel('Abnormality measurments','interpreter','latex')
            title('Abnormality measurement (Control)','interpreter','latex')
            grid minor
            pause(0.005)
            frame = getframe(gcf);
            writeVideo(v1, frame);
        end
    end
end
estimationAbn.errorLocal = errorLocal;
estimationAbn.abnormdb1Odm_Probabilistic = abnormdb1Odm_Probabilistic;
estimationAbn.abnormdb1Cont_Probabilistic = abnormdb1Cont_Probabilistic;
estimationAbn.abnormdb1Odm_ProbabilisticSmooth = abnormdb1Odm_ProbabilisticSmooth;
estimationAbn.abnormdb1Cont_ProbabilisticSmooth = abnormdb1Cont_ProbabilisticSmooth;
estimationAbn.dataTestOdometric = dataTestOdometric;
estimationAbn.abnormdb2Odm_Probabilistic = abnormdb2Odm_Probabilistic;
estimationAbn.abnormdb2Cont_Probabilistic = abnormdb2Cont_Probabilistic;
estimationAbn.abnormdb2Odm_ProbabilisticSmooth = abnormdb2Odm_ProbabilisticSmooth;
estimationAbn.abnormdb2Cont_ProbabilisticSmooth = abnormdb2Cont_ProbabilisticSmooth;
estimationAbn.meanCountOutside = meanCountOutside;
close(v1)
end