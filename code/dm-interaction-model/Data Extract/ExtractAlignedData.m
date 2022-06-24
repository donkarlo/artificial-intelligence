%   FIRST ALGORITHM, Alignment between modalities and extract Aligned data
%   Also some Filtering and Smoothing of data

set(0,'defaultfigurecolor',[1 1 1])
clc
clear
close all
curDir = pwd;

%%
MB = true;
scenario = 2; % % 1(PM): training data; 2(PA), 3(ES), 4(Uturn) testing data

%% Directories
if MB == true
    dirMat = ['C:\Users\mohamad.baydoun\Desktop\Interaction journal (2019)',...
        '\Generated Data\' num2str(scenario)];
else
    dirMat = ['C:\Users\damian.campo\Desktop\Interaction journal (2019)',...
        '\Generated Data\' num2str(scenario)];
end
cd(dirMat)

%% Load modalities and do alignment
load('msgspos__icab1.mat')
load('msgscontrol__icab1.mat')

allTimeStamps{1,1} =  timeStampsPos_icab1;
allTimeStamps{2,1} =  timeStampsvel_icab1;
allTimeStamps{3,1} =  timeStampsste_icab1;

minTimeStamp = min(cell2mat(allTimeStamps));
allTimeStamps{1,1} = allTimeStamps{1,1} - minTimeStamp;
allTimeStamps{2,1} = allTimeStamps{2,1} - minTimeStamp;
allTimeStamps{3,1} = allTimeStamps{3,1} - minTimeStamp;
%

min1 = min(allTimeStamps{1,1});
min2 = min(allTimeStamps{2,1});
min3 = min(allTimeStamps{3,1});
newmin = max([min1,min2,min3]);

ind1 =  find(allTimeStamps{1,1}>= newmin);
ind2 = find(allTimeStamps{2,1}>= newmin);
ind3 = find(allTimeStamps{3,1}>= newmin);

allTimeStamps{1,1} =  allTimeStamps{1,1}(ind1,1) - newmin;
allTimeStamps{2,1} =  allTimeStamps{2,1}(ind2,1)- newmin;
allTimeStamps{3,1} =  allTimeStamps{3,1}(ind3,1)- newmin;
%

[Val,Id] = min([length(allTimeStamps{1,1}),length(allTimeStamps{2,1}),...
    length(allTimeStamps{3,1})]);

allTimeStampsAllign = cell(3,1);
allTimeStampsID = cell(3,1);

for i = 1:Val
    for j = 1:length(allTimeStamps)
        if j ~= Id
            [~, minID] = min(pdist2(allTimeStamps{Id,1}(i,1), allTimeStamps{j,1}));
            allTimeStampsAllign{j,1} = [allTimeStampsAllign{j,1}; allTimeStamps{j,1}(minID,1)];
            allTimeStampsID{j,1} = [allTimeStampsID{j,1}; minID];
        else
            allTimeStampsAllign{j,1} = allTimeStamps{j,1};
            allTimeStampsID{j,1} = [1:Val]';
        end
    end
end

%% Extract Aligned data
allData{1,1} =  msgspos__icab1;
allData{2,1} =  msgsvel__icab1;
allData{3,1} =  msgsste__icab1;

allData{1,1} =  allData{1,1}(ind1,1);
allData{2,1} =  allData{2,1}(ind2,1);
allData{3,1} =  allData{3,1}(ind3,1);

allData{1,1} =  allData{1,1}(allTimeStampsID{1, 1},1);
allData{2,1} =  allData{2,1}(allTimeStampsID{2, 1},1);
allData{3,1} =  allData{3,1}(allTimeStampsID{3, 1},1);

pos = [];
control = [];

msgspos = allData{1,1};
msgsvel = allData{2,1};
msgsste = allData{3,1};

for i = 1:size(msgspos,1)
    % Odometric
    pos1 = [msgspos{i,1}.Pose.Pose.Position.X,...
        msgspos{i,1}.Pose.Pose.Position.Y];
    pos = [pos ; pos1];
    % Control
    control1 = [msgsvel{i, 1}.Data, msgsste{i, 1}.Data];
    control = [control ;control1];
end

%% Filtering data
Vfiltered = filloutliers(control(:,1),'nearest','mean');
Sfiltered =  filloutliers(control(:,2),'nearest','mean');

%%  Smoothing data
V = smooth(Vfiltered,55,'moving');
S = smooth(Sfiltered,15,'moving');

%% calculation of derivatives of data
divV = diff(V);
divS = diff(S);
divV = smooth(divV,15,'moving');
xPos = pos(:,1);
yPos = pos(:,2);
if scenario == 1
    angRot = 86.5;
    posFx = cosd(angRot)*xPos - sind(angRot)*yPos;
    posFy = sind(angRot)*xPos + cosd(angRot)*yPos;
    xPos = posFx + 20;
    yPos = posFy + 20;
elseif scenario == 2
    angRot = 2;
    posFx = cosd(angRot)*xPos - sind(angRot)*yPos;
    posFy = sind(angRot)*xPos + cosd(angRot)*yPos;
    xPos = posFx + 0.6;
    yPos = posFy + 2.3;
elseif scenario == 3
    angRot = -1;
    posFx = cosd(angRot)*xPos - sind(angRot)*yPos;
    posFy = sind(angRot)*xPos + cosd(angRot)*yPos;
    xPos = posFx + 0.6;
    yPos = posFy + 2.3;
elseif scenario == 4
    angRot = -3;
    posFx = cosd(angRot)*xPos - sind(angRot)*yPos;
    posFy = sind(angRot)*xPos + cosd(angRot)*yPos;
    xPos = posFx + 0.6;
    yPos = posFy + 2.3;
end
divxPos = diff(xPos);
divyPos = diff(yPos);

%% Last cut due to derivatives
V = V(1:end-1);
S = S(1:end-1);
xPos = xPos(1:end-1);
yPos = yPos(1:end-1);

cd(curDir)
%% Normalization Filtered data
[normS,minSNorm,maxSNorm] = Normalization(S);
[normV,minVNorm,maxVNorm] = Normalization(V);
[normxPos,minxPosNorm,maxxPosNorm] = Normalization(xPos);
[normyPos,minyPosNorm,maxyPosNorm] = Normalization(yPos);
[normdivS,mindivSNorm,maxdivSNorm] = Normalization(divS);
[normdivV,mindivVNorm,maxdivVNorm] = Normalization(divV);
[normdivxPos,mindivxPosNorm,maxdivxPosNorm] = Normalization(divxPos);
[normdivyPos,mindivyPosNorm,maxdivyPosNorm] = Normalization(divyPos);

%% saving data

structSyncData.Filtered.S = S;
structSyncData.Filtered.V = V;
structSyncData.Filtered.xPos = xPos;
structSyncData.Filtered.yPos = yPos;
structSyncData.Filtered.divS = divS;
structSyncData.Filtered.divV = divV;
structSyncData.Filtered.divxPos = divxPos;
structSyncData.Filtered.divyPos = divyPos;

structSyncData.FilteredNorm.S = normS;
structSyncData.FilteredNorm.V = normV;
structSyncData.FilteredNorm.xPos = normxPos;
structSyncData.FilteredNorm.yPos = normyPos;
structSyncData.FilteredNorm.divS = normdivS;
structSyncData.FilteredNorm.divV = normdivV;
structSyncData.FilteredNorm.divxPos = normdivxPos;
structSyncData.FilteredNorm.divyPos = normdivyPos;

structSyncData.NormMinMax.S = [minSNorm maxSNorm];
structSyncData.NormMinMax.V = [minVNorm maxVNorm];
structSyncData.NormMinMax.xPos = [minxPosNorm maxxPosNorm];
structSyncData.NormMinMax.yPos = [minyPosNorm maxyPosNorm];

structSyncData.NormMinMax.Sdiv = [mindivSNorm maxdivSNorm];
structSyncData.NormMinMax.Vdiv = [mindivVNorm maxdivVNorm];
structSyncData.NormMinMax.xPosdiv = [mindivxPosNorm maxdivxPosNorm];
structSyncData.NormMinMax.yPosdiv = [mindivyPosNorm maxdivyPosNorm];
cd(dirMat)
save('genSyncdata','structSyncData')
cd(curDir)
