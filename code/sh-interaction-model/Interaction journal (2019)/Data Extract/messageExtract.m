%   ZERO ALGORITHM, EXTRACT AND SAVE MESSAGE INFORMATION FROM THE BAGS
clc
clear
close all
curDir = pwd;

%%
MB = true;
scenario = 1; % 1(PM): training data; 2(PA), 3(ES), 4(Uturn) testing data
%% Modalities
odometric = true;
SandV = false;

%% Directories and read bag
if MB == true
    dirData = ['C:\Users\mohamad.baydoun\Desktop\Interaction journal (2019)',...
        '\Raw Data\' num2str(scenario)];    % Directory where dataset is located
    dirMat = ['C:\Users\mohamad.baydoun\Desktop\Interaction journal (2019)',...
        '\Generated Data\' num2str(scenario)]; % Directory where generated information is saved
    bag = bagDataRead(1, dirData, scenario);
else
    dirData = ['C:\Users\damian.campo\Desktop\Interaction journal (2019)',...
        '\Raw Data\' num2str(scenario)];    % Directory where dataset is located
    dirMat = ['C:\Users\damian.campo\Desktop\Interaction journal (2019)',...
        '\Generated Data\' num2str(scenario)];
        bag = bagDataRead(1, dirData, scenario);
end

%% Messges Extract
string = bag.MessageList.Topic;
if odometric == true
    [posID_icab1,~] = find(string =='/icab1/velodyne_odometry');
    timeStampsPos_icab1 = bag.MessageList.Time(posID_icab1);
    msgspos__icab1 = readMessages(bag,posID_icab1);
    filename1 = 'msgspos__icab1.mat';
    cd(dirMat)
    save(filename1, 'msgspos__icab1', 'timeStampsPos_icab1')
    cd(curDir)
end

if SandV == true
    [VID_icab1,~] = find(string =='/icab1/movement_manager/current_velocity');
    [SID_icab1,~] = find(string =='/icab1/movement_manager/current_steering');
    timeStampsvel_icab1 = bag.MessageList.Time(VID_icab1);
    timeStampsste_icab1 = bag.MessageList.Time(SID_icab1);
    msgsvel__icab1 = readMessages(bag,VID_icab1);
    msgsste__icab1 = readMessages(bag,SID_icab1);
    filename1 = 'msgscontrol__icab1.mat';
    cd(dirMat)
    save(filename1, 'msgsvel__icab1','msgsste__icab1','timeStampsvel_icab1','timeStampsste_icab1')
    cd(curDir)
end
