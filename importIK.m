% this script wll read in IK .mot file then filter it for specified
% purposes with the specified butter filter. Then it will creat a smoothed
% .mat, .sto and .mot file. The .sto file can be used to specify the inital
% states during FD simulations

clear;
clc ;

for i=1:5
    fileNo=num2str(i);
    ikFile=['C:\Users\ps605\Documents\PhD\Study_3\Optimisation\Out\IK\IK_drop_1_0' fileNo '.mot'];
    
    [path,ikName,ext]=fileparts(ikFile);
    
    saveFile=['C:\Users\ps605\Documents\PhD\Study_3\Optimisation\In\IK\' ikName '_states_500Hz'];
    
    ikData=importdata(ikFile);
    
    colHeaders=ikData.colheaders;
    %% Filter IK coords
    
    fS=3000; % Sampling f
    fC=500; % Cut off f
    
    [b,a]=butter(4,fC/(fS/2));
    
    filtQ=filtfilt(b,a,ikData.data(:,2:end)); % not include time (first col) - only coords
    numQ=size(filtQ,2);
    
    %% Create time and coords data matrix
    
    ikDataFinal(:,1)=ikData.data(:,1);
    ikDataFinal(:,2:numQ+1)=filtQ;
    
    %% Save data
    
%     save(saveFile,'ikDataFinal', 'colHeaders'); % .mat file
%     
%     generateMotFile(ikDataFinal,colHeaders,[saveFile '.mot']); % .mot
%     generateMotFile(ikDataFinal,colHeaders,[saveFile '.sto']); % .sto
    
    figure(1)
    plot(ikData.data(:,1),ikDataFinal(:,24));
    title('aux4 tx');
    xlabel('time (s)');
    ylabel ('displacement (m)');
    
    hold on
    
    figure(2)
    plot(ikData.data(:,1),ikDataFinal(:,25));
    title('aux4 ty');
    xlabel('time (s)');
    ylabel ('displacement (m)');
    
    hold on
    
    clear filtQ numQ ikDataFinal
end
