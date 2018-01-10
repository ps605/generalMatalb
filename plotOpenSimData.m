% Script to plot data from .sto and .mot files

% Data Directory
mainDir='C:\Users\ps605\Documents\PhD\Study_3\Modelling\Outputs\API\Results\OPC2C4_DR1_12\';

% Output plot directory 
plotDir='C:\Users\ps605\Documents\PhD\Study_3\Modelling\Outputs\API\Results\OPC2C4_DR1_12\Plots\';

% Move to data directory
cd C:\Users\ps605\Documents\PhD\Study_3\Modelling\Outputs\API\Results\OPC2C4_DR1_12\

% List folders of interest
folderWildCard='*HA*DH100*';
folders=dir(folderWildCard);

% OpenSim files of interest
forceRepStr='ForceReporter';
jointReactStr='JointReaction';
statesStr='states';

% Variables to plot
frVarsY={'c6Bush.cerv6.force.Y',...
    'c5Bush.cerv5.force.Y',...
    'c4Bush.cerv4.force.Y'};

frVarsX={'c6Bush.cerv6.force.X',...
    'c5Bush.cerv5.force.X',...
    'c4Bush.cerv4.force.X'};
     
% jrVars={'aux6jnt_on_cerv6_in_cerv6_fy',...
%     'aux6jnt_on_cerv6_in_cerv6_fy',...
%     'aux6jnt_on_cerv6_in_cerv6_fy'};

q={'aux6jnt_t2',...
    'aux5jnt_t2',...
    'aux4jnt_t2',...
    'aux6jnt_t1',...
    'aux5jnt_t1',...
    'aux4jnt_t1'};

qDot={'aux6jnt_t2_u',...
    'aux5jnt_t2_u',...
    'aux4jnt_t2_u',...
    'aux6jnt_t1_u',...
    'aux5jnt_t1_u',...
    'aux4jnt_t1_u'};

% Loop through folders containing data
for iFold=1:numel(folders)
    
    % Move into folder
    folder=[mainDir folders(iFold).name];
    cd (folder)
    [~,condition,~]=fileparts(folder);
    
    % List each of the data files
    jrFile=dir(['*' jointReactStr '*']);
    frFile=dir(['*' forceRepStr '*']);
    stFile=dir(['*' statesStr '*']);
    
    % Get Optimisation/Simulation Name ***** Need check
    simName=jrFile(1).name(1:13);
    
    % Get Data from files
    bushingData=importdata([folder '\' frFile(1).name]);
    stateData=importdata([folder '\' stFile(1).name]);
    
    
    % Plot BushingForce data
    plotDataFromSto(bushingData,...
        frVarsY,...
        [condition '_Y'],forceRepStr, 'Time(s)','BushingForce (N)',plotDir)
    
    plotDataFromSto(bushingData,...
        frVarsX,...
        [condition '_X' ],forceRepStr, 'Time(s)','BushingForce (N)',plotDir)
    
    % Plot position (q) data
     plotDataFromSto(stateData,...
        q,...
        [condition '_pos'],statesStr, 'Time(s)','Translation (m)',plotDir)
    
    % Plot velocity (qDot) data
     plotDataFromSto(stateData,...
        qDot,...
        [condition '_vel'],statesStr, 'Time(s)','Translation (m)',plotDir)
   
end
