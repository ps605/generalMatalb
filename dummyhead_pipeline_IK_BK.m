function data = dummyhead_pipeline_IK_BK(file)

% This function will create the appropriate files to run an OpenSim
% simulation sequence for the stumble recovery data
% Input - file - c3d file to process
% Originally written by Glen Lichtwark (The University of Queensland)
% Copyright 2011
% Contributions by David Graham and Chris Carty (Griffith University) 2012
%
% First, import the classes from the jar file so that these can be called
% directly

import org.opensim.modeling.*

%% Define path and model

osim_path = '\\Mac\Home\Desktop\Work\Research\OpenSim\Projects\DummyHead\ContactModel\DATA\SCALED_Model\';

model = 'DummyHead_SCALED';

ModelFile = [osim_path model '.osim'];

global subdir_current;
global trialdir_current;


if nargin < 1
    [fname, pname] = uigetfile('*.c3d', 'Select C3D file');
else
    if isempty(fileparts(file))
        pname = cd;
        pname = [pname '\'];
        fname = file;
    else [pname, name, ext] = fileparts(file);
        fname = [name ext];
    end
end

cd(pname);

%% Load the c3dfile

% load the c3d file using BTK
data = btk_loadc3d([pname, fname], 5);

%% Naming in case of Static cal files
if ~isempty(strfind(lower(fname),'cal'))
    data.sub_info.Name = fname(1:end-14);
    data.Name = fname(1:end-14);
else
    data.sub_info.Name = fname(7:end-4);
    data.Name = fname(7:end-4);
end

%% BTK uses First/Last so convert to fit existing routine so Start/End
data.marker_data.Last_Frame = data.marker_data.Last_Frame - data.marker_data.First_Frame;
data.marker_data.First_Frame = 1;

data.Start_Frame = data.marker_data.First_Frame;
data.End_Frame = data.marker_data.Last_Frame;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sort the C3D file so we know what is Marker data and what is calculated


if ~isempty(strfind(lower(fname),'cal'))
    marker_names = {'HD_OTL';'HD_OTR';'HD_OBL';'HD_OBR';'HD_O';'HD_LT';'HD_LB';'HD_LF';'HD_RT';'HD_RB';'HD_RF';'HD_C';'HD_N';'HD_T';... %Head
        'NK1_L';'NK1_C';'NK1_R';'NK2_L';'NK2_C';'NK2_R';'NK3_L';'NK3_C';'NK3_R';'TK1_L';'TK2_L';'TK3_L';'TK1_R';'TK2_R';... % Spine
        'PB_T1';'PB_T2';'PB_T3';'PB_T4';'PB_B1';'PB_B2';'PB_B3';'PB_B4';'PB_M1';'PB_M2';'PB_M3';'PB_M4';'PB_M5'}; % punch bag
    
    
else
    marker_names = {'HD_OTL';'HD_OTR';'HD_OBL';'HD_O';'HD_LT';'HD_RT';'HD_RB';'HD_RF';'HD_T';... %Head
        %         'HD_C';'HD_OBR';'HD_N';'HD_LB';'HD_LF';... %low visibility head markers
        'NK1_L';'NK1_C';'NK1_R';'NK2_L';'NK2_C';'NK2_R';'NK3_L';'NK3_C';'NK3_R';'TK1_L';'TK2_L';'TK3_L';'TK1_R';'TK2_R';... % Spine
        'PB_T1';'PB_T2';'PB_T3';'PB_T4';'PB_B1';'PB_B2';'PB_B3';'PB_B4';'PB_M1';'PB_M2';'PB_M3';'PB_M4';'PB_M5'}; % punch bag
end



data.marker_data = btk_sortc3d(data.marker_data,marker_names);

if isempty(strfind(lower(fname),'cal'))
    %% DEFINE EVENTS E1 E2 E3 E4
    
    
    % Load Force and Moment from Dummy Head data to calculate force peak to be selected as E(1)
    LoadCellC1 = importdata([pname(1:78) 'DH_OUTPUT\' trialdir_current 'MRF500.txt']);
    
    % Processing Dummy Head force
    Dummy_force=[];
    Dummy_moment=[];
    %Offset - check why z has high offset - Check 0.8 time offset from
    %Labview
    Dummy_force(:,1)=[ones(0.8*500,1)*mean(LoadCellC1(1:250,1));LoadCellC1(1:(5000-1*500),1)];
    Dummy_force(:,2)=[ones(0.8*500,1)*mean(LoadCellC1(1:250,2));LoadCellC1(1:(5000-1*500),2)];
    Dummy_force(:,3)=[zeros(0.8*500,1);LoadCellC1(1:(5000-1*500),3)-mean(LoadCellC1(1:250,3))];
    
    Dummy_moment(:,1)=[ones(0.8*500,1)*mean(LoadCellC1(1:250,4));LoadCellC1(1:(5000-1*500),4)];
    Dummy_moment(:,2)=[ones(0.8*500,1)*mean(LoadCellC1(1:250,5));LoadCellC1(1:(5000-1*500),5)];
    Dummy_moment(:,3)=[ones(0.8*500,1)*mean(LoadCellC1(1:250,6));LoadCellC1(1:(5000-1*500),6)];
    
    if strcmp(trialdir_current(1:4),'STFF')~=1
        % vpeak is the value, and fpeak is the frame value
        [TOP_M_vpeak,TOP_M_fpeak]=FINDMAX(abs(Dummy_force(:,1)-mean(LoadCellC1(1:250,1))));
        fTekscan=500;
        
        E(1) =  TOP_M_fpeak/fTekscan ; %the time of "impact" call
        
        E(2) = E(1)-0.2; %Start just before set call
        
        E(3) = E(1)+0.2; % Change 2.00 if you want a longer trial
    else
         E(1) =  1 ; % just a couple of sec afetr start as it is a static trial
        
        E(2) = E(1)-0.2; %Start just before 
        
        E(3) = E(1)+0.2; % Change if you want a longer trial
    end
    
    %% Define start and end frames
    % Define start and end frame from the events to write the appropriate TRC
    % and MOT files for the OpenSim simulations
    
    data.Start_Frame = round(E(2)/(1/data.marker_data.Info.frequency));
    data.End_Frame = round(E(3)/(1/data.marker_data.Info.frequency));
    
    %% Define file name prefix structure
    % Added a nameing structure to avoid confusion throughout the script,
    % participant code refers to the actual participant number (eg '066')
    % where as trial code refers to the trial (eg '20_4')
    
% % %     %** Hard coding mass and height at the moment  **
% % %     data.Mass=4.5+1.5;
% % %     data.Height=100;
    
    participant_code = data.Name;
    
    trial_code = fname(1:end-4);
    %% Open and read Tekscan centre of pressure values - TEKSCAN file needs to be in the same folder of .c3D files (xx_cal and trial_xx)
    % ResamplePressure sensors signals from 500 hz up to 2500 hz
    
    
    
    HD_Force=resample(Dummy_force,2500,500);
    HD_Moment=resample(Dummy_moment,2500,500);
    
    if strcmp(trialdir_current(1:4),'STFF')~=1
    %Load Markers for local system on the punch bag ifnot a STIFFNESS trial
    PB_B1=resample(data.marker_data.Markers.PB_B1,2500,250);
    PB_B2=resample(data.marker_data.Markers.PB_B2,2500,250);
    PB_B3=resample(data.marker_data.Markers.PB_B3,2500,250);
    PB_B4=resample(data.marker_data.Markers.PB_B4,2500,250);
    
    
    P1=(PB_B1'+PB_B2'+PB_B3'+PB_B4')/4;
    P2=PB_B3';
    P3=PB_B2';
    end
    
    %% Insert DATA in data
    
% % % % %     % DATA OK for force expressed in ground: 
% % % % %     %(3) Z compression (X in ground OpenSim) - (-)
% % % % %     %(2) Y lat bend(Z in ground OpenSim) --- (-) 
% % % % %     %(1) X FlexExt (Y in ground OpenSim)-- (+)
% % % % %     
% % % % %     data.fp_data.FP_data(3,1).TOP_M=[-HD_Force(:,2) -HD_Force(:,3) HD_Force(:,1)];
% % % % %     
% % % % %     data.fp_data.GRF_data(3,1).F=[-HD_Force(:,2) -HD_Force(:,3) HD_Force(:,1)];
% % % % %     
% % % % %     data.fp_data.GRF_data(3,1).M=[-HD_Moment(:,2) -HD_Moment(:,3) HD_Moment(:,1)];
% % % % %     
% % % % %     data.fp_data.GRF_data(3,1).P=zeros(length(HD_Force),3);
    
    % DATA OK for force expressed in skull
    data.fp_data.FP_data(3,1).TOP_M=[-HD_Force(:,2) HD_Force(:,1) -HD_Force(:,3)];
    
    data.fp_data.GRF_data(3,1).F=[-HD_Force(:,2) HD_Force(:,1) -HD_Force(:,3)];
    
    data.fp_data.GRF_data(3,1).M=[-HD_Moment(:,2) HD_Moment(:,1) -HD_Moment(:,3)];
    
    data.fp_data.GRF_data(3,1).P=zeros(length(HD_Force),3);
    
    
    %% C3D TO TRC
    %'v5' of this script is modified from previous ones since person is now
    %walking in the 'forward' direction of the treadmill (towards the wall) and
    %so coordinate system transformations have changed.
    
    
    data = btk_c3d2trc_v5_creationGRF(data,'off'); % animation off
    
    
    
    %% INVERSE KINEMATIC ANALYSIS 1 (to calculate Load Cell position)
    % Define the standard OpenSim files for IK
    ModelFile = [osim_path model '.osim'];
    ResultsDirectory = [pname 'IK\'];
    InputDirectory=pname;
    ConstraintWeight='Inf';
    MarkerFile=[trial_code '.trc'];
    IKTasksFile =[osim_path 'Dummy_Head_FULL_IK_Tasks.xml'];
    OutputFile= [ResultsDirectory trial_code '_ik.mot'];
    
% % % % %     % Set up the XML file
% % % % %     setup_InverseKinematics('data',data,...
% % % % %         'InputDirectory',InputDirectory,...
% % % % %         'ResultsDirectory',ResultsDirectory,...
% % % % %         'ResultsDirectory',ResultsDirectory,...
% % % % %         'ModelFile',ModelFile,...
% % % % %         'IKTasksFile',IKTasksFile,...
% % % % %         'ConstraintWeight',ConstraintWeight,...
% % % % %         'MarkerFile',MarkerFile,...
% % % % %         'Accuracy',0.00001,...
% % % % %         'OutputFile',OutputFile);%% DC & GT Changed from 0.000005 to 0.00005
% % % % %     
% % % % %     % Call the IK tool
% % % % %     com = ['ik -S ' participant_code '_Setup_InverseKinematics.xml'];
% % % % %     system(com)
    
    
% % % % %         %% Load the data from MOT file and tie kinematic data to data structure
% % % % %         D = load_sto_file([pname trial_code '_ik.mot']);
% % % % %         fnames = fieldnames(D);
    
    
    % % % % %     clear a
    % % % % %     a = [strfind(lower(fnames),'tx') strfind(lower(fnames),'ty') strfind(lower(fnames),'tz') strfind(lower(fnames),'pelvis')];
    % % % % %     for i = 1:size(a,1)
    % % % % %         if isempty([a{i,1:3}])
    % % % % %             data.Kinematics.Angles.(fnames{i}) = D.(fnames{i});
    % % % % %             [D.(fnames{i}),~]= ADF_BUTTERWORTH(D.(fnames{i}),250,[6 10],1,'N','N');
    % % % % %
    % % % % %         elseif ~isempty([a{i,4}])
    % % % % %             data.Kinematics.Angles.(fnames{i}) = D.(fnames{i});
    % % % % %                         [D.(fnames{i}),~]= ADF_BUTTERWORTH(D.(fnames{i}),250,[6 10],1,'N','N');
    % % % % %
    % % % % %         else data.Kinematics.Markers.(fnames{i}) = D.(fnames{i});
    % % % % %                         [D.(fnames{i}),~]= ADF_BUTTERWORTH(D.(fnames{i}),250,[6 10],1,'N','N');
    % % % % %         end
    % % % % %     end
    
    %% BODY KINEMATIC ANALYSIS
    % SetUp BodyKinematics
    % Define the standard OpenSim files for BodyKinematics
    
    % DEFINE STANDARD OPENSIM FILES AND PATHS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CoordinatesFile =[ResultsDirectory trial_code '_ik.mot'];
    MOTFile = [ResultsDirectory trial_code '_ik.mot'];
    % % % % ForceSetFile        = [osim_path model 'GT_CMC_Actuators.xml'];
    % % % % ExternalLoadsFile   = GRFFileXML;
    MADirectoryName     = [pname 'BK_Results'];
    %     ModelFile = [participant_code '_SCALED.osim'];
    
    % Set up the XML file
    setup_kinematicsDG( 'data',data,...
        'ModelFile',ModelFile,...
        'CoordinatesFile',CoordinatesFile,...
        'MOTFILE',MOTFile,...
        'AnalysisOn', true,...
        'UseFLV', false,...
        'StartTime',data.Start_Frame,...
        'SolveForEquilibrium',false,...
        'UseModelForceSet',false,...
        'EndTime',data.End_Frame,...
        'LowPassFilterFreq', 6,...
        'DirectoryName', MADirectoryName);
    
    % Call the Body Kinematics Tool
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    com=['analyze -S ', participant_code '_Setup_BodyKinematics.xml'];
    system(com);
    % SAVE THE WORKSPACE AND PRINT A LOG FILE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % fid = fopen([pname,'\BK_',data.Name,'.log'],'w+');
    % % % % fprintf(fid,'%s\n', log_mes);
    % % % % fclose(fid);
    
    
    
    % Load the data from Body Kinematics MOT file and tie body kinematic data to data structure for
    % later use
    B = load_sto_file([pname 'BK_Results\' participant_code '_' trial_code '_BK_BodyKinematics_pos_global.sto']);
    data.BodyPos = B;
    
    % Add Dummy head ensor position at C1-skull joint (C1 origin at he
    % momeent) - In the future we want to measure exactly the position
    % form DummyHead specs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    onesx_pre=ones(E(2)*250-1,1)*B.cerv1_X(1);
    onesy_pre=ones(E(2)*250-1,1)*B.cerv1_Y(1);
    onesz_pre=ones(E(2)*250-1,1)*B.cerv1_Z(1);
    
    onesx_post=ones(data.marker_data.Info.NumFrames-(E(3)*250),1)*B.cerv1_X(end);
    onesy_post=ones(data.marker_data.Info.NumFrames-(E(3)*250),1)*B.cerv1_Y(end);
    onesz_post=ones(data.marker_data.Info.NumFrames-(E(3)*250),1)*B.cerv1_Z(end);
    
    SOrx_event=[ onesx_pre;B.cerv1_X;onesx_post];
    SOry_event=[ onesy_pre;B.cerv1_Y;onesy_post];
    SOrz_event=[ onesz_pre;B.cerv1_Z;onesz_post];
    
    %  Point of application of the force needs to be synch with force data
    %  as it has been calculated from an IK already cut.
    SOrx=resample(SOrx_event,2500,250);
    SOry=resample(SOry_event,2500,250);
    SOrz=resample(SOrz_event,2500,250);
    
    if length(SOrx)== data.marker_data.Info.NumFrames*10
        if length(HD_Force)>length(SOrx)
            data.fp_data.GRF_data(3,1).P=[SOrx(1: length(SOrx)) SOry(1:length(SOrx)) SOrz(1:length(SOrx))];
        else
            data.fp_data.GRF_data(3,1).P=[SOrx(1: length(HD_Force)) SOry(1:length(HD_Force)) SOrz(1:length(HD_Force))];
        end
    else
        fprintf('Check your FORCE TRACES!!!');
        data.fp_data.GRF_data(3,1).P=[SOrx SOry SOrz];
        
    end
    clear B
    %% C3D TO TRC
    
    data = btk_c3d2trc_v6_no_coord_reordering(data,'off'); % animation off
    
    
    %% INVERSE DYNAMIC ANALYSIS
    
    % % Define the standard OpenSim files for ID
    
    MOTFile = [ResultsDirectory trial_code '_ik.mot'];
    GRFFile = [pname trial_code '_grf.mot'];
    GRFFileXML = [pname trial_code '_grf.xml'];
    ResultsDirectoryID = [pname 'ID\'];
    OutputFile = [trial_code '_ID.sto'];
    
    % Write GRF XML file
    grf2xml(data,'ExternalLoadNames',{'ExternalForce_1','ExternalForce_2','DH_LoadCell'},...
        'AppliedToBodies',{'ground','ground','skull'},...
        'ForceExpressedinBody',{'ground','ground','skull'},...
        'PointExpressedinBody',{'ground','ground','ground'},...
        'GRFFile',GRFFile,'MOTFile',...
        MOTFile,'LowPassFilterForKinematics',6,'OutputFile',GRFFileXML);
    
    % Set up the XML file
    setup_InverseDynamics('data',data, 'ModelFile',ModelFile,'MOTFile',MOTFile,...
        'GRFFile',GRFFileXML,'LowPassFilterForKinematics',6,'OutputFile', OutputFile, 'ResultsDirectory',  ResultsDirectoryID);
    
    % Call the ID tool
    com = ['id -S ' participant_code '_Setup_InverseDynamics.xml'];
    
    system(com);
    
    clear D
    
    % Load the data from ID MOT file and tie kinetic data to data structure
    % %     D = load_sto_file([ResultsDirectoryID trial_code '_ID.sto']);
    % %     data.Kinetics = D;
    
    % % % % % save([(subdir_current) '.mat']');
    % % % % % [outD]= analyse_Dynamics_rugby_scrum(data.Kinetics,[pname 'outKIN\'],[fname '_outKIN'],t_rENG,t_peakSMF);
    % % %     clear D
    
    
    %% Analysis Tool
    %% FORCE REPORTER ANALYSIS
    % SetUp BodyKinematics
    % Define the standard OpenSim files for BodyKinematics
    
    % DEFINE STANDARD OPENSIM FILES AND PATHS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CoordinatesFile =[pname 'IK\' trial_code '_ik.mot'];
    MOTFile = [pname trial_code 'IK\' '_ik.mot'];
    % % % % ForceSetFile        = [osim_path model 'GT_CMC_Actuators.xml'];
    ExternalLoadsFile   = GRFFileXML;
    mkdir(pname,'ForceRep_Results\');
    MADirectoryName     = [pname 'ForceRep_Results'];
    
    
    % Set up the XML file
    setup_ForceReporter( 'data',data,...
        'ModelFile',ModelFile,...
        'CoordinatesFile',CoordinatesFile,...
        'MOTFILE',MOTFile,...
        'AnalysisOn', true,...
        'UseFLV', false,...
        'StartTime',data.Start_Frame,...
        'SolveForEquilibrium',false,...
        'UseModelForceSet',false,...
        'EndTime',data.End_Frame,...
        'LowPassFilterFreq', 6,...
        'OutputPrecision',20,...
        'ExternalLoadsFile',ExternalLoadsFile,...
        'DirectoryName', MADirectoryName);
    
    setup_JointReaction( 'data',data,...
        'ModelFile',ModelFile,...
        'CoordinatesFile',CoordinatesFile,...
        'MOTFILE',MOTFile,...
        'AnalysisOn', true,...
        'UseFLV', false,...
        'StartTime',data.Start_Frame,...
        'SolveForEquilibrium',true,...
        'UseModelForceSet',false,...
        'EndTime',data.End_Frame,...
        'LowPassFilterFreq', 6,...
        'OutputPrecision',20,...
        'ExternalLoadsFile',ExternalLoadsFile,...
        'DirectoryName', MADirectoryName);
    
    % Call the Force Reporter and Joint Reaction Tools
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    com=['analyze -S ', participant_code '_Setup_ForceReporter.xml'];
    system(com);
    
    com=['analyze -S ', participant_code '_Setup_JointReaction.xml'];
    system(com);
    
    % Load the data from joint reaction to data structure
    D = load_sto_file([MADirectoryName '\' trial_code '_JR' '_JointReaction_ReactionLoads.sto']);
    data.JointReaction = D;
    
    clear D
end



