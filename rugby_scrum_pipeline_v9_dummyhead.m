function data = rugby_scrum_pipeline_v9_dummyhead(file)

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

osim_path = 'C:\Users\User\Documents\Pavlos\University of Bath\Year 4\FYP\Models\Modified\'; %Swap

model = 'RModel_CutDown_061216.xml';%swap


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
%%%%%[data.marker_data, data.analog_data, data.fp_data, data.sub_info] = btk_loadc3d([pname, fname], 5);
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
    marker_names = {'HD_OTL';'HD_OTR';'HD_OBL';'HD_OBR';'HD_O';'HD_LT';'HD_LB';'HD_LF';'HD_RT';'HD_RB';'HD_RF';'HD_C';'HD_N';'HD_T';... %Head
        'NK1_L';'NK1_C';'NK1_R';'NK2_L';'NK2_C';'NK2_R';'NK3_L';'NK3_C';'NK3_R';'TK1_L';'TK2_L';'TK3_L';'TK1_R';'TK2_R';... % Spine
        'PB_T1';'PB_T2';'PB_T3';'PB_T4';'PB_B1';'PB_B2';'PB_B3';'PB_B4';'PB_M1';'PB_M2';'PB_M3';'PB_M4';'PB_M5'}; % punch bag
end



data.marker_data = btk_sortc3d(data.marker_data,marker_names);

if isempty(strfind(lower(fname),'cal'))
    %% DEFINE EVENTS E1 E2 E3 E4
    
    
    % Load Pressure data to calculate force peak to be selected as E(1)
    
    load ([pname 'LoadCellC1.mat']);
    
    % Processing Dummy Head force
    Dummy_force=[];
    %Offset - check why z has high offset - Check 0.8 time offset from
    %Labview
    Dummy_force(:,1)=[ones(0.8*500,1)*mean(LoadCellC1(1:250,1));LoadCellC1(1:(5000-0.8*500),1)];
    Dummy_force(:,2)=[ones(0.8*500,1)*mean(LoadCellC1(1:250,2));LoadCellC1(1:(5000-0.8*500),2)];
    Dummy_force(:,3)=[zeros(0.8*500,1);LoadCellC1(1:(5000-0.8*500),3)-mean(LoadCellC1(1:250,3))];
    
    % vpeak is the value, and fpeak is the frame value
    [TOP_M_vpeak,TOP_M_fpeak]=FINDMAX(abs(Dummy_force(:,1)-mean(LoadCellC1(1:250,1))));
    fTekscan=500;
    
    E(1) =  TOP_M_fpeak/fTekscan ; %the time of "impact" call
    
    E(2) = E(1)-0.2; %Start just before set call
    
    E(3) = E(1)+0.2; % Change 2.00 if you want a longer trial
    
    
    %% Define start and end frames
    % Define start and end frame from the events to write the appropriate TRC
    % and MOT files for the OpenSim simulations
    
    data.Start_Frame = round(E(2)/(1/data.marker_data.Info.frequency));
    data.End_Frame = round(E(3)/(1/data.marker_data.Info.frequency));
    
    %% Define file name prefix structure
    % Added a nameing structure to avoid confusion throughout the script,
    % participant code refers to the actual participant number (eg '066')
    % where as trial code refers to the trial (eg '20_4')
    
    %** Hard coding mass and height at the moment  **
    data.Mass=8;
    data.Height=100;
    
    participant_code = data.Name;
    
    trial_code = fname(1:end-4);
    %% Open and read Tekscan centre of pressure values - TEKSCAN file needs to be in the same folder of .c3D files (xx_cal and trial_xx)
    % ResamplePressure sensors signals from 500 hz up to 2500 hz
    
    
    
    HD_Force=[resample(Dummy_force,2500,500)];
    
    
    %Load Markers for local system on the punch bag
    PB_B1=resample(data.marker_data.Markers.PB_B1,2500,250);
    PB_B2=resample(data.marker_data.Markers.PB_B2,2500,250);
    PB_B3=resample(data.marker_data.Markers.PB_B3,2500,250);
    PB_B4=resample(data.marker_data.Markers.PB_B4,2500,250);
    
    
    P1=(PB_B1'+PB_B2'+PB_B3'+PB_B4')/4;
    P2=PB_B3';
    P3=PB_B2';
    
    
    %%Insert DATA in data
    
    %    data.fp_data.FP_data(3,1).TOP_M=[HD_Force(:,2) HD_Force(:,3) HD_Force(:,1)];
    data.fp_data.FP_data(3,1).TOP_M=[-HD_Force(:,2) -HD_Force(:,3) HD_Force(:,1)];
    
    data.fp_data.GRF_data(3,1).F=[-HD_Force(:,2) -HD_Force(:,3) HD_Force(:,1)];
    
    data.fp_data.GRF_data(3,1).M=zeros(length(HD_Force),3);
    
    data.fp_data.GRF_data(3,1).P=zeros(length(HD_Force),3);
    
    
    
    %% C3D TO TRC
    %'v5' of this script is modified from previous ones since person is now
    %walking in the 'forward' direction of the treadmill (towards the wall) and
    %so coordinate system transformations have changed.
    
    
    data = btk_c3d2trc_v5_creationGRF(data,'off'); % animation off
    
    
    
    %% INVERSE KINEMATIC ANALYSIS 1 (to calculate Load Cell position)
    % Define the standard OpenSim files for IK
    ModelFile = [participant_code '_SCALED.osim'];
    ResultsDirectory = pname;
    InputDirectory=pname;
    ConstraintWeight='Inf';
    MarkerFile=[trial_code '.trc'];
    IKTasksFile ='Dummy_Head_FULL_IK_Tasks.xml';
    OutputFile= [trial_code '_ik.mot'];
    
    % Set up the XML file
    setup_InverseKinematics('data',data,...
        'InputDirectory',InputDirectory,...
        'ResultsDirectory',ResultsDirectory,...
        'ResultsDirectory',ResultsDirectory,...
        'ModelFile',ModelFile,...
        'IKTasksFile',IKTasksFile,...
        'ConstraintWeight',ConstraintWeight,...
        'MarkerFile',MarkerFile,...
        'Accuracy',0.00001,...
        'OutputFile',OutputFile);%% DC & GT Changed from 0.000005 to 0.00005
    
    % Call the IK tool
    com = ['ik -S ' participant_code '_Setup_InverseKinematics.xml'];
    system(com)
    
    
    % % % % %     % Load the data from MOT file and tie kinematic data to data structure
    % % % % %     D = load_sto_file([pname trial_code '_ik.mot']);
    % % % % %     fnames = fieldnames(D);
    
    
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
    CoordinatesFile =[pname trial_code '_ik.mot'];
    MOTFile = [pname trial_code '_ik.mot'];
    % % % % ForceSetFile        = [osim_path model 'GT_CMC_Actuators.xml'];
    % % % % ExternalLoadsFile   = GRFFileXML;
    MADirectoryName     = [pname 'BK_Results'];
    ModelFile = [participant_code '_SCALED.osim'];
    
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
    zeros_pre=zeros(E(2)*250-1,1);
    zeros_post=zeros(data.marker_data.Info.NumFrames-(E(3)*250),1);
    
    SOrx_event=[ zeros_pre;B.cerv1_X;zeros_post];
    SOry_event=[ zeros_pre;B.cerv1_Y;zeros_post];
    SOrz_event=[ zeros_pre;B.cerv1_Z;zeros_post];
    
    %  Point of application of the force needs to be synch with force data
    %  as it has been calculated from an IK already cut.
    SOrx=resample(SOrx_event,2500,250);
    SOry=resample(SOry_event,2500,250);
    SOrz=resample(SOrz_event,2500,250);
    
    if length(SOrx)== data.marker_data.Info.NumFrames*10
    
    data.fp_data.GRF_data(3,1).P=[SOrx(1: length(HD_Force)) SOry(1:length(HD_Force)) SOrz(1:length(HD_Force))];
    else
        fprintf('Check your FORCE TRACES!!!');
        data.fp_data.GRF_data(3,1).P=[SOrx SOry SOrz];

    end
    clear B
    %% C3D TO TRC
    %'v5' of this script is modified from previous ones since person is now
    %walking in the 'forward' direction of the treadmill (towards the wall) and
    %so coordinate system transformations have changed.
    
    
    data = btk_c3d2trc_v5(data,'off'); % animation off
    
    
    %% INVERSE DYNAMIC ANALYSIS
    
    % % Define the standard OpenSim files for ID
    MOTFile = 'C:\Users\User\Documents\Pavlos\University of Bath\Year 4\FYP\Models\EXHS01\IK\trial_EXHS_01_ik.mot';
    GRFFile = 'C:\Users\User\Documents\Pavlos\University of Bath\Year 4\FYP\Models\EXHS01\trial_EXHS_01_grf.mot';
    GRFFileXML = 'C:\Users\User\Documents\Pavlos\University of Bath\Year 4\FYP\Models\EXHS01\trial_EXHS_01_grf.xml';
    ModelFile = 'C:\Users\User\Documents\Pavlos\University of Bath\Year 4\FYP\Models\Modified\RModel_CutDown_061216.xml';
    OutputFile = 'C:\Users\User\Documents\Pavlos\University of Bath\Year 4\FYP\Models\EXHS01\trial_EXHS_ID.sto';
    ResultsDirectory = 'C:\Users\User\Documents\Pavlos\University of Bath\Year 4\FYP\Models\EXHS01\ID\';
    
    
    % Write GRF XML file
    grf2xml(data,'ExternalLoadNames',{'ExternalForce_1','ExternalForce_2','DH_LoadCell'},...
        'AppliedToBodies',{'ground','ground','skull'},...
        'ForceExpressedinBody',{'skull','skull','skull'},...
        'PointExpressedinBody',{'skull','skull','skull'},...
        'GRFFile',GRFFile,'MOTFile',...
        MOTFile,'LowPassFilterForKinematics',6,'OutputFile',GRFFileXML);
    
    % Set up the XML file
    setup_InverseDynamics('data',data, 'ModelFile',ModelFile,'MOTFile',MOTFile,...
        'GRFFile',GRFFileXML,'LowPassFilterForKinematics',6,'OutputFile', OutputFile, 'ResultsDirectory',  ResultsDirectory);
    
    % Call the ID tool
    com = ['id -S ' participant_code '_Setup_InverseDynamics.xml'];
    
    system(com);
    
    clear D
    
    % Load the data from ID MOT file and tie kinetic data to data structure
    D = load_sto_file([ResultsDirectory trial_code '_ID.sto']);
    data.Kinetics = D;
    
    % % % % % save([(subdir_current) '.mat']');
    % % % % % [outD]= analyse_Dynamics_rugby_scrum(data.Kinetics,[pname 'outKIN\'],[fname '_outKIN'],t_rENG,t_peakSMF);
    clear D
    %
    
    % % % %     %% Forward Dynamics (inputs from ID needs to be refined!!!)
    % % % %
    % % % %     % % Define the standard OpenSim files for FD
    % % % %
    % % % %     ModelFile = [pname participant_code '_SCALED.osim'];
    % % % %     OutputFile = [trial_code '_FD.sto'];
    % % % %     ResultsDirectory = [pname 'FD\'];
    % % % %     FDControlsFile=[pname 'ID\' 'trial_DHNT02_ID.sto'];
    % % % %
    % % % %     % Set up the XML file
    % % % %     setup_ForwardDynamics('data',data, 'ModelFile',ModelFile,...
    % % % %         'InitialTime',0.1,'FinalTime',0.3,'AuxillaryStates',true,'FDControlsFile',FDControlsFile,...
    % % % %         'MinStepSize',1e-10,'MaxStepSize',0.0001,'ErrorTol',1e-8,...
    % % % %         'Analyses',ForceReporter,...
    % % % %         'OutputPrecision',20,'LowPassFilterFreq',6,'OutputFile', OutputFile,...
    % % % %         'ResultsDirectory',  ResultsDirectory, 'InputDirectory', pname);
    % % % %
    % % % %     % Call the ID tool
    % % % %     com = ['forward -S ' participant_code '_Setup_ForwardDynamics.xml'];
    % % % %
    % % % %     system(com);
    % % % %
    % % % %     clear D
    % % % %
    
    %% Analysis Tool
    %% FORCE REPORTER ANALYSIS
    % SetUp BodyKinematics
    % Define the standard OpenSim files for BodyKinematics
    
    % DEFINE STANDARD OPENSIM FILES AND PATHS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CoordinatesFile =[pname trial_code '_ik.mot'];
    MOTFile = [pname trial_code '_ik.mot'];
    % % % % ForceSetFile        = [osim_path model 'GT_CMC_Actuators.xml'];
    % % % % ExternalLoadsFile   = GRFFileXML;
    MADirectoryName     = [pname 'ForceRep_Results'];
    ModelFile = [pname participant_code '_SCALED.osim'];
    
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
        'DirectoryName', MADirectoryName);
    
    % Call the Body Kinematics Tool
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    com=['analyze -S ', participant_code '_Setup_ForceReporter.xml'];
    system(com);
    
    
    % % % % % % %     %% RESIDUAL REDUCTION ANALYSIS
    % % % % % % %     % Set up RRA - define the standard OpenSim files for RRA
    % % % % % % %     %(MOT and GRF_MOT defined above along with Model file)
    % % % % % % %
    % % % % % % %     MOTFile = [pname trial_code '_ik.mot'];
    % % % % % % %     GRFFile = [pname trial_code '_grf.mot'];
    % % % % % % %     GRFFileXML = [pname trial_code '_grf.xml'];
    % % % % % % %
    % % % % % % %     % % % % % % RRATaskFile = [osim_path model '_FULL_RRA_Tasks.xml'];%[data.Name '_RRA_Tasks.xml'];
    % % % % % % %     % % % % % % RRAForceFile = [osim_path model '_FULL_RRA_Actuators.xml'];%[data.Name '_RRA_Actuators.xml'];
    % % % % % % %     % %  Needed to change 'osim_path' to  '../../../RugbyModel/' to make it
    % % % % % % %     % working. Don't know the reason why it is needed.
    % % % % % % %
    % % % % % % %     RRATaskFile = [' ../../../RugbyModel/'  model '_FULL_RRA_Tasks.xml'];%[data.Name '_RRA_Tasks.xml'];
    % % % % % % %     RRAForceFile = [' ../../../RugbyModel/' model '_FULL_RRA_Actuators.xml'];%[data.Name '_RRA_Actuators.xml'];
    % % % % % % %     %RRAConstraintsFile = [osim_path model '_RRA_ControlConstraints.xml'];%[data.Name '_RRA_ControlConstraints.xml'];
    % % % % % % %     RRADirectoryName = [pname 'RRAResults\'];
    % % % % % % %     ModelFile = [pname participant_code '_SCALED.osim'];
    % % % % % % %
    % % % % % % %     % Write GRF XML file
    % % % % % % %     grf2xml(data,'ExternalLoadNames',{'ExternalForce_1','ExternalForce_2','Scrum_Machine_1','Scrum_Machine_2'},...
    % % % % % % %         'AppliedToBodies',{'toes_r','toes_l','rscapula','lscapula'},...
    % % % % % % %         'ForceExpressedinBody',{'ground','ground','ground','ground'},...
    % % % % % % %         'PointExpressedinBody',{'ground','ground','ground','ground'},...
    % % % % % % %         'GRFFile',GRFFile,'MOTFile',...
    % % % % % % %         MOTFile,'LowPassFilterForKinematics',10,'OutputFile',GRFFileXML);
    % % % % % % %
    % % % % % % %     % % % grf2xml(data,'ExternalLoadNames',{'ExternalForce_1','ExternalForce_2'},...
    % % % % % % %     % % %     'AppliedToBodies',{'calcn_r','calcn_l'},'GRFFile',GRFFile,'MOTFile',...
    % % % % % % %     % % %     MOTFile,'LowPassFilterForKinematics',6,'OutputFile',GRFFileXML);
    % % % % % % %
    % % % % % % %     % Set up the XML file
    % % % % % % %     RRA = setup_ReduceResiduals('data',data, 'ModelFile',ModelFile,...
    % % % % % % %         'MOTFile',MOTFile,'GRFFile',GRFFileXML,'RRATaskFile',RRATaskFile,...
    % % % % % % %         'RRAForceFile',RRAForceFile,...
    % % % % % % %         'AdjCOMRes','true','LowPassFilterFreq',12,'DirectoryName',RRADirectoryName,...
    % % % % % % %         'InitialTime',data.time(1),'FinalTime',data.time(end));
    % % % % % % %
    % % % % % % %     % Call the RRA tool
    % % % % % % %     clear com
    % % % % % % %     com = ['RRA -S ' participant_code '_setup_ReduceResiduals.xml'];
    % % % % % % %     system(com);
    % % % % % % %
    % % % % % % %
    % % % % % % %     clear D
    % % % % % % %     % load the data from RRA
    % % % % % % %     D = load_sto_file([RRADirectoryName  participant_code '_' trial_code '_RRA_Actuation_force.sto']);
    % % % % % % %     data.RRA.Kinetics = D;
    % % % % % % %
    % % % % % % %     clear D
    % % % % % % %     % load the data from RRA
    % % % % % % %     D = load_sto_file([RRADirectoryName  participant_code '_' trial_code '_RRA_states.sto']);
    % % % % % % %     data.RRA.States = D;
    % % % % % % %     % % % % save([(subdir_current) '.mat']);
    % % % % % % %
    % % % % % % %     % Make plots for quick data check of residuals
    % % % % % % %     subplot(2,1,1), plot(data.Kinetics.time,[data.Kinetics.pelvis_tx_force data.Kinetics.pelvis_ty_force data.Kinetics.pelvis_tz_force]); hold on
    % % % % % % %     subplot(2,1,2), plot(data.Kinetics.time,[data.Kinetics.spine_list_moment data.Kinetics.spine_rotation_moment data.Kinetics.spine_tilt_moment]); hold on
    % % % % % % %     subplot(2,1,1), plot(data.RRA.Kinetics.time,[data.RRA.Kinetics.FX data.RRA.Kinetics.FY data.RRA.Kinetics.FZ],'LineWidth',2); hold off
    % % % % % % %     subplot(2,1,2), plot(data.RRA.Kinetics.time,[data.RRA.Kinetics.MX data.RRA.Kinetics.MY data.RRA.Kinetics.MZ],'LineWidth',2); hold off
    % % % % % % %
    % % % % % % %
    % % % % % % %     %% ADJUST THE MASS OF THE MODEL - superseded by strengthScaler routine
    % % % % % % %     % At this stage this is done manually, record the adjusted mass from RAA
    % % % % % % %     % and write to the model for future analysis
    % % % % % % %     %adjust_model_mass;
    % % % % % % %
    % % % % % % %     %% Define file name prefix structure
    % % % % % % %     % Added a nameing structure to avoid confusion throughout the script,
    % % % % % % %     % particpant code refers to the actual participant number (eg '066')
    % % % % % % %     % where as trial code refers to the trial (eg '20_4')
    % % % % % % %
    % % % % % % %     participant_code = data.Name;
    % % % % % % %
    % % % % % % %     trial_code = fname(1:end-4);
    % % % % % % %
    % % % % % % %
    % % % % % % %     %% ADJUST MODEL PROPERTIES ROUTINE
    % % % % % % %     % Using new API functions
    % % % % % % %     InputModel = [pname participant_code '_RRA_adjusted.osim'];
    % % % % % % %     ModelFile = [pname participant_code '_RRA_adjusted_mass.osim'];
    % % % % % % %
    % % % % % % %     % Make sure the strength scaling factor is set correctly
    % % % % % % %     % CAREFUL - withing strengthScalre routine there is currently hard coded in
    % % % % % % %     % the recommended segment mass changes which come out of the RRA routine so
    % % % % % % %     % need to make sure these segment masses are the desired one for the
    % % % % % % %     % participants being analysed.
    % % % % % % %
    % % % % % % %     strengthScalerGT(3.0,InputModel,ModelFile);
    % % % % % % %
    % % % % % % %
    % % % % % % %     %% INVERSE DYNAMIC ANALYSIS 2
    % % % % % % %     % Define the standard OpenSim files for ID
    % % % % % % %     MOTFile = [pname participant_code '_states.sto'];
    % % % % % % %     GRFFile = [pname trial_code '_grf.mot'];
    % % % % % % %     GRFFileXML = [pname trial_code '_grf.xml'];
    % % % % % % %     OutputFile = [trial_code '_ID2.sto'];
    % % % % % % %     ResultsDirectory = [pname 'ID\'];
    % % % % % % %
    % % % % % % %     % Write GRF XML file
    % % % % % % %     grf2xml(data,'ExternalLoadNames',{'ExternalForce_1','ExternalForce_2','Scrum_Machine_1','Scrum_Machine_2'},...
    % % % % % % %         'AppliedToBodies',{'toes_r','toes_l','rscapula','lscapula'},...
    % % % % % % %         'ForceExpressedinBody',{'ground','ground','ground','ground'},...
    % % % % % % %         'PointExpressedinBody',{'ground','ground','ground','ground'},...
    % % % % % % %         'GRFFile',GRFFile,'MOTFile',...
    % % % % % % %         MOTFile,'LowPassFilterForKinematics',12,'OutputFile',GRFFileXML);
    % % % % % % %     %
    % % % % % % %     % grf2xml(data,'ExternalLoadNames',{'ExternalForce_1','ExternalForce_2'},...
    % % % % % % %     %     'AppliedToBodies',{'calcn_r','calcn_l'},'GRFFile',GRFFile,'MOTFile',...
    % % % % % % %     %     MOTFile,'LowPassFilterForKinematics',6,'OutputFile',GRFFileXML);
    % % % % % % %
    % % % % % % %     % Set up the XML file
    % % % % % % %     setup_InverseDynamics('data',data, 'ModelFile',ModelFile,'MOTFile',MOTFile,...
    % % % % % % %         'GRFFile',GRFFileXML,'LowPassFilterFreq',12, 'OutputFile', OutputFile, 'ResultsDirectory', ResultsDirectory);
    % % % % % % %
    % % % % % % %     % Call the ID tool
    % % % % % % %     com = ['id -S ' participant_code '_Setup_InverseDynamics.xml'];
    % % % % % % %     system(com);
    % % % % % % %
    % % % % % % %     clear D
    % % % % % % %
    % % % % % % %     % % % save([(subdir_current) '.mat']');
    % % % % % % %
    % % % % % % %     clear D
    % % % % % % %
    % % % % % % %     %% MUSCLE ANALYSIS 2
    % % % % % % %     % Re-Run to determine if correction has been effective
    % % % % % % %
    % % % % % % %     CoordinatesFile = [pname participant_code '_states.sto'];
    % % % % % % %     ForceSetFile = ['../../../RugbyModel/' 'Rugby_Model_FULL_SO_Actuators.xml'];
    % % % % % % %     ExternalLoadsFile = GRFFileXML;
    % % % % % % %     MADirectoryName = [pname 'MA_Results2\'];
    % % % % % % %     MOTFile = CoordinatesFile;
    % % % % % % %     ModelFile = ModelFile;
    % % % % % % %
    % % % % % % %     % Set up the XML file
    % % % % % % %     setup_muscle_analysisDG('data',data, 'ModelFile',ModelFile,...
    % % % % % % %         'ExternalLoadsFile',ExternalLoadsFile,'CoordinatesFile',CoordinatesFile,...
    % % % % % % %         'MOTFILE',MOTFile,'ForceSetFile',ForceSetFile,'UseFLV', true,...
    % % % % % % %         'StartTime',data.time(1),'EndTime',data.time(end),...
    % % % % % % %         'LowPassFilterFreq', 12, 'DirectoryName', MADirectoryName);
    % % % % % % %
    % % % % % % %     % Call the Analyse tool
    % % % % % % %     clear com
    % % % % % % %     com = ['analyze  -S '  participant_code '_Setup_Muscle_Analysis.xml'];
    % % % % % % %     system(com);
    % % % % % % %
    % % % % % % %     % % % % % save([(subdir_current) '.mat']');
    % % % % % % %
    % % % % % % %     % STATIC OPTIMIZATION ANALYSIS
    % % % % % % %     % SetUp Static Optimazation
    % % % % % % %     % Define the standard OpenSim files for Static Optimization (Setup file
    % % % % % % %     % has been altered to call the Static Optimization tool rather than the
    % % % % % % %     % muscle analysis, added the suffix "DG" to the setup file)
    % % % % % % %
    % % % % % % %     % Write a new GRF XML file
    % % % % % % %     SO_GRFFileXML = [pname trial_code '_grf_SO.xml'];
    % % % % % % %     MOTFile = [pname trial_code '_ik.mot'];
    % % % % % % %     GRFFile = [pname trial_code '_grf.mot'];
    % % % % % % %
    % % % % % % %     grf2xml(data,'ExternalLoadNames',{'ExternalForce_1','ExternalForce_2','Scrum_Machine_1','Scrum_Machine_2'},...
    % % % % % % %         'AppliedToBodies',{'toes_r','toes_l','rclavicle','lclavicle'},...
    % % % % % % %         'ForceExpressedinBody',{'ground','ground','rclavicle','lclavicle'},...
    % % % % % % %         'PointExpressedinBody',{'ground','ground','rclavicle','lclavicle'},...
    % % % % % % %         'GRFFile',GRFFile,'MOTFile',...
    % % % % % % %         MOTFile,'LowPassFilterForKinematics',8,'OutputFile',SO_GRFFileXML);
    % % % % % % %
    % % % % % % %     % grf2xml(data,'ExternalLoadNames',{'ExternalForce_1',...
    % % % % % % %     %     'ExternalForce_2'},'AppliedToBodies',{'calcn_r','calcn_l'},...
    % % % % % % %     %     'GRFFile',GRFFile,'MOTFile',RRAMOTFile,'LowPassFilterForKinematics',-1,...
    % % % % % % %     %     'OutputFile',SO_GRFFileXML);
    % % % % % % %
    % % % % % % %     %
    % % % % % % %     StatesFile = [pname trial_code '_ik.mot'];
    % % % % % % %     ForceSetFile = ['../../../RugbyModel/' 'Rugby_Model_FULL_SO_Actuators.xml'];
    % % % % % % %     ExternalLoadsFile = SO_GRFFileXML;
    % % % % % % %     SODirectoryName = [pname 'SO_Results\'];
    % % % % % % %     ModelFile = [pname participant_code '_SCALED.osim'];
    % % % % % % %
    % % % % % % %     %Set up the XML file
    % % % % % % %     static_optim = setup_StaticOptim('data',data,...
    % % % % % % %         'ModelFile',ModelFile,'ExternalLoadsFile',ExternalLoadsFile,'StatesFile',StatesFile,...
    % % % % % % %         'ForceSetFile',ForceSetFile,...
    % % % % % % %         'UseFLV', true, 'StartTime',data.time(1),'EndTime',data.time(end),...
    % % % % % % %         'LowPassFilterFreq', -1, 'DirectoryName', SODirectoryName);
    % % % % % % %
    % % % % % % %     %Call the Static Optimisation Tool
    % % % % % % %     clear com
    % % % % % % %     com = ['analyze  -S '  participant_code '_Setup_StaticOptim.xml'];
    % % % % % % %     system(com);
    % % % % % % %
    % % % % % % %     addpath(genpath([pname 'SOResults\']));
    % % % % % % %     % % % % % save([(subdir_current) '.mat']');
    % % % % % % %
    % % % % % % %     % %% INDUCED ACCELERATION ANALYSIS - Static Optimisation (SO)
    % % % % % % %     % % StaticOptimisation Pathway: Setup IndAccPl file
    % % % % % % %     % %% Update GRF and GRFXML for IAA
    % % % % % % %     % %We need to update the GRF file to remove the ground tourqe value
    % % % % % % %     % %and update the xml
    % % % % % % %     %
    % % % % % % %     % %MOTFile = [pname participant_code '_' trial_code '_' 'RRA_Kinematics_q_adjusted.sto'];
    % % % % % % %     % %MOTFile = [pname trial_code '_ik.mot'];
    % % % % % % %     % MOTFile = [pname 'RRAResults\' participant_code '_' trial_code '_' 'RRA_states.sto'];
    % % % % % % %     % GRFFile_IAA = [pname trial_code '_IndAcc_grf.mot'];
    % % % % % % %     % IAAGRFFileXML = [pname trial_code '_grf_IAA.xml']; % new GRF xml file
    % % % % % % %     % grf2xml_IAA(data,'ExternalLoadNames',{'ExternalForce_1','ExternalForce_2'},...
    % % % % % % %     %     'AppliedToBodies',{'calcn_l','calcn_r'},'GRFFile',GRFFile_IAA,...
    % % % % % % %     %     'MOTFile',MOTFile,'LowPassFilterForKinematics',-1,'OutputFile',IAAGRFFileXML);
    % % % % % % %     %
    % % % % % % %     % %%
    % % % % % % %     % %CoordinatesFile = [pname participant_code '_' trial_code '_' 'RRA_Kinematics_q_adjusted.sto'];% DOWN SAMPLED STATES FILE
    % % % % % % %     % %CoordinatesFile = [pname trial_code '_ik.mot'];
    % % % % % % %     % ModelFile = ModelFile;
    % % % % % % %     % ForceSetFiles = [osim_path 'FullBodyModelGT_SO_Actuators.xml'];
    % % % % % % %     % CoordinatesFile = MOTFile;
    % % % % % % %     % ExternalLoadsFile = IAAGRFFileXML;
    % % % % % % %     % KineticsFile = [pname trial_code '_IndAcc_grf.mot'];
    % % % % % % %     % ForcesFile = [pname 'SO_Results\' participant_code '_' trial_code '_SO_StaticOptimization_force.sto'];
    % % % % % % %     % %ForcesFile = [pname  participant_code '_' trial_code '_SO_StaticOptimization_force_adjusted.sto'];
    % % % % % % %     % SOIAResultsDirectory_PI = [pname 'SO_IA_Results_PI\'];
    % % % % % % %     %
    % % % % % % %     % % Set up the XML file
    % % % % % % %     % IndAccPI = setup_InducedAccelerationPlugIn('data',data,...
    % % % % % % %     %     'ModelFile',ModelFile, 'ForceSetFiles',ForceSetFiles,...
    % % % % % % %     %     'CoordinatesFile',CoordinatesFile, 'ExternalLoadsFile',...
    % % % % % % %     %     ExternalLoadsFile,'DirectoryName',SOIAResultsDirectory_PI,...
    % % % % % % %     %     'KineticsFile', KineticsFile, 'ForcesFile', ForcesFile,...
    % % % % % % %     %     'LowPassFilterFreq', -1, 'InitialTime',data.time(1),...
    % % % % % % %     %     'FinalTime',data.time(end));
    % % % % % % %     %
    % % % % % % %     % % Call the IAA Tool
    % % % % % % %     % clear com;
    % % % % % % %     % com = ['Analyze -L IndAccPI -S ' participant_code '_Setup_InducedAccelerationPlugIn.xml'];
    % % % % % % %     % system(com);
    % % % % % % %     %
    % % % % % % %     % %     clear com
    % % % % % % %     % %     com = ['Analyze -S ' pname 'Y_010_Setup_IAA_Roll.xml'];
    % % % % % % %     % %     system(com);
    % % % % % % %     %
    % % % % % % %     % save([(subdir_current) '.mat']');
    
else % Scaling process for the static trial
    
    data = btk_c3d2trc_v5(data, 'off');
    
    %  %** Hard coding mass and height at the moment AND again below in dynamic trials **
    data.Mass=8;
    data.Height=100;
    
    
    % Define the starndard OpenSim files for scaling
    ModelFile = [osim_path model '.osim'];
    MeasurementSetFile = [osim_path model '_FULL_Scale_MeasurementSet.xml'];
    ScaleTasksFile = [osim_path model '_FULL_Scale_Tasks.xml'];
    
    
    setup_scale('data',data,'ModelFile',ModelFile,'ScaleTasksFile',ScaleTasksFile,...
        'MeasurementSetFile',MeasurementSetFile, 'PreserveMass','true');
    
    com = ['scale -S ' data.Name '_Setup_Scale.xml'];
    system(com);
    
    % % % % %     save('Scale.mat');
    
end



