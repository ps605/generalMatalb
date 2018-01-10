function f = ObjectiveFunction(coeffs,params,exp_data,analyzeTool,trialdirv)
% Runs model simulation for gvien coefficients, returns integration
%	coeffs 	= initial set of control values
%	params 	= optimization parameters; simulation parameters and
%			  pointers to instantiated OpenSim objects.

% Import OpenSim modeling classes
import org.opensim.modeling.*

% Get access to step counter and current best velocity value.
% global stepCount bestSoFar;

% Prep correct folders

dir_data='\\Mac\Home\Desktop\Work\Research\OpenSim\Projects\DummyHead\ContactModel\DATA\';
pname=trialdirv(3:8);
participant_code=[trialdirv(3:6) '_' trialdirv(7:8)];
trial_code=['trial_' participant_code];


osim_path = '\\Mac\Home\Desktop\Work\Research\OpenSim\Projects\DummyHead\ContactModel\DATA\SCALED_Model\';
model = 'DummyHead_SCALED';
modelFile = [osim_path model '.osim'];

% Get a reference to the model, as the model doesn't change in this example.
osimModel = params.model;
states = params.state;


contact_model=HuntCrossleyForce.safeDownCast(osimModel.getForceSet().get('shoulder_CM_r'));

contact_model.setStiffness(coeffs(1));
contact_model.setDissipation(coeffs(2));

% Perform FD simulation with updated contact parameter(s)
% % % % %     %% Perform FD simulation with updated contact parameter(s)
% % % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %     % %Forward Dynamics
% % % % %
% % % % %     % % Define the standard OpenSim files for FD
% % % % %
% % % % %     ModelFile = [pname participant_code '_SCALED.osim'];
% % % % %     OutputFile = [trial_code '_FD.sto'];
% % % % %     ResultsDirectory = [pname 'FD\'];
% % % % %     FDControlsFile=[pname 'ID\' 'trial_DHNT02_ID.sto'];
% % % % %
% % % % %     % Set up the XML file
% % % % %     setup_ForwardDynamics('data',data, 'ModelFile',ModelFile,...
% % % % %         'InitialTime',0.1,'FinalTime',0.3,'AuxillaryStates',true,'FDControlsFile',FDControlsFile,...
% % % % %         'MinStepSize',1e-10,'MaxStepSize',0.0001,'ErrorTol',1e-8,...
% % % % %         'Analyses',ForceReporter,...
% % % % %         'OutputPrecision',20,'LowPassFilterFreq',6,'OutputFile', OutputFile,...
% % % % %         'ResultsDirectory',  ResultsDirectory, 'InputDirectory', pname);
% % % % %
% % % % %     % Call the ID tool
% % % % %     com = ['forward -S ' participant_code '_Setup_ForwardDynamics.xml'];
% % % % %
% % % % %     system(com);
% % % % %
% % % % %     clear D

%% FORCE REPORTER ANALYSIS
    % SetUp BodyKinematics
    % Define the standard OpenSim files for BodyKinematics

    
    analyzeTool.setModel(osimModel);  
    
    outfile = ['Setup_Analyze_' 'objf' '.xml'];
    analyzeTool.print([dir_data pname '\' participant_code outfile]);

    analyzeTool.run();
   
% Load the data from ID MOT file and tie kinetic data to data structure
    D = load_sto_file([dir_data pname '\ForceRep_Results\' trial_code '_ForceReporter_forces.sto']);
    data.Kinetics = D;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % % % % Get the contact force vector
% % % % % sim_data = [data.Kinetics.shoulder_CM_r_skull_force_X;...
% % % % %     data.Kinetics.shoulder_CM_r_skull_force_Y;...
% % % % %     data.Kinetics.shoulder_CM_r_skull_force_Z;...
% % % % % % %     data.Kinetics.shoulder_CM_r_skull_torque_X;...
% % % % % % %     data.Kinetics.shoulder_CM_r_skull_torque_Y;...
% % % % % % %     data.Kinetics.shoulder_CM_r_skull_torque_Z
% % % % % ];

% Module
for i=1:length(data.Kinetics.shoulder_CM_r_skull_force_X(1:end-15))
   sim_data(i)= sqrt(data.Kinetics.shoulder_CM_r_skull_force_X(i)^2+...
      data.Kinetics.shoulder_CM_r_skull_force_Y(i)^2+...
      data.Kinetics.shoulder_CM_r_skull_force_Z(i)^2);
end

sim_data=abs(sim_data);

% Compute the Least Square Differences between simulated data and
% experimental data
f = sum((sim_data(1:50) - exp_data(1:50)).^2);

% % % figure
% % % plot(sim_data); hold on
% % % plot((exp_data),'LineWidth',2);


% % Update stepCount
% 	stepCount = stepCount + 1;
% 
% % Check step counter, save results when appropriate
%     if f < bestSoFar
% 
%         bestSoFar = f;
%         disp(strcat('Optimization (',num2str(stepCount),') f = ',num2str(f),', bestSoFar = ',num2str(bestSoFar)));
%         disp(coeffs);
%     end
disp(coeffs);
end
