% Run Forward Dynamic Simulations for FYP 2017

clear all
clc

% Main Directory
dir='C:\Users\User\Documents\Pavlos\UniversityofBath\Year4\FYP\Modelling\';
%Directory for .xml FD set-up files
setUpFileDirectory=[dir 'Simulation_Inputs\FD_SetUpFiles\'];
% Model File
modelFile = [dir 'Master Files\Models\RModel_CutDown_061216_LoadCell.osim'];
[~,model,~]=fileparts(modelFile);

% Import OpenSim modeling classes
import org.opensim.modeling.*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Set Scenario conditions for F and EMG
forceScaleScenarios={[100 100] [80 80] [60 60] [50 50] [40 40]}; % cell array containing the scaling values for F and M in pairs of cells

extGrpScenarios={[20 70 30] [30 90 50] [40 100 60]};
fleGrpScenarios={[20 40 20] [40 60 30] [50 70 60]};

for iForceScen=1:numel(forceScaleScenarios);
    
    createForceScenario(forceScaleScenarios{iForceScen})
    
    for iControlScen=1:numel(extGrpScenarios);
        
        [extString, flexString]=createControlScenario(extGrpScenarios{iControlScen}, fleGrpScenarios{iControlScen});
        
%% Forward Dynamics

% % Define the standard OpenSim files for FD

outputFile = 'trial_1_FD.sto'; %needs counter/ what does it do?
resultsDirectory = [dir 'Simulation_Outputs\Trial_force' num2str(forceScaleScenarios{iForceScen}(1)) '_flexion' flexString ... 
    '_extension' extString '\'];%counter

% In future will call functions that will modify EMG(controller files) and
% Force data files to create scenarios
% for loop (controllers)
% for loop (forces)

% Controller .sto file
FDControlsFile=[dir 'Simulation_Inputs\Controllers\controls_flexion' flexString '_extension' extString '.sto']; % needs counter
% External Forces files (.xml - setup; .mot - Force data)
GRFFile=[dir 'Simulation_Inputs\Forces\forceSU_Scale' num2str(forceScaleScenarios{iForceScen}(1)) '.xml']; %needs counter
grfMOTFile=[dir 'Simulation_Inputs\Forces\forceSU_Scale' num2str(forceScaleScenarios{iForceScen}(1)) '.mot']; %needs counter

[~,scenarioCon,~]=fileparts(FDControlsFile);
[~,scenarioFor,~]=fileparts(GRFFile);


% Set up the XML file
setup_ForwardDynamics('MOTFile', grfMOTFile,'ModelFile', modelFile,'FDControlsFile', FDControlsFile,'GRFFile', GRFFile, ...
    'InitialTime',0.0,'FinalTime',0.2,'AuxillaryStates',true,...
    'MinStepSize',1e-10,'MaxStepSize',0.0001,'ErrorTol',1e-8,...
    'OutputPrecision',20,'LowPassFilterFreq',6,'OutputFile', outputFile,...
    'ResultsDirectory',  resultsDirectory, 'SetUpFileDirectory', setUpFileDirectory);

% Call the FD tool
com = ['forward -S ' setUpFileDirectory 'Setup_ForwardDynamics' scenarioFor(8:end) scenarioCon(9:end) '.xml'];

system(com);
    end
end