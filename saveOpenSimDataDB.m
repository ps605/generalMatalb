% Script to read in data from OpenSim outputs and store them as db for
% easier access and handling
%
% Pavlos Silvestros, University of Bath, 2017
clear
clc

% Pull in the modeling classes
import org.opensim.modeling.*

%% SETUP
% Directory of where data is stored in named simulation files
mainDir='C:\Users\ps605\Documents\PhD\Study_3\Modelling\Outputs\API\Results\OPC2C4_DR1_12\';
modelPath='C:\Users\ps605\Documents\PhD\Study_3\Modelling\Inputs\Models\'; %RugbyModel_AUG2017_OPC2C4_DR1_12_HA00MA200

% Conditions, data and variables of interest

% List condition wildcards
conditions={'MA050DH400'...
    'MA075DH400'...
    'MA100DH400'};

% List output data file names to look for OpenSim files of interest
forceRepStr='ForceReporter';
stateStr='states';

% List variables of interest
frVarsY={'c6Bush.cerv6.force.Y',...
    'c5Bush.cerv5.force.Y',...
    'c4Bush.cerv4.force.Y'};

frVarsX={'c6Bush.cerv6.force.X',...
    'c5Bush.cerv5.force.X',...
    'c4Bush.cerv4.force.X'};

bushings={'c6Bush.cerv6',...
    'c5Bush.cerv5',...
    'c4Bush.cerv4'};

%% LOOP CONDITIONS
%% Adjusted 20171212 for force transform problem.
for iBush=1:3
    
    % create bushing variable strings to pass to transformForce func
    bushingVariables{1}=[bushings{iBush} '.force.X'];
    bushingVariables{2}=[bushings{iBush} '.force.Y'];
    bushingVariables{3}=[bushings{iBush} '.force.Z'];
    
    bushingParent=bushings{iBush}(end-4:end);
    
    
    % Get structures that contain folders with condition wildcards
    for iCond=1:numel(conditions)
        
        % Go to main data directory
        cd([mainDir '\Simulations\']);
        
        % List simulation folder of condition specified from wildcard
        condFolders=dir(['*' conditions{iCond} '*']);
        
        for iFold=1:numel(condFolders)
            
            % Move into folder
            folder=[mainDir '\Simulations\' condFolders(iFold).name];
            cd(folder)
            [~,condition,~]=fileparts(folder);
            
            % Get data of interest
            frFile=dir(['*' forceRepStr '*']);
            stFile=dir(['*' stateStr '*']);
            
            %% OPENSIM STUFF
            
            % Get model of this condition that corresponds to file
            modelFile=[modelPath 'RugbyModel_AUG2017_OPC2C4_DR1_12_' condition(1:(end-5)) '.osim'];
            
            % Get model
            model=Model(modelFile);
            
            % Get bodies of interest
            bodyFrom=model.getGroundBody();
            bodyTo=model.getBodySet.get(bushingParent);
            
%             % Use importdata.m and store them in column cell objects
%             frData{iFold,1}=importdata([folder '\' frFile(1).name]);
            
            % Transform force and store new data 
            frData{iFold,1}=transfromForce(model,[stFile(1).folder '\' stFile(1).name],...
                [frFile(1).folder '\' frFile(1).name], bushingVariables, bodyFrom, bodyTo);
        end
        
        % Save condition trials in .mat file
        save([mainDir 'ForceReporter\' conditions{iCond} '_' bushings{iBush} '_transformed.mat'],'frData', 'condition');
    end
end


