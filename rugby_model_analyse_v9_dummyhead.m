
%% This file is used to analyse an individual subjects tackling data
%  It requires additional mfiles and function files that can be found in
%  the OpenSim Folder.


%% Define some directory names, given global status so the directory can be called by various functions. 

close 
clear ; 
clc; 

global subdir;
global trialdir;
           
subdir = ['DummyHead'];%#ok<NBRAK>           

trialdir = ['NTLS05']; %#ok<NBRAK>

[nrow1 ~] = size(subdir);
[nrow2 ~] = size(trialdir);

for j = 1:nrow1
    for k = 1:nrow2

        %% Define the path and file names  
            pname = ['\\Mac\Home\Desktop\Work\Research\OpenSim\Projects\' subdir(j,:) '\' trialdir(k,:) '\'];    
                   path
            cd(pname)
            
            % Current sub and trial directory defined as global for downstream use in
            % calling in results files (mainly within the re-writing of
            % files in "strengthen_model" and "adjust_model_mass"
            global subdir_current;
            global trialdir_current;
            
            subdir_current = subdir(j,:);
            trialdir_current = trialdir(k,:);

            c3d_files = dir('*.c3d');
% % % % % %     TO BE UNCOMMENTED if you WANT TO SCALE!!!        
% % % % % %             % Process the static file to make the model in Opensim
% % % % % %             % first find out file name with 'cal' in it
% % % % % %             
% % % % % %             I = strfind(lower({c3d_files.name}),'cal');
% % % % % %             
% % % % % %             for i = 1:length(I)
% % % % % %                 if isempty(I{i})
% % % % % %                     I{i} = 0;
% % % % % %                 end
% % % % % %             end
% % % % % %             
% % % % % %             cal_file_num = find([I{:}]>0);
% % % % % %             cal_file = c3d_files(cal_file_num).name;
% % % % % %             
% % % % % %             % now run the opensim pipeline on this file to create the scaled model from
% % % % % %             % the FullBodyModel, typically held in the Osim Model folder on
% % % % % %             % the c drive
% % % % % %             
% % % % % %             disp('Running opensim pipeline to create scaled model ...')
% % % % % %             
% % % % % %             
% % % % % %             rugby_scrum_pipeline_v9_dummyhead(cal_file);
% % % % % %             
% % % % % %             disp('Done (cal file).')
            
            %% Process the Trial file 
            % If the trial directory has multiple magnitudes use the below try catch statements
            
            global dat_file
                            I = strmatch('trial',{c3d_files.name});
                            dat_file = c3d_files(I).name;

            data = rugby_scrum_pipeline_v9_dummyhead(dat_file); 
            

     
end
end


return
    
