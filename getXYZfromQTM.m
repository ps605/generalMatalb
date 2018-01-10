function [ labels ] = getXYZfromQTM( qtmFile )
% Returns all labeled marker trajectories from exported .mat file from QTM
% system. Axis have to be changed at the bottom of this file for respective
% lab and software use
%IN:
%   qtmFile : pass qtm file location

%OUT:
%   2xN cell: Row 1 trajectory label name
%             Row 2 trajectory XYZ    
% Pavlos Silvestros, University of Bath, 2017

 load(qtmFile);

[~,fileName,~]=fileparts(qtmFile);

eval(['labels=qtm_' fileName '.Trajectories.Labeled.Labels;']);
eval(['data=qtm_' fileName '.Trajectories.Labeled.Data;']);

for iLabel=1:numel(labels)
    dataToChange =data(iLabel,1:3,:);
    
    dataToChange=reshape(dataToChange,[3,length(data)])';
    
    
    
%   Change QTM XYZ to Opensim XYZ (X=X, Y=Z Z=-Y) and m to mm

    labels{2,iLabel}(:,1)= dataToChange(:,1)/1000;
    labels{2,iLabel}(:,2)= dataToChange(:,3)/1000;
    labels{2,iLabel}(:,3)= -dataToChange(:,2)/1000;
end
end


