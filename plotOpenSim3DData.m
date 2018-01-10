% Script to create 3D plots from continuous OpenSim data
clear
close all
clc

mainDir='C:\Users\ps605\Documents\PhD\Study_3\Modelling\Outputs\API\Results\OPC2C4_DR1_12\ForceReporter\';
plotPath='C:\Users\ps605\Documents\PhD\Study_3\Modelling\Outputs\API\Results\OPC2C4_DR1_12\Plots\3D\';

% Data is stored previously in .mat databases
conditions={'MA100DH400',...
    'MA075DH400',...
    'MA050DH400'}; %

% List variables of interest
%     frVars={'c6Bush.cerv6.force.Y',...
%         'c5Bush.cerv5.force.Y',...
%         'c4Bush.cerv4.force.Y',...
%         'c6Bush.cerv6.force.X',...
%         'c5Bush.cerv5.force.X',...
%         'c4Bush.cerv4.force.X'};
%
% List variables of interest
frVars={'_c4Bush.cerv4_',...
    '_c5Bush.cerv5_',...
    '_c6Bush.cerv6_'};


for iCond=1:numel(conditions)
    
    for iVar=1:numel(frVars)
        
        var=frVars{iVar};
        
        load([mainDir conditions{iCond} frVars{iVar} 'transformed.mat']);
        
        % Initialise matrix
        %         dataMatrix=[];
        dataMatrix_X=[];
        dataMatrix_Y=[];
        dataMatrix_Z=[];
        
        for iTrial=1:numel(frData)
            
            %             varIndex = find(strcmp(frData{iTrial}.colheaders,var));
            %
            %             % Contatenate data to form Trial matrix *** NOTE - have to
            %             % manually cut down number of rows (ie time points) due to
            %             % variable timesteps in OpenSim
            %             dataMatrix=horzcat(dataMatrix,frData{iTrial}.data(1:500,varIndex));
            
            dataMatrix_X=horzcat(dataMatrix_X, frData{iTrial}(1:500,1));
            dataMatrix_Y=horzcat(dataMatrix_Y, frData{iTrial}(1:500,2));
            dataMatrix_Z=horzcat(dataMatrix_Z, frData{iTrial}(1:500,3));
            
        end
        
        %         %% PLOT
        %         figure(1)
        %
        %
        %
        %         if frVars{iVar}(end)=='Y'
        %
        %             surf(0:30,frData{1}.data(1:500,1),-dataMatrix,'EdgeColor','interp')
        %             view(70,10);
        %             zlabel('Compression (N)');
        %
        %         else
        %
        %             surf(0:30,frData{1}.data(1:500,1),dataMatrix,'EdgeColor','interp')
        %             view(70,10);
        %             zlabel('Shear (N)');
        %
        %         end
        %
        %         title([conditions{iCond} ' ' frVars{iVar}]);
        %         xlabel('Angle');
        %         ylabel('Time (s)')
        %
        %         % Save
        %         saveas(gcf,[plotPath  conditions{iCond} '_' frVars{iVar} '.tif'], 'tif');
        %% PLOTS FOR TRANSFORM
        
        %         %% Compression X
        %         figure (1);
        %         surfc(0:30,0.1:0.1:50,dataMatrix_X,'EdgeColor','interp')
        %         view(90,-90); % flat as color map
        %
        %         zlabel('Shear (N)');
        %         title([conditions{iCond} ' ' frVars{iVar}]);
        %         xlabel('Head Flexion (deg)');
        %         ylabel('Time (ms)')
        %
        %         colormap Jet
        %         x=colorbar;
        %         x.Label.String='Shear (N)';
        %
        %
        %         % Save
        % %         saveas(gcf,[plotPath  conditions{iCond} '_' var '_X.tif'], 'tif');
        % %         saveas(gcf,[plotPath  conditions{iCond} '_' var '_X.fig'], 'fig');
        %
        %         %% Compression Y
        %         figure (2);
        %         surfc(0:30,0.1:0.1:50,-dataMatrix_Y,'EdgeColor','interp')
        %         view(90,-90); % flat as color map
        %
        %         zlabel('Compression (N)');
        %         title([conditions{iCond} ' ' frVars{iVar}]);
        %         xlabel('Head Flexion (deg)');
        %         ylabel('Time (ms)')
        %
        %         colormap Jet
        %         y=colorbar;
        %         y.Label.String='Compression (N)';
        %
        %         % Save
        % %         saveas(gcf,[plotPath  conditions{iCond} '_' var '_Y.tif'], 'tif');
        % %         saveas(gcf,[plotPath  conditions{iCond} '_' var '_Y.fig'], 'fig');
        %
        % %         close all
        
        %% Sub plot of all 3 vertebraes
        % % %         figure(3)
        % % %
        % % %         subplot(2,3,iVar) % X
        % % %         surfc(0:30,0.1:0.1:50,dataMatrix_X,'EdgeColor','interp')
        % % %         view(90,-90); % flat as color map
        % % %         zlim([-10 280])
        % % %
        % % %         colormap Jet
        % % %         caxis([-10 280]);
        % % %
        % % %         if iVar==3
        % % %             col=colorbar;
        % % %             col.Location='east';
        % % %         end
        
        %         % Plot vs Time
        figure (3)
        subplot(2,3,iVar+3) % Y
        surfc(0:30,0.1:0.1:50,-dataMatrix_Y,'EdgeColor','interp')
        view(90,-90); % flat as color map
        zlim([-100 2100]);
        
        colormap Jet
        caxis([0 2100]);
        
        if iVar==3
            col=colorbar;
            col.Location='east';
%             col.TickLabels=[0 200 400 600 800 1000 ...
%                 1200 1400 1600 1800 2000 2200];
        end
        
        % Plot Peaks
        [peaks(:,iVar),peaksT(:,iVar)]=max(-dataMatrix_Y);
        
        figure (4)
        scatter(0:1:30,peaks(:,iVar))
        hold on
        
        % Plot C4, C5 and C6 vertebrea
        figure(5)
        plot(0.1:0.1:50,-dataMatrix_Y(:,1))
        hold on
        scatter(peaksT(1,iVar)*0.1,peaks(1,iVar),'r')
        
        
        
        
        %         zlabel('Shear (N)');
        %         title([conditions{iCond} ' ' frVars{iVar}]);
        %         xlabel('Head Flexion (deg)');
        %         ylabel('Time (ms)')
        %
        %         colormap Jet
        %         x=colorbar;
        %         x.Label.String='Shear (N)';
    end
    % Save
    saveas(figure(3),[plotPath  conditions{iCond} '_subplot.tif'], 'tif');
    saveas(figure(3),[plotPath  conditions{iCond} '_subplot.fig'], 'fig');
end