clear all;
clc;

dataDir='C:\Users\ps605\Documents\PhD\Study_3\BUCKLING DATA\';

spines={'cs1', 'cs2', 'cs3', 'cs4', 'cs5', 'cs6', 'cs7', 'cs8'};
vertebras={'C2', 'C3', 'C4', 'C5'};
XYZ={'dx' 'dy' 'dz'};
col=linspecer(8);

timeVect=-0.05:0.0025:0.3;

load([dataDir 'dicData']);
load([dataDir 'loadData']);
load([dataDir 'iVertData']);

[nRow, nCol] = cellfun(@size, dicData);

for iXYZ=5:7
    for iVert=1:4;
        
        iVertData{iXYZ-4,iVert}=nan(max(nRow(:,iVert)),8);
        iVertVel{iXYZ-4,iVert}=nan(max(nRow(:,iVert))-1,8);
        
        for iSpine=1:numel(spines);
            
            iVertData{iXYZ-4,iVert}(1:nRow(iSpine,iVert),iSpine)=dicData{iSpine,iVert}(:,iXYZ);
            iVertPosDiff=diff(iVertData{iXYZ-4,iVert}(:,iSpine));
            iVertVel{iXYZ-4,iVert}(1:length(iVertPosDiff),iSpine)=iVertPosDiff/4000;
            
        end
        
%         plot(timeVect,iVertData{iXYZ-4,iVert}(2:142,:))
%         legend(spines)
%         saveas(gcf,[dataDir 'Plots\Vertebras\' vertebras{iVert} '_' XYZ{iXYZ-4}], 'tif')
%         close all

    end
end








%     vertDisp{iSpine,1}=dicData{iSpine,1}-dicData{iSpine,2};
%     vertDisp{iSpine,2}=dicData{iSpine,2}-dicData{iSpine,3};
%     vertDisp{iSpine,3}=dicData{iSpine,3}-dicData{iSpine,4};
%     vertDisp{iSpine,4}=dicData{iSpine,1}-dicData{iSpine,4};