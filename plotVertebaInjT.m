clear all
clc

dataDir='C:\Users\User\Documents\Pavlos\UniversityofBath\PhD\Study_3\BUCKLING DATA\';

spines={'cs1', 'cs2', 'cs3', 'cs4', 'cs5', 'cs6', 'cs7', 'cs8'};
vertebras={'C2', 'C3', 'C4', 'C5'};
XYZ={'dx' 'dy' 'dz'};
col=linspecer(8,'qualitative');

timeVect=-5:0.25:30;

load([dataDir 'load_vert_inj_data']);

for iVert=1:4
    vert=vertebras{iVert};
    injFrameNo=round(21+timeInjury*4000); % get frame at 4 kHz
    for iXYZ=1:3
        xyz=XYZ{iXYZ};
        for iPlot=1:numel(spines)
            plot(timeVect,iVertData{iXYZ,iVert}(2:end,iPlot))
            hold on
            vline(timeVect(21),'g', 'Impact!');
            vline(timeVect(injFrameNo(iPlot)),'r', 'Injury');
            title([spines{iPlot} ' ' vertebras{iVert} ' ' XYZ{iXYZ}]);
            xlabel('Time (ms)');
            ylabel('Displacement (mm)');
            saveas (gcf, [dataDir 'Plots\Injury\' spines{iPlot} '_' vertebras{iVert} '_' XYZ{iXYZ} ] , 'tif');
            close
        end
        for iSpine=1:8
        plot(timeVect,iVertData{iXYZ,iVert}(2:end,iSpine),'color',col(iSpine,:))
        hold on
        end
        vline(timeVect(21),'g', 'Impact!');
        vline(timeVect(injFrameNo(1)),'r', '1');
        vline(timeVect(injFrameNo(2)),'r', '2');
        vline(timeVect(injFrameNo(4)),'r', '4');
        vline(timeVect(injFrameNo(6)),'r', '6');
        vline(timeVect(injFrameNo(7)),'r', '7');
        vline(timeVect(injFrameNo(8)),'r', '8');
        
        saveas (gcf, [dataDir 'Plots\Injury\' vertebras{iVert} '_' XYZ{iXYZ} '_AllSpines'] , 'tif');
        close
     end
    
end