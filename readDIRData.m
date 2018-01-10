% Read DIC Data from .csv files like the dataset from Tim

clear all;
clc;

dataDir='C:\Users\User\Documents\Pavlos\UniversityofBath\PhD\Study_3\BUCKLING DATA\';

spines={'cs1', 'cs2', 'cs3', 'cs4', 'cs5', 'cs6', 'cs7', 'cs8'};
vertebras={'C2', 'C3', 'C4', 'C5'};

for iSpine=1:numel(spines);
    spine=spines{iSpine};
%     cs=dir([dataDir spines{iSpine} '*']);
    for iVert=1:numel(vertebras);
        vert=vertebras{iVert};
        dicData{iSpine,iVert}=csvread([dataDir spine '_' vert '.csv'],1,0);
    end
end
