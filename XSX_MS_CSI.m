function [fX]=XSX_MS_CSI(trial,path)

% [fX]=INITIALISE_XSX_50Hz_v1(varargin)
%
% INITIALISE_XSX_50Hz_v1.m ... (to be completed).
%
% Input:
%       - varin1       put the characteristics of input variable #1 and
%                      the characteristics it must have
%       - varin2       put the characteristics of input variable #2 and
%                      the characteristics it must have
%
% Output:
%       - varout1      put the characteristics of output variable #1 and
%                      the characteristics it must have
%       - varout2      put the characteristics of output variable #1 and
%                      the characteristics it must have
%
%
% Recalled functions or nested functions ( >> put all the recalled functions)
%       - GUI_SELECTIMUs.m
%           + GUI_SELECTIMUs.fig
%       - FINDMAX.m
%       - ACC_COMP.m
%           + rot3D.m
%
%
% Useful References:
% 1) Preatoni, E. (2009) How designign a readable Matlab Function. Journal
%           of Wasting Time. 99(9),33-333.
% 2) AA.VV. (1840) Il perfetto tesista Bioingegnere. Casa Editrice
%           Inesistente, Milano. pp.1000-1430.
%
%
% Notes/Observation:
% Write here notes about the function, e.g. what has still to be
% done/changed, or if the function produces a file.
%
% Written by:
% -Ezio Preatoni-   University of Bath, October 2012
%
% Modified by:
% -Dario Cazzola-   University of Bath, November 2012 - changed import file
%                                                       and controls on imported files

%% input check
%check arguments number and content
unitslist= {
    '056','na','na','na','390','na','na','na','na',...
    'na','393','na','na','na','na','na','398','na'
    %'056','390','392','393','396','397','398',...
    %'405','406','407','408','409','410',...
    %'895','896','898','899','900','901','na'
    };

units={
    '056','na','na','na','390','na','na','na','na',...
    'na','393','na','na','na','na','na','398','na'
    %'056','390','392','393','396','397','398',...
    %'405','406','407','408','409','410',...
    %'895','896','898','899','900','901','na'
    };

narginchk(0,2)
if nargin==0 ...
        [file,path]= uigetfile('X:\Health\ResearchProjects\GTrewartha\RE-FH1097 - Cervical\Cervical_injury_mechanisms\Phase2_Continuation\DATA\Xsens\*.txt','please select Xsense text file');          %pressure calib file selection

   
end


%% variables declaration
%column headers
count= {'PacketCounter'}; %#ok<NASGU>
count_header= {'PacketCounter'};

A_header= {'Acc_X','Acc_Y','Acc_Z'};
AXYZ_header= {'Acc_X_XYZ','Acc_Y_XYZ','Acc_Z_XYZ'};
Aalpha_header= {'Acc_Xalpha','Acc_Yalpha'};
ACCnoG_header= {'ACCnoG_X','ACCnoG_Y','ACCnoG_Z'};
ACCG_header= {'ACCG_X','ACCG_Y','ACCG_Z'};
T_header= {'Roll','Pitch','Yaw'};
TO_header= {'RollOFF','PitchOFF','YawOFF'};
TxsensA_header= {'Angle_X_axes','Angle_Y_axes','Angle_Z_axes'};
Eu_rel_ANGLES_header= {'EUAngle_X_axes','EUAngle_Y_axes','EUAngle_Z_axes'};
t_header= {'Time'};
R_header= {'RSSI'};

%sampling frequency read from file (it might be different)

sf_string=importdata(strcat(path,file));
sf_c= char(sf_string.textdata(2));
sf=str2double(sf_c(17:18));

%timeline definition

        max_t= 5.5;
        min_t= -5.0;
        t= (min_t:1/sf:max_t)';    %timeline for condition 3
        t_rPAS= -0.9;                   %timeline for condition 1 (CTS)
  
%% data file selection and loading
%initialise variables
fileok= nan(1,numel(units));
filenames= cell(1,numel(units));
XSX= cell(1,numel(units));
A= nan(length(t),numel(A_header),numel(units));
A_4= nan(length(t),numel(A_header),numel(units));
A_3= nan(length(t),numel(A_header),numel(units));
A_89= nan(length(t),numel(A_header),numel(units));
AXYZ= nan(length(t),numel(AXYZ_header),numel(units));
Eu_rel_ANGLES=nan(length(t),numel(Eu_rel_ANGLES_header),numel(units));
Aalpha= nan(length(t),numel(Aalpha_header),numel(units));
T= nan(length(t),numel(A_header),numel(units));
T_4= nan(length(t),numel(A_header),numel(units));
T_3= nan(length(t),numel(A_header),numel(units));
T_89= nan(length(t),numel(A_header),numel(units));
AccNOg= nan(length(t),numel(ACCnoG_header),numel(units));
Accg= nan(length(t),numel(ACCG_header),numel(units));
T_off= nan(length(t),numel(TO_header),numel(units));
T_off_raw= nan(length(t),numel(TO_header),numel(units));
T_off_new= nan(length(t),numel(TO_header),numel(units));
TxsensA= nan(length(t),numel(TxsensA_header),numel(units));
Ax= nan(length(t),numel(units));
Ay= nan(length(t),numel(units));
Az= nan(length(t),numel(units));
AnEuTOTun= nan(length(t),3,numel(units));
Nx= nan(length(t),3,numel(units));
Nz= nan(length(t),3,numel(units));
Ny= nan(length(t),3,numel(units));
NyHE_rot1=nan(length(t),3,numel(units));
NyHE_rot2=nan(length(t),3,numel(units));
NyHE_rot3=nan(length(t),3,numel(units));
NyHE_rot=nan(length(t),3,numel(units));
NxHE_rot=nan(length(t),3,numel(units));
NxHE_rot1=nan(length(t),3,numel(units));
NxHE_rot2=nan(length(t),3,numel(units));
NzHE_rot=nan(length(t),3,numel(units));

angles= nan(numel(units),1);
alphaS_tACC= nan(length(t),1);
anglesTot= nan(size(t,1),numel(units));
t_ONM=nan(3,1);

Fo=nan(numel(units),1);
Th=nan(numel(units),1);
EiNi=nan(numel(units),1);
h= 1;

%import acceleration and orientation data from files
for i=1:numel(units);
    [fileok(i),temp]= fileattrib([path '*' units{i} '.txt']);    %look for suitable file in the given path
    if fileok(i) && isstruct(temp) && ~strcmpi(units{i},'na'),                  %if file is found...
        filenames{i}= temp.Name;
        XSX{i}= importdata(filenames{i});
        XSX{i}.filename= filenames{i};
        %takes missing values away from imported file using RSSI column
        RSSI= zeros(size(XSX{i}.data(:,ismember(XSX{i}.colheaders,R_header)),1));
        RSSI= XSX{i}.data(:,ismember(XSX{i}.colheaders,R_header));
        for v=1:size(RSSI,1),
            if RSSI(v)~= -128,
                switch units{i}(1)
                    case '4'
                        A_4(h,:,i)= XSX{i}.data(v,ismember(XSX{i}.colheaders,A_header));
                        T_4(h,:,i)= XSX{i}.data(v,ismember(XSX{i}.colheaders,T_header));
                        length_manual=length(A_4);
                        if h==1,
                            Fo(i)= XSX{i}.data(v,ismember(XSX{i}.colheaders,count_header));
                        end
                    case {'3','0'}
                        A_3(h,:,i)= XSX{i}.data(v,ismember(XSX{i}.colheaders,A_header));
                        T_3(h,:,i)= XSX{i}.data(v,ismember(XSX{i}.colheaders,T_header));
                        length_manual=length(A_3);
                        if h==1,
                            Th(i)= XSX{i}.data(v,ismember(XSX{i}.colheaders,count_header));
                        end
                    case {'8','9'}
                        A_89(h,:,i)= XSX{i}.data(v,ismember(XSX{i}.colheaders,A_header));
                        T_89(h,:,i)= XSX{i}.data(v,ismember(XSX{i}.colheaders,T_header));
                        length_manual=length(A_89);
                        if h==1,
                            EiNi(i)= XSX{i}.data(v,ismember(XSX{i}.colheaders,count_header));
                        end
                    otherwise
                        return
                end
                h=h+1;
            else
                if h>2          % if RSSI is in the middle of the file (not in the beginning),it leaves raw value (NaN)
                    switch units{i}(1)
                        case '4'
                            A_4(h,:,i)= nan(1,size(A_4(h,:,i),2));
                            T_4(h,:,i)= nan(1,size(A_4(h,:,i),2));
                            h=h+1;
                            length_manual=length(A_4);
                        case {'3','0'}
                            A_3(h,:,i)= nan(1,size(A_4(h,:,i),2));
                            T_3(h,:,i)= nan(1,size(A_4(h,:,i),2));
                            h=h+1;
                            length_manual=length(A_3);
                        case {'8','9'}
                            A_89(h,:,i)= nan(1,size(A_4(h,:,i),2));
                            T_89(h,:,i)= nan(1,size(A_4(h,:,i),2));
                            h=h+1;
                            length_manual=length(A_89);
                        otherwise
                            return
                    end
                end
            end
        end
        h=1;
    else                                                                        %otherwise...
        filenames{i}= 'na';
    end
    clear temp
end

[Fo_min,I_Fo]=min(Fo);
diff_Fo=Fo-Fo_min;

[Th_min,I_Th]=min(Th);
diff_Th=Th-Th_min;

[EiNi_min,I_EiNi]=min(EiNi);
diff_EiNi=EiNi-EiNi_min;

if length(A)<length(A_3) || length(A)<length(A_4)||  length(A)<length(A_89)

    A= nan(length_manual,numel(A_header),numel(units));
    T= nan(length_manual,numel(A_header),numel(units));
    AccNOg= nan(length_manual,numel(ACCnoG_header),numel(units));
Accg= nan(length_manual,numel(ACCG_header),numel(units));
T_off= nan(length_manual,numel(TO_header),numel(units));
T_off_raw= nan(length_manual,numel(TO_header),numel(units));
T_off_new= nan(length_manual,numel(TO_header),numel(units));
TxsensA= nan(length_manual,numel(TxsensA_header),numel(units));
Ax= nan(length_manual,numel(units));
Ay= nan(length_manual,numel(units));
Az= nan(length_manual,numel(units));
AnEuTOTun= nan(length_manual,3,numel(units));
Nx= nan(length_manual,3,numel(units));
Nz= nan(length_manual,3,numel(units));
Ny= nan(length_manual,3,numel(units));
NyHE_rot1=nan(length_manual,3,numel(units));
NyHE_rot2=nan(length_manual,3,numel(units));
NyHE_rot3=nan(length_manual,3,numel(units));
NyHE_rot=nan(length_manual,3,numel(units));
NxHE_rot=nan(length_manual,3,numel(units));
NxHE_rot1=nan(length_manual,3,numel(units));
NxHE_rot2=nan(length_manual,3,numel(units));
NzHE_rot=nan(length_manual,3,numel(units));
Aalpha= nan(length_manual,numel(Aalpha_header),numel(units));
alphaS_tACC= nan(length_manual,1);
AXYZ= nan(length_manual,numel(AXYZ_header),numel(units));
    t= (min_t:1/sf:(length_manual/sf+min_t-1/sf))';
end

for i=1:numel(units),
    if ~isnan(Fo(i)),
        if diff_Fo(i)~=0,
            Mnan=nan(diff_Fo(i),size(A_4,2));
            Atemp=[Mnan;A_4(:,:,i)];
            A(:,:,i)=Atemp(1:size(A(:,:,i),1),1:size(A(:,:,i),2));
            Ttemp=[Mnan;T_4(:,:,i)];
            T(:,:,i)=Ttemp(1:size(T(:,:,i),1),1:size(T(:,:,i),2));
            [AccNOg(:,:,i),Accg(:,:,i),T_off(:,:,i),Aalpha(:,:,i),AXYZ(:,:,i),Nx(:,:,i),Ny(:,:,i),Nz(:,:,i)]= ACC_COMP(T(:,:,i), A(:,:,i),angles(i),alphaS_tACC);
        else
            A(:,:,i)=A_4(:,:,i);
            T(:,:,i)=T_4(:,:,i);
            [AccNOg(:,:,i),Accg(:,:,i),T_off(:,:,i),Aalpha(:,:,i),AXYZ(:,:,i),Nx(:,:,i),Ny(:,:,i),Nz(:,:,i)]= ACC_COMP(T(:,:,i), A(:,:,i),angles(i),alphaS_tACC);
        end
    else
        if ~isnan(Th(i)),
            if diff_Th(i)~=0,
                Mnan=nan(diff_Th(i),size(A_3,2));
                Atemp=[Mnan;A_3(:,:,i)];
                A(:,:,i)=Atemp(1:size(A(:,:,i),1),1:size(A(:,:,i),2));
                Ttemp=[Mnan;T_3(:,:,i)];
                T(:,:,i)=Ttemp(1:size(T(:,:,i),1),1:size(T(:,:,i),2));
                [AccNOg(:,:,i),Accg(:,:,i),T_off(:,:,i),Aalpha(:,:,i),AXYZ(:,:,i),Nx(:,:,i),Ny(:,:,i),Nz(:,:,i)]= ACC_COMP(T(:,:,i), A(:,:,i),angles(i),alphaS_tACC);
            else
                A(:,:,i)=A_3(:,:,i);
                T(:,:,i)=T_3(:,:,i);
                [AccNOg(:,:,i),Accg(:,:,i),T_off(:,:,i),Aalpha(:,:,i),AXYZ(:,:,i),Nx(:,:,i),Ny(:,:,i),Nz(:,:,i)]= ACC_COMP(T(:,:,i), A(:,:,i),angles(i),alphaS_tACC);
            end   %End if Th
        else
            if ~isnan(EiNi(i)),
                if diff_EiNi(i)~=0,
                    Mnan=nan(diff_EiNi(i),size(A_89,2));
                    Atemp=[Mnan;A_89(:,:,i)];
                    A(:,:,i)=Atemp(1:size(A(:,:,i),1),1:size(A(:,:,i),2));
                    Ttemp=[Mnan;T_89(:,:,i)];
                    T(:,:,i)=Ttemp(1:size(T(:,:,i),1),1:size(T(:,:,i),2));
                    [AccNOg(:,:,i),Accg(:,:,i),T_off(:,:,i),Aalpha(:,:,i),AXYZ(:,:,i),Nx(:,:,i),Ny(:,:,i),Nz(:,:,i)]= ACC_COMP(T(:,:,i), A(:,:,i),angles(i),alphaS_tACC);
                else
                    A(:,:,i)=A_89(:,:,i);
                    T(:,:,i)=T_89(:,:,i);
                    [AccNOg(:,:,i),Accg(:,:,i),T_off(:,:,i),Aalpha(:,:,i),AXYZ(:,:,i),Nx(:,:,i),Ny(:,:,i),Nz(:,:,i)]= ACC_COMP(T(:,:,i), A(:,:,i),angles(i),alphaS_tACC);
                end
            end
        end
    end %end if Fo
end %End for

%calculate module of the acceleration NO g
ACCnoG_header= [ACCnoG_header, '|Acc|'];
AccNOg(:,strcmpi(ACCnoG_header,'|Acc|'),:)= (AccNOg(:,1,:).^2+...
    AccNOg(:,2,:).^2+AccNOg(:,3,:).^2).^0.5;

%calculate module of the acceleration
ACCG_header= [ACCG_header, '|Acc|'];
Accg(:,strcmpi(ACCG_header,'|Acc|'),:)= (Accg(:,1,:).^2+...
    Accg(:,2,:).^2+Accg(:,3,:).^2).^0.5;

%% ORIENTATION from GYROSCOPE
Tkin=nan(size(T_off,1),numel(units));
x=[1 0 0];
y=[0 1 0];
z=[0 0 1];
T_off_raw=T_off;

for i=1:3:numel(units),
    % Neck Angle from Xsens|
    
    TxsensA(:,1,i)=(T_off(:,1,i)-T_off(((t<(t_rPAS+0.0065))&(t>(t_rPAS-0.007))),1,i+1));
    TxsensA(:,2,i)=(T_off(:,2,i)-T_off(((t<(t_rPAS+0.0065))&(t>(t_rPAS-0.007))),2,i+1));
    TxsensA(:,3,i)=(T_off(:,3,i)-T_off(((t<(t_rPAS+0.0065))&(t>(t_rPAS-0.007))),3,i+1));
    
    % Neck Angle from KTR|
    %     Tkin(:,i)=(anglesTot(:,i)-anglesTot(((t_ktr<(t_rPAS+0.0001))&(t_ktr>(t_rPAS-0.0001))),i))-(anglesTot(:,i+1)-anglesTot(((t_ktr<(t_rPAS+0.0001))&(t_ktr>(t_rPAS-0.0001))),i+1));
        
    for j=1:size(Nx,1),
        
        % First rotation of Head from fixed (absolute) ref system to C7 ref
        % system - second axis
        NyHE_rot2(j,:,i)=(rot3D(T_off(j,:,i+1),0)*Ny(j,:,i)')/norm(rot3D(T_off(j,:,i+1),0)*Ny(j,:,i)');
        % We need to rotate the system 90 deg on y and then 90 deg on z
        % to match the functional angle between the anatomical
        % segments(Torssion,flex/ext,abd/add)
        
        NyHE_rot3(j,:,i)=(rot3D([0,90,0],0)*NyHE_rot2(j,:,i)')/norm(rot3D([0,90,0],0)*NyHE_rot2(j,:,i)');
        NyHE_rot1(j,:,i)=(rot3D([0,0,90],0)*NyHE_rot3(j,:,i)')/norm(rot3D([0,0,90],0)*NyHE_rot3(j,:,i)');
        
        % First rotation of Head from fixed (absolute) ref system to C7 ref
        % system - X axis
        NxHE_rot1(j,:,i)=(rot3D(T_off(j,:,i+1),0)*Nx(j,:,i)')/norm(rot3D(T_off(j,:,i+1),0)*Nx(j,:,i)');
        
        % We need to rotate the system 90 deg on y and then 90 deg on z
        % to match the functional angle between the anatomical
        % segments(Torssion,flex/ext,abd/add)
        NxHE_rot2(j,:,i)= (rot3D([0,90,0],0)*NxHE_rot1(j,:,i)')/norm(rot3D([0,90,0],0)*NxHE_rot1(j,:,i)');
        NxHE_rot(j,:,i)= (rot3D([0,0,90],0)*NxHE_rot2(j,:,i)')/norm(rot3D([0,0,90],0)*NxHE_rot2(j,:,i)');
        
        % new Z and Y axes calculation
        NzHE_rot(j,:,i)=cross(NyHE_rot1(j,:,i),NxHE_rot(j,:,i));
        NyHE_rot(j,:,i)=cross(NzHE_rot(j,:,i),NxHE_rot(j,:,i));
        
        
        if dot(NzHE_rot(j,:,i),x)< 1,
            if dot(NxHE_rot(j,:,i),y)> -1,
                Ay(j,i)=asin(dot(NzHE_rot(j,:,i),x));
                Ax(j,i)=atan2(dot(NzHE_rot(j,:,i),-y),dot(NzHE_rot(j,:,i),z));
                Az(j,i)=atan2(dot(NyHE_rot(j,:,i),-x),dot(NxHE_rot(j,:,i),x));
            elseif dot(NzHE_rot(j,:,i),x)== -1,
                Ay(j,i)=-pi/2;
                Ax(j,i)= -atan2(dot(NyHE_rot(j,:,i),-x),dot(NyHE_rot(j,:,i),y));
                Az(j,i)= 0;
            end
        elseif dot(NzHE_rot(j,:,i),x)== 1,
            Ay(j,i)=pi/2;
            Ax(j,i)= atan2(dot(NyHE_rot(j,:,i),-x),dot(NyHE_rot(j,:,i),y));
            Az(j,i)= 0;
        end       
    end
        
    AnEuTOTun(:,1,i)=unwrap(Ax(:,i),pi/4)/pi*180;    
    AnEuTOTun(:,2,i)=unwrap(Ay(:,i),pi/4)/pi*180;   
    AnEuTOTun(:,3,i)=unwrap(Az(:,i),pi/4)/pi*180;
    
    AnEuTOTun(:,1,i)=AnEuTOTun(:,1,i)-AnEuTOTun((t<(t_rPAS+0.0065))&(t>(t_rPAS-0.007)),1,i);
    AnEuTOTun(:,2,i)=AnEuTOTun(:,2,i)-AnEuTOTun((t<(t_rPAS+0.0065))&(t>(t_rPAS-0.007)),2,i);
    AnEuTOTun(:,3,i)=AnEuTOTun(:,3,i)-AnEuTOTun((t<(t_rPAS+0.0065))&(t>(t_rPAS-0.007)),3,i);
    
    % % % % % % % %                for k=1:3,
    % % % % % % % %                        for h=1:size(AnEuTOTun,1)-1,
    % % % % % % % %
    % % % % % % % %                                 if AnEuTOTun(h,k,i)-AnEuTOTun(h+1,k,i)>20,
    % % % % % % % %
    % % % % % % % %                                     AnEuTOTun(h+1,k,i)=AnEuTOTun(h+1,k,i)+20;
    % % % % % % % %                                 elseif  AnEuTOTun(h,k,i)-AnEuTOTun(h+1,k,i)<-20,
    % % % % % % % %                                      AnEuTOTun(h+1,k,i)=AnEuTOTun(h+1,k,i)-20;
    % % % % % % % %                                 end
    % % % % % % % %
    % % % % % % % %                         end
    % % % % % % % %                end
    
    %%     GYR CORRECTION ALGORITHM
    for k=1:3,
        for h=1:size(T_off,1)-1,
            
            for co=2.5:2.5:100,
                if T_off(h,k,i+1)-T_off(h+1,k,i+1)>co && T_off(h,k,i+1)-T_off(h+1,k,i+1)<co+2.5,
                    
                    T_off((h+1):size(T_off,1),k,i+1)=T_off((h+1):size(T_off,1),k,i+1)+(co);
                                                          
                else
                    if  T_off(h,k,i+1)-T_off(h+1,k,i+1)<-co && T_off(h,k,i+1)-T_off(h+1,k,i+1)>-co-2.5,
                                               
                        T_off((h+1):size(T_off,1),k,i+1)=T_off((h+1):size(T_off,1),k,i+1)-(co);
                                               
                    else
                             
                    end
                end
            end           
        end       
    end    
    
    %%
    % GRAPHICS
    % Z
        figure('units','normalized','position',[0.1 0.1 0.8 0.8]); hold on;
        % Angles C7 Sensor
    %     plot(t,T(:,3,i+1),'b.-');
          plot(t,T_off(:,2,i+1)-T_off(((t<(t_rPAS+0.0065))&(t>(t_rPAS-0.007))),2,i+1),'-.ok','linewidth',2);
          plot(t,T_off_raw(:,2,i+1)-T_off_raw(((t<(t_rPAS+0.0065))&(t>(t_rPAS-0.007))),2,i+1),'-.ro','linewidth',1);
    %      plot(t,T_off_new(:,3,i+1),'-.yo','linewidth',3);
    %       plot(t_ktr,anglesTot(:,i+1)-anglesTot(((t_ktr<(t_rPAS+0.0001))&(t_ktr>(t_rPAS-0.0001))),i+1),'-.ob');
    
    %      plot(t,(T_off(:,2,i+1)-T_off(((t<(t_rPAS+0.0001))&(t>(t_rPAS-0.0001))),2,i+1))-(T_off_raw(:,2,i+1)-T_off_raw(((t<(t_rPAS+0.0001))&(t>(t_rPAS-0.0001))),2,i+1)),'-.om','linewidth',1);
        % Angles Head Sensor
%         plot(t,T(:,3,1),'m.-');
%         plot(t,T_off(:,3,i)-T_off(((t<(t_rPAS+0.0065))&(t>(t_rPAS-0.007))),3,i),'m--','linewidth',2);
    
    %     plot(t_ktr,anglesTot(:,i)-anglesTot(((t_ktr<(t_rPAS+0.0001))&(t_ktr>(t_rPAS-0.0001))),i+1),'-.om');
        % Neck Angle from Xsens|
%         plot(t,TxsensA(:,3,i),'-.xr');
    
        % Neck Angle from KTR|
    %     plot(t_ktr,Tkin(:,i),'r.-','linewidth',5);
    
    
        box on;
    
% % % % % %         figure('units','normalized','position',[0.1 0.1 0.8 0.8]); hold on;
% % % % % %         % Angles C7 Sensor
% % % % % %         % plot(t,T(:,3,i+1),'b.-');
% % % % % %         % plot(t,T_off(:,3,i+1)-T_off(((t<(t_rPAS+0.0001))&(t>(t_rPAS-0.0001))),3,i+1),'b--','linewidth',2);
% % % % % %         % plot(t_ktr,anglesTot(:,i+1)-anglesTot(((t_ktr<(t_rPAS+0.0001))&(t_ktr>(t_rPAS-0.0001))),i+1),'-.ob');
% % % % % %     
% % % % % %         % Angles Head Sensor
% % % % % %         % plot(t,T(:,3,1),'m.-');
% % % % % %         % plot(t,T_off(:,3,i)-T_off(((t<(t_rPAS+0.0001))&(t>(t_rPAS-0.0001))),3,i),'m--','linewidth',2);
% % % % % %     
% % % % % %         plot(AnEuTOTun(:,1,i),'b.-','linewidth',3);
% % % % % %         plot(AnEuTOTun(:,2,i),'r--','linewidth',3);
% % % % % %         plot(AnEuTOTun(:,3,i),'m:','linewidth',3);
% % % % % %         plot(find((t<(t_rPAS+0.0065))&(t>(t_rPAS-0.007))),0,'-ko','LineWidth',3,'MarkerEdgeColor','k','MarkerSize',10);
% % % % % %         hleg1=legend('Torsion (Neg - Clockwise->Right)','Flex(Neg)/Ext(Pos)','Abduction (Neg - Clockwise ->Left)');
% % % % % %         set(hleg1,'Location','NorthEast');
% % % % % %         % Neck Angle from Xsens|
% % % % % %         % plot(t,TxsensA(:,3,i),'-.xr');
% % % % % %         %
% % % % % %         % % Neck Angle from KTR|
% % % % % %         % plot(t_ktr,Tkin(:,i),'r.-','linewidth',5);
% % % % % %         box on;
end



% % OFFSET SAGITTAL VIEW ?
% % % % % %  X
% % % figure('units','normalized','position',[0.1 0.1 0.8 0.8]); hraw on;
% % %
% % % plot(T(:,1,i+1),'b.-');
% % % plot(T_off(:,1,i+1),'b--','linewidth',2);
% % % plot(T(:,1,1),'m.-');
% % % plot(T_off(:,1,i),'m--','linewidth',2)
% % % plot(TxsensA(:,1,i),'c--','linewidth',5)
% % % box on;
% % % end
% % %

Eu_rel_ANGLES=AnEuTOTun(:,:,:);

%% identification t of actual engagement

if isnan(t_ONM(3)),
    f_i=1;
else
    f_i= find(t>=t_ONM(3),1);

end

[~,f_rENG_NOg]= FINDMAX(squeeze(AccNOg(:,4,:)),f_i*ones(size(AccNOg,3)));
f_rENG_NOg= f_rENG_NOg+f_i-1;
[~,f_XrENG_NOg]= FINDMAX(squeeze(abs(AccNOg(:,1,:))),f_i*ones(size(AccNOg,3)));
f_XrENG_NOg= f_XrENG_NOg+f_i-1;
[~,f_YrENG_NOg]= FINDMAX(squeeze(abs(AccNOg(:,2,:))),f_i*ones(size(AccNOg,3)));
f_YrENG_NOg= f_YrENG_NOg+f_i-1;
[~,f_ZrENG_NOg]= FINDMAX(squeeze(abs(AccNOg(:,3,:))),f_i*ones(size(AccNOg,3)));
f_ZrENG_NOg= f_ZrENG_NOg+f_i-1;

% f_rENG for verticall acc on fixed system
[~,f_Z_XYZrENG]= FINDMAX(squeeze(abs(AXYZ(:,3,:))),f_i*ones(size(AXYZ,3)));
f_Z_XYZrENG= f_Z_XYZrENG+f_i-1;

%%
[~,f_rENG]= FINDMAX(squeeze(Accg(:,4,:)),f_i*ones(size(Accg,3)));
f_rENG= f_rENG+f_i-1;
[~,f_XrENG]= FINDMAX(squeeze(abs(Accg(:,1,:))),f_i*ones(size(Accg,3)));
f_XrENG= f_XrENG+f_i-1;
[~,f_YrENG]= FINDMAX(squeeze(abs(Accg(:,2,:))),f_i*ones(size(Accg,3)));
f_YrENG= f_YrENG+f_i-1;
[~,f_ZrENG]= FINDMAX(squeeze(abs(Accg(:,3,:))),f_i*ones(size(Accg,3)));
f_ZrENG= f_ZrENG+f_i-1;

%% EULER rel ANGLES max,min,sus, -- aP (after pause) - aE (after Engagement) - pE (pre engagement)
% % % % % [maxTors_aP,~]=FINDMAX(squeeze(Eu_rel_ANGLES(find(t<(t_rPAS+0.0001)&(t>(t_rPAS-0.0001))):size(Eu_rel_ANGLES,1),1,:)));
% % % % % [minTors_aP,~]=FINDMIN(squeeze(Eu_rel_ANGLES(find(t<(t_rPAS+0.0001)&(t>(t_rPAS-0.0001))):size(Eu_rel_ANGLES,1),1,:)));
% % % % % [maxFlex_aP,~]=FINDMAX(squeeze(Eu_rel_ANGLES(find(t<(t_rPAS+0.0001)&(t>(t_rPAS-0.0001))):size(Eu_rel_ANGLES,1),2,:)));
% % % % % [minFlex_aP,~]=FINDMIN(squeeze(Eu_rel_ANGLES(find(t<(t_rPAS+0.0001)&(t>(t_rPAS-0.0001))):size(Eu_rel_ANGLES,1),2,:)));
% % % % % [maxAbd_aP,~]=FINDMAX(squeeze(Eu_rel_ANGLES(find(t<(t_rPAS+0.0001)&(t>(t_rPAS-0.0001))):size(Eu_rel_ANGLES,1),3,:)));
% % % % % [minAbd_aP,~]=FINDMIN(squeeze(Eu_rel_ANGLES(find(t<(t_rPAS+0.0001)&(t>(t_rPAS-0.0001))):size(Eu_rel_ANGLES,1),3,:)));
% % % % % 
% % % % % 
% % % % % [maxTors_aE,~]=FINDMAX(squeeze(Eu_rel_ANGLES(f_rENG+10:size(Eu_rel_ANGLES,1),1,:)));
% % % % % [minTors_aE,~]=FINDMIN(squeeze(Eu_rel_ANGLES(f_rENG+10:size(Eu_rel_ANGLES,1),1,:)));
% % % % % [maxFlex_aE,~]=FINDMAX(squeeze(Eu_rel_ANGLES(f_rENG+10:size(Eu_rel_ANGLES,1),2,:)));
% % % % % [minFlex_aE,~]=FINDMIN(squeeze(Eu_rel_ANGLES(f_rENG+10:size(Eu_rel_ANGLES,1),2,:)));
% % % % % [maxAbd_aE,~]=FINDMAX(squeeze(Eu_rel_ANGLES(f_rENG+10:size(Eu_rel_ANGLES,1),3,:)));
% % % % % [minAbd_aE,~]=FINDMIN(squeeze(Eu_rel_ANGLES(f_rENG+10:size(Eu_rel_ANGLES,1),3,:)));
% % % % % 
% % % % % [maxTors_pE,~]=FINDMAX(squeeze(Eu_rel_ANGLES(find(t<(t_rPAS+0.0001)&(t>(t_rPAS-0.0001))):f_rENG-10,1,:)));
% % % % % [minTors_pE,~]=FINDMIN(squeeze(Eu_rel_ANGLES(find(t<(t_rPAS+0.0001)&(t>(t_rPAS-0.0001))):f_rENG-10,1,:)));
% % % % % [maxFlex_pE,~]=FINDMAX(squeeze(Eu_rel_ANGLES(find(t<(t_rPAS+0.0001)&(t>(t_rPAS-0.0001))):f_rENG-10,2,:)));
% % % % % [minFlex_pE,~]=FINDMIN(squeeze(Eu_rel_ANGLES(find(t<(t_rPAS+0.0001)&(t>(t_rPAS-0.0001))):f_rENG-10,2,:)));
% % % % % [maxAbd_pE,~]=FINDMAX(squeeze(Eu_rel_ANGLES(find(t<(t_rPAS+0.0001)&(t>(t_rPAS-0.0001))):f_rENG-10,3,:)));
% % % % % [minAbd_pE,~]=FINDMIN(squeeze(Eu_rel_ANGLES(find(t<(t_rPAS+0.0001)&(t>(t_rPAS-0.0001))):f_rENG-10,3,:)));

%% GRAPHICS
f1= figure('units','normalized','position',[0.02 0.1 0.96 0.8]); %#ok<NASGU>
for i=1:numel(units),
    subplot(3,ceil(numel(units)/3),i,'FontSize',8);hold on;
    plot(t,Accg(:,1:3,i),'-');
    plot(t,Accg(:,4,i),'k-','linewidth',2);
    axis([min_t max_t -8 +8]);
    if i==numel(units),
        legend(ACCnoG_header,'location','best');
    end
    xlabel('time [s]'); ylabel('a [g]'); title(['rENG: ' '-' 'Unit num:' units{i}]);
    hline(0,'k'); vline(0,'k');
    if ~isnan(f_rENG(i));
        vline(t(f_rENG(i)),'k--');
    end
    box on;
end
%%
%%%%% t_rENG NO g
f_rENG_NOg= [f_rENG_NOg min(f_rENG_NOg,[],2)];
t_rENG_NOg= nan(size(f_rENG_NOg));
t_rENG_NOg(~isnan(f_rENG_NOg))= t(f_rENG_NOg(~isnan(f_rENG_NOg)));

%%%%% t_rENG with g
f_rENG= [f_rENG min(f_rENG,[],2)];
t_rENG= nan(size(f_rENG));
t_rENG(~isnan(f_rENG))= t(f_rENG(~isnan(f_rENG)));

% Max X
f_XrENG= [f_XrENG min(f_XrENG,[],2)];
t_XrENG= nan(size(f_XrENG));
t_XrENG(~isnan(f_XrENG))= t(f_XrENG(~isnan(f_XrENG)));

% Max Y
f_YrENG= [f_YrENG min(f_YrENG,[],2)];
t_YrENG= nan(size(f_YrENG));
t_YrENG(~isnan(f_YrENG))= t(f_YrENG(~isnan(f_YrENG)));

% Max Z
f_ZrENG= [f_ZrENG min(f_ZrENG,[],2)];
t_ZrENG= nan(size(f_ZrENG));
t_ZrENG(~isnan(f_ZrENG))= t(f_ZrENG(~isnan(f_ZrENG)));

% Max Z XYZ
f_Z_XYZrENG= [f_Z_XYZrENG min(f_Z_XYZrENG,[],2)];
t_Z_XYZrENG= nan(size(f_Z_XYZrENG));
t_Z_XYZrENG(~isnan(f_Z_XYZrENG))= t(f_Z_XYZrENG(~isnan(f_Z_XYZrENG)));


%%
% % % % % [fC] = COMPARISON_KLA_KRA_KTR_XSX_v2(trial,path,T_off,T_off_raw,t_rENG,t_rENG_NOg,Eu_rel_ANGLES,t,t_rPAS);
%% OUTPUT
fX.rawdata= XSX;
fX.ACCnoG_header= ACCnoG_header;
fX.ACCG_header= ACCG_header;
fX.AXYZ_header= AXYZ_header;
fX.Aalpha_header= Aalpha_header;
fX.Eu_rel_ANGLES_header=Eu_rel_ANGLES_header;
fX.TO_header= TO_header;
fX.t_header= t_header;
fX.AccNOg= AccNOg;
fX.Accg= Accg;
fX.AXYZ= AXYZ;
fX.Aalpha= Aalpha;
fX.T_off= T_off;
fX.T_off_raw= T_off_raw;
fX.t= t;                              %time
fX.sf= sf;                            %sampling frequency
fX.t_rENG_NOg= t_rENG_NOg;            %t @actual engagement
fX.t_rENG= t_rENG;                    %t @actual engagement
fX.t_XrENG= t_XrENG;                  %t X@actual engagement
fX.t_YrENG= t_YrENG;                  %t Y@actual engagement
fX.t_ZrENG= t_ZrENG;                  %t Z@actual engagement
fX.t_Z_XYZrENG= t_Z_XYZrENG;          %t Z_XYZ@actual engagement
fX.f_rENG= f_rENG;                    % Max acc Module
fX.f_XrENG= f_XrENG;                  % Max X acc Module
fX.f_YrENG= f_YrENG;                  % Max Y acc Module
fX.f_ZrENG= f_ZrENG;                  % Max Z acc Module
fX.f_Z_XYZrENG= f_Z_XYZrENG;          % Max Z_XYZ acc Module
fX.Eu_rel_ANGLES=Eu_rel_ANGLES;       % Euler Angles caluculated in respect to trunk relative system
% % % % % % fX.Angles_gyr=fC.Trang_XSX;
% % % % % % fX.Angles_gyr_2b=fC.Trang_XSX_2b;
% % % % % % fX.fC=fC;

%% GRAPHICS
% % % % % %calibrated
% % % % % figure('units','normalized','position',[0.1 0.1 0.8 0.8]), title('calibrated')
% % % % % subplot(2,2,1),
% % % % % plot(t,Rx,'-','linewidth',2),
% % % % % ylim([-1000 5000]),
% % % % % hline(0,'k');vline(0,'k');vline(f_rENG/sf+t(1),'k--',{'1','2','3','4','tot'},[0 0.2]);
% % % % % xlabel('time [s]'), ylabel('lateral F [N]'),
% % % % % legend('beam1','beam2','beam3','beam4')
% % % % %
% % % % % subplot(2,2,2),
% % % % % plot(t,Rz,'-','linewidth',2),
% % % % % ylim([-1000 5000]),
% % % % % hline(0,'k');vline(0,'k');vline(f_rENG/sf+t(1),'k--',{'1','2','3','4','tot'},[0 0.2]);
% % % % % xlabel('time [s]'), ylabel('vertical F [N]'),
% % % % % legend('beam1','beam2','beam3','beam4')
% % % % %
% % % % % subplot(2,2,3),
% % % % % plot(t,Ry,'-','linewidth',2),
% % % % % ylim([-1000 5000]),
% % % % % hline(0,'k');vline(0,'k');vline(f_rENG/sf+t(1),'k--',{'1','2','3','4','tot'},[0 0.2]);
% % % % % xlabel('time [s]'), ylabel('compression F [N]'),
% % % % % legend('beam1','beam2','beam3','beam4')
% % % % %
% % % % % subplot(2,2,4),
% % % % % plot(t,Ry_lc,'-','linewidth',2),
% % % % % ylim([-1000 5000]),
% % % % % hline(0,'k');vline(0,'k');vline(f_rENG/sf+t(1),'k--',{'1','2','3','4','tot'},[0 0.2]);
% % % % % xlabel('time [s]'), ylabel('compression_LC F [N]'),
% % % % % legend('beam1','beam2','beam3','beam4')

%% FILE SAVING
filename= [file 'XSX_MS_(proc' date ')'];
an= questdlg(['Save into a ' filename '.mat|.xls file?'],'File Saving','Y','N','Y');
% % % % % an= 'Y'; % % % % %
if strcmpi(an,'Y');
    for i=1:numel(units),
        %save .mat file
        save([path filename '.mat'],'fX');
        
        %save .xls file
        %general info
        % % % % %         xlswrite([path filename '.xls'],{filename},'rENG','A1');
        % % % % %         xlswrite([path filename '.xls'],{'time @rENG [s]'},'rENG','A4');
        % % % % %         xlswrite([path filename '.xls'],{'B1' 'B2' 'B3' 'B4' 'tot'},'rENG','B3');
        % % % % %         xlswrite([path filename '.xls'],f_rENG/sf+t(1),'rENG','B4');
        
        %???
        if fileok(i),
            xlswrite([path filename '.xls'],{filename},units{i},'A1');
            xlswrite([path filename '.xls'],{'frame #'},units{i},'A3');
            xlswrite([path filename '.xls'],{'time [s]'},units{i},'B3');
            xlswrite([path filename '.xls'],ACCG_header,units{i},'C3');
            xlswrite([path filename '.xls'],TO_header,units{i},'H3');
% % % % %             xlswrite([path filename '.xls'],Eu_rel_ANGLES_header,units{i},'L3');
% % % % %             xlswrite([path filename '.xls'],AXYZ_header,units{i},'O3');
            
            xlswrite([path filename '.xls'],(1:length(t))',units{i},'A4');
            xlswrite([path filename '.xls'],t,units{i},'B4');
            xlswrite([path filename '.xls'],Accg(:,:,i),units{i},'C4');
            xlswrite([path filename '.xls'],T_off(:,:,i),units{i},'H4');
% % % % %             xlswrite([path filename '.xls'],Eu_rel_ANGLES(:,:,i),units{i},'L4');
% % % % %             xlswrite([path filename '.xls'],AXYZ(:,:,i),units{i},'O4');
        end
    end
end
