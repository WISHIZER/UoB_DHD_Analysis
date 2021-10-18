% Calculate Residual Stress from DHD and iDHD test 
clc;
clear all;
close all;


%% INPUTS

% The directory of the excel files for each increment which includs the diameter before and atter EDM 

filenames{1,1} = 'C:\Users\mm17209\OneDrive - University of Bristol\Projects\iDHD_Altlas+\Atlas+\iDHD evaluation\New analysis\EDM1';
filenames{2,1} = 'C:\Users\mm17209\OneDrive - University of Bristol\Projects\iDHD_Altlas+\Atlas+\iDHD evaluation\New analysis\EDM2';
filenames{3,1} = 'C:\Users\mm17209\OneDrive - University of Bristol\Projects\iDHD_Altlas+\Atlas+\iDHD evaluation\New analysis\EDM3';
filenames{4,1} = 'C:\Users\mm17209\OneDrive - University of Bristol\Projects\iDHD_Altlas+\Atlas+\iDHD evaluation\New analysis\EDM4';
filenames{5,1} = 'C:\Users\mm17209\OneDrive - University of Bristol\Projects\iDHD_Altlas+\Atlas+\iDHD evaluation\New analysis\EDM5';
filenames{6,1} = 'C:\Users\mm17209\OneDrive - University of Bristol\Projects\iDHD_Altlas+\Atlas+\iDHD evaluation\New analysis\EDM6';
filenames{7,1} = 'C:\Users\mm17209\OneDrive - University of Bristol\Projects\iDHD_Altlas+\Atlas+\iDHD evaluation\New analysis\EDM7';
filenames{8,1} = 'C:\Users\mm17209\OneDrive - University of Bristol\Projects\iDHD_Altlas+\Atlas+\iDHD evaluation\New analysis\EDM8';
filenames{9,1} = 'C:\Users\mm17209\OneDrive - University of Bristol\Projects\iDHD_Altlas+\Atlas+\iDHD evaluation\New analysis\EDM9';
filenames{10,1} ='C:\Users\mm17209\OneDrive - University of Bristol\Projects\iDHD_Altlas+\Atlas+\iDHD evaluation\New analysis\EDM10';
filenames{11,1} ='C:\Users\mm17209\OneDrive - University of Bristol\Projects\iDHD_Altlas+\Atlas+\iDHD evaluation\New analysis\EDM11';

n_incrs = size(filenames,1);% total number of increments
% select the sheet and interested data range
for ii = 1:n_incrs;
sheets{ii,1} = 'diameter';
ranges{ii,1} = 'B2:Q471';
end



% through edm (last step) - import before and after edm diameters
file_thru = 'C:\Users\mm17209\OneDrive - University of Bristol\Projects\iDHD_Altlas+\Atlas+\iDHD evaluation\New analysis\EDM11';
sheet_thru = 'Diameter';
range_thru = 'B2:Q471';


an = 8;% number of angles
E=195600;%Young's Modulus

%% IMPORT DATA

% increments
for ii = 1:n_incrs
    dia_incrs{ii,1} = xlsread(filenames{ii,1}, sheets{ii,1}, ranges{ii,1});
end

% through EDM

dia_thru = xlsread(file_thru, sheet_thru, range_thru);
dia_depth=flip(xlsread(file_thru, sheet_thru, 'A2:A471'));  %flip depth data to have the reading from front bush to back bush.                                                            
%ind_LOI=find(d_depth1>15.8 & d_depth1<32.1);% find the length of Interest (LOI)
n_depths = size(dia_thru,1); % number of measurements through the length of measurments
EDM_incrs=1:1:n_incrs;

%% Smoothing imported data, all variable with _sm suffix refers to smoothed data 

dia_incrs_sm=cell(n_incrs,1);
for m=1 :n_incrs 
    for n=1:16
        dia_incrs_sm{m,1}(:,n)=smooth(dia_incrs{m,1}(:,n),11,'sgolay',1); %Savitzky-Golay Filter
    end   
end
dia_thru_sm=dia_incrs{n_incrs(end),1};
dia_incrs=dia_incrs_sm;% replace the data with smoothed data
dia_thru=dia_thru_sm;% replace the data with smoothed data

%% Strain Measurement

e_ang=cell(n_incrs,1);% measured strain for each angle
d_ang=cell(n_incrs,1);% diametral distorsion
for i=1:n_incrs
    e_ang{i,1}(:,1)=dia_depth;
    d_ang{i,1}(:,1)=dia_depth;
end
for i=1:n_incrs
    for j=1:an
        e_ang{i,1}(:,j+1)=((dia_incrs{i,1}(:,2*(j))-dia_incrs{i,1}(:,2*(j)-1)))./(dia_incrs{i,1}(:,2*(j)-1));% strain: (dia_afterEDM - dia_BeforeEDM))/dia_BeforeEDM
        d_ang{i,1}(:,j+1)=(dia_incrs{i,1}(:,2*(j))-dia_incrs{i,1}(:,2*(j)-1)); %diametral distorsion: (dia_afterEDM-dia_BeforeEDM)
        
    end
%     hold on
%     plot(e_ang{1,1}(:,1),e_ang{i,1}(:,2));
%     ylim([-0.005,0.005])
%     legend({'1','2','3','4','5','6','7','8','9','10','11','12','13','14'})
%     pause (0.5)
end
%find strain offset
for k=1:n_incrs
    for kk=1:an
        e_ang_avg(k,kk)=mean(e_ang{k,1}(37:122,kk+1)); % Uses the average strain in back bush (assuming it does not have residual stress)to find the offset strain
    end
    
end
% Recalculate the diameter after EDM according to offset strain
for k=1:n_incrs
    for kk=1:an
        for kkk=1:n_depths
        dia_incrs_offset{k,1}(kkk,kk*2)=dia_incrs{k,1}(kkk,kk*2)-(e_ang_avg(k,kk)*dia_incrs{k,1}(kkk,kk*2));
        dia_incrs_offset{k,1}(kkk,kk*2-1)=dia_incrs{k,1}(kkk,kk*2-1);
        end
    end
    
end
%calculate the offset strain and offset diametric distortion based on measured strain offset:
e_ang_offset=cell(n_incrs,1);
d_ang=cell(n_incrs,1);

for i=1:n_incrs
    e_ang_offset{i,1}(:,1)=dia_depth;
    d_ang_offset{i,1}(:,1)=dia_depth;
end
for i=1:n_incrs
    for j=1:an
       e_ang_offset{i,1}(:,j+1)=((dia_incrs_offset{i,1}(:,2*j)-dia_incrs_offset{i,1}(:,2*j-1)))./(dia_incrs_offset{i,1}(:,2*j-1));
       d_ang_offset{i,1}(:,j+1)=(dia_incrs_offset{i,1}(:,2*j)-dia_incrs_offset{i,1}(:,2*j-1));
        
    end
     
end
% 
% % ploting the diametral distortion
% figure
% for i=1:n_incrs
%     hold on
%     plot(dia_depth,d_ang_offset{i,1}(:,2),'linewidth',rand(1)*2)
%     ylim([-0.005,0.005])
%     legend({'1','2','3','4','5','6','7','8','9','10','11','12','13','14'})
%     
% end
% % plotting the diameter after EDM
% 
% 
% for i=1:n_incrs
%     hold on
%     plot(dia_depth,dia_incrs{i,1}(:,2),'linewidth',0.5)
%     ylim([1.495,1.51])
%     legend({'1','2','3','4','5','6','7','8','9','10','11','12','13','14'})
%     
% end

%% Stress Calculation
load M % M is the conversion matrix
PIM=(inv((M.')*M)*(M.')).';
    %Figuers 1 to 3 plotting the Hoop, Axial and Shear stresses for each
    %increament. Being Axial or Hoop depends on the the orientation of
    %plate or pipe with respect to the inial measurement angle. in this
    %case the hoop stress were measured at 0 angle. 
    figure(1) 
    movegui(figure(1),'northwest')
    c=colormap(flip(winter(11)));
    title ('Hoop Stress')
    xlim([0,42.5])
    ylim([-400,400])
    xlabel('Depth, mm','FontSize',12,'FontWeight','bold')
    ylabel('Hoop Stress, MPa', 'FontSize',12,'FontWeight','bold')
    grid on;grid minor
    
    figure(2)
    movegui(figure(2),'north')
    cc=colormap(flip(cool(11)));
    title ('Axial Stress')
    xlim([0,42.5])
    ylim([-400,400])
    xlabel('Depth, mm','FontSize',12,'FontWeight','bold')
    ylabel('Axial Stress, MPa', 'FontSize',12,'FontWeight','bold')
    grid on;grid minor
  
    figure(3)
    title ('Shear Stress')
    movegui(figure(3),'northeast')
    ccc=colormap(flip(parula(11)));
    xlim([0,42.5])
    ylim([-400,400])
    xlabel('Depth, mm','FontSize',12,'FontWeight','bold')
    ylabel('Shear Stress, MPa', 'FontSize',12,'FontWeight','bold')
    grid on;grid minor
    
for i=1:n_incrs
    
        S{i,1}=-E*((e_ang_offset{i,1}(:,2:end)*PIM));
               
            figure(1)
            hold on
            plot(e_ang_offset{1,1}(:,1),S{i,1}(:,1),'Color',c(i,:),'linewidth',i*0.1+0.2);
            legend({'1','2','3','4','5','6','7','8','9','10','11','12','13'})
            figure(2)
            hold on
            plot(e_ang_offset{1,1}(:,1),S{i,1}(:,2),'Color',cc(i,:),'linewidth',i*0.1+0.2);
            legend({'1','2','3','4','5','6','7','8','9','10','11','12','13'})
            figure(3)
            hold on
            plot(e_ang_offset{1,1}(:,1),S{i,1}(:,3),'Color',ccc(i,:),'linewidth',i*0.1+0.2);
            legend({'1','2','3','4','5','6','7','8','9','10','11','12','13'})

end

%% FIND MAX DIAMETRAL DISTORTION (Method developed by Sam Oliver) 
%this calculates the maximum difference in diametral distortion from all of the incremental measurements
dmax_ind_cell = zeros(n_depths,an); %the index of the maximum diametral distortion
dia_iDHD = zeros(n_depths,2*an);    %output array containing diameters before and after edm which correspond to the maximum diametral displacement
for ii = 1:an;  %over the angles
    ind_col_after = 2*ii;    %column index after edm for this angle
    ind_col_before = 2*ii - 1;   %column index before edm for this angle
    for jj = 1:n_depths    %over the measurement depths
        
        %         Find diametric distortions for each measurement depth for
        %         each angle.
        %         Start with reference measurement (taken as before column
        %         of first increment file).
        d_jj(1) = dia_incrs_offset{1,1}(jj,ind_col_before) - dia_incrs_offset{n_incrs,1}(jj,ind_col_after);
%         
        for kk = 1:n_incrs %d_jj = diametral distortion over the measurement depth
            d_jj(kk+1) = dia_incrs_offset{kk,1}(jj,ind_col_after) - dia_incrs_offset{n_incrs,1}(jj,ind_col_after);
        end
        
        %    location in d_jj of max diametral distortion for this measurement depth and angle
        d_jj_abs = abs(d_jj);
        dmax_ind = find(d_jj_abs==max(d_jj_abs));
        dmax_ind_cell(jj,ii) = dmax_ind;
        %         Replace before measurements with those (either from
        %         before column of reference or from after column of any 
        %         other increment) that give the biggest diametral distortion.
        if dmax_ind >=2
        dia_iDHD(jj,ind_col_before) = dia_incrs_offset{dmax_ind-1,1}(jj,ind_col_after);
        else   
            dia_iDHD(jj,ind_col_before) = dia_incrs_offset{1,1}(jj,ind_col_before);
        end
    end
    % replace after measurements in dia_iDHD
    dia_iDHD(:,ind_col_after) = dia_incrs_offset{n_incrs,1}(:,ind_col_after);

end
%% Strain Measurement for iDHD

    e_ang_iDHD(:,1)=dia_depth;
    for j=1:an
        e_ang_iDHD(:,j+1)=((dia_iDHD(:,2*(j))-dia_iDHD(:,2*(j)-1)))./(dia_iDHD(:,2*(j)-1));
    end
     

%% Stress Calculation iDHD
% plots hoop, axial and shear RS from DHD (red line) and IDHD (blue line) 
load M
PIM=(inv((M.')*M)*(M.')).';
S_iDHD=-E*((e_ang_iDHD(:,2:end)*PIM));

    
    figure(4)
    movegui(figure(4),'southwest')
    plot(dia_depth,S{n_incrs,1}(:,1),'r',dia_depth,S_iDHD(:,1),'b');
    title ('Hoope Stress, iDHD vs DHD')
    xlim([0,42.5])
    ylim([-400,400])
    xlabel('Depth, mm','FontSize',12,'FontWeight','bold')
    ylabel('Hoop Stress, MPa', 'FontSize',12,'FontWeight','bold')
    legend({'DHD','IDHD'})
    grid on;grid minor
    
    figure(5)
    movegui(figure(5),'south')
    plot(dia_depth,S{n_incrs,1}(:,2),'r',dia_depth,S_iDHD(:,2),'b');
    legend({'DHD','IDHD'})
    title ('Axial Stress, iDHD vs DHD')
    xlim([0,42.5])
    ylim([-400,400])
    xlabel('Depth, mm','FontSize',12,'FontWeight','bold')
    ylabel('Axial Stress, MPa', 'FontSize',12,'FontWeight','bold')
    grid on;grid minor
    
    figure(6)
    movegui(figure(6),'southeast')
    plot(dia_depth,S{n_incrs,1}(:,3),'r',dia_depth,S_iDHD(:,3),'b');
    legend({'DHD','IDHD'})
    title ('Shear Stress, iDHD vs DHD')
    xlim([0,42.5])
    ylim([-400,400])
    xlabel('Depth, mm','FontSize',12,'FontWeight','bold')
    ylabel('Shear Stress, MPa', 'FontSize',12,'FontWeight','bold')
    grid on;grid minor
    
   
   


