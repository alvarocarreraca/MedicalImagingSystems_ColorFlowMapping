%% SETUP of script for in vivo data:
close all; clear all; clc;

load('cfm_carotis.mat');

%% VARIABLES:
% FIXED PARAMETERS:
fs = 40*10^6; % 40 MHz (sampling frequency)
Resolution = 16; % 16 bits samples
c = 1540; % 1540 m/s
f_prf = 6*10^3; % 6 kHz
T_prf = 1/f_prf;
f0 = 5*10^6; % 5 MHz (center frequency)
focus_depth = 18*10^(-3); % 18 mm
cycles = 8; % 8 cycles in one pulse
start_depth = 1; % 1 mm
end_depth = 30; % 30 mm
% Dimensions left to right:
x_min = -9.75; % -9.75 mm
x_max = 9.75; % 9.75 mm

% VARIABLE PARAMETERS:
segmentSize = 10; % 10 
numPointsCorr = 38; % 10
columnsOverlaped = 9; % It can be: 1, 3, 7, 9 ,21 and 63
velocityRange = 1; % +-1 [m/s] of velocity range in the Carotid

%% CLACULATIONS:

% CALCULATE VELOCITIES:
velocity_matrix = [];
for j = 1:size(vessel,2)
    data = double(rf_cfm_data(:,:,j)).*vessel(:,j);
    velocity_matrix_j = mainFunction(data,fs,f0,cycles,c,T_prf,segmentSize,numPointsCorr,velocityRange,columnsOverlaped);
    velocity_matrix = [velocity_matrix velocity_matrix_j];
end

% CALCULATE B-MODE:
bmode_data(:,end) = [];
Hilb_bmode=hilbert(bmode_data); % hilbert matrix
env_bmode = abs(Hilb_bmode); % signal envelope through hilbert transformation
env_db_bmode=20*log10(env_bmode/max(max(env_bmode)));
E_dB_60_bmode=env_db_bmode+60; %previous values added 60
% We delete the values under the bottom 60 dB range.
E_dB_60_bmode(E_dB_60_bmode<0) = 0; %all the negative values become 0

%% PREPARE RESULTS FOR PLOTING:

% Interpolate Y-axis from velocity matrix:
interpFactor = floor(1520/size(velocity_matrix,1));
newVelocityMatrix = zeros(size(velocity_matrix,1)*interpFactor,size(velocity_matrix,2));
for i = 1:size(velocity_matrix,2)
    newVelocityMatrix(:,i) = interp(velocity_matrix(:,i),interpFactor);
end

% Interpolate X-axis from velocity matrix (for 9 columns overlaped the factor is: 4/7 -> 112(columns)*(4/7) = 64:
interpFactorX = 4;
dicimateX = 7; %7
velocity_matrix_interp = zeros(size(newVelocityMatrix,1),size(newVelocityMatrix,2)*interpFactorX/dicimateX);
for i = 1:size(newVelocityMatrix,1)
    aux = interp(newVelocityMatrix(i,:),interpFactorX);
    velocity_matrix_interp(i,:) = decimate(aux,dicimateX);
end

if(size(velocity_matrix_interp,1) ~= size(E_dB_60_bmode,1))
    diff_lines = size(E_dB_60_bmode,1) - size(velocity_matrix_interp,1); 
    array = zeros(diff_lines/2,size(velocity_matrix_interp,2));
    velocity_matrix_interp = [array; velocity_matrix_interp; array];
end

%% PREPARE MASK FOR OVERLATING B-MODE and VELOCITIES IN ONE PLOT:

% Interpolate original Mask (vessel mask has 16 columns and we need 64 -> factor of 4):
interpFactorMask = 4;
mask_interpolated = zeros(size(vessel,1),size(vessel,2)*interpFactorMask);
for i = 1:size(vessel,1)
    mask_interpolated(i,:) = round(interp(vessel(i,:),interpFactorMask));
end
mask_interpolated(mask_interpolated > 1) = 1;
mask_interpolated(mask_interpolated < 0) = 0;

%% MERGE B-MODE AND VELOCITIES:
% NORMALIZE VALUES FOR B-MODE:
% NewValue = (((OldValue - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
E_dB_60_bmode_color=(((E_dB_60_bmode - 0)*(0.5-(-0.5)))/(60-0))+(-0.5);

% MERGE IMAGES:
merged_matrix = zeros(size(mask_interpolated,1),size(mask_interpolated,2));
for i = 1:size(merged_matrix,1)
    for j = 1: size(merged_matrix,2)
        if (mask_interpolated(i,j) == 0)
            merged_matrix(i,j) = E_dB_60_bmode_color(i,j);
        else
            merged_matrix(i,j) = velocity_matrix_interp(i,j);
        end
    end
end

%% PLOTS:

% CREATE COLORMAP:
jetLine = jet(256);
hotLine = hot(180);
colorMap = zeros(256,3);
colorMap(1:128,:) = flip(jetLine(1:128,:));
colorMap(129:end,:) = hotLine(1:128,:);

figure;
imagesc([x_min x_max],[start_depth end_depth],E_dB_60_bmode_color,'AlphaData',(-1)*(mask_interpolated-1),[-0.5 0.5]);
colorbar;
ylim([start_depth 28.84])
colormap('gray');
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]');
title('B-mode');

figure;
imagesc([x_min x_max],[start_depth end_depth],velocity_matrix_interp,'AlphaData',mask_interpolated,[-0.5 0.5]);
colorbar;
ylim([start_depth 28.84])
colormap(colorMap);
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]');
title('Velocity map of the carotid artery')



figure;
ax1 = axes;
A = imagesc([x_min x_max],[start_depth end_depth],E_dB_60_bmode_color,[-0.5 0.5]);
ylim([start_depth 28.84])
colorbar
ax2 = axes;
B = imagesc([x_min x_max],[start_depth end_depth],velocity_matrix_interp,'AlphaData',mask_interpolated,[-0.5 0.5]);
ylim([start_depth 28.84])
linkaxes([ax1,ax2]);
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
colormap(ax1,'gray');
colormap(ax2,colorMap);
set(ax2,'color','none','visible','off');
colorbar
title('Color flow map from carotid artery')
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]');

% Velocity profile:
figure;
plot(linspace(start_depth,end_depth,length(velocity_matrix_interp(:,50))),velocity_matrix_interp(:,50))
xlim([start_depth end_depth])
title('Velocity profile at a lateral distance of 5.14 [mm]')
xlabel('Axial distance [mm]');
ylabel('Velocity [m/s]');
