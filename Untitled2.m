clear all
close all
clc
%%

load('cfm_carotis.mat')
whos

% H=hilbert(bmode_data); % hilbert matrix
% env = abs(H); % signal envelope through hilbert transformation
% env_db=20*log10(env/max(max(env))); % compress the envelope (range from -inf to 0)
% 
% figure
% imagesc(env_db) % Display image with scaled colors
% set(gca, 'YDir', 'Normal') % GCA: Get handle to current axis to obtain the fig (without it, the image is reverse)

%% 2) Make the scatterer map with a size of 40x40 mm with cyst hole

figure
imagesc(vessel)
colormap(gray(128))
colorbar

%% 3)Make 2D convolution of psf1 and scatterer map. Find compressed envelope data and display the image with dinamic range of 60dB.


% compressed envelope data
H2=hilbert(bmode_data); % hilbert matrix
env2 = abs(H2); % signal envelope through hilbert transformation
env_db2=20*log10(env2/max(max(env2)));

E_dB_60=env_db2+60; %previous values added 60

% We delete the values under the bottom 60 dB range.
E_dB_60(E_dB_60<0) = 0; %all the negative values become 0

% We scale it to a 0-127 scale in order to see the correct grayscale image
% NewValue = (((OldValue - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
E_dB_60_color=(((E_dB_60 - 0)*(127-0))/(60-0))+0;
E_dB_60_color_int=uint32(E_dB_60_color); % The envelope detected and log-compressed data 
                                        % as an integer array as 32 bits values (uint32) 
                                        % uint32 converts elements of the array into elementi interi positivi

figure
imagesc(E_dB_60_color) % Display image with scaled colors
ylim([0 1400])
colormap(gray(128))
colorbar

newImg = imfuse(vessel,E_dB_60_color,'blend','Scaling','joint');
figure
imagesc(newImg)