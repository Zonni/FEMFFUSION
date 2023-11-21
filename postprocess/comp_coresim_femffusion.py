# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 10:42:44 2019

@author: zonni
"""

import scipy.io

problem = '2D_biblis_noise' 


if (problem == '2D_biblis_noise'):
    folder = '../2D_biblis_noise/output/'
    coresim_file = 'RESULTS.mat'
    femffus_file = 'RESULTS_FEMFFUSION.mat'
    
    
mat = scipy.io.loadmat('file.mat')

#%% LOAD
#% CORESIM Results for 2D pertubation 
#load('output/RESULTS.mat');
#
#
#% FEMFFUSION results 
#load('output/RESULTS_FEMFFUSION.mat'); 
#
#%% POSTPROCESS
#mesh_size = size(dFLX1);
#nx =mesh_size(1);
#ny =mesh_size(2);
#
#LX = 23.1226 * 17;
#LY = 23.1226 * 17;
#
#DX = LX / nx;
#DY = LY / ny;
#
#x = DX/2:DX:nx*DX; 
#y = DY/2:DY:ny*DY; 
#
#% Resize variables 
#% 
#[norm, I] = max(FLX1(:));
#[I1,I2,I3,I4] = ind2sub(size(FLX1),I); 
#FLX1 = FLX1(:, :, I3, 1) / norm; 
#FLX2 = FLX2(:, :, I3, 1) / norm;
#dFLX1 = dFLX1(:, :, I3, 1) * 2 * pi; 
#dFLX2 = dFLX2(:, :, I3, 1) * 2 * pi;
#
#
#
#for i=1:nx
#    for j=1:ny
#        if (materials(i, j) == 0)
#            FLX1(i, j) = NaN;
#            FLX2(i, j) = NaN;
#            dFLX1(i, j) = NaN;
#            dFLX2(i, j) = NaN;
#        end
#    end
#end
#
#%% GET FEMFUSSION LINE
#
#tol = 1e-8;
#y_line =  150.297; % cm
#line_indx = find(y_line == y_py);
#x_line = x_py(line_indx);
#norm = max(static_fs);
#static_fs_line = static_fs(line_indx) / norm;
#static_th_line = static_th(line_indx) / norm;
#dfs_line  = py_fs(line_indx) * 2 * pi ;
#dth_line = py_th(line_indx) *  2 * pi ;
# 
#%% STATIC SOLUTIONS
#
#
#line = (138.7356 + 161.8582) / 2;
#central_line = round(line / DY); 
#
#figure
#hold on
#grid on
#title('Static Solution')
#plot(x, abs(FLX1(:, central_line)), '-', 'LineWidth', 3)
#plot(x, abs(FLX2(:, central_line)), '-', 'LineWidth', 3)
#
#plot(x_line, static_fs_line, '-', 'LineWidth', 3)
#plot(x_line, static_th_line, '-', 'LineWidth', 3)
#
#legend("Fast Flux CORE SIM", "Thermal Flux CORE SIM", "Fast Flux FEMFFUSION", "Thermal Flux FEMFFUSION")
#%legend("Fast Flux CORE SIM", "Thermal Flux CORE SIM")
#zlabel("Neutron Flux")
#xlabel("x (cm)")
#ylabel("y (cm)")
#
#
#
#
#%% NOISE AMPLITUDE PLOTS
#% 
#figure
#hold on
#grid on
#plot(x, abs(dFLX1(:, central_line)), '-', 'LineWidth', 3)
#plot(x_line, abs(dfs_line), '--', 'LineWidth', 3)
#ylabel("Fast Flux Noise Amplitude")
#xlabel("x (cm)")
#legend("CORE SIM", "FEMFFUSION")
#saveas(gcf,'comparative_2Dbiblis_fs','epsc')
#
#figure
#hold on
#grid on
#plot(x, abs(squeeze(dFLX2(:, central_line))), '-', 'LineWidth', 3)
#plot(x_line, abs(dth_line), '--', 'LineWidth', 3)
#ylabel("Thermal Flux Noise Amplitude")
#xlabel("x (cm)")
#legend("CORE SIM", "FEMFFUSION")
#saveas(gcf,'comparative_2Dbiblis_th','epsc')
#
#
#figure
#hold on
#title('Fast Flux Noise Amplitude')
#xlabel(' x (cm)');
#ylabel(' y (xm)');
#plot_matrix(x, y, abs(dFLX1)')
#colorbar();
#saveas(gcf,'amplitude_fs_2Dbiblis','epsc')
#
#figure
#hold on
#title('Thermal Flux Noise Amplitude')
#xlabel(' x (cm)');
#ylabel(' y (xm)');
#plot_matrix(x, y, abs(dFLX2)')
#colorbar();
#saveas(gcf,'amplitude_th_2Dbiblis','epsc')
#
#%% PHASE PLOT
#
#figure
#hold on
#grid on
#
#plot(x, angle(dFLX1(:, central_line)) *180/pi , '-', 'LineWidth', 3)
#plot(x_line, angle(dfs_line) *180/pi, '--', 'LineWidth', 3)
#legend("CORE SIM", "FEMFFUSION",  'Location', 'best')
#ylabel("Thermal Flux Noise Phase (deg)")
#xlabel("x (cm)")
#saveas(gcf,'comparative_2Dbiblis_phase_fs','epsc')
#
#figure
#hold on
#grid on
#mean1 = mean(angle(dFLX2(:, central_line))*180/pi);
#plot(x, angle(dFLX2(:, central_line))*180/pi, '-', 'LineWidth', 3)
#plot(x_line, angle(dth_line) *180/pi, '--', 'LineWidth', 3)
#legend("CORE SIM", "FEMFFUSION", 'Location', 'best')
#ylabel("Thermal Flux Noise Phase (deg)")
#xlabel("x (cm)")
#saveas(gcf,'comparative_2Dbiblis_phase_th','epsc')
