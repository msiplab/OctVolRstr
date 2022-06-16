% MAIN_ALL
%
% Run all of the live scripts
setup

%% Modeling process:
% Evaluate the effect of broadening parameter to interference function
% (Fig. 3)
main_sim_1d_interference

%% Simulation process:
% Conduct 1-D simulation by sweeping parameters η and λ      
main_sim_1d_paramswp_proc 

% Analyze the results obtained by the script MAIN_SIM_1D_PARAMSWP_PROC
% (Tables I & II)
main_sim_1d_paramswp_eval

% Evaluate the effect of introducing refractive index distribution 
% (Figs. 2, 5 & 6)
main_sim_1d_paramswp_graph 

% Sweep the optimal values for η and λ by 3-D simulation
main_sim_3d_paramswp_proc

% Draw graphs showing the optimal values of  η and λ from the results 
% by the script MAIN_SIM_3D_PARAMSWP_PROC  (Fig. 10)
main_sim_3d_paramswp_graph 

% Conduct 3-D simulation of sampling adjustment, interference waveform 
% estimation and volumetric data restoration
main_sim_3d_all_proc
  
% Draw graphs from the results obtained by the script MAIN_SIM_3D_ALL_PROC 
% (Figs. 8, 9, 11 & Tab. III)
main_sim_3d_all_graph
    
% Summerize the results obtained by the script MAIN_SIM_3D_ALL_PROC
% (Table III)
main_sim_3d_all_table 

% Sampling adjustment is performed with real data (Fig. 12)
main_exp_3d_sampling_adjust

%% Experimental process:
% Interference waveform estimation with real data
main_exp_3d_waveform_est 

% Draw graphs from the results obtaind 
% by the script MAIN_EXP_3D_WAVEFORM_EST (Fig. 13)
main_exp_3d_waveform_graph 

% Volumetric data restoration with real data
main_exp_3d_rest_prop

% Draw graphs from the results obtained by the script MAIN_EXP_3D_REST_PROP
% (Figs. 14, 15 & 16)
main_exp_3d_rest_graph
      