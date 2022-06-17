# OctVolRstr Project

Supplemental materials are included for the following paper:

- Ruiki Kobayashi, Genki Fujii, Yuta Yoshida, Takeru Ota, Fumiaki Nin, 
  Hiroshi Hibino, Samuel Choi, Shunsuke Ono, Shogo Muramatsu,
  "Sparsity-Aware OCT Volumetric Data Restoration Using Optical Synthesis 
  Model," IEEE Transactions on Computational Imaging, 
  DOI: 10.1109/TCI.2022.3183396, to appear

The latest version is maintained at the following site:

	https://github.com/msiplab/OctVolRstr

## Description:

   This repository contains MATLAB script, function and class definition
   codes so that readers can reproduce the following materials in the
   paper:

    * Section II.A: 
    
      - Fig. 3: Effect of the broadening parameter bp for pz[mz]

    * Section IV.A: 

      - Fig. 5: Examples of 1-D simulation results 
      - Fig. 6: Comparison of computational complexities
      - Table I: 1-D simulation specifications
      - Table II: Simulation results of 1-D reflectance restoration 

    * Section IV.B: Restoration Simulation

      - Fig. 8: Simulation results of the sampling adjustment with STFT
      - Fig. 9: Simulation results of the interference waveform estimation
      - Fig. 10: MSE validation for sparse regularization parameters λ 
        and η
      - Fig. 11: Example set of artificial volumetric arrays
      - Table III: Simulation results for various noise levels

    * Section V: Restoration Experiment

      - Fig. 12: Sampling adjustment result
      - Fig. 13: Waveform estimation results
      - Fig. 14: Observation of the sensory epithelium of an inner ear 
        of a mouse
      - Fig. 15: Restored result
 
   Note that these materials utilize a subset of SaivDr Package for
   MATLAB, which is developed by the authors and the full latest
   version is available from the following site:

    https://github.com/msiplab/SaivDr

## Platform:

   MATLAB R2021b/R2022a

## Environment:

   The attached MATLAB codes have been tested on the following
   operating systems.

    * Windows 10 (64bits) 
    * Ubuntu 20.04 LTS (64bits)

   In order to execute the attached codes, the following MATLAB
   Toolboxes are required.
   
    * Signal Processing Toolbox
    * Image Processing Toolbox
    * Wavelet Toolbox
    * Curve Fitting Toolbox

   In addition, the codes for the design processes require the
   following options.

    * Optimization Toolbox
    
   It is also recommended to prepare the following Toolboxes for
   accelerating the execution.

    * MATLAB Coder
    * Parallel Computing Toolbox

   Some of the MATLAB codes try to download images, volumetric data and
   additional program codes from the Internet. Unless necessary
   materials are not downloaded, those codes must be executed with the
   Internet connection. 

## Major Component Description:

   The MATLAB scripts of which name begin with 'main_' in the top
   layer reproduce the dictionary learning results, figures and
   tables used in the paper.

   Modeling process:
    * main_sim_1d_interference - Evaluate the effect of broadening parameter
      to interference function(Fig. 3)

   Simulation process:

    * main_sim_1d_paramswp_proc - Conduct 1-D simulation by 
      sweeping parameters η and λ      
    * main_sim_1d_paramswp_eval - Analyze the results obtained by 
      the script MAIN_SIM_1D_PARAMSWP_PROC (Tables I & II)
    * main_sim_1d_paramswp_graph - Evaluate the effect of introducing 
      refractive index distribution (Figs. 2, 5 & 6)
    * main_sim_3d_paramswp_proc - Sweep the optimal values for η and λ 
      by 3-D simulation
    * main_sim_3d_paramswp_graph - Draw graphs showing the optimal values
      of  η and λ from the results by the script MAIN_SIM_3D_PARAMSWP_PROC
      (Fig. 10)
    * main_sim_3d_all_proc - Conduct 3-D simulation of sampling adjustment, 
      interference waveform estimation and volumetric data restoration
    * main_sim_3d_all_graph - Draw graphs from the 
      results obtained by the script MAIN_SIM_3D_ALL_PROC
      (Figs. 8, 9, 11 & Tab. III)
    * main_sim_3d_all_table - Summerize the results obtained by the script 
      MAIN_SIM_3D_ALL_PROC (Table III)

   Experimental process:

    * main_exp_3d_sampling_adjust - Sampling adjustment is performed with 
      real data (Fig. 12)
    * main_exp_3d_waveform_est - Interference waveform estimation with 
      real data
    * main_exp_3d_waveform_graph - Draw graphs from the results obtaind
      by the script MAIN_EXP_3D_WAVEFORM_EST (Fig. 13)
    * main_exp_3d_rest_prop - Volumetric data restoration with real data
    * main_exp_3d_rest_graph - Draw graphs from the results obtained by 
      the script MAIN_EXP_3D_REST_PROP (Figs. 14, 15 & 16)

   Folder '+support' contain the functions called by
   the above scripts during their execution.

## Instructions for setup:
    Download the following ZIP file and put it under the directory 'data'.

    /data/materials.zip at https://codeocean.com/capsule/9148728/
    
    Then, move to directiory 'code' and run setup on MATLAB

    >> cd code
    >> setup

## Instructions for experimental scripts:

   Just input a main script name which begins by 'main_' on the MATLAB
   command window. For example,
 
    >> main_sim_1d_paramswp_graph

## Output files:
    The results is saved in folder "/results."
    For the detail, please see the header comments on each script

## Contact Information:
    Ruiki Kobayashi, Master Cource Student
    Graduate School of Science and Technology, Niigata University,
    Email: f20c050g@ieee.org

    Shogo MURAMATSU, Professor
    Faculty of Engineering, Niigata University,
    Email: shogo@eng.niigata-u.ac.jp
           shogo.muramatsu.jp@ieee.org
