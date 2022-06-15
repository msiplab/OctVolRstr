
% Setup script
%
% Requirements: MATLAB R2020b
%
% Copyright (c) 2021-2022, Ruiki KOBAYASHI and Shogo MURAMATSU
%
% All rights reserved.
%
% Contact address: Shogo MURAMATSU,
%    Faculty of Engineering, Niigata University,
%    8050 2-no-cho Ikarashi, Nishi-ku,
%    Niigata, 950-2181, JAPAN
%
% http://msiplab.eng.niigata-u.ac.jp/

%% Create RESULTS folder
if ~exist("../results","dir")
    mkdir("../results")
else
    disp("../results exists.")
end

%% Download SaivDr package
SAIVDR_VERSION = "4.2.2.1";
SAIVDR_ROOT = "SaivDr-"+SAIVDR_VERSION;
if ~exist(SAIVDR_ROOT,"dir")
    unzip("https://github.com/msiplab/SaivDr/archive/refs/tags/"+ ...
        SAIVDR_VERSION+".zip")
end
addpath(SAIVDR_ROOT)
CURRENT_DIR = pwd;
cd(SAIVDR_ROOT)
setpath
cd(CURRENT_DIR)

%% Expand materials
%
% Please download the following ZIP file and put it in this directory
% before execution.
%
% (Limited to IEEE members)
% https://drive.google.com/file/d/13ZrAmw587vPUopEjsF6Se9sDbzMuwF3q/view?usp=sharing
%
MATERIALS="materials.zip";
CURRENT_DIR = pwd;
DATA_DIR = "../data";
if ~exist("../data/materials","dir")
    cd(DATA_DIR)
    unzip(MATERIALS,".")
    cd(CURRENT_DIR)
else
    disp("../data/materials exists.")
end
