%% MAIN_SIM_1D_PARAMSWP_EVAL: Analyze the results obtained by the script MAIN_SIM_1D_PARAMSWP_PROC
% 
% 
% This script examines the effect of parameter choices through simulation MAIN_SWEEP_PARAMETERS_1D. 
% 
% The results are load from 
%% 
% * ../results/sweepparams1d-yyyy-MM-dd-HH-mm.mat
%% 
% which contains
%% 
% * tables: 'tbPropPds','tbPropAdmm','tbRef1dtv'
% * cells:  'configs_proppds','configs_propadmm','configs_ref1dtv'
%% 
% and some materials are stored in folder "../results".
% 
% 
% 
% Requirements: MATLAB R2022a
% 
% 
% 
% Copyright (c) 2022, Ruiki KOBAYASHI and Shogo MURAMATSU, All rights reserved.
% 
% 
% 
% Contact address: Shogo MURAMATSU,
% 
% Faculty of Engineering, Niigata University,
% 
% 8050 2-no-cho Ikarashi, Nishi-ku,
% 
% Niigata, 950-2181, JAPAN
% 
% <http://msiplab.eng.niigata-u.ac.jp/ http://msiplab.eng.niigata-u.ac.jp/>
% 
% 

%
dt = "2022-03-31-09-31";
filename =  "../data/materials/sweepparams1d-" + dt;
disp(filename)
S = load(filename,...
    'tbPropPds','tbPropAdmm','tbRef1dtv',...
    'configs_proppds','configs_propadmm','configs_ref1dtv') %,...
    %'msrparams');
tbPropPds = S.tbPropPds;
tbPropAdmm = S.tbPropAdmm;
tbRef1dtv = S.tbRef1dtv;
configs_proppds = S.configs_proppds;
configs_propadmm = S.configs_propadmm;
configs_ref1dtv = S.configs_ref1dtv;
%msrparams = S.msrparams;
%% Results by proposal w/ PDS

nTrials_ = num2str(max(tbPropPds.iTrial,[],'all'));
disp(['nTrials = ', nTrials_]);
%% 
% Calculate average of MSE for multiple trials

tbPropPds_wo_mse = removevars(tbPropPds,{'mse','iTrial'});
[G,tbPropPds_avemse] = findgroups(tbPropPds_wo_mse);
tbPropPds_avemse.avemse = splitapply(@mean,tbPropPds.mse,G)
%% 
% Find optimal regularization parameters for each sigmaw, nLevels

tbPropPds_avemse.sigmaw = categorical(tbPropPds_avemse.sigmaw);
catSigmaw = categories(tbPropPds_avemse.sigmaw);
tbPropPds_avemse.nLevels = categorical(tbPropPds_avemse.nLevels);
catNumLevels = categories(tbPropPds_avemse.nLevels);

% eta, lambda
tbMinMse = table('Size',[0 6],...
    'VariableTypes',{'categorical','categorical','double','double','double','double'},...
    'VariableNames',{'method','level','sigmaw','minmse','mineta','minlambda'})
for iSigmaw=1:length(catSigmaw)
    sigmaw_ = catSigmaw{iSigmaw};
    for iLevel=1:length(catNumLevels)
        nLevels_ = catNumLevels{iLevel};
        tbPropPds_params = tbPropPds_avemse(tbPropPds_avemse.sigmaw == sigmaw_ & tbPropPds_avemse.nLevels == nLevels_,...
            {'eta','lambda','avemse'});
        figure
        [minmse_,mineta_,minlambda_] = msevsetalambda(tbPropPds_params); 
        tbMinMse = [tbMinMse; {'Prop. w/ PDS', num2str(iLevel), sigmaw_, minmse_, mineta_, minlambda_} ]; %#ok<AGROW> 
        title(['# of tree levels: ' num2str(nLevels_) ', \sigma_w: ' num2str(sigmaw_)])
    end
end

%% Results by proposal w/ ADMM

disp(['nTrials = ', num2str(max(tbPropAdmm.iTrial,[],'all'))]);
%% 
% Calculate average of MSE for multiple trials

tbPropAdmm_wo_mse = removevars(tbPropAdmm,{'mse','iTrial'});
[G,tbPropAdmm_avemse] = findgroups(tbPropAdmm_wo_mse);
tbPropAdmm_avemse.avemse = splitapply(@mean,tbPropAdmm.mse,G)
%% 
% Find optimal regularization parameters for each sigmaw, nLevels, rho

tbPropAdmm_avemse.sigmaw = categorical(tbPropAdmm_avemse.sigmaw);
catSigmaw = categories(tbPropAdmm_avemse.sigmaw);
tbPropAdmm_avemse.nLevels = categorical(tbPropAdmm_avemse.nLevels);
catNumLevels = categories(tbPropAdmm_avemse.nLevels);
tbPropAdmm_avemse.rho = categorical(tbPropAdmm_avemse.rho);
catRho = categories(tbPropAdmm_avemse.rho);
%% 
% Find the best rho

[G,tbPropAdmm_rho] = findgroups(tbPropAdmm.rho);
[~,idx] = min(splitapply(@mean,tbPropAdmm.mse,G));
rho_ = catRho{idx}

% eta, lambda
for iSigmaw=1:length(catSigmaw)
    sigmaw_ = catSigmaw{iSigmaw};
    for iLevel=1:length(catNumLevels)
        nLevels_ = catNumLevels{iLevel};
        tbPropAdmm_params = tbPropAdmm_avemse(tbPropAdmm_avemse.sigmaw == sigmaw_ & ...
            tbPropAdmm_avemse.nLevels == nLevels_ & tbPropAdmm_avemse.rho == rho_,...
            {'eta','lambda','avemse'});
        figure
        [minmse_,mineta_,minlambda_] = msevsetalambda(tbPropAdmm_params); 
        tbMinMse = [tbMinMse; {'Prop. w/ ADMM', num2str(iLevel), sigmaw_, minmse_, mineta_, minlambda_} ]; %#ok<AGROW>
        title(['# of tree levels: ' num2str(nLevels_) ', \sigma_w: ' num2str(sigmaw_) ', \rho: ' num2str(rho_)])
    end
end
%% Results by reference approach w/ PDS

disp(['nTrials = ', num2str(max(tbRef1dtv.iTrial,[],'all'))]);
%% 
% Calculate average of MSE for multiple trials

tbRef1dtv_wo_mse = removevars(tbRef1dtv,{'mse','iTrial'});
[G,tbRef1dtv_avemse] = findgroups(tbRef1dtv_wo_mse);
tbRef1dtv_avemse.avemse = splitapply(@mean,tbRef1dtv.mse,G)
%% 
% Find optimal regularization parameters for each sigmaw

tbRef1dtv_avemse.sigmaw = categorical(tbRef1dtv_avemse.sigmaw);
catSigmaw = categories(tbRef1dtv_avemse.sigmaw);

% eta
for iSigmaw=1:length(catSigmaw)
    sigmaw_ = catSigmaw{iSigmaw};
    tbRef1dtv_params = tbRef1dtv_avemse(tbRef1dtv_avemse.sigmaw == sigmaw_,...
        {'eta','avemse'});
    figure
    [minmse_,mineta_] = msevseta(tbRef1dtv_params); 
    tbMinMse = [tbMinMse; {'Ref. w/ PDS', num2str(0), sigmaw_, minmse_, mineta_, 0 } ]; %#ok<AGROW>
    title(['\sigma_w: ' num2str(sigmaw_)])
end
%% Materials
% Create table of simulation specifications and output 
%% 
% * "../results/tablespec.tex" 
% * "../results/tablemse.tex"
%% 
% Parameters

formatSpec = '%4.2f';
% Sigmaw
sigmaws_ = string(num2str(str2double(catSigmaw),formatSpec));
sigmawset = "$\{";
for idx=1:length(sigmaws_)
    sigmawset = sigmawset.append(sigmaws_(idx));
    if idx<length(sigmaws_)
        sigmawset = sigmawset.append(",\ ");
    else
        sigmawset = sigmawset.append("\ ");
    end
end
sigmawset = sigmawset.append("\}$");

% Eta
etaset_prop_pds = str2double(categories(categorical(tbPropPds.eta)));
etaset_prop_admm = str2double(categories(categorical(tbPropAdmm.eta)));
etaset_ref_pds = str2double(categories(categorical(tbRef1dtv.eta)));
if all(etaset_prop_pds == etaset_prop_admm) && all(etaset_prop_pds == etaset_ref_pds)
    etas = etaset_prop_pds;
else
    etas = 0;
end
etapoints = length(etas);
etaminmax = "["+num2str(min(etas),'%0.2e')+", "+num2str(max(etas),'%0.2e')+"]";

% Lambda
lambdaset_prop_pds = str2double(categories(categorical(tbPropPds.lambda)));
lambdaset_prop_admm = str2double(categories(categorical(tbPropAdmm.lambda)));
if all(lambdaset_prop_pds == lambdaset_prop_admm)
    lambdas = lambdaset_prop_pds;
else
    lambdas = 0;
end
lambdapoints = length(lambdas);
lambdaminmax = "["+num2str(min(lambdas),'%0.2e')+", "+num2str(max(lambdas),'%0.2e')+"]";

% Tree level
% Lambda
levelset_prop_pds = str2double(categories(categorical(tbPropPds.nLevels)));
levelset_prop_admm = str2double(categories(categorical(tbPropAdmm.nLevels)));
if all(levelset_prop_pds == levelset_prop_admm)
    levels_ = categories(categorical(tbPropPds.nLevels));
else
    levels_ = 0;
end
levelset = "$\{";
for idx=1:length(levels_)
    levelset = levelset.append(levels_(idx));
    if idx<length(levels_)
        levelset = levelset.append(",\ ");
    else
        levelset = levelset.append("");
    end
end
levelset = levelset.append("\}$");
%% 
% Configurations

% Prop. w/ PDS
nItems = height(configs_proppds);
tbl_ = struct2table(configs_proppds{1});
for idx = 1:nItems
    tbl_ = [tbl_; struct2table(configs_proppds{idx})]; %#ok<AGROW> 
end
[~,configs_proppds_nIters] = findgroups(tbl_.nIters);
% Prop. w/ ADMM
nItems = height(configs_propadmm);
tbl_ = struct2table(configs_propadmm{1});
for idx = 1:nItems
    tbl_ = [tbl_; struct2table(configs_propadmm{idx})]; %#ok<AGROW> 
end
[~,configs_propadmm_nIters] = findgroups(tbl_.nIters);

% Ref. w/ PDS
nItems = height(configs_ref1dtv);
tbl_ = struct2table(configs_ref1dtv{1});
for idx = 1:nItems
    tbl_ = [tbl_; struct2table(configs_ref1dtv{idx})]; %#ok<AGROW> 
end
[~,configs_ref1dtv_nIters] = findgroups(tbl_.nIters);
%% 
% Generate LaTeX table

spectbl = "";
spectbl = spectbl.append("\begin{tabular}{|c||c|c|c|}\hline"+newline);
spectbl = spectbl.append(" & Ref. w/ PDS & Prop. w/ PDS & Prop. w/ ADMM  \\ \hline\hline"+newline);
%

% Synthesis process
%spectbl = spectbl.append(" \multicolumn{4}{|l|}{Source $\mathbf{u}$, and noise } \\ \hline"+newline);

% Source
spectbl = spectbl.append(" Source $\mathbf{u}$ & \multicolumn{3}{c|}{ 1-D sequence in $[1.00,1.50]^{256}$ } \\ \hline"+newline);
spectbl = spectbl.append(" Target $\mathbf{r}$ & \multicolumn{3}{c|}{ $\bmphi(\mathbf{u})$ in $(-1,1)^{256}$} \\ \hline"+newline);
% Noise
spectbl = spectbl.append(" $\sigma_\mathrm{w}$ & \multicolumn{3}{c|}{" + sigmawset + " } \\ \hline"+newline);
% Algorithm
%spectbl = spectbl.append(" \multicolumn{4}{|l|}{Other settings} \\ \hline"+newline);
% eta,lambda
spectbl = spectbl.append(" $\eta$ & \multicolumn{3}{c|}{" + etapoints + " points in " + etaminmax + " } \\ \hline"+newline);
spectbl = spectbl.append(" $\lambda$ & - & \multicolumn{2}{c|}{" + lambdapoints + " points in " + lambdaminmax +"}  \\ \hline"+newline);
%
% Step sizes
catmu_proppds = categorical(tbPropPds.mu);
catmu_ref1dtv = categorical(tbRef1dtv.mu);
if length(categories(catmu_proppds))==1 && length(categories(catmu_ref1dtv))==1
    spectbl = spectbl.append(" $\mu$ &  "+ num2str(tbRef1dtv.mu(1),'%.4g') + ...
        " & "+num2str(tbPropPds.mu(1),'%.4g') +"& - \\ \hline"+newline);
end
spectbl = spectbl.append(" $\sigma_\mathrm{max}(\mathbf{L})^2$ & 4  & 2 & - \\ \hline"+newline);
spectbl = spectbl.append(" $\gamma_1$ & \multicolumn{2}{c|}{$2/1.05\mu$ } & - \\ \hline"+newline);
spectbl = spectbl.append(" $\gamma_2$ & \multicolumn{2}{c|}{$(1/\gamma_1-\mu/2)/1.05\sigma_\mathrm{max}(\mathbf{L})^2$} & - \\ \hline"+newline);
%gamma1_pds = 2/(1.05*mu);
%gamma2_pds = 1/(1.05*tau2_pds)*(1/gamma1_pds-mu/2);
% rho
spectbl = spectbl.append(" $\rho$ & - & - & " + num2str(str2double(rho_),'%6.3f')+ " \\ \hline"+newline);
%
if configs_proppds_nIters ==  configs_propadmm_nIters && ...
        configs_proppds_nIters ==  configs_ref1dtv_nIters
    spectbl = spectbl.append(" \# of iterations & \multicolumn{3}{c|}{" + num2str(configs_proppds_nIters) + " } \\ \hline"+newline);
end
spectbl = spectbl.append(" \# of trials & \multicolumn{3}{c|}{" + nTrials_ + " } \\ \hline"+newline);


% Measuremnt proc. P
spectbl = spectbl.append(" \multicolumn{4}{|l|}{Measurement process $\mathbf{P}$} \\ \hline"+newline);
spectbl = spectbl.append(" $\omega_\mathrm{p}$ & \multicolumn{3}{c|}{$0.3\pi$} \\ \hline"+newline);
spectbl = spectbl.append(" $\sigma_\mathrm{z}$ & \multicolumn{3}{c|}{30} \\ \hline"+newline);
spectbl = spectbl.append(" $b_\mathrm{p}$ & \multicolumn{3}{c|}{0.05} \\ \hline"+newline);
% Dictionary D
spectbl = spectbl.append(" \multicolumn{4}{|l|}{Synthesis dictionary $\mathbf{D}$} \\ \hline"+newline);
spectbl = spectbl.append(" Type & - & \multicolumn{2}{c|}{1-D Parseval tight UDHT}  \\ \hline"+newline);
spectbl = spectbl.append(" \# of tree levels & - & \multicolumn{2}{c|}{" + levelset + "}  \\ \hline"+newline);
%
spectbl = spectbl.append("\end{tabular}"+newline);
%% 
% Write

fid = fopen('../results/tab1rev.tex','w');
fwrite(fid,spectbl);
fclose(fid);
%% 
% Search minimum MSE for each $\sigma_\mathrm{w}$

minminmsetbl = [];
for idx=1:length(sigmaws_)
    tmp_ =  tbMinMse(tbMinMse.sigmaw == sigmaws_(idx),:);
    [~,imin] = min(tmp_.minmse);
    minminmsetbl = [minminmsetbl; tmp_(imin,:)]; %#ok<AGROW> 
end
minminmsetbl
%% 
% Generate table of estimation results in MSE

nLv_ = length(levels_);
%
msetbl = "\renewcommand{\arraystretch}{1.5}";
msetbl = msetbl.append("\begin{tabular}{|c|c||c|c|c|c|c|c|c|c|c|c|}\hline"+newline);
msetbl = msetbl.append(" \multicolumn{2}{|c||}{$\sigma_\mathrm{w}$} ");
for idx=1:length(sigmaws_)-1
    msetbl = msetbl.append(" & \multicolumn{2}{c|}{"+sigmaws_(idx) +"} ");
end
idx = length(sigmaws_);
msetbl = msetbl.append(" & \multicolumn{2}{c|}{"+sigmaws_(idx) +"} ");
msetbl = msetbl.append(" \\ \hline " + newline);
%

msetbl = msetbl.append("Method & Lv. ");
for idx=1:length(sigmaws_)
    msetbl = msetbl.append(" & MSE & $\begin{smallmatrix}  \lambda \\ \eta \end{smallmatrix}$   ");
end
msetbl = msetbl.append(" \\ \hline " + newline);
% Ref. w/ PDS
minmse_ref_pds = removevars(tbMinMse(tbMinMse.method == "Ref. w/ PDS",:),'method')
msetbl = msetbl.append(" Ref.  w/ PDS & -  " );
for idx=1:length(sigmaws_)
    tmp_ =  minmse_ref_pds(minmse_ref_pds.sigmaw == sigmaws_(idx),["minmse","mineta","minlambda"]);
    if minminmsetbl(minminmsetbl.sigmaw == sigmaws_(idx),:).minmse == tmp_.minmse
        str_minmse = "{\bf "+num2str(tmp_.minmse,'%0.2e')+ "}";
    else
        str_minmse = num2str(tmp_.minmse,'%0.2e');
    end
    str_minlmd = num2str(tmp_.minlambda,'%0.2e');
    str_mineta = num2str(tmp_.mineta,'%0.2e');
    msetbl = msetbl.append(" & "+str_minmse +" &  $\begin{smallmatrix} "+ " - \\" + str_mineta +" \end{smallmatrix}$  ");
end
msetbl = msetbl.append("  \\ \hline "+newline);
% Prop. w/ PDS
minmse_prop_pds = removevars(tbMinMse(tbMinMse.method == "Prop. w/ PDS",:),'method')
for iLevel = 1:nLv_
    if iLevel==1
        msetbl = msetbl.append("\multirow{"+nLv_+"}{*}{\begin{tabular}{c}Prop. \\ w/ PDS \end{tabular}} &   " + levels_(iLevel));
    else
        msetbl = msetbl.append(" &"  + levels_(iLevel) );
    end
    for idx=1:length(sigmaws_)
        tmp_ =  minmse_prop_pds(minmse_prop_pds.sigmaw == sigmaws_(idx) & minmse_prop_pds.level == levels_(iLevel),["minmse","mineta","minlambda"]);
        if minminmsetbl(minminmsetbl.sigmaw == sigmaws_(idx),:).minmse == tmp_.minmse
            str_minmse = "{\bf "+num2str(tmp_.minmse,'%0.2e')+ "}";
        else
            str_minmse = num2str(tmp_.minmse,'%0.2e');
        end
        str_minlmd = num2str(tmp_.minlambda,'%0.2e');
        str_mineta = num2str(tmp_.mineta,'%0.2e');
        msetbl = msetbl.append(" & "+str_minmse +" &  $\begin{smallmatrix} "+ str_minlmd +" \\ " + str_mineta + "\end{smallmatrix}$  ");
    end
    if iLevel == nLv_
        msetbl = msetbl.append(" \\ \hline "+newline);
    else
        msetbl = msetbl.append(" \\ \cline{2-12} "+newline);
    end
end
% Prop. w/ ADMM
minmse_prop_admm = removevars(tbMinMse(tbMinMse.method == "Prop. w/ ADMM",:),'method')
for iLevel = 1:nLv_
    if iLevel==1
        msetbl = msetbl.append("\multirow{"+nLv_+"}{*}{\begin{tabular}{c}Prop. \\ w/ ADMM \end{tabular}} & " + levels_(iLevel) );
    else
        msetbl = msetbl.append(" & "  + levels_(iLevel) );
    end
    for idx=1:length(sigmaws_)
        tmp_ =  minmse_prop_admm(minmse_prop_admm.sigmaw == sigmaws_(idx) & minmse_prop_admm.level == levels_(iLevel),["minmse","mineta","minlambda"]);
        if minminmsetbl(minminmsetbl.sigmaw == sigmaws_(idx),:).minmse == tmp_.minmse
            str_minmse = "{\bf "+num2str(tmp_.minmse,'%0.2e')+ "}";
        else
            str_minmse = num2str(tmp_.minmse,'%0.2e');
        end
        str_minlmd = num2str(tmp_.minlambda,'%0.2e');
        str_mineta = num2str(tmp_.mineta,'%0.2e');
        msetbl = msetbl.append(" & "+str_minmse +" &  $\begin{smallmatrix} "+ str_minlmd +" \\ " + str_mineta + "\end{smallmatrix}$  ");
    end
    if iLevel == nLv_
        msetbl = msetbl.append(" \\ \hline "+newline);
    else
        msetbl = msetbl.append(" \\ \cline{2-12} "+newline);
    end
end

%
msetbl = msetbl.append("\end{tabular}"+newline)
msetbl = msetbl.append("\renewcommand{\arraystretch}{1.0}");
%% 
% Write

fid = fopen('../results/tab2rev.tex','w');
fwrite(fid,msetbl);
fclose(fid);
%% Functions
% Function for visualization

function [msemin,etamin,lambdamin] = msevsetalambda(tb_params)
x = tb_params.eta;
y = tb_params.lambda;
z = tb_params.avemse;
[zmin,idx] = min(z,[],'all');
scatter3(x,y,z,10,z,'filled')
xlabel('\eta')
ylabel('\lambda')
zlabel('MSE')
hold on
%
msemin = zmin;
etamin = x(idx);
lambdamin = y(idx);
plot3(etamin,lambdamin,msemin,'x')
text(etamin,lambdamin,msemin,['\eta=' num2str(etamin) newline '\lambda=' num2str(lambdamin) newline 'MSE=' num2str(msemin)],...
    'Color','blue')
hold off
%
end
% Function for visualization

function [msemin,etamin] = msevseta(tb_params)
x = tb_params.eta;
y = tb_params.avemse;
[ymin,idx] = min(y,[],'all');
scatter(x,y,10,y,'filled')
xlabel('\eta')
ylabel('MSE')
hold on
%
msemin = ymin;
etamin = x(idx);
plot(etamin,msemin,'x')
text(etamin,msemin,['\eta=' num2str(etamin) newline 'MSE=' num2str(msemin)],...
    'Color','blue')
hold off
end