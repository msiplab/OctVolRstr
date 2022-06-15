%% Summaraize the results obtained by the script MAIN_SIM_3D_ALL_PROC
% This script generates tabless from results obtaind by  MAIN_SIM_3D_ALL_PROC. 
% 
% The results are load from
%% 
% * ./results/sim_results_dd-MM-yyy-HH-mm.mat
%% 
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

% Filenames
dtres = "2022-04-14-19-41";
%% Load simulation configurations

targetdir = "./materials/sim_results_" + dtres;
cfilename = targetdir + "/tables";
S = load(cfilename,"tbEvals","tbSimConfs");
tbEvals = S.tbEvals;
tbSimConfs = S.tbSimConfs;
%%
tball = [vertcat(tbEvals{:}),vertcat(tbSimConfs{:})];
disp(tball)
tball.sigmaw = categorical(tball.sigmaw);

%% Summarize resutls and generate tables

[G,ID] = findgroups(tball(:,{'sigmaw'}));
erSample = splitapply(@(x) mean(abs(x)),tball.erDepth,G)
erTheta = splitapply(@(x) mean(abs(x)),tball.erRatio,G)
erWave = mean(erTheta,2);
mses = splitapply(@mean,tball.mseRest,G)
%
taba = "\begin{tabular}{|c||c|c|c|} \hline" + newline;
taba = taba.append("Noise Lv. $\sigma_\mathrm{w}$ & Sampling (\%) & Wave (\%) & MSE \\ \hline\hline" + newline);
for idx = 1:height(ID)
    taba = taba.append(string(ID.sigmaw(idx)) + " & ");
    taba = taba.append(num2str(100*erSample(idx),'%5.2f')+" & ");
    taba = taba.append(num2str(100*erTheta(idx),'%5.2f')+" & ");
    taba = taba.append(num2str(mses(idx),'%5.2e'));
    taba = taba.append("\\ \hline " + newline);
end
taba = taba.append("\end{tabular}"+newline);
disp(taba)
%
fid = fopen('./results/tab3arev.tex','w');
fwrite(fid,taba);
fclose(fid);
%%
tabb = "\begin{tabular}{|c||c|c|c|c|c|} \hline" + newline;
tabb = tabb.append("$\sigma_\mathrm{w}$ & $\alpha_\mathrm{p}$ (\%) & $\sigma_\mathrm{xy}$ (\%) & $\sigma_\mathrm{z}$ (\%) & $\omega_\mathrm{p}$ (\%) & $b_\mathrm{p}$ (\%)\\ \hline\hline"+newline);
for idx = 1:height(ID)
    tabb = tabb.append(string(ID.sigmaw(idx)));
    for iparam = 1:5
        tabb   = tabb.append(" & " + num2str(100*erTheta(idx,iparam),'%5.2f'));
    end
    tabb = tabb.append("\\ \hline " + newline);
end
tabb = tabb.append("\end{tabular}"+newline);
disp(tabb)
%
fid = fopen('./results/tab3brev.tex','w');
fwrite(fid,tabb);
fclose(fid);