%% MAIN_EXP_3D_REST_GRAPH: Visualize the experimental results 
% This script generates graphs from results obtaind by  MAIN_EXP_3D_REST_PROP 
% 
% The results are load from
%% 
% * ../results/exp_restore_yyyy-MM-dd-HH-mm.mat
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
dtres = "2022-04-16-18-18";
targetdir = "../data/materials/exp_restore_"+dtres;
disp(targetdir)
 
%
if ~exist("+support/vol3d.m","file")
    unzip("https://jp.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/22940/versions/5/download/zip")
    movefile("./vol3d.m","./+support/")
end

agamma = 1.7;
rsgamma = 0.5;
vsgamma = 0.5;
%% Experimental results of restoration 

% Filename for output figures
rfigname = @(x,y) "../results/fig"+num2str(x)+string(char(96+y))+"rev.png";

%%
% Filename for read
targetfile = targetdir+"/rest_exp_mouse";
disp(targetfile)
S = load(targetfile,'uEst','rEst','options','vObs','aObs','sbvObs');
%
vObs = S.sbvObs;
%vObs = S.vObs;
rEst = S.rEst;
options = S.options;
disp(options)
obsz = size(vObs,3);
% Revert sampling adjustment
worker = 16;
depthEst = options.depthEst;
devfunc = linspace(1,depthEst,obsz); %1:(depthFreq-1)/obsz:depthFreq-de;
for i = 2:obsz
    devfunc(i) = devfunc(i-1) + 1/devfunc(i);
end
rEst = support.fcn_resampling(double(rEst),devfunc,worker);
%%
ifig = 14;
isubfig = 1;

% Observation
figure
hvol = fcn_volviewer_(vObs,agamma);
drawnow
exportgraphics(hvol.parent,rfigname(ifig,isubfig),"Resolution",300);
isubfig = isubfig + 1;

figure
subplot(2,1,1)
mval = 50*round(max(abs(vObs(end/2,end/2,:)))/50,0);
hslice = fcn_yzslice_(vObs,[-mval mval],"Y-Z slice",vsgamma);
subplot(2,1,2)
hplt = fcn_zplot_(vObs,[-mval mval],"$v[\mathbf{n}]$");
drawnow
exportgraphics(hslice.Parent,rfigname(ifig,isubfig),"Resolution",300);
isubfig = isubfig+1;
exportgraphics(hplt.Parent,rfigname(ifig,isubfig),"Resolution",300);

ifig = ifig + 1;

%%
isubfig = 1;
% Restoration
figure
hvol = fcn_volviewer_(rEst,agamma);
drawnow
exportgraphics(hvol.parent,rfigname(ifig,isubfig),"Resolution",300);
isubfig = isubfig + 1;
figure
subplot(2,1,1)
mval = 1e-2;
hslice = fcn_yzslice_(rEst,[-mval mval],"Y-Z slice",rsgamma);
subplot(2,1,2)
hplt = fcn_zplot_(rEst,[-mval mval],"$\hat{r}[\mathbf{n}]$");
drawnow
exportgraphics(hslice.Parent,rfigname(ifig,isubfig),"Resolution",300);
isubfig = isubfig+1;
exportgraphics(hplt.Parent,rfigname(ifig,isubfig),"resolution",300);

%
ifig = ifig + 1;
%%

obsz = size(vObs,3);
nz = linspace(obsz/2+1,obsz,obsz/2);
d = 12;

% FFT filter in Z direction
%  NOTE: 
%   The original process below was done manually on an application
%   to select spectra in the Fourier domain. This script uses Gaussian 
%   fitting to mimic the manual work that was done by humans.
%
%  Reference:
%   S. Choi, K. Sato, T. Ota, F. Nin, S. Muramatsu, and H. Hibino,
%   "Multifrequency-swept optical coherence microscopy for highspeed full-
%   field tomographic vibrometry in biological tissues," Biomed. Opt.
%   Express, vol. 8, no. 2, pp. 608–621, Feb. 2017
vObsZ = squeeze(sum(sum(vObs,1),2));
sObsZ = fft(vObsZ);
psObsZ = double(sObsZ.*conj(sObsZ));
gfz = fit(nz.',psObsZ(obsz/2+1:obsz),'gauss1');
fwin = zeros(obsz,1);
fwin(round(gfz.b1-d*sqrt(gfz.c1/2)):round(gfz.b1+d*sqrt(gfz.c1/2))) = 1;
fwin = fwin + flipud(fwin);
fwin = permute(fwin(:),[2 3 1]);
fsObs = fwin.*fftn(vObs);
fObs = ifftn(fsObs);

%%
figure
%
atX = 30;
atY = 30;
isubfig = 1;
% X-Z Slice of previous manual approach
fObsXZ = abs(squeeze(fObs(atY,:,:)));
hslice = fcn_imagesc_(fObsXZ,"X-Z cross-section",["x" "z"]);
exportgraphics(hslice.Parent,rfigname(ifig,isubfig),"Resolution",300);
isubfig = isubfig+1;

% Y-Z Slice of previous manual approach
fObsYZ = abs(squeeze(fObs(:,atX,:)));
hslice = fcn_imagesc_(fObsYZ,"Y-Z cross-section",["y" "z"]);
exportgraphics(hslice.Parent,rfigname(ifig,isubfig),"Resolution",300);
isubfig = isubfig+1;

% X-Z Slice of proposed approach
rEstXZ = abs(squeeze(rEst(atY,:,:)));
hslice = fcn_imagesc_(rEstXZ,"X-Z cross-section",["x" "z"]);
exportgraphics(hslice.Parent,rfigname(ifig,isubfig),"Resolution",300);
isubfig = isubfig+1;
% Y-Z Slice of proposed approach
rEstYZ = abs(squeeze(rEst(:,atX,:)));
hslice = fcn_imagesc_(rEstYZ,"Y-Z cross-section",["y" "z"]);
exportgraphics(hslice.Parent,rfigname(ifig,isubfig),"Resolution",300);
%% Functions for visualization

function himg = fcn_imagesc_(array2d,imtitle,labelsub)
    himg = imagesc(abs(array2d));
    ax = himg.Parent;
    ax.Title.String = imtitle;
    ax.FontSize = 20;
    ax.Colormap = pink;
    ax.View = [90 -90];
    ax.XDir = 'reverse';
    ax.YDir = 'normal';
    ax.XTick = (0:5)*1e3;
    ax.YTick = (0:60)*1e1;
    ax.DataAspectRatio = [5000*3 64*2 1];
    ylabel("$n_\mathrm{"+labelsub(1)+"}$","Interpreter","latex")
    xlabel("depth "+"$n_\mathrm{"+labelsub(2)+"}$","Interpreter","latex")
end

function hvol = fcn_volviewer_(voldata,agamma)
    if nargin < 2
        agamma = 1;
    end
    bvol = zeros(size(voldata));
    gvol = voldata>0;
    rvol = voldata<0;
    cdata = cat(4,rvol,gvol,bvol);
    alpha  = (abs(voldata)/max(abs(voldata(:)))).^agamma;
    hvol = support.vol3d('CData',cdata,'Alpha',alpha);
    view(3)
    ax = hvol.parent;
    ax.XLim = [0 size(voldata,1)];
    ax.YLim = [0 size(voldata,2)];
    ax.ZLim = [0 size(voldata,3)];
    ax.Color = [0 0 0];
    ax.DataAspectRatio = [1 1 100];
    ax.XTick = [0 5]*10;
    ax.YTick = [0 5]*10;
    ax.ZTick = (0:5)*1e3;
    ax.XLabel.String = "$n_\mathrm{x}$";
    ax.XLabel.Interpreter = "latex";
    ax.YLabel.String = "$n_\mathrm{y}$";
    ax.YLabel.Interpreter = "latex";
    ax.ZLabel.String = "$n_\mathrm{z}$";
    ax.ZLabel.Interpreter = "latex";
    ax.FontSize = 18;
end


function hslice = fcn_yzslice_(v,range,tstr,sgamma)
if nargin < 4
    sgamma = 1;
end
[szy,szx,~] = size(v);
yzslice = squeeze(v(:,floor(szx/2),:));
% Y-Z slice
yzslice = (yzslice-range(1))/(range(2)-range(1));
if range(1) >= 0
    cmap = colormap(gray(256)).^sgamma;
else
    cmap = zeros(256,3);
    cmap(:,1) = ([127:-1:0 zeros(1,128)].'/127).^sgamma; % R
    cmap(:,2) = ([zeros(1,128)    0:127].'/127).^sgamma; % G
end
hslice = imshow(yzslice);
ax = hslice.Parent;
ax.Colormap = cmap;
ax.DataAspectRatioMode = "manual";
%
ax.DataAspectRatio = [12 1 1];
%
ax.YLabel.String = "$n_\mathrm{y}$";
ax.YLabel.Interpreter = "latex";
ax.XLabel.String = "$n_\mathrm{z}$";
ax.XLabel.Interpreter = "latex";
ax.XLabel.Position(2) = szy*1.1;
ax.FontSize = 14;
ax.Title.String = tstr;
end


function hplot = fcn_zplot_(v,range,ylb)
[szy,szx,szz] = size(v);
zline = squeeze(v(floor(szy/2),floor(szx/2),:));
% Z plot
hplot = plot(1:szz,zline);
ax = hplot.Parent;
ax.YLabel.String = "Intensity";
ax.XLabel.String = "$n_\mathrm{z}$";
ax.XLabel.Interpreter = "latex";
ax.FontSize = 14;
ax.YLim = range;
ax.YTick = linspace(range(1),range(2),3);
ax.XTick = 0:1000:5000;
ax.DataAspectRatioMode = "manual";
%
ax.DataAspectRatio = [round(1000/(range(2)-range(1))) 1 1];
ax.XLabel.Position(2) = range(1)-(range(2)-range(1))*(5/12);
ax.YLabel.String = ylb;
ax.YLabel.Interpreter = "latex";
end