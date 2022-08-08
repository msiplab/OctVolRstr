%% MAIN_EXP_3D_WAVEFORM_EST: Waveform estimation experiments
% 
% 
% This script is an experiment with real data of interference waveform estimation.
% 
% The results are stored in folder "../results" as follows:
%% 
% * ObsData: the observed data, adjusted for sampling and stripped of DC components
% * EstData: Estimated glass data.
% * lnmse: error between estimated and observed data at each training iteration
% * x: the estimated refractive index distribution of the glass
% * theta: parameters required for the obtained interference waveform
%% 
% Requirements: MATLAB R2020b
% 
% 
% 
% Copyright (c) 2021-2022, Ruiki KOBAYASHI and Shogo MURAMATSU
% 
% 
% 
% All rights reserved.
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
% http://msiplab.eng.niigata-u.ac.jp/
%% Initial setting

S = load('../data/materials/glass1.mat','ObsData');
ObsData = S.ObsData;
%dtsa = "2022-04-14-11-45";
dtsa = "2022-04-16-21-20";
S = load("../data/materials/exp_smpadj_"+dtsa+"/adjustFunc","depthEst");
depthEst = S.depthEst;
disp(depthEst)
pSigmaxy = 0.5; 
pB4Sinc = 1e-2; 
[obsy,obsx,obsz] = size(ObsData);
GLASS_REFLECTION = 0.5;
glsdev = 3;
%
isVisible  = true;
mymse = @(x,y) norm(x(:)-y(:),2)^2/numel(x);

%
dt = char(datetime('now','TimeZone','local','Format','yyyy-MM-dd-HH-mm'));
targetdir = "../results/exp_wave_" + dt;
if exist(targetdir,'dir') ~= 7
    mkdir(targetdir)
end
%% Bias removal

nLen = 21;
lpf = ones(nLen,1)/nLen;
lpf = permute(lpf,[2 3 1]);
gsbObs = ObsData - imfilter(ObsData,lpf,'symmetric');
%% Sampling adjustment

funcEst = linspace(depthEst,1,obsz); %depthEst:-(depthEst-1)/obsz:1+1e-10 depthEst:-(depthEst-1)/obsz:1+1e-10;
for i = 2:obsz
    funcEst(i) = funcEst(i-1) + 1/funcEst(i);
end

% Sampling adjusted
gsaObs = support.fcn_resampling(gsbObs,funcEst,0);

%% Objective function

% Objective function
objfun = @(theta,rg,vg) (1/2)*sum( (support.fcn_coherence3d(rg,theta)-vg ).^2 ,'all');

% Gradient of objective function
grdfun = @(theta,rg,vg) [ ... % pScale pSigmaxy pSigmaz pFreq pB4Sinc
    reshape( support.fcn_coherence3d(rg,theta,'Scale'),1,[]);
    reshape( support.fcn_coherence3d(rg,theta,'Sigmaxy'),1,[]);
    reshape( support.fcn_coherence3d(rg,theta,'Sigmaz'),1,[]);
    reshape( support.fcn_coherence3d(rg,theta,'Freq'),1,[]);
    reshape( support.fcn_coherence3d(rg,theta,'Sinc'),1,[])
    ] * reshape( support.fcn_coherence3d(rg,theta)-vg,[],1 );

%% Search peak frequency

% Search peak frequency
fftp = 2^nextpow2(obsz);
mpspct = squeeze(sum(sum(abs(fft(gsaObs,fftp,3)).^2,1),2))/(obsy*obsx);
hpspct = mpspct(1:fftp/2);
omegap = linspace(0,pi*(1-1/fftp),fftp/2);
fgauss = fit(omegap.',sqrt(hpspct),'gauss2');
amx = fgauss.a1;
bmx = fgauss.b1;
cmx = fgauss.c1;
if isVisible
    figure %#ok<UNRCH>
    plot(omegap,sqrt(hpspct))
    axis([0 pi 0 max(sqrt(hpspct))])
    hold on
    plot(fgauss)
    hold off
end
%%
% Initalization
theta0 = zeros(5,1);
theta0(1) = amx/8;      % Scale
theta0(2) = pSigmaxy;
theta0(3) = (2*pi)/cmx; % Sigmaz
theta0(4) = bmx/(2*pi); % Frequency
theta0(5) = pB4Sinc;
thetaEst = theta0(:);
disp(thetaEst.')
scaleunit = 10^(ceil(log10(theta0(1)))-2);
pScale = 600;
% Monitoring
if isVisible
    figure %#ok<UNRCH>
    hpltvg = plot(squeeze(gsaObs(end/2,end/2,:)));
    hold on
    vhat = support.fcn_coherence3d(zeros(size(gsaObs)),thetaEst);
    hpltvh = plot(squeeze(vhat(end/2,end/2,:)));
    htxti = text(floor(obsz/2),0.9*pScale,"Phase: " + num2str(0) + ", \#iter: " + num2str(0,'% 3d'),'Interpreter','latex');
    htxtp = text(floor(obsz/2),-0.9*pScale,..."$\mathbf{\theta}^\star$ "+num2str(theta0(:).','% 6.3f')+newline+...
        "$\mathbf{\theta}_{0}$ "+num2str(theta0(:).','% 6.3f')+newline...
        +"$\hat{\mathbf{\theta}}$   "+num2str(thetaEst(:).','% 6.3f'),'Interpreter','latex');
    axis([0 obsz -pScale pScale])
    drawnow
end

% Search for glass surfaces & Measurement process update
isAdam = true;
isSgd = true;
%
nItersOuter = 1; %2;
nItersInner = 1000; % 500;
if isAdam
    gdsteps = 1e-3*ones(nItersOuter,5); % [pScale pSigmaxy pSigmaz pFreq pB4Sinc]
    gdsteps(:,1) = 1e-4; % pScale
    %gdsteps(:,2) = 2e-3; % pSigmaxy
    %gdsteps(:,3) = 2e-3; % pSigmaz
    gdsteps(:,4) = 1e-5; % pFreq
    gdsteps(:,5) = 1e-4; % pB4Sinc
else
    gdsteps = 1e-9*ones(nItersOuter,5); %#ok<UNRCH>
    % [pScale pSigmaxy pSigmaz pFreq pB4Sinc]
end

gcost = zeros(nItersOuter,nItersInner);
for iterout = 1:nItersOuter

    % Search glass surfaces
    [~,arhat,xhat] = support.fcn_glasssearch_pre(thetaEst,gsaObs,GLASS_REFLECTION);
    disp(arhat)
    disp(find(xhat))
    [rhat,arhat,xhat] = support.fcn_glasssearch_pst(arhat,xhat,thetaEst,gsaObs,objfun,glsdev);
    disp(arhat)
    disp(find(xhat))

    % Interference wave estimation
    beta1_ = []; beta2_ = []; mt_ = []; vt_ = []; %#ok<NASGU>
    if isAdam
        beta1_ = 0.9;
        beta2_ = 0.999;
        mt_ = zeros(size(thetaEst));
        vt_ = zeros(size(thetaEst));
    end

    gdstep_ = gdsteps(iterout,:);
    for iterin = 1:nItersInner
        if isSgd
            patchsz = [obsy obsx]/4;
            %index_y = randi(obsx-patchsz(1));
            %index_x = randi(obsy-patchsz(2));
            stdevy = (obsx-patchsz(1))/3;
            stdevx = (obsy-patchsz(2))/3;
            index_y = min(max(round(stdevy*randn(1)+floor(obsy/2-patchsz(1)/2)),1),obsy-patchsz(2));
            index_x = min(max(round(stdevx*randn(1)+floor(obsx/2-patchsz(2)/2)),1),obsx-patchsz(1));
            r_patch = rhat(index_y:index_y+patchsz(2),index_x:index_x+patchsz(2),:);
            v_patch = gsaObs(index_y:index_y+patchsz(1),index_x:index_x+patchsz(2),:);
        else
            r_patch = rhat; %#ok<UNRCH>
            v_patch = gsaObs;
        end
        if isAdam
            gt_ = grdfun(thetaEst(:),r_patch,v_patch);
            mt_ = beta1_*mt_ + (1-beta1_)*gt_;
            vt_ = beta2_*vt_ + (1-beta2_)*gt_.^2;
            mt_hat = mt_./(1-beta1_^(iterin));
            vt_hat = vt_./(1-beta2_^(iterin));
            thetaEst = max(thetaEst(:) - gdstep_(:).*mt_hat./(sqrt(vt_hat)+sqrt(eps)),0);
        else
            thetaEst = max(thetaEst(:) - gdstep_(:).*grdfun(thetaEst(:),r_patch,v_patch),0); %#ok<UNRCH>
        end
        if isVisible
            hpltvg.YData = squeeze(v_patch(floor(end/2),floor(end/2),:)); %#ok<UNRCH>
            vhat = support.fcn_coherence3d(r_patch,thetaEst);
            hpltvh.YData = squeeze(vhat(floor(end/2),floor(end/2),:)).';
            htxti.String = "Phase: " + num2str(iterout) + ", \#iter: " + num2str(iterin,'% 3d');
            htxtp.String = ..."$\mathbf{\theta}^\star$ "+num2str(theta0(:).','% 6.3f')+newline+...
                "$\mathbf{\theta}_{0}$ "+num2str(theta0(:).','% 6.3f')+newline...
                +"$\hat{\mathbf{\theta}}$   "+num2str(thetaEst(:).','% 6.3f');
            drawnow
            %
            err = objfun(thetaEst,r_patch,v_patch);
            gcost(iterout,iterin) = err;
            if iterin == 1
                figure
                x = iterin;
                y = err;
                hplterr = plot(x,y);
            else
                x = [x iterin];
                y = [y err];
                hplterr.XData = x;
                hplterr.YData = y;
            end
            drawnow
        end
    end
end
%%

%EstData = msrProc.step(x,'Forward');
EstData = support.fcn_coherence3d(rhat,thetaEst);

%% Display results

figure
plot(squeeze(gsaObs(end/2,end/2,:)))
hold on
plot(squeeze(EstData(end/2,end/2,:)))
hold off
theta = thetaEst;
%% Store results

filename = targetdir + "/glass_exp";
save(filename,'ObsData','gsaObs','rhat','EstData','gcost','xhat','theta')