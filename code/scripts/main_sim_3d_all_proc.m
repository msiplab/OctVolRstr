%% Sampling adjustment, waveform estimation and volumetric data restoration
% 
% 
% This script runs all simulations of the proposed method at once. The results 
% are stored in folder "./results" as follows:
% 
% results_sim
%% 
% * tbSimConfs: Table contains configurations
% * tbEvals: Table contains evaluation results
%% 
% tape_sim
%% 
% * tEst: Estimation result for sampling adjustment
% * tdObs: Observed signal for sampling adjustment
% * tObs: Original signal for sampling adjustment
%% 
% glass_sim
%% 
% * gsaEst: Estimation result for interference waveform estimation
% * gsbObs: Observed signal for interference waveform estimation before sampling 
% adjustment
% * gsaObs: Observed signal for interference waveform estimation after sampling 
% adjustment
% * gcost: Cost variation of interference waveform optimization
%% 
% rest_sim
%% 
% * uEst: Refractive index distribution recovered from volumetric data
% * rEst: Recovered reflectance distribution of the volumetric data
% * uSrc: Refractive index distribution for volumetric data
% * rSrc: Reflectance distribution for volumetric data
% * vObs: Observation signal used for volumetric data recovery
% * options: Configuration file for the restoration
%% 
% Requirements: MATLAB R2022a
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
% Reference:
%% 
% * Ruiki Kobayashi, Genki Fujii, Yuta Yoshida, Takeru Ota, Fumiaki Nin, Hiroshi 
% Hibino, Samuel Choi, Shunsuke Ono, and Shogo Muramatsu, "*Sparsity-Aware OCT 
% Volumetric Data Restoration Using an Optical Synthesis Model*,"  submitted to 
% _IEEE Trans. on Computational Imaging_, 2022
%% 
% Copyright © 2022 Ruiki KOBAYASHI and Shogo MURAMATSU

% Set by the results of MAIN_SIM_3D_PARAMSWP (Fig. 10)
iLv = 2;
nSet = 441; 
dtpswp = "2022-04-13-07-55";
filename = "graph-" + num2str(nSet) + "-level-" + num2str(iLv) + "-" + dtpswp;
disp(filename)
% Load result
S = load("./materials/" + filename,'graph','gammas','config');
graph = S.graph;
gammas = S.gammas;
config = S.config;
[minmse,m ] = min(graph(:,1));
mineta = graph(m,3); % mineta 
minlambda = graph(m,2); % minlambda
%%
%
rng(0)
worker = 16;
isDemo = false;
isVerbose  = false;
isVisible  = false;
mymse = @(x,y) norm(x(:)-y(:),2)^2/numel(x);
%% Observation process setup
% Set of noise levels for evaluation
%% 
% * $\sigma_\mathrm{w}\in\{0.01,0.02,\cdots,0.05\}$

sigmawSet = 0.01*(5:-1:1);
%% 
% Set other configurations

depthfactor = [3 2 0.5]; % For sampling adjustment
slopefactor = [0.5 0.25 0]; % For wave estimation
paramDev = 0.3; % For wave estimation
%% 
% Ground truth parameters for measurement process $\mathbf{P}$
%% 
% * $\alpha_\mathrm{p}$
% * $\sigma_\mathrm{xy}$
% * $\sigma_\mathrm{z}$
% * $b_\mathrm{p}$
% * $f_\mathrm{p}=\frac{\omega_\mathrm{p}}{2\pi}$

pScale = config.pScale; % 1;
pSigmaxy = config.pSigmaxy; % 0.5;
pSigmaz = config.pSigmaz; % 30;
pB4Sinc = config.pB4Sinc; % 5e-2;
pFreq  = config.pFreq; %3e-1;

GLASS_REFLECTION = 0.5;
% Configrations for observations

obsx = 64;
obsy = 64;
obsz = 300;
obsSize = [obsx obsy obsz];
assert(all(obsSize==config.obsSize),"Size missmatch")
tmsrProc = support.Coherence3d(...
    'Scale',pScale,...
    'Sigma',pSigmaz,...
    'Sinc',pB4Sinc,...
    'Frequency',pFreq,...
    'Sigmaxy',pSigmaxy);
pKernel = tmsrProc.Kernel;
kernelSize = size(pKernel);

%
rdummy = randn(obsSize);
dummy1_ = tmsrProc.step(rdummy,'Forward');
dummy2_ = support.fcn_coherence3d(rdummy,[pScale pSigmaxy pSigmaz pFreq pB4Sinc],'Forward');
assert(norm(dummy1_-dummy2_,'fro')<1e-15,'Not met')
dummy1_ = tmsrProc.step(rdummy,'Adjoint');
dummy2_ = support.fcn_coherence3d(rdummy,[pScale pSigmaxy pSigmaz pFreq pB4Sinc],'Adjoint');
assert(norm(dummy1_-dummy2_,'fro')<1e-15,'Not met')
%% Restoration settings

nLevels = 2;
splitfactor = []; % Split factor for OLA/OLS method
padsize = 2^(nLevels-1)*ones(1,3); % Padsize for OLA/OLS methods
isintegritytest = false; % Integrity test for OLA/OLS methods
useparallel = false;
usegpu = false;
method = 'udht';
maxIter = 1000;
%% Regulalization parameters for restoration 

sizecomp = false; 
if strcmp(method,'udht')
    bareta     = mineta; %4.89e-2; % Set from Fig. 10
    barlambda  = minlambda; % 5.17e-3;
else
    error('METHOD')
end
%% 
% Range of refractive index
%% 
% * [a,b] = [1.0,1.5]

vrange  = [1.00 1.50];

%% 
% Definition of 3-D Sobel filter
%% 
% * $\mathbf{\Delta}_\mathrm{z}$

d_ = zeros(obsSize); % Fixed by S. Muramatsu
d_(floor(end/2),floor(end/2),floor(end/2)) = 1; % 3-D impulse signal
dltFcn = support.Sobel3d('KernelMode','Normal');
dltImpz = dltFcn.step(d_);
%% 
% Approxmation mode of $\mathbf{\phi}(\cdot)$

% PHIMODE in { 'Linear', 'Signed-Quadratic', 'Reflection' }
phiMode   = 'Linear';
%phiMode   = 'Signed-Quadratic';
%phiMode   = 'Reflection';
%% Configurations for simulation

nTrials = 1;
nSimConfs = length(sigmawSet)*length(slopefactor)*length(depthfactor)*nTrials;
tbSimConfs = cell(nSimConfs,1);
%
iSimulation = 1;
varNamesSimConfs_ = {'slope','depthFreq','sigmaw','iTrial'};
for depthFreq = depthfactor
    for slope = slopefactor
        for sigmaw = sigmawSet
            for iTrial = 1:nTrials
                tbSimConfs{iSimulation} = table('Size',[1,length(varNamesSimConfs_)],...
                    'VariableTypes',{'double','double','double','uint32'},...
                    'VariableNames',varNamesSimConfs_);
                tbSimConfs{iSimulation}(1,:) = ...
                    {slope depthFreq sigmaw iTrial};
                iSimulation = iSimulation + 1;
            end
        end
    end
end
%% Prepare evaluation table

varNamesEvals_ = {...'iSetting',...
    'mseRest','erDepth','depthEst','erRatio','thetahat'};
tbEvals = cell(nSimConfs,1);
%% Display parameters


monint     = 100;
texture    = '2D';
slicePlane = 'YZ';
daspect    = [2 1 1];
obsScale =  20;
estScale = 200;
vdpScale = [1 10];
%% Simulation

% Time stamp
timestamp = char(datetime('now','TimeZone','local','Format','y-MM-d-HH-mm'));

if isDemo
    nSims = 1; %#ok<UNRCH> 
    nWorkers_ = 0;
else
    nSims = height(tbSimConfs);
    nWorkers_ = worker;
end
parfor(iSimulation = 1:nSims, nWorkers_)
%for iSimulation = 1:nSims

    %tic
    tbEvals{iSimulation} = table('Size',[1,5],...
    'VariableTypes',{'double','double','double','double','double'},...
    'VariableNames',varNamesEvals_);
    sigmaw = tbSimConfs{iSimulation}.sigmaw;
    depthFreq = tbSimConfs{iSimulation}.depthFreq;
    slope = tbSimConfs{iSimulation}.slope;
    disp("iSimulation sigmaw depthFreq slope"+newline...
        +num2str(iSimulation)+" "+num2str(sigmaw)+" " ...
        +num2str(depthFreq)+" "+num2str(slope))

%% Sampling adjustment

    % Settings for signal source generation
    tr = zeros(obsSize);
    tr(:,:,round(obsz/4)) = 0.3; % 1st reflection
    tr(:,:,round(obsz/4)*3) = -0.3; % 2nd reflection
    tr = tr + randn(obsSize)*0.1; % Scattering

    % Observation
    tmsrProc = support.Coherence3d(...
        'Scale',pScale,...
        'Sigma',pSigmaz,...
        'Sinc',pB4Sinc,...
        'Frequency',pFreq,...
        'Sigmaxy',pSigmaxy);
    tObs = tmsrProc.step(tr,'Forward') + randn(obsSize)*sigmaw;
    %de = 1e-10*(depthFreq-1);
    devfunc = linspace(1,depthFreq,obsz); %1:(depthFreq-1)/obsz:depthFreq-de;
    for i = 2:obsz
        devfunc(i) = devfunc(i-1) + 1/devfunc(i);
    end
    tdObs = support.fcn_resampling(tObs,devfunc,worker);
%%
    % Sampling adjustment
    len = 100;
    div = 150;
    stsize = div-round(len/obsz*div)+1;
    win = hamming(len,'periodic');
    freq = zeros(stsize,1);
    mval = zeros(stsize,1);
    % Peak search
    %parfor(i = 1:obsx,worker)
    s = [];
    for i=1:obsx
        for j = 1:obsy
            s = abs(stft(squeeze(tdObs(i,j,:)),'Window',win,'OverlapLength',len-obsz/div,'FFTLength',2048));
            [val,index] = max(s(end/2+1:end,:));
            freq = freq+index(:);
            mval = mval+val(:);
        end
    end

    % Generate resampling function


    %{%
    idcs = find(mval>0.5*max(mval,[],'all'));
    p = polyfit(idcs,freq(idcs),1);
    imin = min(idcs);
    imax = max(idcs);
    dc = polyval(p,imin:imin+div-1);
    depthEst = dc(end)/dc(1);
    %}%
    %{
    sizeFreq = size(freq,1);
    t = 1:sizeFreq;
    p = polyfit(t,freq,1);
    dc = polyval(p,round((sizeFreq-div)/4):round((sizeFreq-div)/4)+div-1);
    depthEst = dc(div)/dc(1);
    %}
    %
    %de = 1e-5;
    %devfuncEst = depthEst:-(depthEst-1)/obsz:1+de;
    devfuncEst = linspace(depthEst,1,obsz);
    for i = 2:obsz
        devfuncEst(i) = devfuncEst(i-1) + 1/devfuncEst(i);
    end

    % Resampling
    tEst = support.fcn_resampling(tdObs,devfuncEst,worker);
%%
    % Evaluation
    erDepth = depthEst/depthFreq-1;
    %tbEvals.depthEst(iSetting) = depthEst;
    %tbEvals.erDepth(iSetting) = erDepth;
    tbEvals{iSimulation}.depthEst = depthEst;
    tbEvals{iSimulation}.erDepth = erDepth;
%%
    if isVisible
        disp(num2str(erDepth)) %#ok<UNRCH> 
        stftconfig = [len obsz div];
        %stftconfig = win;
        %stftconfig.len = len;
        %stftconfig.obsz = obsz;
        %stftconfig.div = div;
        %
        hfig = figure;
        fcn_stftview_(hfig,tEst,'Estimation',win,stftconfig)
        drawnow
        %
        hfig = figure;
        fcn_stftview_(hfig,tdObs,'Observation',win,stftconfig)
        drawnow
        %
        hfig = figure;
        fcn_stftview_(hfig,tObs,'Target',win,stftconfig)
        drawnow
    end
%% Estimation of mesurement process

    % Settings for signal source generation
    thick = round(obsz/4);
    p1st = round(obsz/4);
    p2nd = p1st + thick;
    xstar = zeros(obsz,1);
    xstar(p1st) = GLASS_REFLECTION;
    xstar(p2nd) = -GLASS_REFLECTION;

    arstar = slope*ones(2,1); % [ary arx]
    % Simulate glass reflection
    rstar = support.fcn_glasssubstrate(xstar,obsSize,arstar);

    % Observation
    gObs = tmsrProc.step(rstar,'Forward'); % FIXED by S. Muramatsu
    gObs = gObs+randn(size(gObs))*sigmaw;
    gsbObs = support.fcn_resampling(gObs,devfunc,worker);

%% 
% Problem setting    
% 
% $$\{\hat{\mathbf{\theta}},\hat{\mathbf{a}_\mathrm{r}},\hat{\mathbf{x}}\}=\arg\min_{\{\mathbf{\theta},\mathbf{a}_\mathrm{r},\mathbf{x}\}}\frac{1}{2}\|\mathbf{P}_{\mathbf{\theta}}\mathbf{G}_{\mathbf{a}_\mathrm{r}}\mathbf{x}-\mathbf{v}_\mathrm{g}\|_2^2\ 
% \mathrm{s.t.} \|\mathbf{x}\|_0=2\\\hat{\mathbf{r}}_\mathrm{g}=\mathbf{G}_{\hat{\mathbf{a}}_\mathrm{r}}\hat{\mathbf{x}}$$
% 
% Search for glass surfaces  (Peak seach & polynomial fitting)
% 
% $$\{\hat{\mathbf{a}}_\mathrm{r},\hat{\mathbf{x}}\}=\arg\min_{\mathbf{a}_\mathrm{r},\mathbf{x}}\frac{1}{2}\|\mathbf{P}_{\hat{\mathbf{\theta}}}\mathbf{G}_{{\mathbf{a}}_\mathrm{r}}\mathbf{x}-\mathbf{v}_\mathrm{g}\|_2^2\ 
% \mathrm{s.t.} \|\mathbf{x}\|_0\leq 2\\\hat{\mathbf{r}}_\mathrm{g}=\mathbf{G}_{\hat{\mathbf{a}}_\mathrm{r}}\hat{\mathbf{x}}$$
% 
% Measurement process update 
% 
% $$\hat{\mathbf{\theta}}=\arg\min_{\mathbf{\theta}}\frac{1}{2}\|\mathbf{P}_{\mathbf{\theta}}{\hat{\mathbf{r}}_\mathrm{g}}-\mathbf{v}_\mathrm{g}\|_2^2\ 
% $$
% 
% $\left[\nabla_{\mathbf{\theta}}  E(\mathbf{\theta})\right]_i     =\left(\frac{\partial 
% \mathbf{P}_{\mathbf{\theta}}}{\partial \theta_i}\mathbf{r}_\mathrm{g}\right)^T(\mathbf{P}_{\mathbf{\theta}}\mathbf{r}_\mathrm{g}-\mathbf{v}_\mathrm{g})$,
% 
% where
% 
% $E(\mathbf{\theta}):=\frac{1}{2}\|\mathbf{P}_{\mathbf{\theta}}{{\mathbf{r}}_\mathrm{g}}-\mathbf{v}_\mathrm{g}\|_2^2$.

    % Sampling adjusted
    gsaObs = support.fcn_resampling(gsbObs,devfuncEst,worker);

    % Parameter
    thetastar = [pScale pSigmaxy pSigmaz pFreq pB4Sinc];
%%
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

    % Check gradient
    gana = grdfun(thetastar(:),rstar,gObs);
    gnum = zeros(5,1);
    dlttheta = 1e-6;
    for idx = 1:length(thetastar)
        extparam = zeros(5,1);
        extparam(idx) = 1;
        gnum(idx) = (objfun(thetastar(:)+extparam(:)*(dlttheta/2),rstar,gObs)...
            -objfun(thetastar(:)-extparam(:)*(dlttheta/2),rstar,gObs))./dlttheta;
    end
    relerrgrd = max(abs(gana-gnum)./gana);
    assert(relerrgrd<1e-2,'%g: Invalid gradient',relerrgrd);
%%
    % Initalization
    theta0 = thetastar.*(1 + paramDev*( 2*(rand(size(thetastar))-0.5)));
    % Search peak frequency
    fftp = 2^nextpow2(obsz);
    mpspct = squeeze(sum(sum(abs(fft( gsaObs,fftp,3)).^2,1),2))/(obsy*obsx);
    hpspct = mpspct(1:fftp/2);
    omegap = linspace(0,pi*(1-1/fftp),fftp/2);
    fgauss = fit(omegap.',sqrt(hpspct),'gauss1');
    if isVisible
        figure %#ok<UNRCH> 
        plot(omegap,sqrt(hpspct))
        axis([0 pi 0 max(sqrt(hpspct))])
        hold on
        plot(fgauss)
        hold off
    end
    theta0(1) = fgauss.a1/8;      % Scale
    theta0(3) = (2*pi)/fgauss.c1; % Sigmaz NOTE: Too large σz fails glass surface detection 
    theta0(4) = fgauss.b1/(2*pi); % Frequency
    thetaEst = theta0(:);
    disp(thetaEst.')

    % Monitoring
    hpltvg = []; hpltvh = []; htxti = []; htxtp = [];
    if isVisible
        figure %#ok<UNRCH> 
        hpltvg = plot(squeeze(gsaObs(end/2,end/2,:)));
        hold on
        vhat = support.fcn_coherence3d(rstar,thetaEst);
        hpltvh = plot(squeeze(vhat(end/2,end/2,:)));
        htxti = text(210,3,"Phase: " + num2str(0) + ", \#iter: " + num2str(0,'% 3d'),'Interpreter','latex');
        htxtp = text(170,-3,"$\mathbf{\theta}^\star$ "+num2str(thetastar(:).','% 6.3f')+newline...
            +"$\mathbf{\theta}_{0}$ "+num2str(theta0(:).','% 6.3f')+newline...
            +"$\hat{\mathbf{\theta}}$   "+num2str(thetaEst(:).','% 6.3f'),'Interpreter','latex');
        axis([0 obsz -1.1*pScale 1.1*pScale])
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
    rhat = []; hplterr = []; x_ = []; y_ = [];
    for iterout = 1:nItersOuter

        % Search glass surfaces
        [~,arhat,xhat] = support.fcn_glasssearch_pre(thetaEst,gsaObs,GLASS_REFLECTION);
        disp(arhat)
        disp(find(xhat))
        [rhat,arhat,xhat] = support.fcn_glasssearch_pst(arhat,xhat,thetaEst,gsaObs,objfun);
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
                htxtp.String = "$\mathbf{\theta}^\star$ "+num2str(thetastar(:).','% 6.3f')+newline...
                    +"$\mathbf{\theta}_{0}$ "+num2str(theta0(:).','% 6.3f')+newline...
                    +"$\hat{\mathbf{\theta}}$   "+num2str(thetaEst(:).','% 6.3f');
                drawnow
                %
                err = objfun(thetaEst,r_patch,v_patch);
                gcost(iterout,iterin) = err;
                if iterin == 1
                    figure
                    x_ = iterin;
                    y_ = err;
                    hplterr = plot(x_,y_);
                else
                    x_ = [x_ iterin];
                    y_ = [y_ err];
                    hplterr.XData = x_;
                    hplterr.YData = y_;
                end
                drawnow
            end
        end
    end
%%
    % Evaluation
    erRatio = thetaEst(:).'./thetastar(:).'-1;
    %tbEvals.thetahat(iSetting,1:length(thetaEst)) = thetaEst(:).';
    %tbEvals.erRatio(iSetting,1:length(erRatio)) = erRatio;
    tbEvals{iSimulation}.thetahat(1,1:length(thetaEst)) = thetaEst(:).';
    tbEvals{iSimulation}.erRatio(1,1:length(erRatio)) = erRatio;
    % Volumetric data
    gsaEst = support.fcn_coherence3d(rhat,thetaEst);

%%
    %tbEvals
%% Restoration

    fwdDic = []; adjDic = []; gdnFcnG = []; gdnFcnH = [];
    try
        if strcmp(method,'udht')
            import saivdr.dictionary.udhaar.*
            import saivdr.restoration.denoiser.*
            fwdDic = UdHaarSynthesis3dSystem();
            adjDic = UdHaarAnalysis3dSystem('NumberOfLevels',nLevels);
            gdnFcnG = GaussianDenoiserSfth();
            gdnFcnH = GaussianDenoiserSfth();
        elseif strcmp(method,'bm4d')
            import saivdr.dictionary.utility.*
            import saivdr.restoration.denoiser.*
            fwdDic = IdentitySynthesisSystem();
            adjDic = IdentityAnalysisSystem('IsVectorize',false);
            gdnFcnG = GaussianDenoiserBm4d();
            gdnFcnH = GaussianDenoiserSfth();
        else
            error(method+": Invalid METHOD")
        end
    catch
        error('Try execute SETUP script')
    end

    % Target signal
    uSrc = ones(obsSize);
    phtm = phantom('Modified Shepp-Logan',obsz);
    phtm = imresize(phtm,[obsx obsx*2]);
    sliceYZ = permute(phtm,[1 3 2]);
    src = 0.5*repmat(sliceYZ,[1 obsy 1]);
    uSrc(:,:,obsz/2-size(src,3)/2+1:obsz/2+size(src,3)/2) = uSrc(:,:,obsz/2-size(src,3)/2+1:obsz/2+size(src,3)/2)+ src((size(src,1)-obsx)/2+1:(obsx+size(src,1))/2,:,:);

    phi = support.RefractIdx2Reflect();
    % Ideal reflectance distribution
    rSrc = phi.step(uSrc);
    hImg = []; 
    if isVisible
        import saivdr.utility.* %#ok<UNRCH> 
        hImg = figure;
        %
        vdvsrc = VolumetricDataVisualizer(...
            'Texture',texture,...
            'VRange',[0 2]);
        if strcmp(texture,'2D')
            vdvsrc.SlicePlane = slicePlane;
        else
            vdvsrc.DAspect = daspect;
        end

        subplot(2,3,1)
        vdvsrc.step(uSrc);
        title(['Refract. Idx ' slicePlane ' slice'])
        set(gca,'DataAspectRatio',daspect)
        subplot(2,3,4)
        vdpsrc = VolumetricDataPlot(...
            'Direction','Z',...
            'NumPlots',2,...
            'Scales',vdpScale);
        vdpsrc.step(uSrc,rSrc);
        axis([0 size(uSrc,3) -1 2])
        legend('Refraction Idx',['Reflectance x' num2str(vdpScale(2))],'Location','best')
        title('Source')

    end

    % Observation
    vObs = tmsrProc.step(rSrc,'Forward');
    vObs = vObs + sigmaw*randn(size(uSrc)); % FIXED by S. Muramatsu
    vObs = support.fcn_resampling(vObs,devfunc,worker);
           
    vdpobs = [];
    if isVisible
        import saivdr.utility.* %#ok<UNRCH> 
        figure(hImg)
        subplot(2,3,2)
        vdvobs = VolumetricDataVisualizer(...
            'Texture',texture,....
            'VRange',[-1 1],...
            'Scale',obsScale);
        if strcmp(texture,'2D')
            vdvobs.SlicePlane = slicePlane;
        else
            vdvobs.DAspect = daspect;
        end
        vdvobs.step(vObs);

        title(sprintf('Obs %s slice: MSE = %6.4e',slicePlane,mymse(vObs,rSrc)))
        set(gca,'DataAspectRatio',daspect)
        %
        subplot(2,3,5)
        vdpobs = VolumetricDataPlot(...
            'Direction','Z',...
            'NumPlots',2);
        vdpobs.step(vObs,zeros(size(vObs))); % Observation, Initial estimation
        axis([0 size(vObs,3) -1 2])
        legend('Observation','Estimation','Location','best')
        title('Observation')
    end
    % Display settings
    if isVerbose
        disp('-------------------------------') %#ok<UNRCH> 
        disp('Data size')
        disp('-------------------------------')
        disp(['Width × Height × Depth      : ' num2str(obsy) 'x' num2str(obsx) 'x' num2str(obsz)])
        disp('-------------------------------')
        disp('observation model')
        disp('-------------------------------')
        disp(['Scale                       : ' num2str(pScale)])
        disp(['Sigma for depth             : ' num2str(pSigmaz)])
        disp(['Sigma for width and height  : ' num2str(pSigmaxy)])
        disp(['Frequency                   : ' num2str(pFreq)])
        disp(['Frequency fluctuation       : ' num2str(depthFreq)])
        disp(['Noise level                 : ' num2str(sigmaw)])
    end

    % Data monitoring
    phiapx = support.RefractIdx2Reflect(...
        'PhiMode',phiMode,...
        'VRange', vrange);
    hTitle3 = []; vdv = []; vdp = [];
    if isVisible
        import saivdr.utility.* %#ok<UNRCH> 
        vdv = VolumetricDataVisualizer(...
            'Texture',texture,...
            'VRange',[-1 1],...
            'Scale',estScale);
        if strcmp(texture,'2D')
            vdv.SlicePlane = slicePlane;
        else
            vdv.DAspect = daspect;
        end
        r = phiapx.step(vObs);
        %
        figure(hImg)
        subplot(2,3,3)
        vdv.step(r);
        hTitle3 = title(sprintf('Rfl Est(  0): MSE = %6.4e',...
            mymse(vObs,rSrc)));
        set(gca,'DataAspectRatio',daspect)
        %
        subplot(2,3,6)
        vdp = VolumetricDataPlot(...
            'Direction','Z',...
            'NumPlots',2,...
            'Scales',vdpScale);
        vdp.step(vObs,r);
        axis([0 size(vObs,3) -1 2])
        legend('Refraction Idx','Reflectance','Location','best')
        title('Restoration')
    end

    % Restoration
    % Step size calculation
    % Lipschitz constant for the gradient of the fidelity term
    % assuming a Parseval tight dictionary
    msrProc = support.Coherence3d(...
        'Scale',thetaEst(1),...
        'Sigmaxy',thetaEst(2),...
        'Sigma',thetaEst(3),...
        'Frequency',thetaEst(4),...
        'Sinc',thetaEst(5));
    h_ = msrProc.step(dltImpz,'Forward');
    beta1_ = @(x) 2*abs(x(2)-x(1))/(x(2)+x(1))^2;
    mu = beta1_(vrange).^2*max(abs(fftn(h_)).^2,[],'all');
    % Gamma1
    gamma1 = 2/(1.05*mu);
    disp("mu = "+num2str(mu));
    disp("gamma1 = " + num2str(gamma1));
   % disp("gammas = " + num2str(gammas));

    % Instantiation of restoration system
    vtObs = support.fcn_resampling(vObs,devfuncEst,worker);
    pdshshc = support.PdsHsHcOct3(...
        'Observation',    vtObs,...
        'Lambda',         barlambda,...
        'Eta',            bareta,...
        'IsSizeCompensation', sizecomp,...
        'Beta',           mu,...
        'Gamma1',         gamma1,...
        'PhiMode',        phiMode,...
        'VRange',         vrange,...
        'MeasureProcess', msrProc,...
        'Dictionary',     { fwdDic, adjDic },...
        'GaussianDenoiser', { gdnFcnG, gdnFcnH },...
        'SplitFactor',    splitfactor,...
        'PadSize',        padsize,...
        'UseParallel',    useparallel,...
        'UseGpu',         usegpu,...
        'IsIntegrityTest', isintegritytest);


%%
    uEst = []; rEst = [];
    try
        for itr = 1:maxIter
            uEst = pdshshc.step();
            rEst = phiapx.step(uEst);
            if isVerbose && itr==1
                disp(pdshshc)
            end
            % Monitoring
            if isVisible && (itr==1 || mod(itr,monint)==0)
                vdv.step(rEst);
                vdpobs.step(vObs,msrProc((rEst),'Forward'));
                vdp.step(uEst,rEst);
                set(hTitle3,'String',sprintf('Rfl Est (%3d): MSE = %6.4e',itr,mymse(rEst,rSrc)));
                drawnow
            end
        end
    catch
        error('Try execute SETUP script')
    end
    mseRest = mymse(rEst,rSrc);
    tbEvals{iSimulation}.mseRest(1,1) = mseRest;
    disp("MSE: " + num2str(mseRest))
%%
    % Store results
    %dt = char(datetime('now','TimeZone','local','Format','d-MM-y-HH-mm'));
    %targetdir = ['./results/sim_clg_' dt 'lv-' char(string(sigmaw)) '-slope-' char(string(slope)) '-dfreq-' char(string(depthFreq))];
    targetdir = "./results/sim_clg_" +timestamp+"-sgmw-"+num2str(sigmaw)+"-slope-"+num2str(slope)+"-dfreq-"+num2str(depthFreq);
    if exist(targetdir,'dir') ~= 7
        mkdir(targetdir)
    end
    %save(targetdir + "/tape_sim",'tEst','tdObs','tObs')
    support.parsave_tape_sim(targetdir + "/tape_sim"+num2str(iSimulation,'%03d'),tEst,tdObs,tObs)
    %save(targetdir + "/glass_sim",'gsaEst','gsaObs','gcost')
    support.parsave_glass_sim(targetdir + "/glass_sim" +num2str(iSimulation,'%03d'),gsaEst,gsaObs,gsbObs,gcost)

    % Store options
    lambda = pdshshc.LambdaCompensated;
    eta   = pdshshc.EtaCompensated;
    optnames = { 'method', 'barlambda', 'bareta', 'lambda', 'eta', 'gamma1', 'phiMode', 'vrange', 'msrProc', 'fwdDic', 'adjDic', 'gdnFcnG', 'paramDev', 'pScale', 'pSigmaz', 'pB4Sinc', 'pFreq', 'pSigamxy', 'sigmaw'};
   
    options = { method, barlambda, bareta, lambda, eta, gamma1, phiMode, vrange, msrProc, fwdDic, adjDic, gdnFcnG, paramDev, pScale, pSigmaz, pB4Sinc, pFreq, pSigmaxy, sigmaw};
    %{
    options.method = method;
    options.barlambda = barlambda;
    options.bareta    = bareta;
    options.lambda = lambda;
    options.eta    = eta;
    options.gamma1 = gamma1;
    options.phiMode = phiMode;
    options.vrange = vrange;
    options.msrProc = msrProc;
    options.fwdDic = fwdDic;
    options.adjDic = adjDic;
    options.gdnFcnG = gdnFcnG;
    options.paramDev = paramDev;
    options.pScale = pScale;
    options.pSigma = pSigmaz;
    options.pSinc = pB4Sinc;
    options.pFreq  = pFreq;
    options.pSigmaxy = pSigmaxy;
    options.noiseLevel = sigmaw;
    %}
    %save(targetdir + "/rest_sim",'uEst','rEst','uSrc','rSrc','vObs','options')
    support.parsave_rest_sim(targetdir + "/rest_sim"+num2str(iSimulation,'%03d'),uEst,rEst,uSrc,rSrc,vObs,options,optnames)

    %toc
end

%% Store summary of simulation

timestamp = char(datetime('now','TimeZone','local','Format','yyyy-MM-dd-HH-mm'));
targetdir = "./results/sim_results_" +timestamp;
if exist(targetdir,'dir') ~= 7
    mkdir(targetdir)
end
save(targetdir+"/tables",'tbSimConfs','tbEvals')
%%
function fcn_stftview_(hfig,array3d,signalname,win,config) %#ok<DEFNU> 
%win = config.win;
len = config(1); %.len;
obsz = config(2); %.obsz;
div = config(3); %.div;
figure(hfig)
stft(squeeze(array3d(1,1,:)),'Window',win,'OverlapLength',len-obsz/div,'FFTLength',2048);
ylabel('Normalized Frequency (\times\pi radians/sample)')
xlabel('Z')
title("Short-Time Fourier Transform ("+signalname+")")
caxis([-30 30])
yticks(-1:0.5:1)
%
hcb = hfig.Children(1);
hcb.Label.String = 'Magnitude (dB)';
hcb.FontSize = 12;
hax = hfig.Children(2);
hax.FontSize = 12;
end
%% 
% 


% Generate rectangular window
%{
    if isRecWin
        ary_ = mean(arhat(:,1));
        arx_ = mean(arhat(:,2));
        Sigmaz_ = thetahat(3);
        B4Sinc_ = thetahat(5);
        hlenz = round(1.5*Sigmaz_/(B4Sinc_)^0.3/2);
        plenz = max(obsy*ary_,obsx*arx_);
        wrec = zeros(obsz,1);
        for isf = 1:2
            idxz_ = ridcs(isf);
            wrec(max(idxz_-hlenz,1):min(idxz_+hlenz+plenz,obsz)) = 1;
        end
        if isVisibule
            figure
            plot(wrec)
        end
        wrec = permute(wrec,[2 3 1]);
    else
        wrec = 1;
    end
%}