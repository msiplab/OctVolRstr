%% MAIN_SIM_PARAMSWP_3D: Search for optimal values of λ and η
% 
% 
% This script searchs for the best lambda and eta for restoration
% 
% The results are stored in folder "./results" as follows:
%% 
% * graph: Error evaluation of lambda and eta for each tree level
%% 
% Requirements: MATLAB R2020b
% 
% 
% 
% Copyright (c) 2021-2022, Ruiki KOBAYASHI and Shogo MURAMATSU, All rights reserved.
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

minLevel = 2;
maxLevel = 2;
%% Search setting

sizecomp = false;
startlambda = 1e-1;
starteta = 2e0;
div = 5;
range = 4;
nLambdaSet = div*range+1;
nEtaSet = div*range+1;
xlambda = zeros(nLambdaSet,1);
yeta = zeros(nEtaSet,1);
nSet = nLambdaSet*nEtaSet;
graph = zeros(nSet,5);
gammas = zeros(nSet,2);
%% Initial setting

import support.Coherence3d

sigmaw = 0.05; % noiseLevel FIXED by S. Muramatsu on 30 Mar. 2022
pScale = 1; % 4; 10; FIXED by S. Muramatsu on 29 Mar. 2022, Revised 10 Apr. 2022
pSigmaz = 30;
pB4Sinc = 5e-2;
pFreq  = 3e-1; % omegap = 0.6pi
pSigmaxy = 0.5;

obsx = 64;
obsy = 64;
obsz = 300;
obsSize = [obsx obsy obsz];

tmsrProc = Coherence3d(...
    'Scale',pScale,...
    'Sigma',pSigmaz,...
    'Sinc',pB4Sinc,...
    'Frequency',pFreq,...
    'Sigmaxy',pSigmaxy);
pKernel = tmsrProc.Kernel;
kernelSize = size(pKernel);

splitfactor = []; % Split factor for OLA/OLS methods
padsize = 2^0*ones(1,3); % Pad size for OLA/OLS methods
isintegritytest = false; % Integrity test for OLA/OLS methods
useparallel = false;
usegpu = false;

maxIter = 1000;
%% Display parameters

isVerbose  = false;
isVisible  = false;
monint     = 100;
texture    = '2D';
slicePlane = 'YZ';
daspect    = [2 1 1];
obsScale =  20;
estScale = 200;
vdpScale = [1 10];
phiMode   = 'Linear';
vrange  = [1.00 1.50];
%% Gamma1 calculation

import support.Sobel3d
% Lipschitz constant for the gradient of the fidelity term 
% assuming a Parseval tight dictionary
d_ = zeros(obsSize); %size(vObs));
d_(floor(end/2),floor(end/2),floor(end/2)) = 1; % 3-D impulse signal
dltFcn = Sobel3d('KernelMode','Normal');
h_ = tmsrProc.step(dltFcn.step(d_),'Forward');
beta1_ = @(x) 2*abs(x(2)-x(1))/(x(2)+x(1))^2;
mu = beta1_(vrange).^2*max(abs(fftn(h_)).^2,[],'all');
% Gamma1 
gamma1 = 2/(1.05*mu);
%% Target signal

uSrc = ones(obsSize);
phtm = phantom('Modified Shepp-Logan',obsz);
phtm = imresize(phtm,[obsx obsx*2]);
sliceYZ = permute(phtm,[1 3 2]);
src = 0.5*repmat(sliceYZ,[1 obsy 1]);
uSrc(:,:,obsz/2-size(src,3)/2+1:obsz/2+size(src,3)/2) = uSrc(:,:,obsz/2-size(src,3)/2+1:obsz/2+size(src,3)/2)+ src((size(src,1)-obsx)/2+1:(obsx+size(src,1))/2,:,:);
%% Observation

import support.RefractIdx2Reflect

rng(0)
phi = RefractIdx2Reflect();
vObs = tmsrProc.step(phi.step(uSrc),'Forward');
vObs = vObs + sigmaw*randn(size(uSrc)); % FIXED by S. Muramatsu 
mymse = @(x,y) norm(x(:)-y(:),2)^2/numel(x);
phiapx = RefractIdx2Reflect(...
    'PhiMode',phiMode,...
    'VRange', vrange);

% Copy for parallel process
phi_ = cell(nSet,1);
phiapx_ = cell(nSet,1);
for iSetting = 1:nSet
    phi_{iSetting} = phi;
    phiapx_{iSetting} = phiapx;
end
%% Restoration

import support.PdsHsHcOct3
import saivdr.dictionary.udhaar.*
import saivdr.restoration.denoiser.*
for iLv = minLevel:maxLevel
    nLevels = iLv;
    tic
    parfor iSetting = 1:nSet
        uEst = [];
        % Dictionary setting
        fwdDic = UdHaarSynthesis3dSystem();
        adjDic = UdHaarAnalysis3dSystem('NumberOfLevels',nLevels);
        gdnFcnG = GaussianDenoiserSfth();
        gdnFcnH = GaussianDenoiserSfth();

        % Parameters
        barlambda = startlambda * 10^(-floor((iSetting-1)/nEtaSet)/div);
        bareta    = starteta * 10^(-mod(iSetting-1,nEtaSet)/div);

        % Restoration setting
        pdshshc = PdsHsHcOct3(...
            'Observation',    vObs,...
            'Lambda',         barlambda,...
            'Eta',            bareta,...
            'IsSizeCompensation', sizecomp,...
            'Beta',           mu,...
            'Gamma1',         gamma1,...
            'PhiMode',        phiMode,...
            'VRange',         vrange,...
            'MeasureProcess', tmsrProc,...
            'Dictionary',     { fwdDic, adjDic },...
            'GaussianDenoiser', { gdnFcnG, gdnFcnH },...
            'SplitFactor',    splitfactor,...
            'PadSize',        padsize,...
            'UseParallel',    useparallel,...
            'UseGpu',         usegpu,...
            'IsIntegrityTest', isintegritytest);
        gamma1_ = pdshshc.Gamma1;
        gamma2_ = pdshshc.Gamma2;

        % Restoration
        for itr = 1:maxIter
            uEst = pdshshc.step();
        end
        mse_ = mymse(phiapx_{iSetting}.step(uEst),phi_{iSetting}.step(uSrc));
        eta_ = pdshshc.EtaCompensated;
        lambda_ = pdshshc.LambdaCompensated;
        graph(iSetting,:) = [mse_ lambda_ eta_ barlambda bareta];
        gammas(iSetting,:) = [gamma1_ gamma2_];
        disp(['-Eta-' num2str(eta_) '-Lambda-' num2str(lambda_) '-MSE-' num2str(mse_)])
    end
    toc

    % Store results
    config.sigmaw = sigmaw; % noiseLevel
    config.pScale = pScale;
    config.pSigmaz = pSigmaz;
    config.pB4Sinc = pB4Sinc;
    config.pFreq  = pFreq;
    config.pSigmaxy = pSigmaxy;
    config.obsSize = obsSize;

    t = char(datetime('now','Format','yyyy-MM-dd-HH-mm'));
    save(['./results/graph-' num2str(nSet) '-level-' num2str(iLv) '-' t],...
        'graph','gammas','nSet','nLambdaSet','nEtaSet','mu','gamma1','config')
    
end