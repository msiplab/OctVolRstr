function [rhat, mu, gammas] = fcn_refpds1dtv(v,config)
%FCN_REFPDS1DTV PDS realization for solving the proposed restoration model
%
nSamples = length(v);
nIters_1dtv = config.nIters;
eta_1dtv = config.eta;
measureproc = config.measurement{1};
measureadjp = config.measurement{2};

%% Impluse response of the first-order difference filter
fdf = [1 -1];
% Measurement process
diffproc = @(x) circshift(cconv(fdf,x,nSamples),-floor(length(fdf)/2));
diffadjp = @(x) circshift(cconv(fdf(end:-1:1),x,nSamples),mod(length(fdf)-1,2)-floor(length(fdf)/2)); 

%% Initialization
yp = zeros(size(v),'like',v);
rp = v;

%% Stepsize parameters 
d_ = zeros(1,nSamples); d_(1) = 1; % Impluse signal
irp = measureproc(d_);
mu = max(abs(fft(irp,nSamples)).^2,[],'all'); % (σmax(P))^2
tau2_tv = max(abs(fft(fdf,nSamples)).^2,[],'all');  % (σmax(G))^2
gamma1_tv = 2/(1.05*mu);
gamma2_tv = 1/(1.05*tau2_tv)*(1/gamma1_tv-mu/2);
assert((1/gamma1_tv - gamma2_tv*tau2_tv) > mu/2,...
    ['Step size condition is violated. γ1 must be less than ' num2str(2/mu)])
gammas = [gamma1_tv,gamma2_tv];

%% Primal-dual splitting method for 1D TV
rmax = 1;
rmin = -1;
for itr=1:nIters_1dtv
    % Primal step
    rg = measureadjp(measureproc(rp)-v);
    rc = metricproj_(rp-gamma1_tv*(rg + diffadjp(yp)),rmin,rmax);
    % Dual step
    yt = yp + gamma2_tv*diffproc(2*rc - rp);
    yc = yt - softthresh_(yt,eta_1dtv); 
    % Update 
    rp = rc;
    yp = yc;
end
rhat = rc;
end

%% Local functions
function y = softthresh_(x,t)
y = sign(x).*max(abs(x)-t,0);
end

%
function y = metricproj_(x,ymin,ymax)
y= min(max(x,ymin),ymax);
end
