function [rhat, mu, gammas] = fcn_proppds1d(v,config)
%FCN_PROPPDS1D PDS realization for solving the proposed restoration model
%
nSamples = length(v);
nIters_pds = config.nIters;
eta_prop = config.eta;
lambda_prop = config.lambda;
measureproc = config.measurement{1};
measureadjp = config.measurement{2};
analyzer = config.dictionaries{1};
synthesizer = config.dictionaries{2};
ri_lb = config.vrange(1);
ri_ub = config.vrange(2);

%% Difference operators
% Impluse response of the first-order difference filter
dltz = [1 0 -1]/2;
% Convolutional operations
dltzproc = @(x) circshift(cconv(dltz,x,nSamples),-floor(length(dltz)/2));
dltzadjp = @(x) -dltzproc(x);

%% Initialization
beta1_ = fcn_beta1_(ri_lb,ri_ub); % β1(a,b)
u0 =  beta1_*dltzproc(measureadjp(v)); % β1Δz.'P.'v
%u0 = (ri_ub-ri_lb)*rand(size(v))+ri_lb; % Uniform distribution [a,b]
%u0 = ri_lb*ones(size(v),'like',v); 
xp = analyzer(u0); % x0 = D.'s0
y1p = zeros('like',v); % y10 = 0;
y2p = zeros('like',v); % y20 = 0;
qp = u0; 

%% Stepsize parameters
%irpd = cconv(irp,dltz,nSamples); % Impulse response of convolutional operator PΔz
irpd = measureproc(dltz);
mu = beta1_.^2*max(abs(fft(irpd,nSamples)).^2,[],'all'); % (beta1(a,b)*σmax(PΔz))^2
dltzLambdaMax = max(abs(fft(dltz,nSamples)).^2,[],'all'); % λmax(Δz.'*Δz)
tau2_pds = dltzLambdaMax + 1; % (σmax(L))^2 = λmax(Δz.'*Δz) + 1
gamma1_pds = 2/(1.05*mu);
gamma2_pds = 1/(1.05*tau2_pds)*(1/gamma1_pds-mu/2);
assert((1/gamma1_pds - gamma2_pds*tau2_pds) > mu/2,...
    ['Step size condition is violated. γ1 must be less than ' num2str(2/mu)])
gammas = [gamma1_pds,gamma2_pds];

%% Primal-dual splitting for proposal
for itr=1:nIters_pds
    % Primal step
    nablaF = beta1_*dltzproc( measureadjp( beta1_*measureproc(dltzadjp(qp)) - v ) );
    t_ = analyzer(nablaF + dltzadjp(y1p) + y2p);
    xc = softthresh_( xp - gamma1_pds*t_, lambda_prop*gamma1_pds );
    qc = synthesizer(xc);
    % Dual step
    u_ = 2*qc - qp;
    y1p = y1p + gamma2_pds * dltzproc(u_);
    y2p = y2p + gamma2_pds * u_;
    y1c = y1p - softthresh_(y1p,eta_prop);
    y2c = y2p - gamma2_pds * metricproj_(y2p/gamma2_pds,ri_lb,ri_ub);
    % Update
    xp = xc;
    y1p = y1c;
    y2p = y2c;
    qp = qc;
end
rhat = phi1d_linear_(qc,[ri_lb ri_ub]);
end

%% Local functions
function y = softthresh_(x,t)
y = sign(x).*max(abs(x)-t,0);
end

%
function y = metricproj_(x,ymin,ymax)
y= min(max(x,ymin),ymax);
end

%
function y = phi1d_linear_(u,ribound)
dz = [1 0 -1];
dzu = circshift(cconv(dz,u,length(u)),-1);
if nargin < 2
    ribound = [1.0 1.5];
end
b = ribound(2);
a = ribound(1);
y = -fcn_beta1_(a,b)*dzu;
end

function b = fcn_beta1_(a,b)
b = 2*abs(b-a)/(b+a)^2;
end
