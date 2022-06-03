function rhat = fcn_propadmm1d(v,config)
%FCN_PROPADMM1D ADMM realization for solving the proposed restoration model
%   
nSamples = length(v);
nIters_admm = config.nIters;
rho_admm = config.rho;
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

%% Preperation to generate process by matrix Q
d_ = zeros(1,nSamples); d_(1) = 1; % Impluse signal
t_ = dltzproc(d_); % Impluse response of Δz
beta1_ = fcn_beta1_(ri_lb,ri_ub); % β1(a,b)
irh_ = dltzadjp( beta1_^2*measureadjp(measureproc(t_)) + rho_admm*(eta_prop/lambda_prop)^2*t_ ); % Impluse response of H
h_ = fft(irh_,nSamples); % DFT of H
%H_ = diag(h_);
diagH_ = h_(:);
% Preperation to generate process by matrix inv(Q)
idxnz_ = find(abs(h_)>1e-15); % Extract indices of non-zero elements of h
invCheckH_ = diag(1./h_(idxnz_));  % Inverse of Non-zero DFT Coefs. of H
%invCheckH_ = inv(checkH_);
In_ = eye(nSamples); % NxN identity matrix
S_ = In_(idxnz_,:); % Downsampling matrix S to extract non-zero elements in h
Ik_ = eye(length(idxnz_)); % KxK identity matrix, where K is # of non-zero elements in h
R_ = S_.'*((invCheckH_+Ik_/(2*rho_admm))\S_); % Matrix R
diagR_ = diag(R_);
% Check if invpQ is inverse of procQ
nCoefs = length(analyzer(v)); %(nLevels+1)*nSamples;
y_ = randn(nCoefs,1);
q_ = fcn_procQ_(y_,diagH_,rho_admm,{analyzer, synthesizer},nSamples);
t_ = fcn_invpQ_(q_,diagR_,rho_admm,{analyzer, synthesizer},nSamples);
assert(norm(t_-y_)<1e-3,'Process by matrix Q is not inverted',norm(t_-y_))

%% Initialization
y_ = analyzer(beta1_*dltzproc(measureadjp(v)));
%procG1 = @(x) [ x(:) ; (eta_prop/lambda_prop)*reshape(dltzproc(synthesizer(x)),[],1) ];
%procG2 = @(x) reshape(synthesizer(x),[],1);
%adjpG1 = @(x) reshape(x(1:nCoefs),[],1) ...
%    + (eta_prop/lambda_prop)*reshape(analyzer(dltzadjp(x(nCoefs+1:end))),[],1);
%adjpG2 = @(x) reshape(analyzer(x),[],1);
procG11 = @(x) x(:);
procG12_wo_syn = @(x) (eta_prop/lambda_prop)*reshape(dltzproc(x),[],1);
procG2_wo_syn = @(x) x(:);
%
%procG1 = @(x) [ procG11(x) ; procG12_wo_syn(synthesizer(x)) ]; 
%procG2 = @(x) reshape(synthesizer(procG2_wo_syn(x)),[],1); 
%
adjpG11 = @(x) x(:);
adjpG12_wo_ana = @(x) reshape((eta_prop/lambda_prop)*dltzadjp(x),[],1);
adjpG2_wo_ana = @(x) x(:);
%
%adjpG1 = @(x) adjpG11(x(1:nCoefs)) ... 
%    + reshape(analyzer(adjpG12_wo_ana(x(nCoefs+1:end))),[],1); 
%adjpG2 = @(x) reshape(analyzer(adjpG2_wo_ana(x)),[],1);
%
z1p = zeros(nCoefs+nSamples,1);
z2p = zeros(nSamples,1);
d1p = zeros(nCoefs+nSamples,1);
d2p = zeros(nSamples,1);

%% ADMM steps for proposal
invpQ_ = @(x) fcn_invpQ_(x,diagR_,rho_admm,{analyzer, synthesizer},nSamples);
for itr=1:nIters_admm
    % Step 1
    %t_ = y_ + rho_admm*(adjpG1(z1p-d1p) + adjpG2(z2p-d2p));
    z11p = z1p(1:nCoefs);
    d11p = d1p(1:nCoefs);
    z12p = z1p(nCoefs+1:end);
    d12p = d1p(nCoefs+1:end);
    adjrev = adjpG11(z11p-d11p) ...
        + reshape(analyzer(adjpG12_wo_ana(z12p-d12p)+adjpG2_wo_ana(z2p-d2p)),[],1);
    t_ = y_ + rho_admm*(adjrev);
    xc = invpQ_(t_);
    % Step 2
    %procG1xc = procG1(xc);
    %procG2xc = procG2(xc);
    sc = synthesizer(xc);
    procG1xc = [ procG11(xc); procG12_wo_syn(sc) ];
    procG2xc = procG2_wo_syn(sc);
    z1c = softthresh_(procG1xc + d1p, lambda_prop/rho_admm); %eta_prop);
    z2c = metricproj_(procG2xc + d2p, ri_lb, ri_ub);
    % Step 3
    d1c = d1p + procG1xc - z1c;
    d2c = d2p + procG2xc - z2c;
    % Update
    z1p = z1c;
    z2p = z2c;
    d1p = d1c;
    d2p = d2c;
end
qc = synthesizer(xc);
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

%
function y = fcn_procQ_(x,diagH,gamma,dictionaries,nSamples)
analyzer = dictionaries{1};
synthesizer = dictionaries{2};
% Definition of process by matrix Q
s_ = synthesizer(x);
y = gamma*x(:) + reshape(analyzer(ifft(diagH.*reshape(fft(s_,nSamples),[],1),nSamples)+gamma*s_(:)),[],1);
end

%
function y = fcn_invpQ_(x,diagR,gamma,dictionaries,nSamples)
analyzer = dictionaries{1};
synthesizer = dictionaries{2};
% Definition of process by matrix inv(Q)
s_ = synthesizer(x);
y = x(:)/gamma - (1/gamma)^2*reshape(analyzer((1/4)*ifft(diagR.*reshape(fft(s_,nSamples),[],1),nSamples)+(gamma/2)*s_(:)),[],1);
end