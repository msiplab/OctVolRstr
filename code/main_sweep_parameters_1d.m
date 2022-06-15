%% Local functions
disp("Dummy script for consistency with old files")
%---
%%
function y = fcn_phi1d_(u,mode,ribound)
%FCN_PHI1D Function to convert refractive index to reflection
%
dz = [1 0 -1];
az = [1 0  1];
dzu = circshift(cconv(dz,u,length(u)),-1);
azu = circshift(cconv(az,u,length(u)),-1);
if nargin < 3
    ribound = [1.0 1.5];
end
if strcmp(mode,'reflection')
    y = abs(dzu).*(-dzu)./(azu.^2);
elseif strcmp(mode,'linear')
    b = ribound(2);
    a = ribound(1);
    y = -beta1(a,b)*dzu;
else
    error('Not supported')
end
end

%---
function [measureproc,measureadjp] = fcn_gencoh1d_(nSamples,params)
%FCN_GENCOH1D Generation of measurement process
alpha_p = params.alpha_p; % Amplitude
omega_p = params.omega_p; % Frequency of fluctuation
sigma_z = params.sigma_z; % Standard deviation of Gaussian function
b_p = params.b_p; % Broadening parameter

% Impluse response of measurement process
p = @(m) alpha_p.*exp(-m.^2/(2*sigma_z.^2)).*cos(omega_p*m).*sinc(b_p*m);
% Impulse response of P
hlen = round(3*sigma_z);
zsmpl = -hlen:hlen;
irp = p(zsmpl); % impluse response of p

% Measurement process
measureproc = @(x) circshift(cconv(irp,x,nSamples),-floor(length(irp)/2));
measureadjp = @(x) circshift(cconv(irp(end:-1:1),x,nSamples),mod(length(irp)-1,2)-floor(length(irp)/2));

% Check if <Px,y> = <x,P.'y> (is adjoint?)
fcn_isinadjrel_(measureproc,measureadjp,nSamples)
end

%---
function fcn_isinadjrel_(fwdp,adjp,nSamples1,nSamples2)
%FCN_ISINADJREL Function to check the adjoint relationship
%
if nargin < 4
    nSamples2 = nSamples1;
end
x1 = randn(nSamples1,1);
y1 = fwdp(x1);
y2 = randn(nSamples2,1);
x2 = adjp(y2);
err = abs(x1(:).'*x2(:)-y1(:).'*y2(:));
assert(err<1e-9,'Adjoint relation is not met. (abs err: %f)',err)
end

%---
function u = fcn_rigen1d_(nSamples,vrange)
%FCN_RIGEN1D Generation of refractive index distribution
% Generation of a refractive index distribution for simulation
u2d = phantom('Modified Shepp-Logan',nSamples);
minu2d = min(u2d(:));
maxu2d = max(u2d(:));
ri_lb = vrange(1);
ri_ub = vrange(2);
u2d = (ri_ub-ri_lb)*(u2d-minu2d)/(maxu2d-minu2d)+ri_lb; % Normalize to [a,b]
% Extraction of 1-D sequence
u = u2d(floor(end/2),:);
end

%---
function plotur_(u, r, nSamples)
yyaxis right
plot(u)
ylabel('Refractive index')
xlabel('z')
axis([1 nSamples 0 2])
yyaxis left
plot(r)
ylabel('Reflection')
axis([1 nSamples -0.1 0.1])
end

%---
function ploturv_(u, r, v, nSamples)
subplot(3,1,1)
plot(u)
title('Refractive index {\bf u}')
xlabel('z')
axis([1 nSamples 0 2])
subplot(3,1,2)
plot(r)
title('Refletion {\bf r}')
xlabel('z')
axis([1 nSamples -0.1 0.1])
subplot(3,1,3)
plot(v)
title('Observation {\bf v}')
xlabel('z')
axis([1 nSamples -0.2 0.2])
hold off
end

%---
% Undecimated Haar decomposition
function y = udhaardec_(x,nLevels)
    swc = swtscale_(nLevels)*swt(x,nLevels,'haar');
    y = swc(:);
end

% Undecimated Haar reconstruction
function x = udhaarrec_(y,nLevels)
    swc = swtscale_(nLevels)\reshape(y,nLevels+1,[]);
    x = iswt(swc,'haar');
end

% Scaling factor for Parseval tight property 
function s = swtscale_(nLevels)
    s = diag([2.^(-(1:nLevels)/2) 2^(-nLevels/2)]);
end