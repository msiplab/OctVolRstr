function y = fcn_coherence3d(x,theta,mode)
% COHERENCD3D: 3-D coherence function
%
% Requirements: MATLAB R2022a
%
% Copyright (c) 2022, Ruiki KOBAYASHI and Shogo MURAMATSU
%
% All rights reserved.
%
% Contact address: Shogo MURAMATSU,
%    Faculty of Engineering, Niigata University,
%    8050 2-no-cho Ikarashi, Nishi-ku,
%    Niigata, 950-2181, JAPAN
%
% http://msiplab.eng.niigata-u.ac.jp/
%
if nargin < 3
    mode = 'Forward';
end

% Parameters
Scale     = theta(1);
Sigmaxy   = theta(2);
Sigmaz    = theta(3);
Frequency = theta(4);
B4Sinc    = theta(5);

len = [round(Sigmaxy*12)+1, 2*round(Sigmaz/(B4Sinc)^0.3/2)+1 ];
lenz = len(2); %2*round(Sigmaz/(B4Sinc)^0.3/2)+1; % Set an odd number
lenxy = len(1); %round(Sigmaxy*12)+1
%
nz = -floor(lenz/2):floor(lenz/2);
nx = -floor(lenxy/2):floor(lenxy/2);
ny = -floor(lenxy/2):floor(lenxy/2);
%
gcx = exp(-nx.^2./(2*Sigmaxy.^2));
gcy = exp(-ny.^2./(2*Sigmaxy.^2));
gcz = Scale ...
    *exp(-nz.^2./(2*Sigmaz.^2)).*cos(2*pi*Frequency*nz).*sinc(nz*B4Sinc);
gck = permute(gcz(:),[2 3 1]);
%
KernelZ = gck;
KernelX = gcx(:).';
KernelY = gcy(:);

if strcmp(mode,'Forward')
    v = imfilter(x,KernelY,'conv','circ');
    w = imfilter(v,KernelX,'conv','circ');
    y = imfilter(w,KernelZ,'conv','circ');
elseif strcmp(mode,'Adjoint')
    v = imfilter(x,KernelY,'corr','circ');
    w = imfilter(v,KernelX,'corr','circ');
    y = imfilter(w,KernelZ,'corr','circ');
elseif strcmp(mode,'Sinc')
    gc = Scale*exp(-nz.^2./(2*Sigmaz.^2)).*cos(2*pi*Frequency*nz).* ...
        (cos(pi*nz*B4Sinc)/B4Sinc - sinc(nz*B4Sinc)/B4Sinc);
    gck = permute(gc(:),[2 3 1]);
    SinckernelZ = gck;
    %
    v = imfilter(x,KernelY,'corr','circ');
    w = imfilter(v,KernelX,'corr','circ');
    y = imfilter(w,SinckernelZ,'corr','circ');
elseif strcmp(mode,'Freq')
    gc = Scale*exp(-nz.^2./(2*Sigmaz.^2)).*sinc(nz*B4Sinc)...
        .*(-2*pi*nz.*sin(2*pi*Frequency*nz));
    gck = permute(gc(:),[2 3 1]);
    FreqkernelZ = gck;
    %
    v = imfilter(x,KernelY,'corr','circ');
    w = imfilter(v,KernelX,'corr','circ');
    y = imfilter(w,FreqkernelZ,'corr','circ');
elseif strcmp(mode,'Sigmaz')
    gc = Scale*cos(2*pi*Frequency*nz).*sinc(nz*B4Sinc)...
        .*(nz.^2./Sigmaz.^3).*exp(-nz.^2./(2*Sigmaz.^2));
    gck = permute(gc(:),[2 3 1]);
    SigmazkernelZ = gck;
    %
    v = imfilter(x,KernelY,'corr','circ');
    w = imfilter(v,KernelX,'corr','circ');
    y = imfilter(w,SigmazkernelZ,'corr','circ');
elseif strcmp(mode,'Sigmaxy')
    gc = Scale...
        *exp(-nz.^2./(2*Sigmaz.^2)).*cos(2*pi*Frequency*nz).*sinc(nz*B4Sinc); 
    gck = permute(gc(:),[2 3 1]);
    nxx = nx(:).';
    nyy = ny(:);
    gcDiffxy = ((nxx.^2+nyy.^2)./(Sigmaxy.^3)).*exp(-(nxx.^2+nyy.^2)./(2*Sigmaxy.^2));
    SigmaxykernelXY = gcDiffxy;
    SigmaxykernelZ = gck;
    %
    v = imfilter(x,SigmaxykernelXY,'corr','circ');
    y = imfilter(v,SigmaxykernelZ,'corr','circ');
elseif strcmp(mode,'Scale')
    gc = exp(-nz.^2./(2*Sigmaz.^2)).*cos(2*pi*Frequency*nz).*sinc(nz*B4Sinc);
    gck = permute(gc(:),[2 3 1]);
    ScalekernelZ = gck;
    %
    v = imfilter(x,KernelY,'corr','circ');
    w = imfilter(v,KernelX,'corr','circ');
    y = imfilter(w,ScalekernelZ,'corr','circ');
else
    error('Invalid mode')
end
end
