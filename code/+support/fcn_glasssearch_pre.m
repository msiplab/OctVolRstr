function [rhat,arhat,xhat] = fcn_glasssearch_pre(theta,vg,glsref,isEnv,isBackProj)
%FCN_GLASSSERCH_PRE Search surface by polyfit
%   Function for searachng glass surfaces
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
isVisible = false;
if nargin < 5 
    isBackProj = false;
end
if nargin < 4 || isempty(isEnv)
    isEnv = true;
end
if nargin < 3 || isempty(glsref)
    glsref = 0.5;
end
[szy,szx,szz] = size(vg);
sigmaz = theta(3);

% Find peaks in Z direction for each X-Y position
K = 2;
pz1 = zeros(szy,szx);
pz2 = zeros(szy,szx);
if isBackProj % Back projection
    vgg = support.fcn_coherence3d(vg,theta,'Adjoint');
else
    vgg = vg;
end
for ix = 1:szx
    for iy = 1:szy
        zline = squeeze(vgg(iy,ix,:));
        if isEnv % Take envelope
            zline = envelope(zline);
        end
        [~,idzs] = findpeaks(abs(zline),...
            'NPeaks',K,...
            'SortStr','descend',...
            'MinPeakDistance',sigmaz);    
        idzs = sort(idzs);
        pz1(iy,ix) = idzs(1);
        pz2(iy,ix) = idzs(2);
    end
end
if isVisible
    figure
    zline = squeeze(vgg(end/2,end/2,:));
    if isEnv
        zline = envelope(zline);
    end
    findpeaks(zline,...
        'NPeaks',K,'SortStr','descend','MinPeakDistance',sigmaz)
    drawnow
end
[px,py] = meshgrid((1:szx)-floor(szx/2),(1:szy)-floor(szy/2));
sf1 = fit([px(:),py(:)],pz1(:),'poly11');
sf2 = fit([px(:),py(:)],pz2(:),'poly11');
if isVisible
    figure
    plot(sf1,[px(:),py(:)],pz1(:))
    drawnow
    figure
    plot(sf2,[px(:),py(:)],pz2(:))
    drawnow
end
xhat = zeros(szz,1);
xhat(round(sf1.p00)) = glsref;
xhat(round(sf2.p00)) = -glsref;
arhat(1,:) = [sf1.p01 sf1.p10]; % [ary arx]
arhat(2,:) = [sf2.p01 sf2.p10];
% Estimation of r
rhat = support.fcn_glasssubstrate(max(xhat,0),[szy szx szz],arhat(1,:));
rhat = rhat + support.fcn_glasssubstrate(min(xhat,0),[szy szx szz],arhat(2,:));