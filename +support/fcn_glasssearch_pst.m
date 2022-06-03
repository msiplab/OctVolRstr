function [rhat,arhat,xhat] = fcn_glasssearch_pst(ar,x,theta,vg,objfun,dev,niters)
%FCN_GLASSSERCH_PST Exhaustive search for glass surface
%   Function for searaching glass surfaces exhaustively
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
if nargin < 7
    niters = 1;
end
if nargin < 6
    dev = 2;
end
[szy,szx,szz] = size(vg);
%
xhat = x;
arhat = ar;
% Compensation
maxErr = Inf;
for iter=1:niters
    [indcs,~,vals] = find(xhat);
    %
    for dz1st = -dev:dev
        idz1st = indcs(1);
        for dz2nd = -dev:dev
            idz2nd = indcs(2);
            % 1st surface
            x1 = zeros(size(x));
            x1(idz1st+dz1st) = vals(1);
            r_ = support.fcn_glasssubstrate(x1,[szy,szx,szz],arhat(1,:));
            % 2nd surface
            x2 = zeros(size(x));
            x2(idz2nd+dz2nd) = vals(2);
            r_ = r_ + support.fcn_glasssubstrate(x2,[szy,szx,szz],arhat(2,:));
            err = objfun(theta,r_,vg);
            %
            %rc_ = r_(round(end/4):round(end*(3/4)),round(end/4):round(end*(3/4)),:);
            %vc_ = vg(round(end/4):round(end*(3/4)),round(end/4):round(end*(3/4)),:);
            %err = objfun(theta,rc_,vc_);
            if err < maxErr
                maxErr = err;
                x1hat = x1;
                x2hat = x2;
                xhat = x1 + x2;
                rhat = r_;
            end
        end
    end
    %{%
    arhatpre = arhat;
    for rar1st = 0.95:0.0125:1.05
        ar1st = rar1st*arhatpre(1,:);
        for rar2nd = 0.95:0.0125:1.05
            ar2nd = rar2nd*arhatpre(2,:);
            % 1st surface
            r_ = support.fcn_glasssubstrate(x1hat,[szy,szx,szz],ar1st);
            % 2nd surface
            r_ = r_ + support.fcn_glasssubstrate(x2hat,[szy,szx,szz],ar2nd);
            %
            err = objfun(theta,r_,vg);
            if err < maxErr
                maxErr = err;
                arhat(1,:) = ar1st;
                arhat(2,:) = ar2nd;
                rhat = r_;
            end
        end
        %}%
    end
end
end

%%
%{
function y = dicFwd(x,obsSize,theta,ar)
y = support.fcn_coherence3d(...
    support.fcn_glasssubstrate(x,obsSize,ar,'Forward'),...
    theta,'Forward');
end

function y = dicAdj(x,obsSize,theta,ar)
y = support.fcn_glasssubstrate(...
    support.fcn_coherence3d(x,theta,'Adjoint'),...
    obsSize, ar,'Adjoint');
end

function y = hardthresh(x,glsref)
[~,idx1] = max(x);
[~,idx2] = min(x);
y = zeros(size(x));
y(idx1) = glsref;
y(idx2) = -glsref;
end

function isadjcheck()
% Check if adjoint
arr = randn(2,1);
% <y,Dx> = <D.'y,x> % x0
v = dicFwd(xhat,obsSize,theta,arr); % Dx0
% <vg,Dx0> == <D.'vg,x0>
x = dicAdj(vg,obsSize,theta,arr);
%
a = sum(vg(:).*v(:),'all');
b = sum(x(:).*xhat(:),'all');
err = abs(a-b);
assert(err<1e-9,'%g: Adjoint relation is violated!',err);
end
%}