function r = fcn_glasssubstrate(s,sz,ar,mode)
%FCN_GLASSSUBSTRATE Simulate glass substrate
%   Function for generating artifical glass substrate
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
if nargin < 4
    mode = 'Forward';
end

if strcmp(mode,'Forward')
    % Zero-order hold in XY
    u_ = zohxy_(s,[sz(1) sz(2)]);
    % Shift in Z
    r = shiftz_(u_,ar);
%elseif strcmp(mode,'Adjoint')
    % Shift in Z
%    u_ = shiftz_(s,-ar);
    % Adjoint of zero-order hold in XY
%    r = sum(sum(u_,1),2);
else
    error('Invalid mode')
end

%{
[inz,~,snz] = find(s);
thick = inz(2)-inz(1);
r_ = zeros(sz);
t0x = 1:sz(2);
t0y = (1:sz(1)).';
xind = ones(sz(1))*round(sz(3)/4);
fi = polyval([ar(1) 0],t0x-1);
fj = polyval([ar(2) 0],t0y-1);
for i = 1:sz(2)
    xind(i,:) = xind(i,:)+fi;
    xind(:,i) = xind(:,i)+fj;
end
xtind = round(xind);
for i = 1:sz(2)
    for j = 1:sz(1)
        r_(i,j,xtind(i,j)) = 0.5;
        r_(i,j,xtind(i,j)+thick) = -0.5;
    end
end
err = norm(r_(:)-r(:))
assert(err<1e-6)

end
%}
end
function y = zohxy_(s, ufactor)
y = repmat(permute(s(:),[2 3 1]),[ufactor 1]);
end

function y = shiftz_(u, ar)
[szy,szx,~] = size(u);

y = u;
for ix = 1:szx
    px = ix-floor(szx/2);
    for iy = 1:szy
        py = iy-floor(szy/2);
        iz = round(ar(2)*px+ ar(1)*py);
        y(iy,ix,:) = circshift(y(iy,ix,:),[0 0 iz]);
    end
end
end

%{
function y = impres_(sz, inz_, val_, ar)
y = zeros(sz);
for ix = 1:sz(2)
    for iy = 1:sz(1)
        iz = mod(round(inz_ + ar(2)*(ix-1) + ar(1)*(iy-1))-1,sz(3))+1;
        y(iy,ix,iz) = val_;
    end
end
end
%}