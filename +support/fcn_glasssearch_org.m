function [xEstIndex,thickEst] = fcn_glasssearch_org(vg)
%FCN_GLASSSEARCH Search surface by polyfit
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
[obsy,obsx,~] = size(vg);

% 
thickEst = 0;
indexEst = zeros(obsy,obsx);
for i = 1:obsx
    for j = 1:obsy
        [x,y] = support.fcn_sparseone(vg(i,j,:));
        thickEst = thickEst+abs(y-x);
        indexEst(i,j) = x;
    end
end
thickEst = round(thickEst/(obsx*obsy));

t = 1:obsy;
xEstIndex = zeros(obsx,obsy);
px = [0 0];
py = [0 0];
for i = 1:obsx
    px = px + polyfit(t,indexEst(:,i),1);
end
for j = 1:obsy
    py = py + polyfit(t,indexEst(j,:),1);
end
px = px/obsx;
py = py/obsy;
for i = 1:obsx
    fx = polyval([px(1) px(2)/2-px(1)*obsx/4],t);
    fy = polyval([py(1) py(2)/2-py(1)*obsy/4],t);
    xEstIndex(i,:) = xEstIndex(i,:)+fy;
    xEstIndex(:,i) = xEstIndex(:,i)+fx(:);
end
xEstIndex = round(xEstIndex);
