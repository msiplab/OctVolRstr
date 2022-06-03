function y = fcn_resampling(y,tx,worker,method)
% RESAMPLING:
%
% Requirements: MATLAB R2022a
%
% Copyright (c) 2021-2022, Ruiki KOBAYASHI and Shogo MURAMATSU
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
if nargin < 4
    method = 'spline';
end
[obsx, obsy , ~] = size(y);
%if worker > 0
parfor(i = 1:obsx,worker)
    for j = 1:obsy
        y(i,j,:) = resample(squeeze(y(i,j,:)),tx,method);
    end
end
%{
else
    for i = 1:obsx
        for j = 1:obsy
            y(i,j,:) = resample(squeeze(y(i,j,:)),func);
        end
    end
end
%}
end