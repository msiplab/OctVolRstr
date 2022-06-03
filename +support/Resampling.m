function y = Resampling(y,func,worker)
% RESAMPLING:
%
% Requirements: MATLAB R2020b
%
% Copyright (c) 2021, Ruiki KOBAYASHI
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
[obsx, obsy , ~] = size(y);
if worker > 0
    parfor(i = 1:obsx,worker)
        for j = 1:obsy
            y(i,j,:) = resample(squeeze(y(i,j,:)),func);
        end
    end
else
    for i = 1:obsx
        for j = 1:obsy
            y(i,j,:) = resample(squeeze(y(i,j,:)),func);
        end
    end
end
end