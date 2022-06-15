function [i,j] = fcn_sparseone(y)
% SPARSEONE:
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
y(abs(y/max(y,[],'all'))<0.3) = 0; % hard thresholding
[~,i] = max(abs(y));
k = find(y(max(i-20,1):min(i+20,size(y,3))));
i = k(round(end/2))+max(i-20,1)-1;
y(max(i-20,1):i+20) = 0;
[~,j] = max(abs(y));
k = find(y(max(j-20,1):min(j+20,size(y,3))));
j = k(round(end/2))+max(j-20,1)-1;
if i > j
    temp = j;
    j = i;
    i = temp;
end
end
