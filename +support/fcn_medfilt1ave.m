function y = fcn_medfilt1ave(x,n)
%FCN_STAGE1MEDIAN Stage 1 median filtering
%   
% Reference:
%
% - M. Shamouilian and I. Selesnick, "Total Variation Denoising for Optical
% Coherence Tomography," 2019 IEEE Signal Processing in Medicine and Biology
% Symposium (SPMB), 2019, pp. 1-5, doi: 10.1109/SPMB47826.2019.9037832.
%
if nargin < 2
    n=3;
end
d = ndims(x);
y = 0;
for idim=1:d
    y = y + medfilt1(x,n,[],idim);
end
y = y/d;

