function mymlx2m(varargin)
%% Live script to MATLAB script conversion
%
% Copyright (c) Shogo MURAMATSU, 2020
% All rights resereved

srcDir = fullfile(pwd,'./');
dstDir = './scripts/';
isVerbose = true;

%% Fetch files
if nargin < 1
    srcFiles = [srcDir '*.mlx'];
else
    srcFiles = sprintf('%s/%s*.mlx',srcDir,varargin{1});
end
if isunix
    list = dir(srcFiles);
else
    list = ls(srcFiles);
end

%% File conversion
for idx = 1:size(list,1)
    % Extract filenames
    if isunix
        [~,fname,~] = fileparts(list(idx).name);
    else
        [~,fname,~] = fileparts(list(idx,:));
    end
    % Convert to script
    support.fcn_mlx2m(srcDir,fname,dstDir,isVerbose)
    % Contents of converted script
    %open(fullfile(dstDir,[fname '.m']))
end