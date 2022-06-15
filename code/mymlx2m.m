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
    list = ls([srcDir '*.mlx']);
else
    list = ls(sprintf('%s/%s*.mlx',srcDir,varargin{1}));    
end

%% File conversion
for idx = 1:size(list,1)
    % Extract filenames
    [~,fname,~] = fileparts(list(idx,:));
    % Convert to script
    support.fcn_mlx2m(srcDir,fname,dstDir,isVerbose)
    % Contents of converted script
    %open(fullfile(dstDir,[fname '.m']))
end