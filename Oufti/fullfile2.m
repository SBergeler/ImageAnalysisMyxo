% This function was obtained from the Oufti software package from the
% Jacobs-Wagner lab, Yale University. The package is published under the
% GNU General Public License as published by the Free Software Foundation
% (version 3).
%
% Copyright (c) 2015, the Jacobs-Wagner lab, Yale University. 

function res = fullfile2(varargin)
    % This function replaces standard fullfile function in order to correct
    % a MATLAB bug that appears under Mac OS X
    % It produces results identical to fullfile under any other OS
    arg = '';
    for i=1:length(varargin)
        if ~strcmp(varargin{i},'\') && ~strcmp(varargin{i},'/')
            if i>1, arg = [arg ',']; end
            arg = [arg '''' varargin{i} ''''];
        end
    end
    eval(['res = fullfile(' arg ');']);
end