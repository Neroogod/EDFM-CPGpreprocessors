function root = webRoot()
% Repository root, derived from this file's own location.
% web/matlab/webRoot.m  ->  two levels up.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
