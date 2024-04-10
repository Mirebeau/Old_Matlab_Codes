% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

% getoptions - retrieve options with optional default value v0.
% Adapted from: 2007 Gabriel Peyre.

function v = GetOptions(options, name, v0)
if nargin >3
    error('Too many arguments');
elseif nargin<2
    error('Not enough arguments.');
elseif isfield(options, name)
    v = options.(name);
elseif nargin==3
    v=v0;  
    if isfield(options,'verbose') && options.verbose
        disp(['Using default value for options.' name ' : ']);
        disp(v)
    end
else
    error(['You have to provide options.' name '.']);
end 