% TimeStep for balls with friction.

% Input : 
%  - p : ball positions
%  - v : ball velocities
%  - a : ball angular velocities
%  - options : 
%       - r : radius (default : 1)
%       - rTol : penetration tolerance
%       - maxMove : max motion in a timestep (default : r)
%       - maxTau : max timestep (default : +infty)
%       - alpha : friction coefficient (Ball/Ball, and Ball/Obs for each obstacle)
%       - obs : TODO obstacles (polygons ? Fixed balls ?)

% Output : 
% - p : advected positions
% - v : projected velocities
% - a : projected angular velocities
% - tau : timestep

% Notes :
% - angular position not stored, only angular velocity.
% - You may want to apply external forces v->v+tau g to the output.

% - Make model should be solved ... ?


% Copyright ...

function [p,v,a,tau] = TimeStep(p,v,a,options)

Dimension = size(p,1);
nPoints = size(p,2);
ASymDimension = Dimension*(Dimension-1)/2;
assert(all(size(p)==[Dimension,nPoints]));
assert(all(size(v)==[Dimension,nPoints]));
assert(all(size(a)==[ASymDimension,nPoints]));

r = GetOptions(options,'r',1);
rTol = GetOptions(options,'rTol',r/10);    
maxMove = GetOptions(options,'maxMove',r);
maxTau = GetOptions(options,'maxTau',inf);
alpha = GetOptions(options,'alpha');

% Get all pairs of balls which interact

% Get all pais of balls/obstacle which interact

% Project velocities onto admissible ones

% Get all pairs of balls which may interact

% Adjust timestep to meet tolerance

% Advect
p = p+tau*v;

end