function [alpha, delta, sc] = car2AlfaDelta(y)
 
% car2AlfaDelta     This function returns Rigth Ascension and declintion given the state
%                   vector y (matrix containing position and velocity
%                   vectors) (function vectorised)
%
% PROTOTYPE:
%   [alpha, delta, sc] = car2AlfaDelta(y)
% 
% INPUT:
%   Y[6x1] State of the body ( rx , ry , rz , vx , vy , vz ) [ L, L/T ]
%                       OR
%   R[3x1] Position of the body ( rx , ry , rz ) [ L ]
% 
% OUTPUT:
%   alpha:    Rigth Acension (celestial longitude) [rad]
%   delta:    Declination (celestial latitude)  [rad]
%   sc = [l; m; n]:    spherical coordinates 
% 
% CALLED FUNCTIONS:
%   [N.A]
%
% CONTRIBUTORS:
%   Hariharan Venkatesh Vitaladevuni
%
% VERSIONS:
%   2021-11-13: First version
%
%--------------------------------------------------------------------------

r = y(1:3,:);
r_vec = vecnorm(r);

l = r(1,:)./r_vec;     
m = r(2,:)./r_vec;    
n = r(3,:)./r_vec;
sc = [l;m;n];

delta = asin(n);

alpha = acos(l./cos(delta));
alpha = (m<=0).*(2*pi-alpha) + (m>0).*alpha;

end

