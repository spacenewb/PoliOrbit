function y = kep2car(kep,mu) 

% kep2car   Returns cartesian position and velocity given the s/c keplerian
%           elements and the planetary constant mu (Function Vectorised)
%
% PROTOTYPE
%   y = kep2car(kep,mu)
%
% INPUT:
%       kep[6xn]        keplerian elements
%                       kep = [a; e; i; Om; om; theta] [km, rad]
%       mu[1]           Gravitational costants of the gravitational body [km^3/s^2]
% 
% OUTPUT: 
%       y[6xn]          State of the body ( rx , ry , rz , vx , vy , vz ) [ L, L/T ]
%
% CALLED FUNCTIONS:
%   [N.A]
%
% CONTRIBUTORS:
%   Hariharan Venkatesh Vitaladevuni
%
% VERSIONS:
%   2021-10-03: First version
%
%--------------------------------------------------------------------------

a=kep(1,:); 
e=kep(2,:);
i=kep(3,:); 
Om=kep(4,:); 
om=kep(5,:); 
theta=kep(6,:);

p=a.*(1-e.^2);
r=p./(1+e.*cos(theta)); 
p = (p<0).*(-p) + (p>=0).*p;

% Calculate R --> Perifocal
rx = r.*cos(theta); 
ry = r.*sin(theta); 
rz = zeros(1,length(i)); 
rv = [rx;ry;rz];

% Calculate V --> Perifocal
vx = -sqrt(mu./p).*sin(theta);
vy = sqrt(mu./p).*(e+cos(theta));
vz = zeros(1,length(i)); 
vv = [vx;vy;vz]; 

% Perifocal --> Celestial
r_vec = zeros(size(rv));
v_vec = zeros(size(rv));
for idx=1:length(i)
    R313=rotz(rad2deg(Om(idx)))*rotx(rad2deg(i(idx)))*rotz(rad2deg(om(idx)));
    r_vec(:,idx)=R313*rv(:,idx); 
    v_vec(:,idx)=R313*vv(:,idx);
end

y = [r_vec;v_vec];

end