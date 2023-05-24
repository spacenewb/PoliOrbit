function kep = car2kep(y,mu)

% car2kep   Returns Keplerian orbital elements given the s/c position and velocity
%           vectors in cartesian coordinates and the planetary constant mu (Function Vectorised)
%
% PROTOTYPE
%   kep = car2kep(y,mu)
%
% INPUT:
%   y[6xn]  State of the body ( rx , ry , rz , vx , vy , vz ) [ L, L/T ]
%   mu[1]   Gravitational parameter of the primary [ L^3/T^2 ]
%
% OUTPUT:
%       kep[6xn]        Keplerian elements
%                       kep = [a; e; i; Om; om; theta] [km, rad]
%
%       a       Semi-major axis [km]
%       e       Eccentricity [-]
%       i       Inclination [rad]
%       Om      RAAN (Right Ascension of the Ascending Node) [rad]
%       om      Periapsis true anomaly [rad]
%       theta   True anomaly [rad]
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

r_vec = y(1:3,:);
v_vec = y(4:6,:);

r = vecnorm(r_vec);
v = vecnorm(v_vec);
k = [0; 0; 1];

a = -mu./(v.^2 - 2.*mu./r);

h_vec = cross(r_vec,v_vec,1);
h = vecnorm(h_vec);

e_vec = cross(v_vec,h_vec,1)./mu - r_vec./r;
e = vecnorm(e_vec);

i = acos(h_vec(3)./h);

N = cross(repmat(k,1,length(r)),h_vec,1)./vecnorm(cross(h_vec,repmat(k,1,length(r)),1));

Om = acos(N(1,:));
om = acos(dot(N,e_vec,1)./e);
theta = acos(dot(r_vec,e_vec,1)./(r.*e));

if length(r)>1
    Om = (N(2,:) < 0).*(2.*pi - Om) + (N(2,:) >= 0).*Om;
    om = (e_vec(3,:) < 0).*(2.*pi - om) + (e_vec(3,:) >= 0).*om;
    theta = (dot(r_vec,v_vec,1) < 0).*(2.*pi - theta) + (dot(r_vec,v_vec,1) >= 0).*theta;
else
    
    if N(2) < 0
        Om = 2*pi - Om;
    end

    if e_vec(3) < 0
        om = 2*pi - om;
    end

    if dot(r_vec,v_vec,1) < 0
        theta = 2*pi - theta;
    end

end

kep = real([a; e; i; Om; om; theta]);

end

