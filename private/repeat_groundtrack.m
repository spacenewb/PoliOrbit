function a_rep = repeat_groundtrack(k, m, planetID)

% groundtrack     This function returns the groundtrack plot and values of topocentric coordinates
%
% PROTOTYPE:
%   [lon, lat] = groundtrack(y ,greenwich0, t_span, omega_planet, planetID)
% 
% INPUT:
%   y[6x1]          State of the body ( rx , ry , rz , vx , vy , vz ) [ L, L/T ]
%                                   OR
%   R[3x1]          Position of the body ( rx , ry , rz ) [ L ]
%   greenwich0      Longitude of Greenwich meridian [ rad ]
%   t_span:         Time span for ground track [ s ]
%   omega_Planet:   Rotational speed of the planet [ rad/s ]
%   planetID:       Integer number identifying the celestial body (< 10)
%                       1:   Mercury
%                       2:   Venus
%                       3:   Earth
%                       4:   Mars
%                       5:   Jupiter
%                       6:   Saturn
%                       7:   Uranus
%                       8:   Neptune
%                       9:   Pluto
%                       10:  Sun
%
% OUTPUT:
%   lon:            longitude [deg]
%   lat:            latitude [deg]
% 
% CALLED FUNCTIONS:
%   car2AlfaDelta 
%
% CONTRIBUTORS:
%   Hariharan Venkatesh Vitaladevuni
%
% VERSIONS:
%   2021-11-13: First version
%
%--------------------------------------------------------------------------

%% Calculating Ground Track

if planetID == 3
    omega_planet = 2*pi/(23*3600+56*60+4.09); %angular velocity of Earth
end

mu = astroConstants(10+planetID);

n = omega_planet*k/m;
a_rep = (mu/n^2)^double(1/3);

end