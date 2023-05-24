function [ r, v ] = uplanet_rv( time_mjd2000, planetNum )
% Calculate analytical ephemeris for planets
% Function is the same as uplanet.m, but output is already in cartesian coords
%
% PROTOTYPE:
%   [ r, v ] = uplanet_rv( time_mjd2000, planetNum )
%
% INPUT:
%   time_mjd2000[1]     MJD2000 time
%   planetNum[1]        ID of planet for ephemeris
%                           1:   Mercury
%                           2:   Venus
%                           3:   Earth
%                           4:   Mars
%                           5:   Jupiter
%                           6:   Saturn
%                           7:   Uranus
%                           8:   Neptune
%                           9:   Pluto
%                           10:  Sun
% OUTPUT:
%   r[3]                Position vector of planet in HECI frame (heliocentric)
%   v[3]                Velocity of planet in HECI frame (heliocentric)

[kep, mu_sun] = uplanet ( time_mjd2000, planetNum );

[ r, v ] = kp2rv( kep(1), kep(2), kep(3), kep(4), kep(5), kep(6), mu_sun );

end