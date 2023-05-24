function [ alpha, delta, lon, lat, plotHandle ] = groundTrack( a0, e0, i0, OM0, om0, theta0,...
                                                               t0, RA0_G, orbits,...
                                                               perturbationModel, relTol,...
                                                               absTol, points, figureNum,...
                                                               linecolor, fontsize,...
                                                               markersizeEndStart, markersizeLine )
% groundTrack plots groundTrack of satellite on GEOCENTRIC orbit
%
% PROTOTYPE:
%     [ alpha, delta, lon, lat, plotHandle ] = groundTrack( a, e, i, OM, om, theta0,...
%                                                                    t0, RA0_G, orbits,...
%                                                                    perturbationModel, relTol,...
%                                                                    absTol, points, figureNum,...
%                                                                    linecolor, fontsize,...
%                                                                    markersizeEndStart, markersizeLine )
%
% INPUT:
%   a[1]                    [km]        Semi-major axis
%   e[1]                    [-]         Eccentricity
%   i[1]                    [rad]       Inclination
%   OM[1]                   [rad]       RAAN
%   om[1]                   [rad]       Periapsis anomaly
%   theta0[1]               [rad]       Initial true anomaly on orbit
%   t0[1]                   [s]         Initial reference time on orbit
%   RA0_G[1]                [rad]       Initial right ascension of Greenwhich meridian with respect to earth centered inertial equatorial frame
%   orbits[1]               [-]         Number of orbits to represent (can be decimal)
%   perturbationModel[1]    [-]         'unperturbed' or 'J2-perturbed'
%   relTol[1]               [-]         Relative tolerance of integration for orbit propagation
%   absTol[1]               [-]         Absolute tolerance of integration for orbit propagation
%   points[1]               [-]         Number of points to represent on the plot
%   figureNum[1]            [-]         'new' for new figure or integer for pre-existing figure
%   linecolor[1]            [-]         RGB color of plot line
%   fontsize[1]             [-]         Fontsize in plots
%   markersizeEndStart[1]   [-]         Markersize of starting and ending points
%   markersizeLine[1]       [-]         Markersize of points on orbit
%
% OUTPUT:
%   alpha[1xtime]           [rad]       Right ascension in earth centered inertial equatorial frame
%   delta[1xtime]           [rad]       Declination in earth centered inertial equatorial frame
%   lon[1xtime]             [rad]       Longitude of ground track with respect to rotating earth (0 deg at Greenwhich meridian)
%   lat[1xtime]             [rad]       Latitude of ground track with respect to rotating earth

%% Constants/Parameters

% Earth parameters

mu_E = astroConstants(13);
w_e = 15.04*pi/180/3600; % rad/s

% Initial orbit parameters
[ T0, ~, ~, ~, ~ ] = orbitParameters( a0, e0, i0, OM0, om0, theta0, mu_E);

%% Orbit propagation from t0

t_end = T0 .* orbits;
TSPAN = linspace( t0, t_end, points );
[ r0, v0 ] = kp2rv( a0, e0, i0, OM0, om0, theta0, mu_E );


[ r, v, TSPAN ] = orbit_propagation( r0, v0, TSPAN, mu_E, perturbationModel, relTol, absTol );

% Find theta(t)
theta = zeros( 1, length( r(:,1) ) );

for k = 1:size( r(:,1) )

    [ ~, ~, ~, ~, ~, theta_temp ] = rv2kp( r(k,:), v(k,:) , mu_E );
    theta(k) = theta_temp;

end

%% Calculate right ascension and declination

delta = zeros(size(theta));
alpha = zeros(size(theta));

for k = 1:length(theta)
    
    r_norm = norm( r(k,:) );

    delta(k) = asin( r(k,3) / r_norm );

    if r(k,2) > 0

        alpha(k) = acos( ( r(k,1) / r_norm ) / cos( delta(k) ) );

    else

        alpha(k) = 2*pi - acos( ( r(k,1) / r_norm ) / cos( delta(k) ) );

    end

end

%% Calculate longitude and latitude

RA_G = wrapTo2Pi(RA0_G + w_e .* ( TSPAN - t0 )); % Right ascension of Greenwhich in time
lon = wrapTo180(rad2deg(alpha - RA_G'));
lat = rad2deg(delta);

%% Plot results

% Figure and image

switch figureNum

    case 'new'

        figNew = figure;
        image = imread('https://www.solarsystemscope.com/textures/download/8k_earth_daymap.jpg');
        imagesc( [-180,180], [-90,90], flipud(image) );
        set(gca,'ydir','normal');
        hold on
        xlabel('\bf{Longitude [deg]}','interpreter', 'latex','fontsize', fontsize)
        ylabel('\bf{Latitude [deg]}','interpreter', 'latex','fontsize', fontsize)
        title('\bf{Ground Track}','interpreter', 'latex','fontsize', fontsize)
        axis equal
        ylim([-90,90])
        xticks(-180:30:180)
        yticks(-90:30:90)
        plotHandle = figNew.Number;

    otherwise

        figOld = figure(figureNum);
        plotHandle = figOld.Number;
        hold on
        grid on

end

plot( lon, lat, '.', 'color', linecolor, 'MarkerSize', markersizeLine, 'handlevisibility', 'off' );

switch figureNum

    case 'new'

        plot(lon(1),lat(1),'color', linecolor)
        plot(lon(1), lat(1), 'x', 'color', linecolor, 'markersize', markersizeEndStart, 'linewidth', 15);
        plot(lon(end), lat(end), 'o', 'color', linecolor, 'markersize', markersizeEndStart, 'linewidth', 15);
        legend([num2str(orbits) ,' orbits'], 'Start', 'End', 'interpreter', 'latex','fontsize', fontsize);

    otherwise

        plot(lon(1),lat(1),'color', linecolor, 'linewidth', 2.5, 'Displayname', [num2str(orbits) ,' orbits'])
        plot(lon(1), lat(1), 'x', 'color', linecolor, 'markersize', markersizeEndStart, 'linewidth', 15, 'Displayname', 'Start');
        plot(lon(end), lat(end), 'o', 'color', linecolor, 'markersize', markersizeEndStart, 'linewidth', 15, 'Displayname', 'End');

end

drawnow;

end