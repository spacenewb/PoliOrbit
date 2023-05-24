function [ orbit_object ] = orbitPlotrv( r , linestyle, figureNum, color, linewidth, LegendVisibility)
% orbitPlot creates a plot of the orbit on a new or existing figure
%
% PROTOTYPE:
%   [ orbit_object ] = orbitPlotrv( r , linestyle, figureNum, color, linewidth, LegendVisibility)
%
% INPUT:
%   r[timex3]           Position vector IN TIME    [km]
%   linestyle[string]   Linestyle '-', ''
%   figureNum           Plots the orbit on pre-existing figure numbered 'figureNum'
%   color               RGB color for plot
%   linewidth           Linewidth of orbit plot  
%   LegendVisibility    1: Visible in legend
%                       0: Not visible in legend
%
% OUTPUT:
%   orbit_object[1]     Object of the orbit plot
%
% NOTE:
%   Plot is done with axis equal turned on. This way all versors
%   represented will be in scale, and you can see if they are perpendicular
%   to eachother.

figure(figureNum);
hold on
grid on

fontsize = 18;

if LegendVisibility == 1
    orbit_object = plot3( r(:,1), r(:,2), r(:,3), linestyle, 'linewidth', linewidth, 'color', color );
else
    orbit_object = plot3( r(:,1), r(:,2), r(:,3), linestyle, 'linewidth', linewidth, 'color', color, 'handlevisibility', 'off' );
end

hold on
grid on
xlabel('\bf{x [km]}', 'interpreter', 'latex', 'fontsize', fontsize)
ylabel('\bf{y [km]}', 'interpreter', 'latex', 'fontsize', fontsize)
zlabel('\bf{z [km]}', 'interpreter', 'latex', 'fontsize', fontsize)
axis equal
end