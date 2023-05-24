function [lon, lat] = groundtrack(y ,greenwich0, t_span, planetID)

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

%% Setup Plotting Environment
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

%% Calculating Ground Track

if planetID == 3
    omega_planet = 2*pi/(23*3600+56*60+4.09); %angular velocity of Earth
end

a = 1./(2./vecnorm(y(1:3,1)) - dot(y(4:6,1),y(4:6,1),1)./astroConstants(10+planetID)); % Semi major axis [km]
T = 2.*pi.*(a.^3./astroConstants(10+planetID)).^0.5; % Time period of Orbit [s]

% Repeating track ccalculations
delta_lambda = T*omega_planet;

thetaG = wrapTo2Pi(greenwich0 - omega_planet*t_span);

[~, delta, sc] = car2AlfaDelta(y);

lat = rad2deg(delta);
lon = rad2deg(wrapTo2Pi(atan2(sc(2,:),sc(1,:))+thetaG)-pi); % vector of longitudes

tol = 270; % distance > tol indicates discontinuity
dl = diff([lon;lat],1,2); % look up what that command does if you don't know it
euler_dist = sqrt((dl(1,:)+dl(2,:)).^2); % distance between data points
jumpind = [0 euler_dist>tol]; % now if jumpind(i) = true, we know that the 
                 %   point [lat(i) lon(i)] is the first after a jump
blocks = cumsum(jumpind); % points that belong to the same continuous part
                           % have the same value in blocks

%% Plot ground track

planetNames = {'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto', 'Sun'};
N = {dir(fullfile(strcat(pwd, '\textures'),'*.jpg')).name};
cdata = flip(imread(fullfile(strcat(pwd, '\textures'),N{contains(N,planetNames{planetID})})));

figure();
%set(OrbitEvolution, 'Units', 'Normalized', 'OuterPosition', [.15 .25 .7 .7]);
pbaspect([16 9 1])
hold on;
xlim([-180,180]); 
ylim([-90 90]);
set(gca,'XTick',-180:15:180,'XTickMode','manual');
set(gca,'YTick',-90:15:90,'YTickMode','manual');
xlabel('Longitude $\lambda$ [deg]');
ylabel('Latitude $\phi$ [deg]');
title('Ground Track Plot');

start_pt = plot(lon(1),lat(1), 'ok', 'LineWidth',2, 'MarkerSize',10);
uistack(start_pt,"top")
end_pt = plot(lon(end),lat(end), 'sk', 'LineWidth',2, 'MarkerSize',10);
uistack(end_pt,"top")
bg=imagesc([-180,180],[-90, 90],cdata);
uistack(bg,"bottom")

%plot(lon, lat, '.r');
% Now just loop over the continuous blocks to draw a separate line for each one
for i=0:blocks(end)
    track_pt = plot(lon(blocks==i),lat(blocks==i),'LineWidth',2,'Color','r');
    uistack(track_pt,"down",2)
    hold on;
end
legend([start_pt,end_pt],'Start','End','Location','best')

%% Reset Plotting Environment
set(groot, 'defaultLegendInterpreter','tex');
set(groot, 'defaultAxesTickLabelInterpreter','tex');
set(groot, 'defaultTextInterpreter','tex');

end