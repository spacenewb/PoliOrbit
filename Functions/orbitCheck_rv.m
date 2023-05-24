function [ eps, h, h_norm, e, e_norm ] = orbitCheck_rv( r, v, mu, graph_parameters, linestyle, TSPAN )
% orbitCheck_rv Used to see if input orbit is correct using conservation laws
%
% PROTOTYPE:
%   [ eps, h, h_norm, e, e_norm ] = orbitCheck_rv( r, v, mu, graph_parameters, linestyle, TSPAN )
%
% INPUT:
%   r[timex3]           Position vector         IN TIME    [km]
%   v[timex3]           Velocity vector         IN TIME    [km/s]
%   mu[1x1]             Planetary constant                 [km^3/s^2]
%   graph_parameters    If specified, creates time-dependent graphs of parameters eps, ||h||, ||e||
%                           0 : Does not graph parameters
%                           1 : Graphs parameters
%   linestyle of lines  '-', '--', ...
%   TSPAN               Time vector used for graphs        [s]
%  
% OUTPUT:
%   eps[timex1]         Specific energy         IN TIME    [km^2/s^2]
%   h[timex3]           Ang. momentum vector    IN TIME    [km^2/s]
%   h_norm[timex1]      Ang. momentum magnitude IN TIME    [km^2/s]
%   e[timex3]           Eccentricity vector     IN TIME    [adimensional]
%   e_norm[timex1]      Eccentricity magnityde  IN TIME    [adimensional]
%

%% Check inputs

if ~isequal(size(r),size(v))
    error('Input vectors (arrays) must be of same size')
end

if ~isequal(length(TSPAN), size(r,1))
    error('TSPAN and dimensions of r,v must be of same size')
end


%% Calculate orbit parameters

% h, e
h = cross( r , v );
h_norm = norm_Timedependent( h );
r_norm = norm_Timedependent( r );
eps = 1/2 .* norm_Timedependent(v).^2 - mu ./ r_norm;
e = cross( v, h ) ./ mu - r ./ r_norm;
e_norm = norm_Timedependent( e );

r_unitvect = r ./ r_norm;
h_unitvect = h ./ h_norm;
v_r = [];
v_theta = [];
for i = 1:size(v,1)
    theta_unitvect = cross( h_unitvect(i,:), r_unitvect(i,:) );
    v_r = [v_r; dot(v(i,:), r_unitvect(i,:))];
    v_theta = [v_theta; dot(v(i,:), theta_unitvect)];
end

%% Graph

fontsize = 18;

if graph_parameters == 1
    % Create graphs

    % Specific energy
    figure();
    hold on
    grid on
    plot(TSPAN, eps, '-', 'linewidth', 2.5)
    title('\bf{Specific energy \boldmath$\varepsilon$}', 'interpreter', 'latex', 'fontsize', fontsize)
    ylabel('\boldmath$\bf{\varepsilon\ [km^2/s^2]}$', 'interpreter', 'latex', 'fontsize', fontsize)
    xlabel('\bf{t [s]}', 'interpreter', 'latex', 'fontsize', fontsize)
    legend('\boldmath$\bf{\varepsilon}$', 'interpreter', 'latex', 'fontsize', fontsize)
    xlim([TSPAN(1), TSPAN(end)])

    % Angular momentum
    figure();
    hold on
    grid on
    plot(TSPAN, h(:,1), linestyle, 'linewidth', 2.5, 'color', [0,0,1] );
    plot(TSPAN, h(:,2), linestyle, 'linewidth', 2.5, 'color', [1,0,0] );
    plot(TSPAN, h(:,3), linestyle, 'linewidth', 2.5, 'color', [0,1,0] );
    plot(TSPAN, h_norm, linestyle, 'linewidth', 2.5, 'color', [0,0,0] );
    title('\bf{Angular momentum \boldmath$\underline{h}$}', 'interpreter', 'latex', 'fontsize', fontsize)
    ylabel('$\bf{h_x,\ h_y,\ h_z,\ ||\underline{h}||\ [km^2/s]}$', 'interpreter', 'latex', 'fontsize', fontsize)
    xlabel('\bf{t [s]}', 'interpreter', 'latex', 'fontsize', fontsize)
    legend( '\boldmath{$h_x$}', '\boldmath{$h_y$}', '\boldmath{$h_z$}', '$\bf{||\underline{h}||}$', 'interpreter', 'latex', 'fontsize', fontsize)
    xlim([TSPAN(1), TSPAN(end)])

    % Eccentricity
    figure();
    hold on
    grid on
    plot(TSPAN, e(:,1), linestyle, 'linewidth', 2.5, 'color', [0,0,1])
    plot(TSPAN, e(:,2), linestyle, 'linewidth', 2.5, 'color', [1,0,0])
    plot(TSPAN, e(:,3), linestyle, 'linewidth', 2.5, 'color', [0,1,0])
    plot(TSPAN, e_norm, linestyle, 'linewidth', 2.5, 'color', [0,0,0] )
    title('\bf{Eccentricity \boldmath$\underline{e}$}', 'interpreter', 'latex', 'fontsize', fontsize)
    ylabel('$\bf{||\underline{e}||}$', 'interpreter', 'latex', 'fontsize', fontsize)
    xlabel('\bf{t [s]}', 'interpreter', 'latex', 'fontsize', fontsize)
    legend( '\boldmath{$e_x$}', '\boldmath{$e_y$}', '\boldmath{$e_z$}', '\boldmath$\bf{||\underline{e}||}$', 'interpreter', 'latex', 'fontsize', fontsize)
    xlim([TSPAN(1), TSPAN(end)])

    % e*h
    figure();
    hold on
    grid on
    plot(TSPAN, dot_Timedependent(e,h), linestyle, 'linewidth', 2.5, 'color', [0,0,1])
    title('\boldmath$\underline{e}\cdot\underline{h}$', 'interpreter', 'latex', 'fontsize', fontsize)
    ylabel('\boldmath$\underline{e}\cdot\underline{h}\ [km^2/s]$', 'interpreter', 'latex', 'fontsize', fontsize)
    xlabel('\bf{t [s]}', 'interpreter', 'latex', 'fontsize', fontsize)
    legend('\boldmath$\underline{e}\cdot\underline{h}$', 'interpreter', 'latex', 'fontsize', fontsize)
    xlim([TSPAN(1), TSPAN(end)])
    

    % v_r, v_theta
    figure();
    hold on
    grid on
    plot(TSPAN, v_r, linestyle,'linewidth', 2.5, 'color', [0,0,1] )
    plot(TSPAN, v_theta, linestyle, 'linewidth', 2.5, 'color', [1,0,0] )
    title('Radial and transversal velocities', 'interpreter', 'latex', 'fontsize', fontsize)
    xlabel('\bf{t [s]}', 'interpreter', 'latex', 'fontsize', fontsize)
    ylabel('\boldmath$v_r,\ v_{\theta}\ [km/s]$', 'interpreter', 'latex', 'fontsize', fontsize)
    legend('\boldmath{$v_r$}', '\boldmath{$v_{\theta} $}', 'interpreter', 'latex', 'fontsize', fontsize)
    
    % Plot remaining orbital elements a, i, OM, om, theta

    for k = 1:length(r(:,1))
        [ a(k), ~, i(k), OM(k), om(k), theta(k) ] = rv2kp( r(k,:), v(k,:) , mu );
    end

    % a - semi-major axis
    figure();
    hold on
    grid on
    plot( TSPAN, a, linestyle, 'linewidth', 2.5,'color', [0,0,1] );
    title( 'Semi-major axis', 'interpreter', 'latex', 'fontsize', fontsize )
    xlabel( '\bf{t [s]}', 'interpreter', 'latex', 'fontsize', fontsize )
    ylabel( '\boldmath$a\ [km]$', 'interpreter', 'latex', 'fontsize', fontsize )
    legend( '\boldmath{$a$}', 'interpreter', 'latex', 'fontsize', fontsize )

    % i - inclination
    figure();
    hold on
    grid on
    plot(TSPAN, rad2deg(i), linestyle, 'linewidth', 2.5,'color', [0,0,1]);
    title('Inclination', 'interpreter', 'latex', 'fontsize', fontsize)
    xlabel('\bf{t [s]}', 'interpreter', 'latex', 'fontsize', fontsize)
    ylabel('\boldmath$i\ [deg]$', 'interpreter', 'latex', 'fontsize', fontsize)
    legend('\boldmath{$i$}', 'interpreter', 'latex', 'fontsize', fontsize)

    % OM - RAAN
    figure();
    hold on
    grid on
    plot(TSPAN, rad2deg(OM), linestyle, 'linewidth', 2.5,'color', [0,0,1]);
    title('RAAN', 'interpreter', 'latex', 'fontsize', fontsize)
    xlabel('\bf{t [s]}', 'interpreter', 'latex', 'fontsize', fontsize)
    ylabel('\boldmath$\Omega\ [deg]$', 'interpreter', 'latex', 'fontsize', fontsize)
    legend(['\boldmath{$\Omega$}'], 'interpreter', 'latex', 'fontsize', fontsize)

    % om - Argument of pericenter
    figure();
    hold on
    grid on
    plot(TSPAN, rad2deg(om), linestyle, 'linewidth', 2.5,'color', [0,0,1]);
    title('Argument of pericenter', 'interpreter', 'latex', 'fontsize', fontsize)
    xlabel('\bf{t [s]}', 'interpreter', 'latex', 'fontsize', fontsize)
    ylabel('\boldmath$\omega\ [deg]$', 'interpreter', 'latex', 'fontsize', fontsize)
    legend('\boldmath{$\omega$}', 'interpreter', 'latex', 'fontsize', fontsize)

     % theta - True anomaly
    figure();
    hold on
    grid on
    plot(TSPAN, rad2deg(theta), linestyle, 'linewidth', 2.5,'color', [0,0,1]);
    title('True anomaly', 'interpreter', 'latex', 'fontsize', fontsize)
    xlabel('\bf{t [s]}', 'interpreter', 'latex', 'fontsize', fontsize)
    ylabel('\boldmath$\theta\ [deg]$', 'interpreter', 'latex', 'fontsize', fontsize)
    legend('\boldmath{$\theta$}', 'interpreter', 'latex', 'fontsize', fontsize)

end

end

