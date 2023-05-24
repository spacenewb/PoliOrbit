function [] = plot_orb_mot_const(y0,mu,n_orbs,n_eval,pert,options)

%plot_orb_mot_const Calculate and plot evolution of orbital constants 
%               for the two body problem ( Keplerian motion). Function is
%               vectorised, avoid for-loops for performance.
%
% PROTOTYPE:
%   [] = plot_orb_mot_const(y0,mu,tspan,pert,options)
%
% INPUT:
%   y0[6x1]         State of the body at initial time( rx , ry , rz , vx , vy , vz ) [ L, L/T ]
%   mu[1]           Gravitational parameter of the primary [ L^3/T^2 ]
%   tspan[nx1]      Array of times at which to evaluate integral [ T ]
%   pert[-]         Type of perturbation to include [ - ]
%   options[-]      Options to setup ode solver [ - ]
%
% OUTPUT:
%   Graphs
%
% CALLED FUNCTIONS:
%   astroConstants.m
%   ode_2bp
%   orb_mot_const
%   plotPlanet
%   orbitPropagator
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

%% Identifying applied perturbations
pert_logical = zeros([1,2]);
pert_logical(1) = contains(pert,'0','IgnoreCase',true);
pert_logical(2) = contains(pert,'J2','IgnoreCase',true);

%% Either Perturbed or Unperturbed

if isempty(pert)
    pert = '0';
end

% Single Computation Case - No Comparission in Graph
if xor(pert_logical(1),pert_logical(2))
    % Determine which case for perturbation
    pert_switch = find(pert_logical);
    switch pert_switch
    case 1
        pert = '0';
        pert_text = '';
    case 2
        pert = 'J2';
        pert_text = ' - J2';
    end
    
    % Orbit Propagation
    % Perform the integration
    %[ ~, Y ] = ode113( @(t,y) ode_2bp(t,y,mu,pert), tspan , y0, options );
    [tspan,Y] = orbitPropagator(y0, mu, pert, n_orbs, n_eval);

    % Calculate orbital motion constants
    [eps,h,e,~,v_r,v_t,e_dot_h] = orb_mot_const(Y,mu);

    % Plot the graphs
    % Orbit
    figure()
    hold on
    plotPlanet(3, [0,0,0]);
    plot3( Y(1,:), Y(2,:), Y(3,:), 'b-', 'LineWidth', 2);
    scatter3(y0(1),y0(2),y0(3),40,'filled','ok');
    xlabel('$X$ [km]'); ylabel ('$Y$ [km]'); zlabel ('$Z$ [km]');
    title('Two body problem - Orbit');
    legend('Orbit','$r_o$','Location','best');
    axis equal;
    grid on;
    view(3);
    hold off

    % Epsilon - Specific Energy
    figure();
    hold on
    plot( tspan, eps, 'b-', 'LineWidth', 2);
    xlabel('Time [s]'); ylabel ('Specific Orbital Energy - $\epsilon$ [km$^2$/s$^2$]');
    ytickformat('%.8f')
    title('Two body problem - Specific Energy');
    legend(strcat('Specific Energy',pert_text),'Location','best');
    grid on;
    hold off
    
    % h - Angular Momentum
    figure();
    hold on
    plot( tspan, h(1,:), 'r-', 'LineWidth', 2); plot( tspan, h(2,:), 'g-', 'LineWidth', 2); plot( tspan, h(3,:), 'b-', 'LineWidth', 2)
    plot( tspan, vecnorm(h), 'k-', 'LineWidth', 2);
    xlabel('Time [s]'); ylabel ('Specific Angular Momentum - $\vec{h}$ [km$^2$/s]');
    title('Two body problem - Specific Angular Momentum');
    legend(strcat('$h_x$',pert_text),strcat('$h_y$',pert_text),strcat('$h_z$',pert_text),strcat('$||h||$',pert_text)','Location','best');
    grid on;
    hold off
    
    % e - Eccentricity
    figure();
    hold on
    plot( tspan, e(1,:), 'r-', 'LineWidth', 2); plot( tspan, e(2,:), 'g-', 'LineWidth', 2); plot( tspan, e(3,:), 'b-', 'LineWidth', 2)
    plot( tspan, vecnorm(e), 'k-', 'LineWidth', 2);
    xlabel('Time [s]'); ylabel ('Orbital Eccentricity - $\vec{e}$ [-]');
    title('Two body problem - Orbital Eccentricity');
    legend(strcat('$e_x$',pert_text),strcat('$e_y$',pert_text),strcat('$e_z$',pert_text),strcat('$||e||$',pert_text)','Location','best');
    grid on;
    hold off
    
    % Velocity - V_r, V_t
    figure();
    hold on
    plot( tspan, v_r, 'r-', 'LineWidth', 2); plot( tspan, v_t, 'g-', 'LineWidth', 2);
    xlabel('Time [s]'); ylabel ('Velocity - $V_r$, $V_t$ [km/s]');
    title('Two body problem - Velocity');
    legend(strcat('$V_r$',pert_text),strcat('$V_t$',pert_text),'Location','best');
    grid on;
    hold off
    
    % e dot h - Perpendicularity of e and h vectors
    figure();
    hold on
    plot( tspan, e_dot_h, 'r-', 'LineWidth', 2);
    xlabel('Time [s]'); ylabel ('$\vec{e}\cdot\vec{h}$ [km/s]');
    title('Two body problem - Perpendicularity');
    legend(strcat('$\vec{e}\cdot\vec{h}$',pert_text),'Location','best');
    grid on;
    hold off
end

%% Perturbed and Unperturbed
% Twin Computation Case - Comparission in Graph
if pert_logical(1) && pert_logical(2)

    % Orbit Propagation
    % Perform the integration
    %[ ~, Y ] = ode113( @(t,y) ode_2bp(t,y,mu,'0'), tspan , y0, options );
    [~,Y] = orbitPropagator(y0, mu, '0', n_orbs, n_eval);
    %[ ~, Y_pert ] = ode113( @(t,y) ode_2bp(t,y,mu,'J2'), tspan , y0, options );
    [tspan,Y_pert] = orbitPropagator(y0, mu, 'J2', n_orbs, n_eval);
    size(Y)
    
    % Calculate orbital motion constants
    [eps,h,e,~,v_r,v_t,e_dot_h] = orb_mot_const(Y,mu);
    [eps_pert,h_pert,e_pert,~,v_r_pert,v_t_pert,e_dot_h_pert] = orb_mot_const(Y_pert,mu);

    % Plot the graphs
    % Orbit
    figure();
    hold on
    plotPlanet(3, [0,0,0]);
    plot3( Y(1,:), Y(2,:), Y(3,:), 'b-', 'LineWidth', 2);
    plot3( Y_pert(1,:), Y_pert(2,:), Y_pert(3,:), 'r--', 'LineWidth', 2);
    scatter3(y0(1),y0(2),y0(3),40,'filled','ok');
    xlabel('$X$ [km]'); ylabel ('$Y$ [km]'); zlabel ('$Z$ [km]');
    title('Two body problem - Orbit');
    legend('Orbit','Orbit - J2','$r_o$','Location','best');
    axis equal;
    grid on;
    view(3);
    hold off

    % Epsilon - Specific Energy
    figure();
    hold on
    plot( tspan, eps, 'b-', 'LineWidth', 2);
    plot( tspan, eps_pert, 'r--', 'LineWidth', 2);
    xlabel('Time [s]'); ylabel ('Specific Orbital Energy - $\epsilon$ [km$^2$/s$^2$]');
    ytickformat('%.8f')
    title('Two body problem - Specific Energy');
    legend('$\epsilon$','$\epsilon$ - J2','Location','best');
    grid on;
    hold off
    
    % h - Angular Momentum
    figure();
    hold on
    plot( tspan, h(1,:), 'r-', 'LineWidth', 2); plot( tspan, h(2,:), 'g-', 'LineWidth', 2); plot( tspan, h(3,:), 'b-', 'LineWidth', 2)
    plot( tspan, vecnorm(h), 'k-', 'LineWidth', 1);
    plot( tspan, h_pert(1,:), 'm--', 'LineWidth', 2); plot( tspan, h_pert(2,:), 'y--', 'LineWidth', 2); plot( tspan, h_pert(3,:), 'c--', 'LineWidth', 2)
    plot( tspan, vecnorm(h_pert), 'k--', 'LineWidth', 1);
    xlabel('Time [s]'); ylabel ('Specific Angular Momentum - $\vec{h}$ [km$^2$/s]');
    title('Two body problem - Specific Angular Momentum');
    legend('$h_x$','$h_y$','$h_z$','$||h||$','$h_x$ - J2','$h_y$ - J2','$h_z$ - J2','$||h||$ - J2','Location','best');
    grid on;
    hold off
    
    % e - Eccentricity
    figure();
    hold on
    plot( tspan, e(1,:), 'r-', 'LineWidth', 2); plot( tspan, e(2,:), 'g-', 'LineWidth', 2); plot( tspan, e(3,:), 'b-', 'LineWidth', 2)
    plot( tspan, vecnorm(e), 'k-', 'LineWidth', 1);
    plot( tspan, e_pert(1,:), 'm--', 'LineWidth', 2); plot( tspan, e_pert(2,:), 'y--', 'LineWidth', 2); plot( tspan, e_pert(3,:), 'c--', 'LineWidth', 2)
    plot( tspan, vecnorm(e_pert), 'k--', 'LineWidth', 1);
    xlabel('Time [s]'); ylabel ('Orbital Eccentricity - $\vec{e}$ [-]');
    title('Two body problem - Orbital Eccentricity');
    legend('$e_x$','$e_y$','$e_z$','$||e||$','$e_x$ - J2','$e_y$ - J2','$e_z$ - J2','$||e||$ - J2','Location','best');
    grid on;
    hold off
    
    % Velocity - V_r, V_t
    figure();
    hold on
    plot( tspan, v_r, 'r-', 'LineWidth', 2); plot( tspan, v_t, 'g-', 'LineWidth', 2);
    plot( tspan, v_r_pert, 'm--', 'LineWidth', 2); plot( tspan, v_t_pert, 'y--', 'LineWidth', 2);
    xlabel('Time [s]'); ylabel ('Velocity - $V_r$, $V_t$ [km/s]');
    title('Two body problem - Velocity');
    legend('$V_r$','$V_t$','$V_r$ - J2','$V_t$ - J2','Location','best');
    grid on;
    hold off
    
    % e dot h - Perpendicularity of e and h vectors
    figure();
    hold on
    plot( tspan, e_dot_h, 'b-', 'LineWidth', 2);
    plot( tspan, e_dot_h_pert, 'r--', 'LineWidth', 2);
    xlabel('Time [s]'); ylabel ('$\vec{e}\cdot\vec{h}$ [km/s]');
    title('Two body problem - Perpendicularity');
    legend('$\vec{e}\cdot\vec{h}$','$\vec{e}\cdot\vec{h}$ - J2','Location','best');
    grid on;
    hold off

end

%% Reset Plotting Environment
set(groot, 'defaultLegendInterpreter','tex');
set(groot, 'defaultAxesTickLabelInterpreter','tex');
set(groot, 'defaultTextInterpreter','tex');

end
