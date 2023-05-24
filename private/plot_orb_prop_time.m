function [] = plot_orb_prop_time(Mu_Central, y0, tspan, n_pts, color, line_width, line_alfa, line_style)

% Calculate the orbit trajectory
[~,Y] = orbitPropagator_time(y0, Mu_Central, [], tspan);
y2 = Y(:,end);

%% Plot

% Setup Plotting Environment
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

% Orbit
hold on
if isempty(line_width)
    line_width = 2;
end
if isempty(line_alfa)
    line_alfa = 1;
end
if isempty(line_style)
    line_style = '-';
end

pp = plot3( Y(1,:), Y(2,:), Y(3,:), line_style, 'LineWidth', line_width, 'Color',color);
pp.Color(4) = line_alfa;
if n_pts == 1
    scatter3(y0(1),y0(2),y0(3),80,'filled','ob','MarkerEdgeColor','k');
elseif n_pts == 2
    scatter3(y0(1),y0(2),y0(3),80,'filled','ob','MarkerEdgeColor','k');
    scatter3(y2(1),y2(2),y2(3),80,'filled','oy','MarkerEdgeColor','k');
end
hold off

% Reset Plotting Environment
set(groot, 'defaultLegendInterpreter','tex');
set(groot, 'defaultAxesTickLabelInterpreter','tex');
set(groot, 'defaultTextInterpreter','tex');

end
