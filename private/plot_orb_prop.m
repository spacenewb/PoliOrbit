function [] = plot_orb_prop(Mu_Central, y0, n_pts, color)

% Calculate the orbit trajectory
[~,Y] = orbitPropagator(y0, Mu_Central, [], 1, 1000);
y2 = Y(:,end);

%% Plot

% Setup Plotting Environment
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

% Orbit
hold on
plot3( Y(1,:), Y(2,:), Y(3,:), '-', 'LineWidth', 2, 'Color',color);
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
