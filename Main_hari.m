clc, clearvars, close all;

%% Initialisation
Mu_Sun = astroConstants(4);     % Sun's gravitational parameter [km^3/s^2];

Dep.planetID = 5;
Fly.planetID = 3;
Arr.planetID = 2;

Window.Earliest = [2031 03 01 00 00 00];
Window.Latest = [2071 03 01 00 00 00];

Leg2.Synodic = Synodic_period( Fly.planetID, Arr.planetID ) / 24 / 3600; % [d]


%% Porkchop Leg - 1

Leg1.Synodic = Synodic_period( Dep.planetID, Fly.planetID ) / 24 / 3600; % [d]

Leg1.DepartureEarliest = Window.Earliest;
Leg1.DepartureLatest = mjd20002date( date2mjd2000( Window.Earliest ) + 5*Leg1.Synodic );

Leg1.ArrivalEarliest = mjd20002date( date2mjd2000( Window.Earliest ) + 0*Leg1.Synodic );
Leg1.ArrivalLatest = mjd20002date( date2mjd2000( Leg1.DepartureLatest ) + 0*Leg1.Synodic );

Leg1.C3Launcher = inf;
Leg1.points = 200;
Leg1.DV_range_plot = 7;
Leg1.plot_params = [1 0 1 inf 10];

[Leg1.DV_req, Leg1.min_DV_req, Leg1.best_dep_mjd2k, Leg1.best_arr_mjd2k,...
    Leg1.Best_times_opt, Leg1.best_DV_opt, Leg1.PorkchopPlot] = Pork_plot(Mu_Sun, Dep.planetID,...
    Leg1.DepartureEarliest, Leg1.DepartureLatest, Fly.planetID, Leg1.ArrivalEarliest,...
    Leg1.ArrivalLatest, Leg1.points, Leg1.DV_range_plot, Leg1.plot_params);

%% New Porkchop Leg - 1

Leg1n.DV_cutoff = 8.5;
Leg1n.maxToF_estimate = 4000;

Leg1n.DepartureEarliest = mjd20002date( Leg1.Best_times_opt(1) - 4.5*Leg1.Synodic );
Leg1n.DepartureLatest = mjd20002date( Leg1.Best_times_opt(1) + 2*Leg1.Synodic );

Leg1n.ArrivalEarliest = mjd20002date( Leg1.Best_times_opt(2) - 0.5*Leg1.Synodic );
Leg1n.ArrivalLatest = mjd20002date( Leg1.Best_times_opt(2) + 0.5*Leg1.Synodic );

[Leg1n.DV_req, Leg1n.min_DV_req, Leg1n.best_dep_mjd2k, Leg1n.best_arr_mjd2k,...
    Leg1n.Best_times_opt, Leg1n.best_DV_opt, Leg1n.PorkchopPlot] = Pork_plot(Mu_Sun, Dep.planetID,...
    Leg1n.DepartureEarliest, Leg1n.DepartureLatest, Fly.planetID, Leg1n.ArrivalEarliest,...
    Leg1n.ArrivalLatest, Leg1.points, Leg1.DV_range_plot, Leg1.plot_params);

[Leg1n.ToF_limits, Leg1n.bounding_rect, Leg1n.rect_ratio] = pork_window_estimator(Leg1.Synodic, Leg1n.DepartureEarliest, Leg1n.DepartureLatest, Leg1n.ArrivalEarliest, Leg1n.ArrivalLatest, Leg1.points, Leg1n.DV_req, Leg1n.Best_times_opt, Leg1n.best_DV_opt, Leg1.DV_range_plot, Leg1n.DV_cutoff, Leg1n.maxToF_estimate, Leg1n.PorkchopPlot);


% Handling figures
fig_handle = Leg1.PorkchopPlot;
allaxes = findall(fig_handle, 'type', 'axes');
all_legend = findall(fig_handle, 'type', 'axes', 'Tag', 'legend');
all_colorbar = findall(fig_handle, 'type', 'axes', 'Tag', 'Colorbar');
all_plotting_axes = setdiff(allaxes, union(all_legend, all_colorbar));
optims_xx = Leg1n.Best_times_opt(1):Leg1.Synodic:date2mjd2000(Leg1.DepartureLatest);
optims_yy = Leg1n.Best_times_opt(2):Leg1.Synodic:date2mjd2000(Leg1.ArrivalLatest);
hold on
for i = 1:min([length(optims_xx), length(optims_yy)])
    rectangle(all_plotting_axes, 'Position', [optims_xx(i)-(Leg1n.bounding_rect(3)*Leg1n.rect_ratio(1)), optims_yy(i)-(Leg1n.bounding_rect(4)*Leg1n.rect_ratio(2)), Leg1n.bounding_rect(3:4)], 'EdgeColor','y', 'LineWidth',2);    
end

Dep_mjd2000_vec = [date2mjd2000(Leg1.DepartureEarliest) date2mjd2000(Leg1.DepartureLatest)];
ytt1 = tand(45)*Dep_mjd2000_vec + Leg1n.ToF_limits(1);
ytt2 = tand(45)*Dep_mjd2000_vec + Leg1n.ToF_limits(2);
% plot(all_plotting_axes, Dep_mjd2000_vec, ytt1, '-.r', 'LineWidth', 2);
% plot(all_plotting_axes, Dep_mjd2000_vec, ytt2, '-.r', 'LineWidth', 2);
hold off


%% Porkchop Leg - 2

Leg2.Synodic = Synodic_period( Fly.planetID, Arr.planetID ) / 24 / 3600; % [d]

Leg2.DepartureEarliest = Window.Earliest;
Leg2.DepartureLatest = mjd20002date( date2mjd2000( Window.Earliest ) + 5*Leg2.Synodic );

Leg2.ArrivalEarliest = mjd20002date( date2mjd2000( Window.Earliest ) + 0*Leg2.Synodic );
Leg2.ArrivalLatest = mjd20002date( date2mjd2000( Leg2.DepartureLatest ) + 0*Leg2.Synodic );

Leg2.C3Launcher = inf;
Leg2.points = 200;
Leg2.DV_range_plot = 7;
Leg2.plot_params = [1 0 2 inf 25];

[Leg2.DV_req, Leg2.min_DV_req, Leg2.best_dep_mjd2k, Leg2.best_arr_mjd2k,...
    Leg2.Best_times_opt, Leg2.best_DV_opt, Leg2.PorkchopPlot] = Pork_plot(Mu_Sun, Fly.planetID,...
    Leg2.DepartureEarliest, Leg2.DepartureLatest, Arr.planetID, Leg2.ArrivalEarliest,...
    Leg2.ArrivalLatest, Leg2.points, Leg2.DV_range_plot, Leg2.plot_params);

%% New Porkchop Leg - 2

Leg2n.DV_cutoff = 5;
Leg2n.maxToF_estimate = 500;

Leg2n.DepartureEarliest = mjd20002date( Leg2.Best_times_opt(1) - 0.5*Leg2.Synodic );
Leg2n.DepartureLatest = mjd20002date( Leg2.Best_times_opt(1) + 0.5*Leg2.Synodic );

Leg2n.ArrivalEarliest = mjd20002date( Leg2.Best_times_opt(2) - 0.5*Leg2.Synodic );
Leg2n.ArrivalLatest = mjd20002date( Leg2.Best_times_opt(2) + 0.5*Leg2.Synodic );

[Leg2n.DV_req, Leg2n.min
    _DV_req, Leg2n.best_dep_mjd2k, Leg2n.best_arr_mjd2k,...
    Leg2n.Best_times_opt, Leg2n.best_DV_opt, Leg2n.PorkchopPlot] = Pork_plot(Mu_Sun, Fly.planetID,...
    Leg2n.DepartureEarliest, Leg2n.DepartureLatest, Arr.planetID, Leg2n.ArrivalEarliest,...
    Leg2n.ArrivalLatest, Leg2.points, Leg2.DV_range_plot, Leg2.plot_params);

[Leg2n.ToF_limits, Leg2n.bounding_rect, Leg2n.rect_ratio] = pork_window_estimator(Leg2.Synodic, Leg2n.DepartureEarliest, Leg2n.DepartureLatest, Leg2n.ArrivalEarliest, Leg2n.ArrivalLatest, Leg2.points, Leg2n.DV_req, Leg2n.Best_times_opt, Leg2n.best_DV_opt, Leg2.DV_range_plot, Leg2n.DV_cutoff, Leg2n.maxToF_estimate, Leg2n.PorkchopPlot);
Leg2n.ToF_limits = [50, 600];

% Handling figures
fig_handle = Leg2.PorkchopPlot;
allaxes = findall(fig_handle, 'type', 'axes');
all_legend = findall(fig_handle, 'type', 'axes', 'Tag', 'legend');
all_colorbar = findall(fig_handle, 'type', 'axes', 'Tag', 'Colorbar');
all_plotting_axes = setdiff(allaxes, union(all_legend, all_colorbar));
optims_xx = Leg2n.Best_times_opt(1):Leg2.Synodic:date2mjd2000(Leg2.DepartureLatest);
optims_yy = Leg2n.Best_times_opt(2):Leg2.Synodic:date2mjd2000(Leg2.ArrivalLatest);
hold on
for i = 1:length(optims_xx)
    rectangle(all_plotting_axes, 'Position', [optims_xx(i)-(Leg2n.bounding_rect(3)*Leg2n.rect_ratio(1)), optims_yy(i)-(Leg2n.bounding_rect(4)*Leg2n.rect_ratio(2)), Leg2n.bounding_rect(3:4)], 'EdgeColor','y', 'LineWidth',2);    
end

Dep_mjd2000_vec = [date2mjd2000(Leg2.DepartureEarliest) date2mjd2000(Leg2.DepartureLatest)];
ytt1 = tand(45)*Dep_mjd2000_vec + Leg2n.ToF_limits(1);
ytt2 = tand(45)*Dep_mjd2000_vec + Leg2n.ToF_limits(2);
% plot(all_plotting_axes, Dep_mjd2000_vec, ytt1, '-.r', 'LineWidth', 2);
% plot(all_plotting_axes, Dep_mjd2000_vec, ytt2, '-.r', 'LineWidth', 2);
hold off



%% Global Time Window Creation
% resolution = 5; % Days
% Global_dep_times = date2mjd2000(Window.Earliest):5:date2mjd2000(Window.Latest);
% Global_arr_times = date2mjd2000(Window.Earliest):5:date2mjd2000(Window.Latest);
% 
% num_syn_leg1_win = ceil( ( date2mjd2000(Window.Latest) - date2mjd2000(Window.Earliest) )  / Leg1.Synodic );
% num_syn_leg2_win = ceil( ( date2mjd2000(Window.Latest) - date2mjd2000(Window.Earliest) )  / Leg2.Synodic );
% 
% length_moving_window_x = ceil( Leg1n.bounding_rect(3) / resolution );
% length_moving_window_x_preceeding = floor( length_moving_window_x*Leg1n.rect_ratio(1) );
% 
% length_moving_window_y = ceil( Leg2n.bounding_rect(4) / resolution );
% length_moving_window_y_preceeding = floor( length_moving_window_y*Leg2n.rect_ratio(2) );
% 
% x_vals_mov_window = [-length_moving_window_x_preceeding:1:length_moving_window_x].*resolution;
% y_vals_mov_window = [-length_moving_window_y_preceeding:1:length_moving_window_y].*resolution;
% 
% [Mov_Win_X, Mov_Win_Y] = ndgrid(x_vals_mov_window, y_vals_mov_window);
% 
% figure()
% plot(Mov_Win_X, Mov_Win_Y, 'ok')



%% Bounding Box Minimiser


% [ dV_sat, dV_planet1, dV_planet3, dV_poweredGA, TOF1, TOF2 ] = interplanetaryTransfer_3Body( departureTime_mjd, flybyTime_mjd, arrivalTime_mjd, planet1, planet2, planet3 )

resolution = 2; % days

Dep_opt_Epoch = Leg1n.Best_times_opt(1);
Arr_opt_Epoch = Leg2n.Best_times_opt(2);

Dep_dates_optimums = Dep_opt_Epoch:Leg1.Synodic:date2mjd2000(Window.Latest);
Arr_dates_optimums = Arr_opt_Epoch:Leg2.Synodic:date2mjd2000(Window.Latest);

DV_real_box_min = NaN.*ones(length(Dep_dates_optimums), length(Arr_dates_optimums));
Box_min_dates = NaN.*ones(length(Dep_dates_optimums), length(Arr_dates_optimums), 3);

wb = waitbar(0,'Please wait...Minimising Bounding Boxes');

for i = 1:length(Dep_dates_optimums)

    waitbar(i/length(Dep_dates_optimums),wb,'Please wait...Minimising Bounding Boxes');

    for j = 1:length(Arr_dates_optimums)

        Dist_boxes_Dep_Arr = Arr_dates_optimums(j) - Dep_dates_optimums(i);

        if Dist_boxes_Dep_Arr < 6*365.25 && Dist_boxes_Dep_Arr > 1*365.25

            bb1_x = ( Dep_dates_optimums(i)-(Leg1n.bounding_rect(3)*Leg1n.rect_ratio(1)) ) + [ 0, Leg1n.bounding_rect(3) ];
            bb2_y = ( Arr_dates_optimums(j)-(Leg2n.bounding_rect(4)*Leg2n.rect_ratio(2)) ) + [ 0, Leg2n.bounding_rect(4) ];
            
            GA_leg1_lims = [ ( min(bb1_x) + min(Leg1n.ToF_limits) ) ( max(bb1_x) + max(Leg1n.ToF_limits) ) ];
            GA_leg2_lims = [ ( min(bb2_y) - max(Leg2n.ToF_limits) ) ( max(bb2_y) - min(Leg2n.ToF_limits) ) ];
            GA_leg1_win = GA_leg1_lims(1) - mod(GA_leg1_lims(1),resolution):resolution:GA_leg1_lims(2) - mod(GA_leg1_lims(2),resolution);
            GA_leg2_win = GA_leg2_lims(1) - mod(GA_leg2_lims(1),resolution):resolution:GA_leg2_lims(2) - mod(GA_leg2_lims(2),resolution);
            final_GA_win = intersect(GA_leg1_win, GA_leg2_win);
            
            lb = [min(bb1_x), min(final_GA_win), min(bb2_y)];
            ub = [max(bb1_x), max(final_GA_win), max(bb2_y)];
%             lb = [min(bb1_x), min(bb1_x), min(bb2_y)];
%             ub = [max(bb1_x), max(bb2_y), max(bb2_y)];
            
%             ga_opts = optimoptions('ga', 'Display', 'off', 'PopulationSize', 50, 'FunctionTolerance', 1e-5,...
%                 'CrossoverFraction', 0.8, 'MaxGenerations', 100, 'UseParallel', false);
            
            part_opts = optimoptions('particleswarm', 'SwarmSize', 100, 'FunctionTolerance', 1e-5,...
                'MaxIterations', 300, 'UseParallel', true, 'Display', 'off'); %, 'HybridFcn', @fmincon);

            % [Box_min_dates(i,j,:), fval] = ga(@(a) interplanetaryTransfer_3Body(a(1), a(2), a(3), Dep.planetID, Fly.planetID, Arr.planetID) , 3, [],[],[],[],lb,ub,[],ga_opts);
            
            %rng default  % For reproducibility
            %[Box_min_dates(i,j,:), fval, ~,~] = particleswarm(@(a) interplanetaryTransfer_3Body(a(1), a(2), a(3), Dep.planetID, Fly.planetID, Arr.planetID),3,lb,ub,part_opts);
            
            %[Box_min_dates(i,j,:), fval] = fmincon(@(a) interplanetaryTransfer_3Bodys(a(1), a(2), a(3), Dep.planetID, Fly.planetID, Arr.planetID), (lb+ub)'/2,[],[],[],[],lb,ub);

            %[Box_min_dates(i,j,:), fval] = patternsearch(@(a) interplanetaryTransfer_3Bodys(a(1), a(2), a(3), Dep.planetID, Fly.planetID, Arr.planetID), (lb+ub)'/2,[],[],[],[],lb,ub);
            
            fmins_opts = optimset('MaxFunEvals',2000,'MaxIter',2000);
            simanneal_opts = optimoptions('simulannealbnd');
            
            %[Box_min_dates(i,j,:), fval, ~,~] = fminsearchbnd(@(a) interplanetaryTransfer_3Bodys(a(1), a(2), a(3), Dep.planetID, Fly.planetID, Arr.planetID), [Dep_dates_optimums(i),(lb(2)+ub(2))/2,Arr_dates_optimums(j)],lb,ub,fmins_opts);
            [Box_min_dates(i,j,:), fval] = fminsearch(@(a) interplanetaryTransfer_3Bodys(a(1), a(2), a(3), Dep.planetID, Fly.planetID, Arr.planetID), [Dep_dates_optimums(i),(lb(2)+ub(2))/2,Arr_dates_optimums(j)],fmins_opts);
            %[Box_min_dates(i,j,:), fval] = simulannealbnd(@(a) interplanetaryTransfer_3Bodys(a(1), a(2), a(3), Dep.planetID, Fly.planetID, Arr.planetID), [Dep_dates_optimums(i),(lb(2)+ub(2))/2,Arr_dates_optimums(j)],lb,ub);

            DV_real_box_min(i,j) = fval(1);
            

        else
            Box_min_dates(i,j,:) = [NaN, NaN, NaN];
            DV_real_box_min(i,j) = NaN;
        end

    end
end

delete(wb);

%% Finding Minima among Bounding Box Evaluations

no_candidates = 10;
viable = numel(DV_real_box_min(DV_real_box_min<25));
no_candidates = min([no_candidates, viable]);
DV_real_box_min_temp = DV_real_box_min;
min_global_dv = zeros(no_candidates,1);
best_global_dates = zeros(no_candidates,3);
Best_dep = zeros(no_candidates,6);
Best_fly = zeros(no_candidates,6);
Best_arr = zeros(no_candidates,6);

for i=1:no_candidates
    [min_global_dv(i), min_global_dv_idx] = min(DV_real_box_min_temp(:));
    [Dep_idx_global, Arr_idx_global] = ind2sub(size(DV_real_box_min_temp), min_global_dv_idx);
    best_global_dates(i,:) = Box_min_dates(Dep_idx_global, Arr_idx_global, :);
    Best_dep(i,:) = mjd20002date(best_global_dates(i,1));
    Best_fly(i,:) = mjd20002date(best_global_dates(i,2));
    Best_arr(i,:) = mjd20002date(best_global_dates(i,3));
    DV_real_box_min_temp(Dep_idx_global, Arr_idx_global) = NaN;
end
format long g
[Best_dep(:,1:3), Best_fly(:,1:3), Best_arr(:,1:3)]
min_global_dv

%% Penultimate Search

for i=1:no_candidates
    curr_dates = best_global_dates(i,:);
    %[curr_min_dates(i,:), curr_min_val(i)] = fmincon(@(a) interplanetaryTransfer_3Bodys(a(1), a(2), a(3), Dep.planetID, Fly.planetID, Arr.planetID), curr_dates,[],[],[],[],lb,ub);
    %[curr_min_dates(i,:), curr_min_val(i)] = fminunc(@(a) interplanetaryTransfer_3Bodys(a(1), a(2), a(3), Dep.planetID, Fly.planetID, Arr.planetID), curr_dates);
    [curr_min_dates(i,:), curr_min_val(i)] = fminsearch(@(a) interplanetaryTransfer_3Body(a(1), a(2), a(3), Dep.planetID, Fly.planetID, Arr.planetID), curr_dates);
    Best_curr_dep(i,:) = mjd20002date(curr_min_dates(i,1));
    Best_curr_fly(i,:) = mjd20002date(curr_min_dates(i,2));
    Best_curr_arr(i,:) = mjd20002date(curr_min_dates(i,3));
end
format long g
[Best_curr_dep(:,1:3), Best_curr_fly(:,1:3), Best_curr_arr(:,1:3)]
curr_min_val

%% Final Search
% Varying GA only
Final_sol_temp = curr_min_dates(1,:);
[fin_min_dates, fin_min_val] = fminunc(@(a) interplanetaryTransfer_3Body(Final_sol_temp(1), a, Final_sol_temp(3), Dep.planetID, Fly.planetID, Arr.planetID), Final_sol_temp(2));
%[fin_min_dates, fin_min_val] = fminbnd(@(a) interplanetaryTransfer_3Bodys(Final_sol_temp(1), a, Final_sol_temp(3), Dep.planetID, Fly.planetID, Arr.planetID), Final_sol_temp(1), Final_sol_temp(2));
Best_fin_dep = mjd20002date(Final_sol_temp(1));
Best_fin_fly = mjd20002date(fin_min_dates);
Best_fin_arr = mjd20002date(Final_sol_temp(3));
format long g
[Best_fin_dep(1:3), Best_fin_fly(1:3), Best_fin_arr(1:3)]
fin_min_val

%% Final Values
Best_MJD2k_Global = [Final_sol_temp(1),fin_min_dates,Final_sol_temp(3)]
Best_Dates_Global = [Best_fin_dep; Best_fin_fly; Best_fin_arr]
Min_DV_Global = fin_min_val

[ Verif_dV_sat, Verif_dV_planet1, Verif_dV_planet3, Verif_dV_poweredGA, Verif_TOF1, Verif_TOF2 ] = interplanetaryTransfer_3Body(Best_MJD2k_Global(1), Best_MJD2k_Global(2), Best_MJD2k_Global(3), Dep.planetID, Fly.planetID, Arr.planetID);
[ Verif_dV_sat, Verif_dV_planet1, Verif_dV_planet3, Verif_dV_poweredGA, Verif_TOF1/86400, Verif_TOF2/86400 ]

save('filename.mat');








%% Local Functions


function [ToF_limits, bounding_rect, rect_ratio] = pork_window_estimator(Synodic, DepartureEarliest, DepartureLatest, ArrivalEarliest, ArrivalLatest, points, DV_req, Best_times_opt, best_DV_opt, DV_range_plot, DV_cutoff, guess_max_tof, fig_handle)

% Handling figures
allaxes = findall(fig_handle, 'type', 'axes');
all_legend = findall(fig_handle, 'type', 'axes', 'Tag', 'legend');
all_colorbar = findall(fig_handle, 'type', 'axes', 'Tag', 'Colorbar');
all_plotting_axes = setdiff(allaxes, union(all_legend, all_colorbar));

%figure()
Dep_mjd2000_vec = linspace(date2mjd2000(DepartureEarliest), date2mjd2000(DepartureLatest), points);
Arr_mjd2000_vec = linspace(date2mjd2000(ArrivalEarliest), date2mjd2000(ArrivalLatest), points);
[Dp_grid_mjd2k, Ar_grid_mjd2k] = ndgrid(Dep_mjd2000_vec, Arr_mjd2000_vec);

tempX = Dp_grid_mjd2k(DV_req <= DV_cutoff);
tempY = Ar_grid_mjd2k(DV_req <= DV_cutoff);
ch = convhull(tempX, tempY,'Simplify',true); % Convex Hull Boundary of All the points with DV <= DV_cutoff

pgon = simplify(polyshape(tempX(ch(1:end-1)), tempY(ch(1:end-1))), 'KeepCollinearPoints', false);
[xclim,yclim] = boundingbox(pgon); % Bounding Box X & Y Limits

xcclim = Best_times_opt(1) + [-Synodic/2 Synodic/2]; % Synodic period Bounding Box X-Limits
ycclim = Best_times_opt(2) + [-Synodic/2 Synodic/2]; % Synodic period Bounding Box Y-Limits
%Bounding_Box = [ [xcclim(1) xcclim(1) xcclim(2) xcclim(2) xcclim(1)]' ; [ycclim(1) ycclim(2) ycclim(2) ycclim(1) ycclim(1)]' ];
bounding_rect = [ min(xclim), min(yclim), abs(diff(xclim)), abs(diff(yclim)) ];
rect_x_ratio = abs(Best_times_opt(1) - min(xclim))/abs(diff(xclim));
rect_y_ratio = abs(Best_times_opt(2) - min(yclim))/abs(diff(yclim));
rect_ratio = [rect_x_ratio rect_y_ratio];

Y_Intercept_1 = tangenter(tempX(ch), tempY(ch), 45, Dep_mjd2000_vec, 0); % Y-Intercept of Min. ToF Limit line
Y_Intercept_2 = tangenter(tempX(ch), tempY(ch), 45, Dep_mjd2000_vec, guess_max_tof); % Y-Intercept of Max. ToF Limit line
ToF_limits = [Y_Intercept_1 Y_Intercept_2];

ytt1 = tand(45).*Dep_mjd2000_vec + Y_Intercept_1; % Y-points of Min. TOF Limit Line
ytt2 = tand(45).*Dep_mjd2000_vec + Y_Intercept_2; % Y-points of Max. TOF Limit Line

% Plotting
%contour_ticks = floor(best_DV_opt):0.5:floor(best_DV_opt) + DV_range_plot;
%contour(Dp_grid_mjd2k, Ar_grid_mjd2k, DV_req, contour_ticks, 'ShowText', 'on', 'LineWidth', 2);
hold on
boundary_plot = plot(all_plotting_axes, tempX(ch), tempY(ch), '--y', 'LineWidth', 2); % Boundary
bounding_box_plot = plot(all_plotting_axes, [xclim(1) xclim(1) xclim(2) xclim(2) xclim(1)],[yclim(1) yclim(2) yclim(2) yclim(1) yclim(1)], '-r', 'LineWidth', 2);
%plot(all_plotting_axes, [xcclim(1) xcclim(1) xcclim(2) xcclim(2)
%xcclim(1)],[ycclim(1) ycclim(2) ycclim(2) ycclim(1) ycclim(1)], '-g', 'LineWidth', 2);
ToF_lines_min_plot = plot(all_plotting_axes, Dep_mjd2000_vec, ytt1, '-.r', 'LineWidth', 2);
ToF_lines_max_plot = plot(all_plotting_axes, Dep_mjd2000_vec, ytt2, '-.r', 'LineWidth', 2);
hold off

boundary_plot.Annotation.LegendInformation.IconDisplayStyle = 'off';
bounding_box_plot.Annotation.LegendInformation.IconDisplayStyle = 'off';
ToF_lines_min_plot.Annotation.LegendInformation.IconDisplayStyle = 'off';
ToF_lines_max_plot.Annotation.LegendInformation.IconDisplayStyle = 'off';

xlim(all_plotting_axes, [Dep_mjd2000_vec(1) Dep_mjd2000_vec(end)]);
ylim(all_plotting_axes, [Arr_mjd2000_vec(1) Arr_mjd2000_vec(end)]);

end




function [New_Y_Intercept] = tangenter(shape_X, shape_Y, inclination_deg, X_vals, initial_y_intercept)

ly = tand(inclination_deg).*X_vals + initial_y_intercept;

line = [ reshape(X_vals,[length(X_vals), 1]), reshape(ly,[length(ly), 1]) ];

shape = [ reshape(shape_X,[length(shape_X), 1]), reshape(shape_Y,[length(shape_Y), 1]) ];

for i = 1:length(shape)
    x = shape(i,:);
    for j = 1:length(line)
        y = line(j,:);
        dists(i,j) = norm(x-y);
    end
end

[~, min_dist_id] = min(dists, [], "all");
[min_dist_idx_shape, ~] = ind2sub(size(dists), min_dist_id);
shape_tangent_pt = shape(min_dist_idx_shape,:);

line_Y_at_shape_tangent_pt_X = tand(inclination_deg).*shape_tangent_pt(1) + initial_y_intercept;

New_Y_Intercept = ( shape_tangent_pt(2) - line_Y_at_shape_tangent_pt_X ) + initial_y_intercept;

end


