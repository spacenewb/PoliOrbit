function [DV_req, min_DV_req, best_dep_mjd2k, best_arr_mjd2k, Best_times_opt, best_DV_opt, porkchop] = Pork_plot(Mu_Central, Dep_bodyID, Dep_date_1, Dep_date_2, Arr_bodyID, Arr_date_1, Arr_date_2, n_samples_axis, DV_contour_range, plot_params)

% Pork_plot     This function returns the minimum delta-V required and the 
%               porkchop plot for an interplanetary transfer between two 
%               bodies orbiting the same central body. This function assumes
%               that the initial and final orbits are the equal to the body
%               heliocentric orbits
%
% PROTOTYPE
%   [DV_tot, min_DV_tot, best_dep_mjd2k, best_arr_mjd2k, Best_times_opt, best_DV_opt] = Pork_plot(Mu_Central, Dep_bodyID, Dep_date_1, Dep_date_2, Arr_bodyID, Arr_date_1, Arr_date_2, n_samples_axis, DV_contour_range, plot_params)
%
% INPUT:
%   Mu_Central[1]           Gravitational parameter of the primary [ L^3/T^2 ]
%   Dep_bodyID[1]           Departure body planetID number [ - ]
%   Dep_date_1[-]           Date: Start of Dep window [yyyy mm dd HH MM SS]
%   Dep_date_2[-]           Date: End of Dep window [yyyy mm dd HH MM SS]
%   Arr_bodyID[1]           Arrival body planetID number   [ - ]
%   Arr_date_1[-]           Date: Start of Arr window [yyyy mm dd HH MM SS]
%   Arr_date_2[-]           Date: End of Arr window [yyyy mm dd HH MM SS]
%   n_samples_axis[1]       Number of points in time axis
%   DV_contour_range[1]     The max. DV in contour plot rel. to DV min. [km/s]
%   plot_params:
%       plot_flag[1]        Flag for plotting porkchop? [1/0]
%       plot_grid_flag[1]   Flag for plotting porkchop? [1/0]
%       DV_req_flag[1]      Which Delta-V is required? [1:Dep,2:Arr,3:Tot]
%       LauncherC3Lim[1]    Launcher C3 Limit [Inf, any number]
%       NumTOFlines[1]      Number of TOF contour lines
%
%   Reference: planetID:    Integer number identifying the celestial body (< 10)
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
%   DV_req[1]               Delta-V requested for the manouvre/s [ L/T ]
%   min_DV_req[1]           Min. Delta-V requested from grid search [ L/T ]
%   best_dep_mjd2k[1]       Best Dep MJD2000 from grid search corrrsponding to min_DV_req
%   best_arr_mjd2k[1]       Best Arr MJD2000 from grid search corrrsponding to min_DV_req
%   Best_times_opt[1x2]     Best Dep & Arr MJD2000 from minimisation around neighbourhood of min_DV_req
%   best_DV_opt[1]          Min. Delta-V requested from minimisation around neighbourhood of min_DV_req [ L/T ]
%   porkchop[]              Figure Handle of the porkchop plot
%
% CALLED FUNCTIONS:
%   astroConstants
%   uplanet
%   kep2car
%   lambertMR
%
% CONTRIBUTORS:
%   Hariharan Venkatesh Vitaladevuni
%
% VERSIONS:
%   2021-10-03: First version
%
%--------------------------------------------------------------------------

%% Initialisation
plot_flag = plot_params(1);
plot_grid_flag = plot_params(2);
DV_req_flag = plot_params(3);
LauncherC3Lim = plot_params(4);
NumTOFlines = plot_params(5);

% Departure Data stored in a structure
Dep.bodyID = Dep_bodyID; 
Dep.date_1 = Dep_date_1;
Dep.date_2 = Dep_date_2;

% Arrival Data stored in a structure
Arr.bodyID = Arr_bodyID; 
Arr.date_1 = Arr_date_1; 
Arr.date_2 = Arr_date_2;

%% Dep & ToF combinations
% Converting Inputs to relevant format
Dep.mjd2000_1 = date2mjd2000(Dep.date_1); 
Dep.mjd2000_2 = date2mjd2000(Dep.date_2); 
Dep.mjd2000_vec = linspace(Dep.mjd2000_1, Dep.mjd2000_2, n_samples_axis);
Arr.mjd2000_1 = date2mjd2000(Arr.date_1);
Arr.mjd2000_2 = date2mjd2000(Arr.date_2);
Arr.mjd2000_vec = linspace(Arr.mjd2000_1, Arr.mjd2000_2, n_samples_axis);

Dep.WinLength_S = (Dep.mjd2000_2 - Dep.mjd2000_1)*86400;
Arr.WinLength_S = (Arr.mjd2000_2 - Arr.mjd2000_1)*86400;

% Generating times grid
[Dp_grid_mjd2k, Ar_grid_mjd2k] = ndgrid(Dep.mjd2000_vec, Arr.mjd2000_vec);
ToF_grid = (Ar_grid_mjd2k - Dp_grid_mjd2k); % Time of Flight in days [s];
% ToF_s = ToF_grid.*86400; % Time of Flight in [s];

%% Computation
DV_req = NaN.*ones(size(Dp_grid_mjd2k));

for i=1:length(Dep.mjd2000_vec)

    Dep_mjd2000 = Dep.mjd2000_vec(i);

    for j=1:length(Arr.mjd2000_vec)

        Arr_mjd2000 = Arr.mjd2000_vec(j);

        DV_req(i,j) = LambertDV([Dep_mjd2000; Arr_mjd2000], Dep.bodyID, Arr.bodyID, Mu_Central, DV_req_flag, LauncherC3Lim);

    end
end

%% Finding Optimum
% Best grid solution
[min_DV_req, min_DV_tot_idx] = min(DV_req, [], 'all');
[min_DV_tot_idx_dep, min_DV_tot_idx_arr] = ind2sub(size(DV_req), min_DV_tot_idx);
best_dep_mjd2k = Dp_grid_mjd2k(min_DV_tot_idx_dep, min_DV_tot_idx_arr); % Minima X Index (Grid Eval)
best_arr_mjd2k = Ar_grid_mjd2k(min_DV_tot_idx_dep, min_DV_tot_idx_arr); % Minima Y Index (Grid Eval)
best_dep_date = mjd20002date(best_dep_mjd2k);
best_arr_date = mjd20002date(best_arr_mjd2k);

% Finding unconstrained minima based on grid evaluations
fminunc_options = optimoptions('fminunc', 'MaxIterations', 200, 'OptimalityTolerance',1e-13,...
    'FiniteDifferenceType', 'central', 'StepTolerance', 1e-13, 'Display', 'off');

[Best_times_opt, best_DV_opt] = fminunc(@(Times_optim) LambertDV(Times_optim, Dep.bodyID, Arr.bodyID, Mu_Central, DV_req_flag, LauncherC3Lim),...
                                [best_dep_mjd2k;best_arr_mjd2k], fminunc_options);

best_dep_date_opt = mjd20002date(Best_times_opt(1));
best_arr_date_opt = mjd20002date(Best_times_opt(2));

%%  Plotting
if plot_flag
% datenum, datetick, xtickangle;

% Setup Plotting Environment
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

contour_ticks = floor(best_DV_opt):1:floor(best_DV_opt)+DV_contour_range;
planetNames = {'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto', 'Sun'};

% clf
porkchop = figure();
set(gca,'Color','w')
hold on
scatter(best_dep_mjd2k, best_arr_mjd2k, 50, 'Xg','LineWidth', 2);
scatter(Best_times_opt(1), Best_times_opt(2), 50, 'Or','LineWidth', 2);
if plot_grid_flag == 1
    plot(Dp_grid_mjd2k, Ar_grid_mjd2k, '.k');
end
[C,h] = contour(Dp_grid_mjd2k, Ar_grid_mjd2k, DV_req, contour_ticks, 'ShowText', 'on', 'LineWidth', 2);
cb = colorbar('Ticks',contour_ticks);
caxis([contour_ticks(1),contour_ticks(end)])
cb.Label.String = 'Delta-V [km/s]';
cb.Label.Interpreter = 'latex';
contour(Dp_grid_mjd2k, Ar_grid_mjd2k, ToF_grid, floor(linspace(ceil(min(ToF_grid, [], 'all')),floor(max(ToF_grid, [], 'all')),NumTOFlines)),'--k','ShowText', 'on')
hold off
colormap("winter")
grid on;
clabel(C,h,'Color','k');

if plot_grid_flag == 1
    lgnd_txt = [strcat('Grid Optimum:', {''}, num2str(min_DV_req), ' km/s'),...
                strcat('Solved Optimum:', {''}, num2str(best_DV_opt), ' km/s'),...
                'Grid Points'];
else
    lgnd_txt = [strcat('Grid Optimum:', {''}, num2str(min_DV_req), ' km/s'),...
                strcat('Solved Optimum:', {''}, num2str(best_DV_opt), ' km/s')];
end

legend(lgnd_txt,'Location','best');

text_offset = 0.05;

text_ann_1 = ['Dep:', char(datetime(best_dep_date_opt(1:3),'InputFormat','yyyy-MM-dd-HH-mm-ss' ))
    'Arr:', char(datetime(best_arr_date_opt(1:3),'InputFormat','yyyy-MM-dd-HH-mm-ss' ))];
text(Best_times_opt(1)+(Dep.mjd2000_vec(end)-Dep.mjd2000_vec(1))*text_offset,Best_times_opt(2),text_ann_1,...
                        'Color','red','FontSize',10,'Interpreter','latex','HorizontalAlignment', 'left');

text_ann_2 = ['Dep:', char(datetime(best_dep_date(1:3),'InputFormat','yyyy-MM-dd-HH-mm-ss' ))
    'Arr:', char(datetime(best_arr_date(1:3),'InputFormat','yyyy-MM-dd-HH-mm-ss' ))];
text(best_dep_mjd2k-(Dep.mjd2000_vec(end)-Dep.mjd2000_vec(1))*text_offset,best_arr_mjd2k,text_ann_2,...
                        'Color','g','FontSize',10,'Interpreter','latex','HorizontalAlignment', 'right');

title('Porkchop Plot - Interplanetary Transfer');
subtitle(strcat('Lambert Solver - [', planetNames{Dep.bodyID}, ' to',{' '}, planetNames{Arr.bodyID},']'));
axis equal;
n_ticks = 6;
xrt = linspace(Dep.mjd2000_vec(1),Dep.mjd2000_vec(end),n_ticks); xticks(xrt)
n_xt = cellfun(@mjd20002date, num2cell(xrt), 'UniformOutput',false);
n_xt = cellfun(@(x) strcat(num2str(x(1)),'-',num2str(x(2)),'-',num2str(x(3))), n_xt, 'UniformOutput',false);

yrt = linspace(Arr.mjd2000_vec(1),Arr.mjd2000_vec(end),n_ticks); yticks(yrt)
n_yt = cellfun(@mjd20002date, num2cell(yrt), 'UniformOutput',false);
n_yt = cellfun(@(y) strcat(num2str(y(1)),'-',num2str(y(2)),'-',num2str(y(3))), n_yt, 'UniformOutput',false);

xticklabels(cellfun(@(a) char(datetime(a,'InputFormat','yyyy-MM-dd'), 'yyyy-MMM-dd'), n_xt, 'UniformOutput', false))
yticklabels(cellfun(@(a) char(datetime(a,'InputFormat','yyyy-MM-dd'), 'yyyy-MMM-dd'), n_yt, 'UniformOutput', false))

% Reset Plotting Environment
set(groot, 'defaultLegendInterpreter','tex');
set(groot, 'defaultAxesTickLabelInterpreter','tex');
set(groot, 'defaultTextInterpreter','tex');

else
    porkchop = '';
end

end


%% Functions
function DV_tot = LambertDV(Times_mjd2000, Dep_bodyID, Arr_bodyID, Mu_Central, DV_req_flag, C3)

% LambertDV     This function returns the delta-V required for an interplanetary transfer
%               between two bodies orbiting the same central body. This function assumes
%               that the initial and final parking orbits are the equal to the body
%               surface altitude
%
% PROTOTYPE
%   DV_tot = LambertDV(Times_mjd2000, Dep_bodyID, Arr_bodyID, Mu_Central)
%
% INPUT:
%   Times_mjd2000[2x1]      [Dep_mjd2000; Arr_mjd2000] - Vector   [MJD2000]
%                           Dep_mjd2000[1]:     Departure time    [MJD2000]
%                           Arr_mjd2000[1]:     Arrival time      [MJD2000]
%   Mu_Central[1]           Gravitational parameter of the primary [ L^3/T^2 ]
%   Dep_bodyID[-]           Departure body planetID number [ - ]
%   Arr_bodyID[-]           Arrival body planetID number   [ - ]
%   DV_req_flag[-]          Which Delta-V is required? [1:Dep,2:Arr,3:Tot]
%
%
%   Reference: planetID:    Integer number identifying the celestial body (< 10)
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
%   DV_tot[1]           Total Delta-V for the transfer manouvre [ L/T ]
%
% CALLED FUNCTIONS:
%   astroConstants
%   uplanet
%   kep2car
%   lambertMR
%
% CONTRIBUTORS:
%   Hariharan Venkatesh Vitaladevuni
%
% VERSIONS:
%   2021-10-03: First version
%
%--------------------------------------------------------------------------
Dep_mjd2000 = Times_mjd2000(1);
Arr_mjd2000 = Times_mjd2000(2);

ToF = (Arr_mjd2000 - Dep_mjd2000)*86400;

if ToF > 0

    [Dep_kep, ~] = uplanet(Dep_mjd2000, Dep_bodyID);
    Dep_cart = kep2car(Dep_kep', Mu_Central);
    
    [Arr_kep, ~] = uplanet(Arr_mjd2000, Arr_bodyID);
    Arr_cart = kep2car(Arr_kep', Mu_Central);
    
    [~,~,~,~,VI,VF,~,~] = lambertMR( Dep_cart(1:3), Arr_cart(1:3), ToF, Mu_Central, 0, 0, 0 );
    
    Vs_p_d = VI' - Dep_cart(4:6);
    Vinf_d = norm(Vs_p_d); 
    DV_d = Vinf_d;
    
    Vs_p_a = VF' - Arr_cart(4:6);
    Vinf_a = norm(Vs_p_a); 
    DV_a = Vinf_a;
    
    switch DV_req_flag
        case 1
            DV_tot = DV_d;
        case 2
            DV_tot = DV_a;
        case 3
            DV_tot = DV_d + DV_a;
    end

    if sqrt(C3) < DV_d
    DV_tot = NaN;
    end

else
    DV_tot = NaN;
end

end