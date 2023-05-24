
%% FMINCON
ideal_Dep_MJD = 17179;
ideal_Arr_MJD = 18402;
ideal_GA_MJD_vec = [17629, 18037, 18838, 19239.1];

num_possible_global_opt = length(ideal_GA_MJD_vec);

%DV_trip = DV_MP_PGA(ideal_Dep_MJD,ideal_Arr_MJD,ideal_GA_MJD);
fmin_DV_trip_idx = NaN.*ones(num_possible_global_opt,3);
fmin_DV_trip = NaN.*ones(num_possible_global_opt,1);

for i = 1:num_possible_global_opt

    ideal_MJDs = [ideal_Dep_MJD, ideal_GA_MJD, ideal_Arr_MJD];
    search_bounding_window_tolerance = 15.*[1, 1, 1]; % days
    
    lb = ideal_MJDs - search_bounding_window_tolerance;
    ub = ideal_MJDs + search_bounding_window_tolerance;
    
    options = optimoptions('fmincon','MaxIterations', 10000, 'Display', 'off');
    [fmin_DV_trip_idx(i,:), fmin_DV_trip(i)] = fmincon(@DV_MP_PGA,ideal_MJDs,[],[],[],[],lb,ub,[],options);

end

[Global_optimum_DV, Global_opt_idx] = min(fmin_DV_trip(:));
Global_optimum_MJDs = fmin_DV_trip_idx(Global_opt_idx,:)

Global_optimum_ToF_1 = (Global_optimum_MJDs(2)-Global_optimum_MJDs(1));
Global_optimum_ToF_2 = (Global_optimum_MJDs(3)-Global_optimum_MJDs(2));

Global_optimum_ToF = (Global_optimum_ToF_1 + Global_optimum_ToF_2)/365.25
Global_optimum_DV


%% Functions

function DV_trip = DV_MP_PGA(ideal_MJDs)

ideal_Dep_MJD = ideal_MJDs(1);
ideal_GA_MJD = ideal_MJDs(2);
ideal_Arr_MJD = ideal_MJDs(3);

Mu_Central = astroConstants(4);

DV_trip = NaN;

[GA_kep, ~] = uplanet(ideal_GA_MJD, 3);
GA_cart = kep2car(GA_kep', Mu_Central);

[Dep_kep, ~] = uplanet(ideal_Dep_MJD, 5);
Dep_cart = kep2car(Dep_kep', Mu_Central);
ToF_1 = (ideal_GA_MJD - ideal_Dep_MJD)*86400;

[Arr_kep, ~] = uplanet(ideal_Arr_MJD, 2);
Arr_cart = kep2car(Arr_kep', Mu_Central);
ToF_2 = (ideal_Arr_MJD - ideal_GA_MJD)*86400;

[~,~,~,~,VI_1,VF_1,~,~] = lambertMR( Dep_cart(1:3), GA_cart(1:3), ToF_1, Mu_Central, 0, 0, 0, 0 );
DV_d_1 = norm(VI_1' - Dep_cart(4:6));

[~,~,~,~,VI_2,VF_2,~,~] = lambertMR( GA_cart(1:3), Arr_cart(1:3), ToF_2, Mu_Central, 0, 0, 0, 0 );
DV_a_2 = norm(VF_2' - Arr_cart(4:6));

[DV_ga,~,~,~,~,~,~,~] = PGAflyby(VF_1, VI_2, astroConstants(13), GA_cart(4:6)',100,astroConstants(23));

if ~isnan(DV_ga)
    DV_trip = DV_d_1 + DV_ga + DV_a_2;
end

end