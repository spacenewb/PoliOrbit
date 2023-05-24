clc
close all
%%

%ideal_MJDs = [17178.9947, 18034.45499, 18401.8563];
ideal_MJDss = [15398.6, 16435.1, 16735.4];

search_bounding_window_tolerance = 15.*[1, 1, 1]; % days

lb = ideal_MJDss - search_bounding_window_tolerance;
ub = ideal_MJDss + search_bounding_window_tolerance;

options = optimoptions('fmincon','OptimalityTolerance', 1e-7, 'Display', 'off');
[fmin_DV_trip_idx fmin_DV_trip] = fmincon(@(x) DV_MP_PGA(x) ,ideal_MJDss,[],[],[],[],lb,ub,[],options);

%%
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

[DV_ga,r_p,~,~,~,~,~,~] = PGAflyby(VF_1, VI_2, astroConstants(13), GA_cart(4:6)',100,astroConstants(23));
% DV_ga
% r_p
if ~isnan(DV_ga)
    DV_trip = DV_d_1 + DV_ga + DV_a_2;
end

if isnan(DV_ga)
    DV_trip = 4e6;
end

end