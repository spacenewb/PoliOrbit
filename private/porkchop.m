function [DV_tot] = porkchop(Dep_mjd2000_vec, Arr_mjd2000_vec, Dep_bodyID, Arr_bodyID)

Mu_Sun = astroConstants(4);

%% Computation
DV_tot = NaN*ones(Dep_mjd2000_vec, Arr_mjd2000_vec);
DV_departure = NaN*ones(Dep_mjd2000_vec, Arr_mjd2000_vec);
DV_arrival = NaN*ones(Dep_mjd2000_vec, Arr_mjd2000_vec);

for i=1:length(Dep_mjd2000_vec)

    Dep_mjd2000 = Dep_mjd2000_vec(i);

    for j=1:length(Arr_mjd2000_vec)

        Arr_mjd2000 = Arr_mjd2000_vec(j);

        [DV_t, DV_dep, DV_arr]= LambertDV([Dep_mjd2000; Arr_mjd2000], Dep_bodyID, Arr_bodyID, Mu_Sun);

        DV_tot(i,j) = DV_t;
        DV_departure(i,j) = DV_dep;
        DV_arrival(i,j) = DV_arr;

    end
end