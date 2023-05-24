function S_ref_2 = Synodic_period(planet_ref_ID, planet_2_ID)

ref_mjd2000 = date2mjd2000([2003, 4, 1, 0, 0, 0]); % Any random date

[kep_ref, Mu_Sun]  = uplanet(ref_mjd2000, planet_ref_ID);
[kep_2, ~]       = uplanet(ref_mjd2000, planet_2_ID);

s2yr = 1/(astroConstants(32)*24*60*60);

T_ref = 2*pi*sqrt(kep_ref(1)^3/Mu_Sun)*s2yr;
T_2 = 2*pi*sqrt(kep_2(1)^3/Mu_Sun)*s2yr;

S_ref_2 = ((1/T_ref - 1/T_2)^-1)/s2yr;
S_ref_2 = abs(S_ref_2);

end

