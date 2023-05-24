function DV_tot = LambertDV(Times_mjd2000, Dep_bodyID, Arr_bodyID, Mu_Central)

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
[Dep_kep, ~] = uplanet(Dep_mjd2000, Dep_bodyID);
Dep_cart = kep2car(Dep_kep', Mu_Central);

Arr_mjd2000 = Times_mjd2000(2);
[Arr_kep, ~] = uplanet(Arr_mjd2000, Arr_bodyID);
Arr_cart = kep2car(Arr_kep', Mu_Central);

ToF = (Arr_mjd2000 - Dep_mjd2000)*86400;

[~,~,~,~,VI,VF,~,~] = lambertMR( Dep_cart(1:3), Arr_cart(1:3), ToF, Mu_Central, 0, 0, 0 );

Vs_p_d = VI' - Dep_cart(4:6);
Vinf_d = norm(Vs_p_d); 
% V0_d = sqrt(Vinf_d^2 + 2*astroConstants(10+Dep_bodyID)/astroConstants(20+Dep_bodyID));
% DV_d = V0_d - sqrt(astroConstants(10+Dep_bodyID)/astroConstants(20+Dep_bodyID));
DV_d = Vinf_d;

Vs_p_a = VF' - Arr_cart(4:6);
Vinf_a = norm(Vs_p_a); 
% V0_a = sqrt(Vinf_a^2 + 2*astroConstants(10+Arr_bodyID)/astroConstants(20+Arr_bodyID));
% DV_a = V0_a - sqrt(astroConstants(10+Arr_bodyID)/astroConstants(20+Arr_bodyID));
DV_a = Vinf_a;

DV_tot = DV_d + DV_a;

end