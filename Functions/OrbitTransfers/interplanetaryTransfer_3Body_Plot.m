function [] = interplanetaryTransfer_3Body_Plot( departureTime_mjd, flybyTime_mjd, arrivalTime_mjd, planet1, planet1Multiplier, planet2, planet2Multiplier, planet3, planet3Multiplier )

[ r_P1_planet1, ~ ] = uplanet_rv( departureTime_mjd, planet1 ); % Planet1 ephemeris @ P1
[ r_P2_planet1, ~ ] = uplanet_rv( flybyTime_mjd, planet1 );     % Planet1 ephemeris @ P2
[ r_P3_planet1, ~ ] = uplanet_rv( arrivalTime_mjd, planet1 );   % Planet1 ephemeris @ P3

[ r_P1_planet2, ~ ] = uplanet_rv( departureTime_mjd, planet2 ); % Planet2 ephemeris @ P1
[ r_P2_planet2, ~ ] = uplanet_rv( flybyTime_mjd, planet2 );     % Planet2 ephemeris @ P2
[ r_P3_planet2, ~ ] = uplanet_rv( arrivalTime_mjd, planet2 );   % Planet2 ephemeris @ P3

[ r_P1_planet3, ~ ] = uplanet_rv( departureTime_mjd, planet3 ); % Planet3 ephemeris @ P1
[ r_P2_planet3, ~ ] = uplanet_rv( flybyTime_mjd, planet3 );     % Planet3 ephemeris @ P2
[ r_P3_planet3, ~ ] = uplanet_rv( arrivalTime_mjd, planet3 );   % Planet3 ephemeris @ P3

% Plot planets
planet3D(10, [0,0,0], 50); % Sun
planet3D_ADD(planet1, r_P1_planet1', gcf, planet1Multiplier);
planet3D_ADD(planet1, r_P2_planet1', gcf, planet1Multiplier);
planet3D_ADD(planet1, r_P3_planet1', gcf, planet1Multiplier);
planet3D_ADD(planet2, r_P1_planet2', gcf, planet2Multiplier);
planet3D_ADD(planet2, r_P2_planet2', gcf, planet2Multiplier);
planet3D_ADD(planet2, r_P3_planet2', gcf, planet2Multiplier);
planet3D_ADD(planet3, r_P1_planet3', gcf, planet3Multiplier);
planet3D_ADD(planet3, r_P2_planet3', gcf, planet3Multiplier);
planet3D_ADD(planet3, r_P3_planet3', gcf, planet3Multiplier);

% Calculate planet orbits from ephemerides
Mjd_vect = linspace(departureTime_mjd,arrivalTime_mjd, ceil(arrivalTime_mjd-departureTime_mjd));

r_planet1 = [];
r_planet2 = [];
r_planet3 = [];

for i = 1:length(Mjd_vect)

    [ r, v ] = uplanet_rv( Mjd_vect(i), planet1 );
    r_planet1 = [r_planet1; r'];

    [ r, v ] = uplanet_rv( Mjd_vect(i), planet2 );
    r_planet2 = [r_planet2; r'];

    [ r, v ] = uplanet_rv( Mjd_vect(i), planet3 );
    r_planet3 = [r_planet3; r']; 

end

r_planet2_flyby = uplanet_rv( flybyTime_mjd, planet2 );

% Plot planet orbits
planet(1).name = 'Mercury';
planet(2).name = 'Venus';
planet(3).name = 'Earth';
planet(4).name = 'Mars';
planet(5).name = 'Jupiter';
planet(6).name = 'Saturn';
planet(7).name = 'Uranus';
planet(8).name = 'Neptune';
planet(9).name = 'Pluto';

orbitPlotrv( r_planet1 , '-', gcf, [1,0,0], 2.5, 1);
orbitPlotrv( r_planet2 , '-', gcf, [0,1,0], 2.5, 1);
orbitPlotrv( r_planet3 , '-', gcf, [0,0,1], 2.5, 1);


% Plot departure, flyby, arrival points

rx = [r_planet1(1,1);r_planet2_flyby(1);r_planet3(end,1)];
ry = [r_planet1(1,2);r_planet2_flyby(2);r_planet3(end,2)];
rz = [r_planet1(1,3);r_planet2_flyby(3);r_planet3(end,3)];
plot3(rx(1),ry(1),rz(1),'o', 'markeredgecolor', [1,0,0], 'markersize', 40, 'linewidth', 2);
plot3(rx(2),ry(2),rz(2),'o', 'markeredgecolor', [0,1,0], 'markersize', 18, 'linewidth', 2);
plot3(rx(3),ry(3),rz(3),'o', 'markeredgecolor', [0,0,1], 'markersize', 18, 'linewidth', 2);

% Calculate lambert arcs
mu_sun = astroConstants(4);
TOF1 = (flybyTime_mjd - departureTime_mjd) * 24 * 3600;
TOF2 = (arrivalTime_mjd - flybyTime_mjd) * 24 * 3600;
[ ~, ~, ~, ~, V_sc_i_1, V_sc_f_1, ~, ~ ] = lambertMR( r_P1_planet1, r_P2_planet2, TOF1, mu_sun, 0, 0, 0, 1 );
[ ~, ~, ~, ~, V_sc_i_2, V_sc_f_2, ~, ~ ] = lambertMR( r_P2_planet2, r_P3_planet3, TOF2, mu_sun, 0, 0, 0, 1 );

relTol = 1e-16;
absTol = 1e-16;
perturbationModel.type = 'unperturbed';
[ r_arc1, ~, ~ ] = orbit_propagation( r_P1_planet1, V_sc_i_1, linspace(0,TOF1,TOF1/3600), mu_sun, relTol, absTol, perturbationModel );

[ r_arc2, ~, ~ ] = orbit_propagation( r_P2_planet2, V_sc_i_2, linspace(0,TOF2,TOF2/3600), mu_sun, relTol, absTol, perturbationModel );

% Plot lambert arcs
orbitPlotrv( r_arc1 , '--', gcf, [0,128,255]/255, 2.5, 1);
orbitPlotrv( r_arc2 , '--', gcf, [255,179,0]/255, 2.5, 1);

[ dV_sat, ~, ~, ~, ~, ~ ] = interplanetaryTransfer_3Body( departureTime_mjd, flybyTime_mjd, arrivalTime_mjd, planet1, planet2, planet3 );

% Create legend
planet(1).name = 'Mercury';
planet(2).name = 'Venus';
planet(3).name = 'Earth';
planet(4).name = 'Mars';
planet(5).name = 'Jupiter';
planet(6).name = 'Saturn';
planet(7).name = 'Uranus';
planet(8).name = 'Neptune';
planet(9).name = 'Pluto';

plot3(0,0,0,'o', 'markersize', 1, 'color', [1,1,1]); % Add extra legend entry for dV
legend([planet(planet1).name, ' orbit'], [planet(planet2).name, ' orbit'], [planet(planet3).name, ' orbit'], ['Departure: ', mjd20002DateString(departureTime_mjd)], ['Flyby: ', mjd20002DateString(flybyTime_mjd)], ['Arrival: ', mjd20002DateString(arrivalTime_mjd)] , 'First Lambert arc', 'Second Lambert arc', ['$\Delta V = $', num2str(dV_sat), '\ [km/s]'], 'interpreter', 'latex');
hold off
end




