function DKep = gauss_planetary_equations( t, Kep, Date_initial )

mu = astroConstants(13);

a = Kep(1);
e = Kep(2);
i = Kep(3);
OM = Kep(4);
om = Kep(5);
theta = Kep(6);

[ r_sc_ECEI, ~ ] = kp2rv( a, e, i, OM, om, theta, mu );
[ a_r, a_t, a_h ] = a_perturbing( t, a, e, i, OM, om, theta, r_sc_ECEI, Date_initial );

b  = a * sqrt( 1 - e.^2 );
p = b.^2 / a;
h = sqrt( p * mu );
r = p / ( 1 + e * cos( theta ) );

Da = 2 * a.^2 / h * ( e * sin( theta ) * a_r + p / r * a_t );
De = 1 / h * ( p * sin( theta ) * a_r + ( ( p + r ) * cos( theta ) + r * e ) * a_t ); 
Di = r * cos( theta + om ) / h * a_h;
DOM = r * sin( theta + om ) / ( h * sin(i) ) * a_h;
Dom = 1 / ( h * e ) * ( -p * cos( theta ) * a_r + ( p + r ) * sin( theta ) * a_t ) - r * sin( theta + om ) * cos( i ) / ( h * sin(i) ) * a_h;
Dtheta = h / r.^2 + 1 / ( e * h ) * (p * cos( theta ) * a_r - ( p + r ) * sin( theta ) * a_t ); % Wrong Formula Corrected

DKep = [Da;De;Di;DOM;Dom;Dtheta];

end

%%

function [ a_r, a_t, a_h ] = a_perturbing( t, a, e, i, OM, om, theta, r_sc_ECEI, Date_initial )

mu = astroConstants(13);
J2 = astroConstants(9);
R_E = astroConstants(23);
AU = astroConstants(2);
Cr = 1;
A_M = 5 * 1e-6; % Converted from [ m^2/Kg ] to [ Km^2/Kg ]
p_sr = 4.5e-6 * 1e3; % Converted from [ N/m^2 ] to [ (Kg Km) / ( Km^2 s^2 ) ]
Time_actual = date2mjd2000( Date_initial ) + t / ( 24 * 3600 );

r_sun2earth_HECI = uplanet_rv( Time_actual, 3 );
r_sun2earth_ECEI = HECI2ECEI(r_sun2earth_HECI);
r_earth2sun_ECEI = -r_sun2earth_ECEI;
r_sc2sun_ECEI = r_earth2sun_ECEI - r_sc_ECEI; 
r_sc2sun_LVLH = ECEI2LVLH( i, OM, om, theta, r_sc2sun_ECEI );
r = norm(r_sc_ECEI);

aSRP = p_sr * AU.^2 / norm( r_sc2sun_LVLH ).^2 * Cr * A_M;
aSRP_vect = -aSRP * r_sc2sun_LVLH ./ norm( r_sc2sun_LVLH );

aJ2_vect = -3/2 * J2 * mu * R_E.^2 / r.^4 * ...
    [ 1 - 3 * sin(i).^2 * sin( theta + om ).^2; ...
    sin(i).^2 * sin( 2 * (theta + om) ); ...
    sin( 2 * i ) * sin( theta + om ) ];

a_tot = aJ2_vect + aSRP_vect;

a_r = a_tot(1);
a_t = a_tot(2);
a_h = a_tot(3);

end