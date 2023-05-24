function [ r, v, TOUT ] = orbit_propagation( r0, v0, TSPAN, mu, relTol, absTol, perturbationModel )
% orbit_propagation numerically integrates the 2-body problem equation of motion
%--------------------------------------------------------
% NOTE:
%   Unperturbed case works for all 2-body problems.
%   Perturbed cases are ONLY for geocentric orbits.
%   -> J2-Perturbation
%   -> SRP-Perturbation
%      
%--------------------------------------------------------
% PROTOTYPE:
%   [ r, v, TOUT ] = orbit_propagation( r0, v0, TSPAN, mu, relTol, absTol, perturbationModel )
%
% INPUTS:
%   r0[3x1]                         Initial position vector                  [km]
%   v0[3x1]                         Initial velocity vector                  [km/s]
%   TSPAN[time]                     Time vector for integration              [s]
%   mu[1]                           Planetary constant of planet             [km^3/s^2]   
%   relTol[1]                       Relative tolerance for ode45             [-]
%   absTol[1]                       Absolute tolerance for ode45             [-]
%   perturbationModel [struct]      Choose model from:   
%                                       Unperturbed:
%                                           perturbationModel.type = 'unperturbed'
%                                       J2-perturbed:
%                                           perturbationModel.type = 'J2-perturbed'
%                                       SRP-perturbed: HP: s/c-sun distance is the same as sun-earth
%                                           perturbationModel.type = 'SRP-perturbed'
%                                           perturbationModel.AreaToMass = <s/c area to mass ratio wrt sun>
%                                           perturbationModel.Cr = <optical properties of s/c surface>
%                                           perturbationModel.p_sr = <solar radiation pressure>
%                                           perturbationModel.startDate = <[year,month,day,hour,minute,second]> --> when the propagation starts  
%                                       J2_SRP-perturbed: s/c-sun distance is the same as sun-earth
%                                           perturbationModel.type = 'J2_SRP-perturbed'
%                                           perturbationModel.AreaToMass = <s/c area to mass ratio wrt sun>
%                                           perturbationModel.Cr = <optical properties of s/c surface>
%                                           perturbationModel.p_sr = <solar radiation pressure>
%                                           perturbationModel.startDate = <[year,month,day,hour,minute,second]> --> when the propagation starts

%
% OUTPUTS:
%   r[timex3]                      Timeseries of position vectors            [km]
%   v[timex3]                      Timeseries of velocity vectors            [km/s]
%   TSPAN[time]                    Time instants of r, v                     [s]

%% Check inputs

if length(TSPAN) < 10000
    warning('Use lots of points for long TSPANS')
end

if isequal(size(r0), [1 3])

    r0 = r0';

end
if isequal(size(v0), [1 3])

    v0 = v0';

end

% Set state vector of initial conditions [ r0x; r0y; r0z; v0x; v0y; v0z ]
Y0 = [ r0; v0 ];

%% Integration

% Tolerance of integration
options = odeset( 'RelTol', relTol, 'AbsTol', absTol );

% Choose the perturbations to use -> All ODE functions are defined as local function

switch perturbationModel.type

    case 'unperturbed'

        [ TOUT, YOUT ] = ode45( @( t , y ) ode_2BP( t, y , mu ), TSPAN, Y0, options );

    case 'J2-perturbed'

        R_E = astroConstants(23);
        J2 = astroConstants(9);
        [ TOUT, YOUT ] = ode45( @( t , y ) ode_2BP_perturbed_J2( t, y, mu, R_E, J2 ), TSPAN, Y0, options );

    case 'SRP-perturbed'

        [ TOUT, YOUT ] = ode45( @( t , y ) ode_2BP_perturbed_SRP( t, y, mu, perturbationModel.AreaToMass, perturbationModel.Cr, perturbationModel.p_sr, perturbationModel.startDate ), TSPAN, Y0, options );


    case 'J2_SRP-perturbed'

        R_E = astroConstants(23);
        J2 = astroConstants(9);
        [ TOUT, YOUT ] = ode45( @( t , y ) ode_2BP_perturbed_SRP( t, y, mu, R_E, J2, perturbationModel.AreaToMass, perturbationModel.Cr, perturbationModel.p_sr, perturbationModel.startDate ), TSPAN, Y0, options );

    otherwise

        error( 'Incorrect perturbationModel entered' );

end

% Set r and v from output of integration
r = YOUT( :, 1:3 );
v = YOUT( :, 4:6 );

end

%% Local functions

%--------------------- UNPERTURBED ---------------------

function dy = ode_2BP( ~, y, mu )
%
% ode_2bp ODE system for the two-body integration
%
% PROTOTYPE: 
%   dy = ode_2bp( t, y, mu )
%
% INPUT:
%
%   y   State vector
%   mu  Planetary constant of earth
%
% OUTPUT: 
%   dy  Time derivatives of state vector
% 

% Set r and v from input state vector (unnecessary)
r = y(1:3);
v = y(4:6);

% 2BP equation of motion in state-space
dy = [ v; -mu./norm(r).^3 .* r ]; 

end

%--------------------- J2-PERTURBED ---------------------

function dy = ode_2BP_perturbed_J2( ~, y, mu, R_E, J2 )
%
% ode_Perturbed2bp ODE system for the perturbed two-body problem around a
% GEOCENTRIC orbit
%
% PROTOTYPE: 
%   dy = ode_Perturbed2bp( ~, y, mu, Re, J2 )
%   ~ is the time, it is not needed here
%
% INPUT:
%
%   t       Time-step [s]
%   y[6x1]  State vector ( rx, ry, rz, vx, vy, vz )
%   mu      Planetary constant [km^3/s^2]
%   R_E     Earth's mean radius [km]
%   J2      Second zonal harmonic of zonal variations
% 
% OUTPUT: 
%   dy      Time derivatives of state vector

r = y(1:3);
r_norm = norm(r);
v = y(4:6);

aJ2_vect_ECEI = 3/2 * J2 * mu * R_E^2 / r_norm^5 *       ... % Multiplicative constant
    [ r(1) * ( 5 * (r(3) / r_norm)^2 - 1 );   ... % x-component
      r(2) * ( 5 * (r(3) / r_norm)^2 - 1 );   ... % y-component
      r(3) * ( 5 * (r(3) / r_norm)^2 - 3 )];      % z-component

dy = [ v; -mu ./ r_norm.^3 .* r + aJ2_vect_ECEI];

end

%--------------------- SRP-PERTURBED ---------------------

%--------------------- J2_SRP-PERTURBED ---------------------




