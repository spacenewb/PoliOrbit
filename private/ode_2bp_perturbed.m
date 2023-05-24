function dy = ode_2bp_perturbed( t, y, mu, pert)

% ode_2bp   ODE system for the two body problem ( Keplerian motion)
%
% PROTOTYPE
%   dy = ode_2bp( t, y, mu, pert)
%
% INPUT:
%   t[1]        Time (can be omitted, as the system is autonomous) [T]
%   y[6x1]      State of the body ( rx , ry , rz , vx , vy , vz ) [ L, L/T ]
%   mu[1]       Gravitational parameter of the primary [ L^3/T^2 ]
%   pert[-]     String describing what perturbation to implement ['J2'/'0']
%
% OUTPUT:
%   dy[6x1]     Derivative of the state [ L/T^2, L/T^3 ]
%
% CALLED FUNCTIONS:
%   astroConstants.m
%
% CONTRIBUTORS:
%   Hariharan Venkatesh Vitaladevuni
%
% VERSIONS:
%   2021-10-03: First version
%
%--------------------------------------------------------------------------

% Position and velocity
r = y(1:3);
v = y(4:6);

% Distance from the primary
rnorm = norm(r);

% Calculate perturbation acceleration
if isempty(pert)

    acc_pert = [0;0;0];

else

    if isfield(pert,'J2')
        J2 = astroConstants(9); % Radius of earth value compatibility check
        R_E = 6378.1363; % [km] Check documentation of astroConstants
        J2_const_1 = 3/2*J2*mu*R_E^2/rnorm^4;
        J2_const_2 = r/rnorm.*[(5*r(3)^2/rnorm^2 - 1); (5*r(3)^2/rnorm^2 - 1); (5*r(3)^2/rnorm^2 - 3)];
        acc_pert_J2 = J2_const_1*J2_const_2;
    else
        acc_pert_J2 = [0;0;0];
    end

    if isfield(pert,'SRP')

        R_P = astroConstants(20 + pert.SRP(4)); % Radius of Planet
        [Kep_Planet, Mu_Sun] = uplanet(pert.SRP(3)+(t/3600/24),pert.SRP(4));
        y_Planet = kep2car(Kep_Planet', Mu_Sun);
        R_Planet = y_Planet(1:3); % Heliocentric Position of Planet
        R_sat = R_Planet + r; % Heliocentric Position of Satellite

        ShadeHalfAngleSun = atan(R_P/norm(R_Planet)); % Angle Subtended By Planet Shade at Sun
        SatAngleSun = real(acos(dot(R_Planet,R_sat)/(norm(R_Planet)*norm(R_sat)))); % Angle Subtended By Satellite and Planet at Sun
        
        if abs(SatAngleSun) <= abs(ShadeHalfAngleSun) || norm(R_sat) >= norm(R_Planet)
            ShadeFunction = 1; % Satellite Illumination Occluded by Planet
        else
            ShadeFunction = 0; % Satellite Illumination Not Occluded by Planet
        end

        acc_mag_pert_SRP = ShadeFunction*astroConstants(31)/astroConstants(5)/1000*pert.SRP(1)*pert.SRP(2);
        acc_pert_SRP = acc_mag_pert_SRP.*(-R_sat./norm(R_sat));

    else
        acc_pert_SRP = [0;0;0];
    end

    acc_pert = acc_pert_SRP + acc_pert_J2;

end

% Set the derivatives of the state
dy = [ v; (-mu/rnorm^3)*r + acc_pert ];

end