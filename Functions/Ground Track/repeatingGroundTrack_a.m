function a_rgt = repeatingGroundTrack_a( m, k, e, i, perturbationModel )
% repeatingGroundTrack_a Calculates semi-major axis of orbit needed to have
% a repeating ground track (Ground track remains the same)
%
% PROTOTYPE: 
%   a_rgt = repeatingGroundTrack_a( m, k, e, i, perturbationModel )
%
% INPUT:
%   m[1]                [-]     Rotations of the planet after which the ground track repeats itself
%   k[1]                [-]     Revolutions of the satellite after which the ground track repeats itself
%   e[1]                [-]     Eccentricity of orbit
%   i[1]                [rad]   Inclination of orbit
%   perturbationModel   [-]     'unperturbed' or 'J2-perturbed'
%
% OUTPUT: 
%   a      [1x1]    [km]    Required semi-major axis
%

%% Check inputs

if mod(m,1) ~= 0
    
    error('m must be a whole number')
    
end

if mod(k,1) ~= 0
    
    error('k must be a whole number')

end

%% Calculate a

switch perturbationModel

    case 'unperturbed'

        w_earth = 15.04*pi/180/3600; % rad/s

        mu = astroConstants(13);

        T_sat = 2*pi .* m ./ ( w_earth .* k );

        a_rgt = (mu .* T_sat.^2 ./ (4 .* pi.^2)).^(1/3);

    case 'J2-perturbed'

        w_earth = 15.04*pi/180/3600; % rad/s
        mu = astroConstants(13);
        J2 = astroConstants(9);
        Re = astroConstants(23);

        syms a;
        OM_dot = - ( 3/2 .* (sqrt(mu) .* J2 .* Re.^2) ./ ((1-e.^2).^2.*a.^(7/2)) ) .* cos(i);
        om_dot = -( 3/2 .* (sqrt(mu) .* J2 .* Re.^2) ./ ((1-e.^2).^2.*a.^(7/2)) ) .* (5/2 .* sin(i).^2 - 2);
        M0_dot = ( 3/2 .* (sqrt(mu) .* J2 .* Re.^2) ./ ((1-e.^2).^(3/2).*a.^(7/2)) ) .* (1 - 3/2 .* sin(i).^2);
        n = sqrt(mu ./ a.^3);

        FUN = m ./ k -  (w_earth - OM_dot) ./ ( n + om_dot + M0_dot) ;
        FUN = matlabFunction(FUN);

        % Guess is from unperturbed repeating ground track
        T_sat = 2*pi .* m ./ ( w_earth .* k );
        a_guess = (mu .* T_sat.^2 ./ (4 .* pi.^2)).^(1/3);
        
        a_rgt = fzero(FUN, a_guess);

end

end

