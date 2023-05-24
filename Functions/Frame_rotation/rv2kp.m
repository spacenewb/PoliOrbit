function [ a, e_norm, i, OM, om, theta ] = rv2kp( r_ECI, v_ECI , mu )
% rv2kp Calculates Keplerian parameters from input vectors position (r) and velocity (v)
%
% PROTOTYPE:
%   [ a, e, i, OM, om, theta ] = rv2kp( r, v , mu )
%
% INPUT:     SIZE:        DESCRIPTION:                      UNITS:
%   r          [3x1]        Position vector in ECI frame    [km]
%   v          [3x1]        Velocity vector in ECI frame    [km/s]
%   mu         [1x1]        Planetary constant              [km^3/s^2]
%
% OUTPUT:     SIZE:        DESCRIPTION:                     UNITS:
%   a          [1x1]        Semi-axis major                 [km]
%   e          [1x1]        Eccentricity                    [-]
%   i          [1x1]        Inclination                     [rad]
%   OM         [1x1]        RAAN                            [rad]
%   om         [1x1]        Pericenter anomaly              [rad]
%   theta      [1x1]        True anomaly                    [rad]
%
% CONTRIBUTORS:
%   Matteo D'Ambrosio

%% Check inputs

if length(r_ECI)~=3
    error('rv2kp: r must be a vector')
end

if length(v_ECI)~=3
    error('rv2kp: v must be a vector')
end

%% Calculate parameters

% This tolerance is required to have accurate results from MATLAB
tol = .2e-15;

% Define xyz axis
x = [ 1, 0, 0 ];
y = [ 0, 1, 0 ];
z = [ 0, 0, 1 ];

r_norm = norm(r_ECI);

% Calculate semi-major axis
a = 1 ./ ( 2 ./ r_norm - norm(v_ECI).^2 ./ mu );

% Calculate angular momentum vector
h_ECI = cross( r_ECI, v_ECI );

% Calculate eccentricity vector
e = cross( v_ECI, h_ECI ) ./ mu - r_ECI ./ r_norm;
e_norm = norm(e);

% Calculate inclination
i = acos( dot( h_ECI, z ) ./ norm(h_ECI) ); % 0<i<pi

if i - tol < 0

    % RAAN is undefined if i = 0 ( i < tol )
    warning('rv2kp: Setting OM = 0 since i = 0')
    OM = 0;
    AN = [1;0;0];

else

    % RAAN is defined

    % Calculate ascending node (AN)
    AN = cross( z, h_ECI ) ./ norm( cross( z, h_ECI ) ); % Ascending node unit vector

    % Solve OM angle ambiguity
    if dot( AN, y ) > 0 % 0<OM<2pi

        OM = acos( dot( x, AN ) );

    elseif dot( AN, y ) < 0 % 0<OM<2pi

        OM = 2*pi - acos( dot( x, AN ) );

    else % dot( AN, y ) = 0

        % In this case, OM could be 0 or 180 ???? -> Check if it can
        % actually be 180 because idk
        OM = 0;

    end

end

if e_norm - tol > 0

    % om is defined

    % Solve om angular ambiguity
    if dot( z, e ) >= 0

        om = acos( dot(AN, e) ./ e_norm );

    elseif dot( z, e ) < 0

        om = 2*pi - acos( dot(AN, e) ./ e_norm );

    end

else

    % om is undefined if e = 0 ( e < tol )
    warning('rv2kp: Setting om = 0 since e = 0')
    om = 0;

end

if e_norm - tol > 0

    % theta can be calculated from eccentricity vector if e_norm > tol

    % Solve theta angular ambiguity
    if dot( r_ECI, v_ECI ) >= 0

        theta = acos( dot( r_ECI, e ) ./ ( r_norm .* e_norm ) );

    elseif dot( r_ECI, v_ECI ) < 0

        theta = 2*pi - acos( dot( r_ECI, e ) ./ ( r_norm .* e_norm ) );

    end

else

    % theta cannot be calculated from eccentricity vector if e_norm < tol

    % Set theta to start from x-axis and solve for angular ambiguity
    if dot( r_ECI, v_ECI ) >= 0

        theta = acos( dot( r_ECI, [1,0,0] ) ./ ( r_norm ) );

    elseif dot( r_ECI, v_ECI ) < 0

        theta = 2*pi - acos( dot( r_ECI, [1,0,0] ) ./ ( r_norm ) );

    end

    % If you arrive in this if, e_norm < tol -> it was actually 0
    e_norm = 0;

end

end

