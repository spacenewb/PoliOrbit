function [ r_ECEI ] = HECI2ECEI( r_HECI )

eps = astroConstants(8); % Obliquity angle

% Rotation of -eps around vernal equinox
R = [ 1, 0,         0         ; ...
      0, cos(eps), -sin(eps)  ; ...
      0, sin(eps),  cos(eps) ];

r_ECEI = R * r_HECI;

end