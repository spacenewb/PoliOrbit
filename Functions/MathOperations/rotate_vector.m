function v_rotated  = rotate_vector( v, u, delta)
% rotate_vector rotates vector v around axis u of angle delta
% Vectors must all be in the same reference frame, delta is a rotation
% following the right-hand-rule around axis u (counter-clockwise)
%
% PROTOTYPE:
%   v_rotated  = rotate_vector( v, u, delta)
%
% INPUTS:
%   v[3]        [-]        Vector to rotate
%   u[3]        [-]        Rotation axis expressed as vector
%   delta[1]    [rad]      Angle of rotation around axis u
%
% OUTPUTS:
% v_rotated[3]  [-]        Vector after rotation

if size(u) == [1,3]

    u = u';

end

u = u ./ norm(u); % Normalize rotation axis
v_rotated = v .* cos(delta) + cross( u, v ) .* sin(delta) + u .* dot( u, v ) .* ( 1 - cos(delta) );


end