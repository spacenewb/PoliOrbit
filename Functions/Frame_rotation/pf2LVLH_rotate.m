function [  R_pf2LVLH, v_LVLH ] = pf2LVLH_rotate( theta, v_pf )
% pf2LVLH_rotate Generates rotation matrix R_pf2LVLH from perifocal frame
% pf{e,p,h} to Local-vertical Local-horizontal LVLH{r,theta,h}
% In v_pf is specified, v_LVLH is also outputted.
%
% PROTOTYPE: 
%   [  R_pf2LVLH, v_LVLH ] = pf2LVLH_rotate( theta, v_pf )
%
% INPUT:
%   theta     [1x1]   True anomaly                                     [rad]
%   v_pf      [3x1]   If specified, vector in the pf frame             []
%
% OUTPUT: 
%   R_pf2LVLH [3x3]   Rotation matrix from pf to LVLH frame such that:
%                           v_LVLH = R_pf2LVLH * v_pf
%   v_LVLH    [3x1]   If specified, vector in the perifocal frame      []
%
% CONTRIBUTORS:
%   Matteo D'Ambrosio

%% Check inputs

if nargin == 2

    if isequal(size(v_pf), [1,3])

        v_pf = v_pf';

    end
end
%% Calculate rotation matrix

R_pf2LVLH = [   cos(theta), sin(theta),  0   ;
               -sin(theta), cos(theta),  0   ;
                0,          0,           1  ];

if nargin == 2

    v_LVLH = R_pf2LVLH * v_pf;

end

end