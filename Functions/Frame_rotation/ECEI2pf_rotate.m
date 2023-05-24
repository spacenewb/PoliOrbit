function [  R_ECEI2pf, v_pf ] = ECEI2pf_rotate( i, OM, om, v_ECEI )
% ECEI2pf_rotate generates rotation matrix R_ECEI2pf from Earth-centered
% equatorial intertial frame ECEI{i,j,k} to perifocal frame pf{e,p,h} such 
% that: v_pf = R_ECEI2pf * v_ECEI
% If vector v_ECEI is specified, v_pf is also outputted.
%
% PROTOTYPE: 
%   [  R_ECEI2pf, v_pf ] = ECEI2pf_rotate( i, OM, om, v_ECEI )
%
% INPUT:
%   i         [1x1]                 Inclination of the orbit                      [rad]
%   OM        [1x1]                 Right ascension of the ascending node         [rad]
%   om        [1x1]                 Argument of pericenter                        [rad]
%   v_ECEI    [3x1] or [timex3]     If specified, vector in the ECEI frame        []
%
% OUTPUT: 
%   R_ECEI2pf [3x3]                 Rotation matrix from ECEI to pf frame such that:
%                                           v_pf = R_ECEI2pf * v_ECEI
%   v_pf      [3x1] or [timex3]     If specified, vector in the perifocal frame   []
%
% CONTRIBUTORS:
%   Matteo D'Ambrosio

%% Check inputs

if nargin == 4
    if isequal(size(v_ECEI), [1,3]) % Make v_ECEI a column vector

        v_ECEI = v_ECEI';

    elseif ~isequal(size(v_ECEI,2), 3)

        error('If you input a timeseries, vector must have dimensions [timex3]')

    end

    if isequal(size(v_ECEI,2),1)

        type = 'single';

    else

        type = 'timeseries';

    end

end



%% Calculate rotation matrix

R_OM  = [  cos(OM), sin(OM), 0       ;
          -sin(OM), cos(OM), 0       ;
           0,       0,       1      ];

R_i  = [   1,       0,       0       ;
           0,       cos(i),  sin(i)  ;
           0,      -sin(i),  cos(i) ];

R_om = [   cos(om), sin(om), 0                     
          -sin(om), cos(om), 0
           0,       0,       1,     ];

R_ECEI2pf = R_om*R_i*R_OM;

if nargin == 4

switch type

    case 'single'

    v_pf = R_ECEI2pf * v_ECEI;

    case 'timeseries'

        v_pf = zeros(size(v_ECEI));

        for i = 1:size(v_ECEI,1)
            
            v_pf(i,:) = (R_ECEI2pf*v_ECEI(i,:)')';

        end

end

end

end