function [ norm_a ] = norm_Timedependent( a )
% norm_Timedependent Calculates norm of vector (given as array) for each timestep
%
% PROTOTYPE: 
% [ norm_a ] = norm_Timedependent( a )
%
% INPUT:
% a[timex3]   Vector in each timestep given as array
%
%        x    y    z
% a = [ a11, a12, a13;       Timestep_1
%       a21, a22, a23;       Timestep_2
%       ...  ...  ...;        ...
%       a31, a32, a33 ]      Timestep_n
%
% OUTPUT: 
% norm_a [timex1]  Column vector containing norm of vector in every timestep
%
% norm_a = [norm([a11, a12, a13]);
%           norm([a21, a22, a23]);
%           ... ]; 
%
% CONTRIBUTORS: 
%   Matteo D'Ambrosio, 2021 (MD)

%% Check input

if ~isequal( size(a), [size(a,1), 3])
    error('Dimension of vector input must be [timex3]')
end

%% Calculate norm

norm_a = sqrt(a(:,1).^2 + a(:,2).^2 + a(:,3).^2);

end