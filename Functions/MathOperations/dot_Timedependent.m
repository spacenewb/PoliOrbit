function [ a_dot_b ] = dot_Timedependent( a, b )
% dot_Timedependent Calculates dot product of vectors (given as array) for each timestep
%
% PROTOTYPE: 
% [ a_dot_b ] = dot_Timedependent( a, b )
%
% INPUT:
% a[timex3]   Vector in each timestep given as array
% b[timex3]   Vector in each timestep given as array
%
% Example:
%        x    y    z
% a = [ a11, a12, a13;       Timestep_1
%       a21, a22, a23;       Timestep_2
%       ...  ...  ...;        ...
%       a31, a32, a33 ]      Timestep_n
%
% OUTPUT: 
% a_dot_b [timex1]  Column vector containing dot product of a and b in each timestep
%
% a_dot_b = [ dot(a(1,:), b(1,:)) ]
%             dot(a(2,:), b(2,:)) 
%             ...                 ];
%
% CONTRIBUTORS: 
%   Matteo D'Ambrosio, 2021 (MD)

%% Check input

if ~isequal( size(a), size(b))
    error('Dimension of vector inputs must be the same')
end

%% Calculate norm

a_dot_b = [];

for i = 1:size(a,1)
    a_dot_b = [a_dot_b; dot(a(i,:), b(i,:))];
end

end