function [ quiver_Handles ] = referenceSystem( origin, A, length, figureNum , color)
% referenceSystem adds reference system to plot figureNum
%       Versors of the system must be in ROWS (direction cosines). If you
%       are rotating a MATRIX, each vector must be a column vector so: 
%       A_ECEI=R_ECEI2pf' * A_pf' (remember the second apex ONLY FOR THE FIRST ROTATION).
%       Once all rotations are complete, the vectors in the matrix need to become row vectors
%       again so: A_ECEI = A_ECEI' is the new director cosine matrix.
%       Remember to only transpose the matrix during first and last
%       rotation, but not for the ones in between.
%
% PROTOTYPE: 
%   [ quiver_Handles ] = referenceSystem( origin, A, length, figureNum , color)
%
% INPUT:
%   origin [3x1]    Origin of reference system [x,y,z]
%   A [3x3]         Cosine director matrix in the ECEI system
%   length            Length of vectors in plot
%   figureNum       Figure handle where to add the reference system
%   color           If specified, normalized RGB color triplet [0.1, .2, .3] 
%                   or string 'different': 'different' generates 3 axis with colors
%                   x:red, y:green, z:blue.
%
% OUTPUTS:
%   quiver_Handles  Returns vector of handles for the quivers
% 
% For x,y,z reference A = [1,0,0;0,1,0;0,0,1];
%

v1 = A(1,:);
v2 = A(2,:);
v3 = A(3,:);

v1 = length * v1 ./ norm(v1);
v2 = length * v2 ./ norm(v2);
v3 = length * v3 ./ norm(v3);

hold on

    quiver_Handles(1) = quiver3(origin(1), origin(2), origin(3), v1(1), v1(2), v1(3), 'color', color, 'linewidth', 2.5 );
    quiver_Handles(2) = quiver3(origin(1), origin(2), origin(3), v2(1), v2(2), v2(3), 'color', color, 'linewidth', 2.5, 'handlevisibility', 'off' );
    quiver_Handles(3) = quiver3(origin(1), origin(2), origin(3), v3(1), v3(2), v3(3), 'color', color, 'linewidth', 2.5, 'handlevisibility', 'off' );

end