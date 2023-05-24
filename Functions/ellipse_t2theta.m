function theta = ellipse_t2theta( e, a, mu, t, t0, theta0 )
% ellipse_t2theta Calculates true anomaly on orbit corresponding to input times t
%
% PROTOTYPE: 
%   [theta] = ellipse_t2theta( e, a, mu, t, t0, theta0_radians )
%
% INPUT:         Units:      Description:
%   e[1x1]         [-]         Eccentricity of orbit
%   a[1x1]         [km]        Semi-axis major of orbit
%   mu[1x1]        [km^3/s^2]  Planetary constant
%   t[timex1]      [s]         Timeseries in which to evaluate theta
%   t0[1x1]        [s]         If specified, initial reference time on orbit (must be specified if theta0 is specified)
%   theta0[1x1]    [rad]       If speficied, initial true anomaly on orbit (only needed if reference time t0 is not given wrt pericenter)
%                              If theta0 is not specified, assumes theta0 = 0 (t0 is reference time at pericenter)
%
% OUTPUT: 
%   theta[timex1]  [rad]
% 

%% Check inputs

if nargin == 5
   theta0 = 0;
end

if nargin == 4
   t0 = 0;
   theta0 = 0;
end

if ~isequal(size(e), [1,1])
   error('e must be a constant') 
end

if ~isequal(size(a), [1,1])
   error('a must be a constant') 
end
if ~isequal(size(mu), [1,1])
   error('mu must be a constant') 
end

if isequal(size(t), [1,1])
    error('t must be a vector - Example: t = 0:T')
end
    
if size(t,1) == 1 % Makes t a row vector
    % do nothing
elseif size(t,2) == 1
    t = t';
end

if t(1) < t0
   error('The first element of t must be such that: t(1) >= t0') 
end

%% Parameters

% Find M0 from theta0
E0 = 2*atan( sqrt( (1-e) ./ (1+e)) .* tan( theta0 ./ 2 ) );
M0 = E0 - e .* sin(E0);

% Calculate M
M = M0 + sqrt( mu ./ a.^3 ) .* (t-t0);
M = wrapTo2Pi(M); % This way the solutions E that we find are 0 < E < 2pi, we have to fix these angles later

% Calculate E with fsolve
E_guess = M + e .* sin(M) ./ ( 1 - sin(M+e) + sin(M) );
FUN = @(E) M - E + e .* sin(E);

options = optimoptions( 'fsolve', 'functiontolerance', 1e-16 );
E = fsolve(FUN, E_guess, options);

theta = rad2deg(2*atan( sqrt( (1+e) ./ (1-e) ) .* tan( E ./ 2 ))); % True anomaly mapped to 0 < theta < 2pi --> Does not use atan2 since 
theta(1) = rad2deg(theta0); % Sets initial true anomaly to initial condition f0, so the rest get fixed accordingly (otherwise for angles pi<f0<2pi the plots start from -180

%% Fix true anomaly in case of n revolutions

for i = 1:length(theta) % Since the true anomaly is a continuous function, this is used to fix the values theta to make it continuous --> without this it would jump from -180 < theta < 180
    
    if i~=length(theta)
        
        while theta(i) > theta(i+1) 
            theta(i+1) = theta(i+1) + 360;
        end
        
    end
    
end

while theta(end) < theta(end-1)
   theta(end) = theta(end)+360; 
end

theta = deg2rad(theta);

end