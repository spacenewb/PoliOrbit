function t = ellipse_theta2t( e, a, mu, theta, tp )

% ellipse_theta2t: 
%       If theta is a 2-element vector [theta_start, theta_end], calculates time passed
%       between theta_start and theta_end. If theta is a scalar, calculates time
%       passed from pericenter to theta.
%       *** Note: Function only considers one revolution. For repeated
%       revolutions, add n periods to the answer: t = t + n*T.
%
% PROTOTYPE: 
%   t = ellipse_theta2t( e, a, mu, theta_radians, tp_optional )
%
% INPUT:
%   e[1x1]             [-]         Eccentricity of orbit
%   a[1x1]             [km]        Semi-axis major of orbit
%   mu[1x1]            [km^3/s^2]  Planetary constant
%   theta[1x2 OR 1x1]  [rad]       If 2 angles are specified (ex: [theta_start, theta_end]), 
%                                  calculates the deltaT passed between the 2 points. If one
%                                  angle is specified, calculates time passed from pericenter
%                                  to that point.
%   tp[1x1]            [s]         If specified, reference time at pericenter
%
% OUTPUT: 
%   t[1x1]             [s]         If 2 angles specified: Time passed from theta_start to theta_end
%                                  If 1 angle specified: Time from pericenter to theta
% 
% CONTRIBUTORS: 
%   Matteo D'Ambrosio, 2021 (MD)

%% Check inputs

if ~isequal(size(e), [1,1])
   error('e must be a scalar') 
end
if ~isequal(size(a), [1,1])
   error('a must be a scalar') 
end
if ~isequal(size(mu), [1,1])
   error('mu must be a scalar') 
end
if ~isequal(size(tp), [1,1])
   error('tp must be a scalar') 
end

if nargin == 4
   tp = 0; 
end

%% Calculate time

T = sqrt( 4 .* pi.^2 .* a.^3 ./ mu );

if length(theta) == 1
    
    while theta < 0
       theta = theta + 2*pi; 
    end
      
    if theta < pi
                  
        E = 2 * atan( sqrt( (1-e) ./ (1+e) ) .* tan( theta ./ 2 ));
        t = tp + sqrt(a.^3 ./ mu) .* (E - e .* sin(E));
    
    elseif theta > pi

        theta_temp = 2*pi - theta; % This angle is now between 0 and pi
        E = 2 * atan( sqrt( (1-e) ./ (1+e) ) .* tan( theta_temp ./ 2 ));
        t = tp + T - sqrt( a.^3 ./ mu ) * (E - e .* sin(E));
        
    else
        
        t = T/2;
        
    end
    
elseif length(theta) == 2
        
    while theta(1) < 0
       theta = theta + 2*pi;
    end
    
    if theta(1) < pi
        
        E1 = 2 * atan( sqrt( (1-e) ./ (1+e) ) .* tan( theta(1) ./ 2 ));
        t1 = sqrt(a.^3 ./ mu) .* (E1 - e .* sin(E1)); % Does not need tp since we are calculating a deltaT
        
    elseif theta(1) > pi
        
        theta_temp = 2*pi - theta(1); % Angle between 0 and pi
        E1 = 2 * atan( sqrt( (1-e) ./ (1+e) ) .* tan( theta_temp ./ 2 ));
        t1 = T - sqrt( a.^3 ./ mu ) * (E1 - e .* sin(E1));
        
    else
        
        t1 = T/2;
        
    end
    
     if theta(2) < pi
        
        E2 = 2 * atan( sqrt( (1-e) ./ (1+e) ) .* tan( theta(2) ./ 2 ));
        t2 = sqrt(a.^3 ./ mu) .* (E2 - e .* sin(E2)); % Does not need tp since we are calculating a deltaT
        
    elseif theta(2) > pi
        
        theta_temp = 2*pi - theta(2); % Angle between 0 and pi
        E2 = 2 * atan( sqrt( (1-e) ./ (1+e) ) .* tan( theta_temp ./ 2 ));
        t2 = T - sqrt( a.^3 ./ mu ) * (E2 - e .* sin(E2));
        
     else
         
        t2 = T/2;
         
     end
    
     if theta(2) > theta(1)
         
         t = t2 - t1;
         
     else
         
         t = t2 + T - t1; % T-t1 is the time from theta(1) to 2pi, t2 is the time from 2pi to theta(2)
         
     end
    
else
    
   error('theta must be a 1x1 or 1x2 vector: theta or [theta_start, theta_end]') 
    
end

end