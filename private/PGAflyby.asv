function [DV,r_p,e_in,e_out,a_in,a_out,V_p_in,V_p_out] = PGAflyby(V_in, V_out,Mu_body, V_body,h_atm,R_planet)

% PGAflyby   Function to evaluate all the data from a powered flyby around a celestial body 
%
% PROTOTYPE:
% [DV,r_p,e_in,e_out,a_in,a_out,V_p_in,V_p_out] = PGAflyby(V_in, V_out,Mu_body, V_body,h_atm,R_planet)
%
% INPUT:
%       VV_in[1x3]      Velocity of incoming arc in heliocentric centered frame [km/s]
%       VV_out[1x3]     Velocity of outgoing arc in heliocentric centered frame [km/s]
%       mu_p[1]         Gravitational costant of the planet [km^3/s^2]
%       vv_p[1x3]       Velocity of the planet at the moment of flyby [km/s]
%       h_atm[1]        Height of the atmosphere of the celestial body [km]
%       R_planet[1]     Radius of the celestial body of the flyby [km]
%
% OUTPUT:
%       DV[1]           Cost of the manoeuvre at perigee in terms of delta velocity [km/s]
%       r_p[1]          Radius of perigee of the two arcs of hyperbola of the flyby [km]   
%       e_in[1]         Eccentricity of the incoming hyperbola [-]
%       e_out[1]        Eccentricity of the outgoing hyperbola [-]
%       a_in[1]         Semi-major axis of the incoming hyperbola [km]
%       a_out[1]        Semi-major axis of the incoming hyperbola [km]
%       v_p_in[1]       Velocity at perigee of incoming hyperbola [km/s]
%       v_p_out[1]      Velocity at perigee of outgoing hyperbola [km/s]
% 

%% Establishing Velocities
% Evaluate infinite velocity in planet center frame 
V_inf_in_vec = V_in - V_body;
V_inf_in = norm(V_inf_in_vec);

V_inf_out_vec = (V_out - V_body);
V_inf_out = norm(V_inf_out_vec);

% Evaluate the delta angle between the asymptote 
Delta = acos((V_inf_in_vec*V_inf_out_vec')./(V_inf_in.*V_inf_out));

%% Finding the r_p of the flyby 
delta_minus  = @(r_p)2.*asin(1./(1+ r_p .* V_inf_in.^2./Mu_body));
delta_plus   = @(r_p)2.*asin(1./(1+ r_p .* V_inf_out.^2./Mu_body));

r_p_SolveFun = @(r_p)(delta_minus(r_p) + delta_plus(r_p))./2 - Delta;

r_p_min = R_planet + h_atm;

options = optimoptions('fsolve','TolFun',1e-14,'Display','off');
%r_p = fsolve(r_p_SolveFun, r_p_min, options);
r_p = fzero(r_p_SolveFun, r_p_min);

%% Check if the flyby is possible 
if r_p < r_p_min
    
    % If the flyby is not possible set the variable to nan
    DV=NaN;
    V_p_in=NaN;
    V_p_out=NaN;
    
else
    
    % Calculate the velocities at perigee and necessary delta V
    V_p_in = sqrt(V_inf_in.^2 + 2.*Mu_body./r_p);
    V_p_out = sqrt(V_inf_out.^2 + 2.*Mu_body./r_p);
    DV = abs(V_p_out - V_p_in);
    
end

%% Compute some keplerian elements of the two hyperbola
e_in = 1/sin(delta_minus(r_p)/2);
e_out = 1/sin(delta_plus(r_p)/2);
a_in = -r_p/(1-e_in);
a_out = -r_p/(1-e_out);




