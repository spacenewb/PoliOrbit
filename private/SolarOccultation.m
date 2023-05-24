function LineOfSight = SolarOccultation(r_sat, planetID, ref_mjd2000)

R_P = astroConstants(20 + planetID); % Radius of Planet

[Kep_Planet, Mu_Sun] = uplanet(ref_mjd2000,planetID);

y_Planet = kep2car(Kep_Planet', Mu_Sun);
R_Planet = y_Planet(1:3); % Heliocentric Position of Planet

R_sat = R_Planet + r_sat; % Heliocentric Position of Satellite

ShadeHalfAngleSun = atan(R_P/norm(R_Planet)); % Angle Subtended By Planet Shade at Sun
SatAngleSun = real(acos(dot(R_Planet,R_sat)/(norm(R_Planet)*norm(R_sat)))); % Angle Subtended By Satellite and Planet at Sun

if abs(SatAngleSun) <= abs(ShadeHalfAngleSun) && norm(R_sat) >= norm(R_Planet)
    LineOfSight = 0; % Satellite Illumination Occluded by Planet
else
    LineOfSight = 1; % Satellite Illumination Not Occluded by Planet
end

end