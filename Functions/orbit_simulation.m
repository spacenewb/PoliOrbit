function [] = orbit_simulation( r, TSPAN, omega_planet, planetNum, figureNum, multiplier, orbit_linestyle, orbit_linewidth, orbit_color, satellite_color, satellite_marker, satellite_marker_size, view_vector )



% Add planet
planet_object = planet3D(planetNum, [0,0,0], figureNum, multiplier);
view( view_vector(1), view_vector(2) );
% Add orbit
LegendVisibility = 0;
orbitPlotrv( r , orbit_linestyle, figureNum, orbit_color, orbit_linewidth, LegendVisibility)

length_r = length(r(:,1));
deltaT = TSPAN(end)/length_r;

rotation_angle = rad2deg(omega_planet*deltaT);

r_unit = zeros(size(r));

for i = 1 : ( length_r - 1 )
    r_unit(i,:) = r(i,:)/norm(r(i,:));
end

direction = [ 0, 0, 1 ];

ORIGIN = [0,0,0];

for i = 1 : ( length_r - 1 )

    sat = plot3( r(i,1), r(i,2), r(i,3), satellite_marker, 'color', satellite_color, 'markerfacecolor', satellite_color, 'MarkerSize', satellite_marker_size );

    % Update plot
    drawnow;
    rotate( planet_object, direction, rotation_angle, ORIGIN);
    delete(sat);

end
