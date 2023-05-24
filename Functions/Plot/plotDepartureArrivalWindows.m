function [] = plotDepartureArrivalWindows( planet1, planet2, dateDepartureEarliest, dateDepartureLatest, dateArrivalEarliest, dateArrivalLatest, points, colorDeparture, colorArrival, linewidth, linealpha, figureNum )

mu_sun = astroConstants(4);

%% Time

% Departure mdj2000 time
mjd_earliestDeparture = date2mjd2000( dateDepartureEarliest ); % [year, month, day, hour, minute, and second]
mjd_latestDeparture = date2mjd2000( dateDepartureLatest );
% Arrival mdj2000 time
mjd_earliestArrival = date2mjd2000( dateArrivalEarliest ); % [year, month, day, hour, minute, and second]
mjd_latestArrival = date2mjd2000( dateArrivalLatest );
% Time vectors
mjd_Departure = linspace( mjd_earliestDeparture, mjd_latestDeparture, points ); % [days]
mjd_Arrival = linspace( mjd_earliestArrival, mjd_latestArrival, points ); % [days]

r_departure = zeros(length(mjd_Departure),3);
r_arrival = zeros(length(mjd_Arrival),3);

for k = 1:length(mjd_Departure)

    departureTime = mjd_Departure(k);

    % Planet 1 ephemeris at departure
    [kep, ~] = uplanet ( departureTime, planet1 );
    a = kep(1);
    e = kep(2);
    i = kep(3);
    OM = kep(4);
    om = kep(5);
    theta = kep(6);
    [ r_i, ~ ] = kp2rv( a, e, i, OM, om, theta, mu_sun ); % Cartesian coords of planet 1 at departure in sun-centerd ecliptic frame

    r_departure(k,:) = r_i';

end

for k = 1:length(mjd_Arrival)

    arrivalTime = mjd_Arrival(k);

    % Planet 2 ephemeris at arrival
    [kep, ~] = uplanet ( arrivalTime, planet2 );
    a = kep(1);
    e = kep(2);
    i = kep(3);
    OM = kep(4);
    om = kep(5);
    theta = kep(6);
    [ r_f, ~ ] = kp2rv( a, e, i, OM, om, theta, mu_sun ); % Cartesian coords of planet 2 at arrival in sun-centerd ecliptic frame

    r_arrival(k,:) = r_f';

end

figure(figureNum);
hold on
grid on

plot3( r_departure(:,1), r_departure(:,2), r_departure(:,3), '-', 'color', [ colorDeparture, linealpha ], 'linewidth', linewidth )

plot3( r_arrival(:,1), r_arrival(:,2), r_arrival(:,3), '-', 'color', [ colorArrival, linealpha ], 'linewidth', linewidth )


end