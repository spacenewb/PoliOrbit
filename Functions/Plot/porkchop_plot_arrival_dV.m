function [ dV_min, dV1, departureDate_minDeltaV, arrivalDate_minDeltaV, a_minDeltaV,...
    P_minDeltaV, e_minDeltaV, v_i_minDeltaV, v_f_minDeltaV, t_parabolic_minDeltaV, ...
    dTheta_minDeltaV, TOF, r_i_planet1, r_f_planet1, r_i_planet2, r_f_planet2  ] = ...
    porkchop_plot_arrival_dV( planet1, planet2, dateDepartureEarliest, dateDepartureLatest, ...
    dateArrivalEarliest, dateArrivalLatest, C3_launcher, points, linewidth, marker, color, ...
    markersize, titlee )
%   THIS FUNCTION IS THE SAME AS PORKCHOP PLOT, BUT THE PLOT IS MADE
%   CONSIDERING ONLY THE FINAL dV GIVEN AT THE SECOND PLANET
%
%   porkchop_plot Calculates porkchop plot for transfer between planets
%   planet1 and planet 2 with specified transfer window. Outputs correspond
%   to optimized conditions for this transfer using function fminunc()
%
% PROTOTYPE:
%     [ dV_min, dV1, departureDate_minDeltaV, arrivalDate_minDeltaV, a_minDeltaV,...
%        P_minDeltaV, e_minDeltaV, v_i_minDeltaV, v_f_minDeltaV, t_parabolic_minDeltaV, ...
%        dTheta_minDeltaV, TOF, r_i_planet1, r_f_planet1, r_i_planet2, r_f_planet2  ] = ...
%        porkchop_plot( planet1, planet2, dateDepartureEarliest, dateDepartureLatest, ...
%        dateArrivalEarliest, dateArrivalLatest, C3_launcher, points, linewidth, marker, color, ...
%        markersize, titlee )
%
% INPUTS:
%   planet1[1]                       Integer number identifying the first
%                                    celestial body around the sun
%   planet2[1]                       Integer number identifying the second
%                                    celestial body around the sun
%                                          1:   Mercury
%                                          2:   Venus
%                                          3:   Earth
%                                          4:   Mars
%                                          5:   Jupiter
%                                          6:   Saturn
%                                          7:   Uranus
%                                          8:   Neptune
%                                          9:   Pluto
%                                          10:  Sun
%   dateDepartureEarliest[6]         [ year, month, day, hour, minute, second ]
%   dateDepartureLatest[6]           [ year, month, day, hour, minute, second ]
%   dateArrivalEarliest[6]           [ year, month, day, hour, minute, second ]
%   dateArrivalLatest[6]             [ year, month, day, hour, minute, second ]
%   C3_launcher[1]                   Characteristic energy of rocker: C3 = v_inf^2 -> Hyperbolic excess speed
%                                       Use C3_launcher = inf for unconstrained solution
%                                       Use C3_launcher = v_inf^2 for
%                                       constrained solution: Minimum found
%                                       in this way has a bigger dV_min
%                                       than unconstrained case.
%                                    deltaV1 < sqrt(C3).
%   points[1]                        Number of points in time vectors
%   linewidth[1]                     Linewidth for plot
%   minimum_marker[string]           Marker for minimum
%   min_markersize[1]                Size of marker for plot
%   title[string]                    Title for plot
%
% OUTPUTS:
%   dV_min[1]                        Minimum deltaV required          [km/s]
%   dV1[1]                           First deltaV required            [km/s]
%   departureDate_minDeltaV[6]       Departure date                   [ year, month, day, hour, minute, second ]
%   arrivalDate_minDeltaV[6]         Arrival date                     [ year, month, day, hour, minute, second ]
%   a_minDeltaV                      Semi-major axis of orbit         [km]
%   P_minDeltaV                      Semi-latus rectum of orbit       [km]
%   e_minDeltaV                      Eccentricity of orbit            [-]
%   v_i_minDeltaV                    Initial v on transfer orbit      [km/s]
%   v_f_minDeltaV                    Final v on transfer orbit        [km/s]
%   t_parabolic_minDeltaV            Parabolic time for transfer      [s]
%   dTheta_minDeltaV                 Angle between r1 and r2          [rad]
%   TOF                              Time of flight on transfer orbit [s]
%   r_i_planet1                      Initial position of planet 1     [km]
%   r_f_planet1                      Final position of planet 1     [km]
%   r_i_planet2                      Initial position of planet 2     [km]
%   r_f_planet2                      Final position of planet 2     [km]


%% Parameters
mu_sun = astroConstants( 4 );

%% Plot settings
set(0,'defaulttextinterpreter','latex')
set(0,'defaultlegendinterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'DefaultTextFontSize',12)
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontName','Times')
set(0,'DefaultAxesFontName','Times')

%% Time window

% Departure mjd2000 time
mjd_earliestDeparture = date2mjd2000( dateDepartureEarliest ); % [year, month, day, hour, minute, and second]
mjd_latestDeparture = date2mjd2000( dateDepartureLatest );

% Arrival mdj2000 time
mjd_earliestArrival = date2mjd2000( dateArrivalEarliest ); % [year, month, day, hour, minute, and second]
mjd_latestArrival = date2mjd2000( dateArrivalLatest );

% Time vectors
mjd_Departure = linspace( mjd_earliestDeparture, mjd_latestDeparture, points ); % [days]
mjd_Arrival = linspace( mjd_earliestArrival, mjd_latestArrival, points+1 ); % [days]

% Time of flight matrix for each condition in transfer window
[ dep, arr ] = ndgrid( mjd_Departure, mjd_Arrival );
TOF_matrix_mjd = arr - dep ; % [days]

%% Calculate dV for each condition in transfer window -> Only to create porkchop plot, optimisation is done afterwards

dV_t = zeros( points );
k = 1;
for departureTime = mjd_Departure
    j = 1;
    for arrivalTime = mjd_Arrival

        % Planet 1 ephemeris at departure
        [kep, ~] = uplanet ( departureTime, planet1 );
        a = kep(1);
        e = kep(2);
        i = kep(3);
        OM = kep(4);
        om = kep(5);
        theta = kep(6);
        [ r_i, ~ ] = kp2rv( a, e, i, OM, om, theta, mu_sun ); % Cartesian coords of planet 1 at departure in sun-centerd ecliptic frame

        % Planet 2 ephemeris at arrival
        [kep, ~] = uplanet ( arrivalTime, planet2 );
        a = kep(1);
        e = kep(2);
        i = kep(3);
        OM = kep(4);
        om = kep(5);
        theta = kep(6);
        [ r_f, v_f ] = kp2rv( a, e, i, OM, om, theta, mu_sun ); % Cartesian coords of planet 2 at arrival in sun-centerd ecliptic frame

        % Time of flight
        deltaT = ( arrivalTime - departureTime ) * 24 * 60 * 60;  % [s]
        
        if deltaT < 0
            
            % If TOF is negative:
            dV_t(k,j) = NaN;
            
        else

        [ ~, ~, ~, ~, ~, v_f_t, ~, ~ ] = lambertMR( r_i, r_f, deltaT, mu_sun, 0, 0, 0, 1 ); % Lambert problem solver

        dV_t(k,j) = norm( v_f - v_f_t' ); % Second deltaV required for transfer

        end

        j = j+1; % Columns
    end
    k = k+1; % Rows
end

%% Problem optimisation
% Minimisation of dV = dV( depTime, arrTime ) = dV(x): x is a 2-element vector
% Intial guess x0 = [ depTime_guess; arrTime_guess ]

% Choose constrained or unconstrained optimisation (based on C3_launcher), and solve
switch C3_launcher

    case inf % Unconstrained optimisation

        % FMINCON is still used to constrain solution to axis limits, but
        % no other constraint is considered
        
        fprintf('\nUsing unconstrained optimisation with genetic algorithm\n')

OPTIONS = optimoptions( 'ga', 'functiontolerance', 1e-16, 'constrainttolerance', 1e-16, 'populationsize', 200 );
        % Find unconstrained optimal x = [depTime, arrTime]
        A = []; % No linear inequalities for constraint
        B = [];
        Aeq = []; % No linear equalities for constraint
        Beq = [];
        LB = [ mjd_earliestDeparture; mjd_earliestArrival ]; % Lower bounds of x sets limits for solution
        UB = [ mjd_latestDeparture; mjd_latestArrival ]; % Upper bounds of x sets limits for solution

        NVARS = 3;
        x_optimal = ga(@(x) lambert_minimization_FUN_arrival_dV( x, mu_sun, planet1, planet2 ), NVARS,A,B,Aeq,Beq,LB,UB, [], OPTIONS);
        
        departureTime_minDeltaV = x_optimal(1);
        arrivalTime_minDeltaV = x_optimal(2);

    otherwise % Constrained optimisation

        fprintf('\nUsing constrained optimisation with genetic algorith\n')

OPTIONS = optimoptions( 'ga', 'functiontolerance', 1e-16, 'constrainttolerance', 1e-16, 'populationsize', 200 );
        % Constraints:
        % See local function NONLCON for non-linear constraint
        %       Non-linear constraint is such that: dV1(x) < sqrt(C3_launcher)

        A = []; % No linear inequalities for constraint
        B = [];
        Aeq = []; % No linear equalities for constraint
        Beq = [];
        LB = [ mjd_earliestDeparture; mjd_earliestArrival ]; % Lower bounds of x sets limits for solution
        UB = [ mjd_latestDeparture; mjd_latestArrival ]; % Upper bounds of x sets limits for solution

        % Find constrained optimal x = [depTime, arrTime]
        NVARS = 3;
        x_optimal = ga(@(x) lambert_minimization_FUN_arrival_dV( x, mu_sun, planet1, planet2 ), NVARS,A,B,Aeq,Beq,LB,UB, @(x) NONLCON( x, C3_launcher, mu_sun, planet1, planet2 ), OPTIONS);

        departureTime_minDeltaV = x_optimal(1);
        arrivalTime_minDeltaV = x_optimal(2);

end

% Value of first deltaV required -> Tied to laucher C3 constraint such that: dV1(x_optimal) < sqrt(C3_launcher)
[ ~, dV1 ] = lambert_minimization_FUN_arrival_dV( x_optimal, mu_sun, planet1, planet2 );

% Optimal dates (convert from mjd2000 time)
departureDate_minDeltaV = mjd20002date( x_optimal(1) );
arrivalDate_minDeltaV = mjd20002date( x_optimal(2) );

% TOF corresponding to optimal condition
TOF = ( arrivalTime_minDeltaV - departureTime_minDeltaV ) * 24 * 60 * 60;

% Planet 1 ephemeris at dV_min departure
[kep, ~] = uplanet ( departureTime_minDeltaV, planet1 );
a = kep(1);
e = kep(2);
i = kep(3);
OM = kep(4);
om = kep(5);
theta = kep(6);
[ r_i_planet1, v_i ] = kp2rv( a, e, i, OM, om, theta, mu_sun ); % Cartesian coords of planet 1 at departure in sun-centerd ecliptic frame

% Planet 1 ephemeris at dV_min arrival
[kep, ~] = uplanet( arrivalTime_minDeltaV, planet1 ); % Planet 1
a = kep(1);
e = kep(2);
i = kep(3);
OM = kep(4);
om = kep(5);
theta = kep(6);
[ r_f_planet1, ~ ] = kp2rv( a, e, i, OM, om, theta, mu_sun );

% Planet 2 ephemeris at dV_min departure
[ kep, ~ ] = uplanet( departureTime_minDeltaV, planet2 ); % Planet 2
a = kep(1);
e = kep(2);
i = kep(3);
OM = kep(4);
om = kep(5);
theta = kep(6);
[ r_i_planet2, ~ ] = kp2rv( a, e, i, OM, om, theta, mu_sun );

% Planet 2 ephemeris at dV_min arrival
[kep, ~] = uplanet ( arrivalTime_minDeltaV, planet2 );
a = kep(1);
e = kep(2);
i = kep(3);
OM = kep(4);
om = kep(5);
theta = kep(6);
[ r_f_planet2, v_f ] = kp2rv( a, e, i, OM, om, theta, mu_sun ); % Cartesian coords of planet 2 at arrival in sun-centerd ecliptic frame

% Solve Lambert once in optimal conditions for optimal parameters
[ a_minDeltaV, P_minDeltaV, e_minDeltaV, ~, v_i_minDeltaV, v_f_minDeltaV, t_parabolic_minDeltaV, dTheta_minDeltaV ] = lambertMR( r_i_planet1, r_f_planet2, TOF, mu_sun, 0, 0, 0, 1 );
dV_min = norm( v_i_minDeltaV' - v_i ) + norm( v_f - v_f_minDeltaV' ); % Total deltaV required for transfer

%% Create porkchop plot

% Axis grid
[ meshgrid_DepartureTime, meshgrid_ArrivalTime ] = meshgrid( mjd_Departure, mjd_Arrival );

% Plot contour lines with numbers
figure();
hold on
grid on
[ c, h ] = contour( meshgrid_DepartureTime, meshgrid_ArrivalTime, dV_t', [ floor( dV_min ) : 0.5 : dV_min + 5 ], 'linewidth', linewidth );
clabel( c, h );
a = colorbar;
a.Label.String = '\Delta V_{tot} [km/s]';
caxis( [dV_min, dV_min + 5 ] );


% Plot extra contour lines without numbers
contour( meshgrid_DepartureTime, meshgrid_ArrivalTime, dV_t', [ floor( dV_min + 5 ) : 0.5 : dV_min + 7.5 ], 'linewidth', linewidth, 'handlevisibility', 'off' );

% Plot constant TOF lines (only with TOF>0)
sx = size(TOF_matrix_mjd,1);
sy = size(TOF_matrix_mjd,2);
for i = 1:sx
    for j = 1:sy

        if TOF_matrix_mjd(i,j) < 0
            TOF_matrix_mjd(i,j) = NaN;
        end

    end
end
contour( mjd_Departure, mjd_Arrival, TOF_matrix_mjd', round ( linspace (min( TOF_matrix_mjd, [], 'all' ), max( TOF_matrix_mjd, [], 'all' )  , 10 ) ), 'k', 'showtext', 'on', 'linewidth', linewidth );

% Plot minimum point
plot( departureTime_minDeltaV, arrivalTime_minDeltaV, marker, 'linewidth', linewidth, 'markersize', markersize, 'color', color )

% Create title and legend
title( titlee )
legend( '$\Delta V$ [km/s]', 'TOF [days]' , ['$\Delta V_{min}$ = ', num2str( dV_min ), ' [km/s]' ] );

% Create x-axis and y-axis ticks and tick labels
N = 10; % Number of x-labels
x_ = linspace( mjd_earliestDeparture, mjd_latestDeparture, N );   % Label position
y_ = linspace( mjd_earliestArrival, mjd_latestArrival, N );   % Label position

tick_label_x = cell( 1, N );
for s = 1:N
    date_array_x = mjd20002date( x_(s) );
    tick_label_x{s} = datestr( date_array_x,'dd mmm yyyy' );
    date_array_x = datevec( tick_label_x{s},'dd mmm yyyy' );
    x_(s) = date2mjd2000( date_array_x );
end

tick_label_y = cell( 1, N );
for s = 1:N
    date_array_y = mjd20002date( y_(s) );
    tick_label_y{s} = datestr( date_array_y, 'dd mmm yyyy' );
    date_array_y = datevec( tick_label_y{s}, 'dd mmm yyyy' );
    y_(s) = date2mjd2000( date_array_y );
end

xticks( x_ );
xticklabels( tick_label_x )
xtickangle( 45 )
yticks( y_ );
yticklabels( tick_label_y )

xlabel( '\bf{Departure Window}' )
ylabel( '\bf{Arrival Window}' )

%% Create 3D surface plot

figure();
surf(meshgrid_DepartureTime,meshgrid_ArrivalTime,dV_t','edgecolor', 'none');
hold on
grid on
title( [ titlee, ' - 3D'] )
a = colorbar;
a.Label.String = '\Delta V_{tot} [km/s]';
caxis( [dV_min, dV_min + 5 ] );
zlim( [ dV_min-1, dV_min + 5] );
[ c, h ] = contour( meshgrid_DepartureTime, meshgrid_ArrivalTime, dV_t', [ floor( dV_min ) : 0.5 : dV_min + 5 ], 'linewidth', linewidth );
clabel( c, h );
xticks( x_ );
xticklabels( tick_label_x )
xtickangle( 45 )
yticks( y_ );
yticklabels( tick_label_y )

xlabel( '\bf{Departure Window}' )
ylabel( '\bf{Arrival Window}' )

end

%% Local functions:

% Function dV = dV(x) for fmincon
function [ FUN, dV1 ] = lambert_minimization_FUN_arrival_dV( x, mu_sun, planet1, planet2 )

departureTime = x(1);
arrivalTime = x(2);

% Planet 1 ephemeris at departure
[kep, ~] = uplanet ( departureTime, planet1 );
a = kep(1);
e = kep(2);
i = kep(3);
OM = kep(4);
om = kep(5);
theta = kep(6);
[ r_i, v_i ] = kp2rv( a, e, i, OM, om, theta, mu_sun ); % Cartesian coords of planet 1 at departure in sun-centerd ecliptic frame

% Planet 2 ephemeris at arrival
[kep, ~] = uplanet ( arrivalTime, planet2 );
a = kep(1);
e = kep(2);
i = kep(3);
OM = kep(4);
om = kep(5);
theta = kep(6);
[ r_f, v_f ] = kp2rv( a, e, i, OM, om, theta, mu_sun ); % Cartesian coords of planet 2 at arrival in sun-centerd ecliptic frame

% Time of flight
deltaT = ( arrivalTime - departureTime ) * 24 * 60 * 60;  % [s]

[ ~, ~, ~, ~, v_i_t, v_f_t, ~, ~ ] = lambertMR( r_i, r_f, deltaT, mu_sun, 0, 0, 0, 1 );

dV1 = norm( v_i_t' - v_i );
dV2 = norm( v_f - v_f_t' );

FUN = dV2; % deltaV given at final planet

end

% Non-linear constraint for constrained optimisation
function [ C, Ceq ] = NONLCON( x, C3_launcher, mu_sun, planet1, planet2 )

[ ~, dV1 ] = lambert_minimization_FUN_arrival_dV( x, mu_sun, planet1, planet2 );

C = dV1 - sqrt(C3_launcher);
Ceq = [];

end
