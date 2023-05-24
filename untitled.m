%% ORBITAL MECHANICS PROJECT 2021/2022
clear all;
close all;
clc;

% GROUP 2164:
%   LUCIA BIANCHI
%   ALESSIA CREMASCO
%   MATTEO D'AMBROSIO
%   HARIHARAN VITALADEVUNI

% !!!!!! DO NOT TOUCH THE TOLLERANCES !!!!!!
% !!!!!! DO NOT TOUCH THE TOLLERANCES !!!!!!
% !!!!!! DO NOT TOUCH THE TOLLERANCES !!!!!!

% --------------- Assignment 2 - Planetary Explorer Mission ---------------

%% ORBITAL DATA
Earth = 3;
mu_Earth = astroConstants(Earth+10) ;

% TLE data
TLE = TLE_read("TLE.txt");
[TLE_kep, TLE_EpochMJD] = TLE2kep(TLE, mu_Earth);                          % [-] Keplerian elements of the real satellite


% Orbit data
a = 22152 ;                                                                % [km] - Semi-major axis
e = 0.6464 ;                                                               % [-] - Eccentricity
i = deg2rad(55.5215) ;                                                     % [rad] - Inclination
OM = 3/2*pi;                                                               % [rad] - Right ascention - change it
om = pi/4;                                                                 % [rad] - Pericenter anomaly - change it
theta = 0;                                                                 % [rad] - True anomaly - change it
T = 2*pi*sqrt(a^3/mu_Earth) ;                                              % [s] - Nominal orbital period

Kep0 = [a e i OM om theta] ;                                               % [-] Keplerian elements of the initial orbit
[r0, v0] = kp2rv( a, e, i, OM, om, theta, mu_Earth ); % [km] [km/s] - Initial position and velocityof the S/C

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GROUND TRACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GROUND TRACK - Unperturbed nominal orbit 
relTol = 1e-15 ;
absTol = 1e-15 ;
t0 = 0;
mjd2000 = date2mjd2000([2017 01 05 00 00 00]) ;
%mjd2000 = mjd2mjd2000(TLE_EpochMJD);
RA0_G = deg2rad(mjd20002gmst(mjd2000)*15); %[rad] - Greenwich right ascension

% perturbationModel.type = 'unperturbed' ;
% 
% % Ground track - 1 orbit
% orbits = 1 ; % Number of orbits
% % fig(1)
% [ alpha_1orb, delta_1orb, lon_1orb, lat_1orb, plotHandle_1orb ] = groundTrack( Kep0,...
%                                                                t0, RA0_G, orbits,...
%                                                                perturbationModel, relTol,...
%                                                                absTol, Earth, 10000, 'new',...
%                                                                [0 1 1], 12,...
%                                                                5, 7 ) ;
% 
% % Ground track - 1 day
% orbits = 1/(T/(24*3600)) ; % Number of orbits in 1 day
% % fig(2)
% [ alpha_1day, delta_1day, lon_1day, lat_1day, plotHandle_1day ] = groundTrack(Kep0,...
%                                                                t0, RA0_G, orbits,...
%                                                                perturbationModel, relTol,...
%                                                                absTol, Earth, 10000, 'new',...
%                                                                [0 1 0], 12,...
%                                                                5, 7 ) ;
% 
% % Ground track - 10 days
% orbits = 10/(T/(24*3600)) ; % Number of orbits in 10 days
% % fig(3)
% [ alpha_10days, delta_10days, lon_10days, lat_10days, plotHandle_10days ] = groundTrack( Kep0,...
%                                                                t0, RA0_G, orbits,...
%                                                                perturbationModel, relTol,...
%                                                                absTol, Earth, 500000, 'new',...
%                                                                [1 1 0], 12,...
%                                                                3, 7 ) ;
% 
% %% REPEATING GROUND TRACK - Unperturbed orbit 
% K = 13 ;
% m = 5 ;
% 
% perturbationModel.type = 'unperturbed' ;
% a_rgt = repeatingGroundTrack_a( m, K, e, i, perturbationModel ) ;          % [km] - Modified semi-major axis
% Kep0_new = [a_rgt e i OM om theta] ;
% T_new = 2*pi*sqrt(a_rgt^3/mu_Earth) ;                                      % [s] - Orbital period of the modified orbit
% orbits = 10/(T_new/(24*3600));                                              % [-] - Number of orbits in 10 days
% 
% % fig(4)
% [ alpha_rgt, delta_rgt, lon_rgt, lat_rgt, plotHandle_rgt ] = groundTrack( Kep0_new,...
%                                                                t0, RA0_G, orbits,...
%                                                                perturbationModel, relTol,...
%                                                                absTol, Earth, 500000, 'new',...
%                                                                [0.34 1 0], 12,...
%                                                                3, 5 ) ;
% 
% %% GROUND TRACK - Perturbed nominal orbit (J2 and SRP)
% relTol = 1e-14 ;
% absTol = 1e-14 ;
% t0 = 0;
% RA0_G = 0;
% perturbationModel.type = 'J2_SRP-perturbed' ;
% perturbationModel.startDate = mjd2000;
% perturbationModel.AreaToMass = 5 ;
% perturbationModel.Cr = 1;
% perturbationModel.p_sr = 4.5e-6;
% 
% % Ground track - 1 orbit
% orbits = 1 ; % Number of orbits
% % fig(1)
% [ alpha_1orb_per, delta_1orb_per, lon_1orb_per, lat_1orb_per, plotHandle_1orb_per ] = groundTrack( Kep0,...
%                                                                t0, RA0_G, orbits,...
%                                                                perturbationModel, relTol,...
%                                                                absTol, Earth, 500000, 1,...
%                                                                [0 0 0], 12,...
%                                                                3, 3 ) ;
% 
% % Ground track - 1 day
% orbits = 1/(T/(24*3600)) ; % Number of orbits in 1 day
% % fig(2)
% [ alpha_1day_per, delta_1day_per, lon_1day_per, lat_1day_per, plotHandle_1day_per ] = groundTrack( Kep0,...
%                                                                t0, RA0_G, orbits,...
%                                                                perturbationModel, relTol,...
%                                                                absTol, Earth, 500000, 2,...
%                                                                [0 0 0], 12,...
%                                                                3, 3 ) ;
% 
% % Ground track - 10 days
% orbits = 10/(T/(24*3600)) ; % Number of orbits in 10 days
% % fig(3)
% [ alpha_10days_per, delta_10days_per, lon_10days_per, lat_10days_per, plotHandle_10days_per ] = groundTrack( Kep0,...
%                                                                t0, RA0_G, orbits,...
%                                                                perturbationModel, relTol,...
%                                                                absTol, Earth, 500000, 3,...
%                                                                [0 0 0], 12,...
%                                                                3, 3 ) ;
% 
% %% GROUND TRACK - Unperturbed modified orbit
% relTol = 1e-14 ;
% absTol = 1e-14 ;
% t0 = 0;
% RA0_G = 0;
% perturbationModel.type = 'unperturbed' ;
% 
% % Ground track - 1 orbit
% orbits = 1 ; % Number of orbits
% % fig(5)
% [ alpha_1orb_per, delta_1orb_per, lon_1orb_per, lat_1orb_per, plotHandle_1orb_per ] = groundTrack( Kep0_new,...
%                                                                t0, RA0_G, orbits,...
%                                                                perturbationModel, relTol,...
%                                                                absTol, Earth, 50000, 'new',...
%                                                                [0 1 1], 12,...
%                                                                3, 7 ) ;
% 
% % Ground track - 1 day
% orbits = 1/(T_new/(24*3600)) ; % Number of orbits in 1 day
% % fig(6)
% [ alpha_1day_per, delta_1day_per, lon_1day_per, lat_1day_per, plotHandle_1day_per ] = groundTrack( Kep0_new,...
%                                                                t0, RA0_G, orbits,...
%                                                                perturbationModel, relTol,...
%                                                                absTol, Earth, 500000, 'new',...
%                                                                [0 1 0], 12,...
%                                                                3, 7 ) ;
% 
% % Ground track - 10 days
% orbits = 10/(T_new/(24*3600)) ; % Number of orbits in 10 days
% % fig(7)
% [ alpha_10days_per, delta_10days_per, lon_10days_per, lat_10days_per, plotHandle_10days_per ] = groundTrack( Kep0_new,...
%                                                                t0, RA0_G, orbits,...
%                                                                perturbationModel, relTol,...
%                                                                absTol, Earth, 500000, 'new',...
%                                                                [1 1 0], 12,...
%                                                                3, 7 ) ;
% 
% %% GROUND TRACK - Perturbed modified orbit (J2 and SRP)
% relTol = 1e-14 ;
% absTol = 1e-14 ;
% t0 = 0;
% RA0_G = 0;
% perturbationModel.type = 'J2_SRP-perturbed' ;
% perturbationModel.startDate = mjd2000;
% perturbationModel.AreatoMass = 5 ;
% perturbationModel.Cr = 1;
% perturbationModel.p_sr = 4.5e-6;
% 
% % Ground track - 1 orbit
% orbits = 1 ; % Number of orbits
% % fig(5)
% [ alpha_1orb_per, delta_1orb_per, lon_1orb_per, lat_1orb_per, plotHandle_1orb_per ] = groundTrack( Kep0_new,...
%                                                                t0, RA0_G, orbits,...
%                                                                perturbationModel, relTol,...
%                                                                absTol, Earth, 50000, 5,...
%                                                                [0 0 0], 12,...
%                                                                3, 3 ) ;
% 
% % Ground track - 1 day
% orbits = 1/(T_new/(24*3600)) ; % Number of orbits in 1 day
% % fig(6)
% [ alpha_1day_per, delta_1day_per, lon_1day_per, lat_1day_per, plotHandle_1day_per ] = groundTrack( Kep0_new,...
%                                                                t0, RA0_G, orbits,...
%                                                                perturbationModel, relTol,...
%                                                                absTol, Earth, 500000, 6,...
%                                                                [0 0 0], 12,...
%                                                                3, 3 ) ;
% 
% % Ground track - 10 days
% orbits = 10/(T_new/(24*3600)) ; % Number of orbits in 10 days
% % fig(7)
% [ alpha_10days_per, delta_10days_per, lon_10days_per, lat_10days_per, plotHandle_10days_per ] = groundTrack( Kep0_new,...
%                                                                t0, RA0_G, orbits,...
%                                                                perturbationModel, relTol,...
%                                                                absTol, Earth, 500000, 7,...
%                                                                [0 0 0], 12,...
%                                                                3, 3 ) ;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         %% GROUND TRACK - Unerturbed nominal orbit 
% %         relTol = 1e-14 ;
% %         absTol = 1e-14 ;
% %         t0 = 0;
% %         RA0_G = 0;
% %         perturbationModel.type = 'unperturbed' ;
% %         perturbationModel.startDate = mjd2000;
% %         perturbationModel.AreaToMass = 0 ;
% %         perturbationModel.Cr = 0;
% %         perturbationModel.p_sr = 0;
% %         
% %        % Ground track - 100 orbits
% %         orbit = 100 ; % Number of orbits
% %         fig(8)
% %         [ alpha_1orb_per, delta_1orb_per, lon_1orb_per, lat_1orb_per, plotHandle_1orb_per ] = groundTrack( Kep0,...
% %                                                                        t0, RA0_G, orbit,...
% %                                                                        perturbationModel, relTol,...
% %                                                                        absTol, Earth, 1000000, 'new',...
% %                                                                        [0 1 1], 12,...
% %                                                                        3, 5 ) ;
% %         
% %         %% GROUND TRACK - Perturbed nominal orbit (J2 and SRP)
% %         relTol = 1e-14 ;
% %         absTol = 1e-14 ;
% %         t0 = 0;
% %         RA0_G = 0;
% %         perturbationModel.type = 'J2_SRP-perturbed' ;
% %         perturbationModel.startDate = mjd2000;
% %         perturbationModel.AreaToMass = 5 ;
% %         perturbationModel.Cr = 1;
% %         perturbationModel.p_sr = 4.5e-6;
% %         
% %        % Ground track - 100 orbits
% %         orbit = 100 ; % Number of orbits
% %         fig(8)
% %         [ alpha_1orb_per, delta_1orb_per, lon_1orb_per, lat_1orb_per, plotHandle_1orb_per ] = groundTrack( Kep0,...
% %                                                                        t0, RA0_G, orbit,...
% %                                                                        perturbationModel, relTol,...
% %                                                                        absTol, Earth, 1000000, 8,...
% %                                                                        [0 0 0], 12,...
% %                                                                        3, 3 ) ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%% ORBIT PROPAGATION - Perturbed orbit %%%%%%%%%%%%%%%%%%%
relTol = 1e-15;
absTol = 1e-15;

perturbationModel.type = 'J2_SRP-perturbed'; % [-] - Orbit type
perturbationModel.Cr = 1;             % [-] - Reflectivity coefficient
perturbationModel.AreaToMass = 5;    % [m^2/kg] - Area to mass ratio relative to the Sun
perturbationModel.p_sr = 4.5e-6;           % [N/m^2] - Solar radiation ressure
perturbationModel.startDate = mjd2000;      % [-] - MJD2000 date


%% CARTESIAN COORDINATES

% 10 orbits propagation
orbits = 10; % [-] Number of orbits 
points = ceil((orbits*T)/600);
[time_car_10orb, Kep_car_10orb] = Car_propagate(orbits, [r0 v0], Kep0, [relTol absTol], points, perturbationModel, Earth);

% 1000 days propagation (approximately)
orbits = 2700; % [-] Number of orbits 
points = ceil((orbits*T)/(3600));
[time_car_2700orb, Kep_car_2700orb] = Car_propagate(orbits, [r0 v0], Kep0, [relTol absTol], points, perturbationModel, Earth);

%% GUASS PLANETARY EQUATIONS

% 10 orbits propagation
orbits = 10; % [-] Number of orbits 
points = ceil((orbits*T)/600);
[time_gauss_10orb, Kep_gauss_10orb] = Gauss_propogate(orbits, Kep0, [relTol absTol], points, perturbationModel, Earth);   

% 1000 days propagation (approximately)
orbits = 2700; % [-] Number of orbits 
points = ceil((orbits*T)/(3600));
[time_gauss_2700orb, Kep_gauss_2700orb] = Gauss_propogate(orbits, Kep0, [relTol absTol], points, perturbationModel, Earth);   

%% PLOT
% a
element.type = 'semi-major axis';
plot_kep_error(Kep0, Kep_gauss_10orb, Kep_gauss_2700orb, ...
                             Kep_car_10orb, Kep_car_2700orb, time_gauss_10orb , ...
                             time_gauss_2700orb, element) ;


% e
element.type = 'eccentricity';
plot_kep_error(Kep0, Kep_gauss_10orb, Kep_gauss_2700orb, ...
                             Kep_car_10orb, Kep_car_2700orb, time_gauss_10orb , ...
                             time_gauss_2700orb, element) ;

% i
element.type = 'inclination';
plot_kep_error(Kep0, Kep_gauss_10orb, Kep_gauss_2700orb, ...
                             Kep_car_10orb, Kep_car_2700orb, time_gauss_10orb , ...
                             time_gauss_2700orb, element) ;

% OM
element.type = 'RAAN';
plot_kep_error(Kep0, Kep_gauss_10orb, Kep_gauss_2700orb, ...
                             Kep_car_10orb, Kep_car_2700orb, time_gauss_10orb , ...
                             time_gauss_2700orb, element) ;

% om
element.type = 'pericenter anomaly';
plot_kep_error(Kep0, Kep_gauss_10orb, Kep_gauss_2700orb, ...
                             Kep_car_10orb, Kep_car_2700orb, time_gauss_10orb , ...
                             time_gauss_2700orb, element) ;
% theta
element.type = 'true anomaly';
plot_kep_error(Kep0, Kep_gauss_10orb, Kep_gauss_2700orb, ...
                             Kep_car_10orb, Kep_car_2700orb, time_gauss_10orb , ...
                             time_gauss_2700orb, element) ;

% %%%%%%%%%%%%%%%%%%%%%%%%%% EVOLUTION OF THE ORBIT %%%%%%%%%%%%%%%%%%%%%%%%%

% 
% %%%%%%%%%%%%%%%%%%%% FILTERING OF HIGH FREQUENCIES %%%%%%%%%%%%%%%%%%%%%%%%

%% FREQUENCY - find the cut-off frequency
perturbationModel.type = 'J2-perturbed';
orbits = 10; % [-] Number of orbits 
points = ceil((orbits*T)/600);
[time_gauss_10orb_J2, Kep_gauss_10orb_J2] = Gauss_propogate(orbits, Kep0, [relTol absTol], points, perturbationModel, Earth);   

perturbationModel.type = 'SRP-perturbed';
[time_gauss_10orb_SRP, Kep_gauss_10orb_SRP] = Gauss_propogate(orbits, Kep0, [relTol absTol], points, perturbationModel, Earth); 

% Semi-major axis
figure
plot(time_car_10orb/T, Kep_gauss_10orb_J2(:,1),'r')
hold on
grid on
plot(time_car_10orb/T, Kep_gauss_10orb_SRP(:,1),'b')
xlabel('Orbits')
ylabel('Semi-major axis [km]');
legend('J2','SRP')

% Eccentricity
figure
plot(time_car_10orb/T, Kep_gauss_10orb_J2(:,2),'r')
hold on
grid on
plot(time_car_10orb/T, Kep_gauss_10orb_SRP(:,2),'b')
xlabel('Orbits')
ylabel('Eccentricity [-]');
legend('J2','SRP')

% Inclination
figure
plot(time_car_10orb/T, rad2deg(Kep_gauss_10orb_J2(:,3)),'r')
hold on
grid on
plot(time_car_10orb/T, rad2deg(Kep_gauss_10orb_SRP(:,3)),'b')
xlabel('Orbits')
ylabel('Inclination [deg]');
legend('J2','SRP')

% RAAN
figure
plot(time_car_10orb/T, rad2deg(Kep_gauss_10orb_J2(:,4)),'r')
hold on
grid on
plot(time_car_10orb/T, rad2deg(Kep_gauss_10orb_SRP(:,4)),'b')
xlabel('Orbits')
ylabel('Right ascension [deg]');
legend('J2','SRP')

% Pericenter anomaly
figure
plot(time_car_10orb/T, rad2deg(Kep_gauss_10orb_J2(:,5)),'r')
hold on
grid on
plot(time_car_10orb/T, rad2deg(Kep_gauss_10orb_SRP(:,5)),'b')
xlabel('Orbits')
ylabel('Pericenter anomaly [deg]');
legend('J2','SRP')
% True anomaly
figure
plot(time_car_10orb/T, rad2deg(Kep_gauss_10orb_J2(:,6)),'r')
hold on
grid on
plot(time_car_10orb/T, rad2deg(Kep_gauss_10orb_SRP(:,6)),'b')
xlabel('Orbits')
ylabel('True anomaly [deg]');
legend('J2','SRP')

%% Filter data
f_sampling = 10000/(orbits*T); % [Hz] - Sampling frequency
cut_off_orbit = orbits/5; 
f_cutoff = 1/(cut_off_orbit*T) ; % [Hz] - Cut off frequency (function of number of revolutions)

% Semi-major axis
a_fil = Kep_gauss_2700orb(:,1); 
filtered_a = filterData(a_fil, f_cutoff, f_sampling);

figure
plot(time_car_2700orb/T, Kep_gauss_2700orb(:,1),'r')
hold on
grid on
plot(time_car_2700orb/T, filtered_a,'b')
xlabel('Orbits')
ylabel('Semi-major axis [km]');
legend('Unfiltered','Filtered')

% Eccentricity
e_fil = Kep_gauss_2700orb(:,2); 
filtered_e = filterData(e_fil, f_cutoff, f_sampling);

figure
plot(time_car_2700orb/T, Kep_gauss_2700orb(:,2),'r')
hold on
grid on
plot(time_car_2700orb/T, filtered_e,'b')
xlabel('Orbits')
ylabel('Eccentricity [-]');
legend('Unfiltered','Filtered')

% Inclination
i_fil = Kep_gauss_2700orb(:,3); 
filtered_i = filterData(i_fil, f_cutoff, f_sampling);

figure
plot(time_car_2700orb/T, rad2deg(Kep_gauss_2700orb(:,3)),'r')
hold on
grid on
plot(time_car_2700orb/T, rad2deg(filtered_i),'b')
xlabel('Orbits')
ylabel('Inclination [deg]');
legend('Unfiltered','Filtered')

% RAAN
OM_fil = Kep_gauss_2700orb(:,4); 
filtered_OM = filterData(OM_fil, f_cutoff, f_sampling);

figure
plot(time_car_2700orb/T, rad2deg(Kep_gauss_2700orb(:,4)),'r')
hold on
grid on
plot(time_car_2700orb/T, rad2deg(filtered_OM),'b')
xlabel('Orbits')
ylabel('Right ascension [deg]');
legend('Unfiltered','Filtered')

% Pericenter anomaly
om_fil = Kep_gauss_2700orb(:,5); 
filtered_om = filterData(om_fil, f_cutoff, f_sampling);

figure
plot(time_car_2700orb/T, rad2deg(Kep_gauss_2700orb(:,5)),'r')
hold on
grid on
plot(time_car_2700orb/T, rad2deg(filtered_om),'b')
xlabel('Orbits')
ylabel('Pericenter anomaly [deg]');
legend('Unfiltered','Filtered')

% True anomaly
theta_fil = Kep_gauss_2700orb(:,6); 
filtered_theta = filterData(theta_fil, f_cutoff, f_sampling);

figure
plot(time_car_2700orb/T, rad2deg(Kep_gauss_2700orb(:,6)),'r')
hold on
grid on
plot(time_car_2700orb/T, rad2deg(filtered_theta),'b')
xlabel('Orbits')
ylabel('True anomaly [deg]');
legend('Unfiltered','Filtered')


%%%%%%%%%%%%%%%%%%%%%%% COMPARINSON WITH REAL ORBIT %%%%%%%%%%%%%%%%%%%%%%
%% Gauss propagation of TLE for 2700 Orbits
[TLE_time_gauss_2700_orb, TLE_Kep_gauss_2700_orb] = Gauss_propogate(2700, TLE_kep_0, perturbationModel, Earth);
[TLE_time_car_2700_orb, TLE_Kep_car_2700_orb] = Car_propogate(2700, TLE_kep_0, perturbationModel, Earth);

%% Plots
plot_compare_real(TLE_time_gauss_2700_orb, Kep_gauss_2700orb, TLE_Kep_gauss_2700_orb, 1, 0);