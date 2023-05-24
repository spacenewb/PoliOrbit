clear, close all, clc;

%%
abstol = 1e-15;
reltol = 1e-14;

Date_initial = [ 2021, 12, 14, 0, 0, 0];
Kep0 = [ 7571; 0.01; deg2rad(87.9); pi; pi; 0 ];
T = sqrt( 4 * pi.^2 .* Kep0(1).^3 / astroConstants(13) );
orbits = 10;
TSPAN10 = linspace(0, orbits*T, 100*orbits/1000*T);

options = odeset( 'reltol', reltol, 'abstol', abstol);
[ timevect10, KepGauss10 ] = ode45( @( t, Kep ) gauss_planetary_equations( t, Kep, Date_initial ), TSPAN10, Kep0, options );

[ r0, v0 ] = kp2rv( Kep0(1), Kep0(2), Kep0(3), Kep0(4), Kep0(5), Kep0(6), astroConstants(13) );

perturbationModel.type = 'J2-perturbed';

[ r_2BP, v_2BP, TSPAN10 ] = orbit_propagation( r0, v0, TSPAN10, astroConstants(13), reltol, abstol, perturbationModel );
Kep2BP10 = zeros(length(TSPAN10),6);

for k = 1:length(TSPAN10)
    [ Kep2BP10(k,1), Kep2BP10(k,2), Kep2BP10(k,3), Kep2BP10(k,4), Kep2BP10(k,5), Kep2BP10(k,6) ] = rv2kp( r_2BP(k,:), v_2BP(k,:) , astroConstants(13) );
end

orbits = 100;
TSPAN100 = linspace(0, orbits*T, orbits*T);

options = odeset( 'reltol', 1e-16, 'abstol', 1e-16);
[ timevect100, KepGauss100 ] = ode45( @( t, Kep ) gauss_planetary_equations( t, Kep, Date_initial ), TSPAN100, Kep0, options );

[ r0, v0 ] = kp2rv( Kep0(1), Kep0(2), Kep0(3), Kep0(4), Kep0(5), Kep0(6), astroConstants(13) );

perturbationModel.type = 'J2-perturbed';
relTol = 1e-16;
absTol = 1e-16;
[ r_2BP, v_2BP, TSPAN100 ] = orbit_propagation( r0, v0, TSPAN100, astroConstants(13), reltol, abstol, perturbationModel );
Kep2BP100 = zeros(length(TSPAN100),6);

for k = 1:length(TSPAN100)
    [ Kep2BP100(k,1), Kep2BP100(k,2), Kep2BP100(k,3), Kep2BP100(k,4), Kep2BP100(k,5), Kep2BP100(k,6) ] = rv2kp( r_2BP(k,:), v_2BP(k,:) , astroConstants(13) );

end

%% Plot elements

% a
figure
subplot(3,1,1)
semilogy(TSPAN100, abs((KepGauss100(:,1)-Kep2BP100(:,1))./Kep0(1)));
hold on
grid on
ylabel('a_{car} - a_{gauss} / a0')
subplot(3,1,2)
hold on
grid on
plot(timevect100, KepGauss100(:,1))
plot(TSPAN100, Kep2BP100(:,1))
subplot(3,1,3)
hold on
grid on
plot(timevect10, KepGauss10(:,1))
plot(TSPAN10, Kep2BP10(:,1))

% e
figure
subplot(3,1,1)
semilogy(TSPAN100, abs((KepGauss100(:,2)-Kep2BP100(:,2))));
hold on
grid on
ylabel('e_{car} - e_{gauss}')
subplot(3,1,2)
hold on
grid on
plot(timevect100, KepGauss100(:,2))
plot(TSPAN100, Kep2BP100(:,2))
subplot(3,1,3)
hold on
grid on
plot(timevect10, KepGauss10(:,2))
plot(TSPAN10, Kep2BP10(:,2))

% i
figure
subplot(3,1,1)
semilogy(TSPAN100, rad2deg(abs(KepGauss100(:,3)-Kep2BP100(:,3))/(2*pi)));
hold on
grid on
ylabel('i_{car} - i_{gauss} / 2*pi')
subplot(3,1,2)
hold on
grid on
plot(timevect100, rad2deg(KepGauss100(:,3)))
plot(TSPAN100, rad2deg(Kep2BP100(:,3)))
subplot(3,1,3)
hold on
grid on
plot(timevect10, rad2deg(KepGauss10(:,3)))
plot(TSPAN10, rad2deg(Kep2BP10(:,3)))

% OM
figure
subplot(3,1,1)
semilogy(TSPAN100, rad2deg(abs(KepGauss100(:,4)-Kep2BP100(:,4))/(2*pi)));
hold on
grid on
ylabel('OM_{car} - OM_{gauss} / 2*pi')
subplot(3,1,2)
hold on
grid on
plot(timevect100, rad2deg(KepGauss100(:,4)))
plot(TSPAN100, rad2deg(Kep2BP100(:,4)))
subplot(3,1,3)
hold on
grid on
plot(timevect10, rad2deg(KepGauss10(:,4)))
plot(TSPAN10, rad2deg(Kep2BP10(:,4)))

% om
figure
subplot(3,1,1)
semilogy(TSPAN100, rad2deg(abs(KepGauss100(:,5)-Kep2BP100(:,5))/(2*pi)));
hold on
grid on
ylabel('om_{car} - om_{gauss} / 2*pi')
subplot(3,1,2)
hold on
grid on
plot(timevect100, rad2deg(KepGauss100(:,5)))
plot(TSPAN100, rad2deg(Kep2BP100(:,5)))
subplot(3,1,3)
hold on
grid on
plot(timevect10, rad2deg(KepGauss10(:,5)))
plot(TSPAN10, rad2deg(Kep2BP10(:,5)))

% theta
figure
subplot(3,1,1)
semilogy(TSPAN100, rad2deg(abs(KepGauss100(:,6)-Kep2BP100(:,6))/(2*pi)));
hold on
grid on
ylabel('theta_{car} - theta_{gauss} / 2*pi')
subplot(3,1,2)
hold on
grid on
plot(timevect100, rad2deg(KepGauss100(:,6)))
plot(TSPAN100, rad2deg(Kep2BP100(:,6)))
subplot(3,1,3)
hold on
grid on
plot(timevect10, rad2deg(KepGauss10(:,6)))
plot(TSPAN10, rad2deg(Kep2BP10(:,6)))

