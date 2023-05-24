function [] = merge_orbitCheckrv( energyFigNum1, angmomentumFigNum1, ...
                                  eccentricityFigNum1, edothFigNum1, ...
                                  vrvthetaFigNum1, aFigNum1,         ...
                                  iFigNum1, OMFigNum1,               ...
                                  omFigNum1, thetaFigNum1,           ...
                                  energyFigNum2, angmomentumFigNum2, ...
                                  eccentricityFigNum2, edothFigNum2, ...
                                  vrvthetaFigNum2, aFigNum2,         ...
                                  iFigNum2, OMFigNum2,               ...
                                  omFigNum2, thetaFigNum2,           ...
                                  energyColor2, edothColor2, legendFontsize )
% merge_orbitCheckrv merges plots coming from orbitCheckrv in 2 cases
% Example: unperturbed and perturbed case

% Specific energy
secondplotObject = findobj( energyFigNum2, 'type', 'line' );
secondplotObject.Color = energyColor2;
copyobj( secondplotObject, findobj( energyFigNum1, 'type', 'axes' ) );
figure( energyFigNum1 );
hold on
grid on
legend('\bf{2BP}', '\bf{J2-perturbed}', 'interpreter', 'latex', 'fontsize', legendFontsize );

% Angular momentum
secondplotObject = findobj( angmomentumFigNum2, 'type', 'line' );
copyobj( secondplotObject, findobj( angmomentumFigNum1, 'type', 'axes' ) );
figure( angmomentumFigNum1 );
hold on
grid on
legend( '\boldmath{$h_x$, 2BP}', '\boldmath{$h_y$, 2BP}', '\boldmath{$h_z$, 2BP}', '$\bf{||\underline{h}||}$, 2BP', '\boldmath{$h_x$, J2}', '\boldmath{$h_y$, J2}', '\boldmath{$h_z$, J2}', '$\bf{||\underline{h}||}$, J2', 'interpreter', 'latex', 'fontsize', legendFontsize);

% Eccentricity
secondplotObject = findobj( eccentricityFigNum2, 'type', 'line' );
copyobj( secondplotObject, findobj( eccentricityFigNum1, 'type', 'axes' ) );
figure( eccentricityFigNum1 );
hold on
grid on
legend( '\boldmath{$e_x$, 2BP}', '\boldmath{$e_y$, 2BP}', '\boldmath{$e_z$, 2BP}', '$\bf{||\underline{e}||}$, 2BP', '\boldmath{$e_x$, J2}', '\boldmath{$e_y$, J2}', '\boldmath{$e_z$, J2}', '$\bf{||\underline{e}||}$, J2', 'interpreter', 'latex', 'fontsize', legendFontsize);

% e*h
secondplotObject = findobj( edothFigNum2, 'type', 'line' );
secondplotObject.Color = edothColor2;
copyobj( secondplotObject, findobj( edothFigNum1, 'type', 'axes' ) );
figure( edothFigNum1 );
hold on
grid on
legend('\bf{2BP}', '\bf{J2-perturbed}', 'interpreter', 'latex', 'fontsize', legendFontsize );

% vr vtheta
secondplotObject = findobj( vrvthetaFigNum2, 'type', 'line' );
copyobj( secondplotObject, findobj( vrvthetaFigNum1, 'type', 'axes' ) );
figure( vrvthetaFigNum1 );
hold on
grid on
legend('\boldmath{$v_r$}, 2BP', '\boldmath{$v_{\theta} $}, 2BP','\boldmath{$v_r$}, J2', '\boldmath{$v_{\theta} $, J2}', 'interpreter', 'latex', 'fontsize', legendFontsize);

% a
secondplotObject = findobj( aFigNum2, 'type', 'line' );
copyobj( secondplotObject, findobj( aFigNum1, 'type', 'axes' ) );
figure( aFigNum1 );
hold on
grid on
legend('\boldmath{$a$}, 2BP', '\boldmath{$a$}, J2', 'interpreter', 'latex', 'fontsize', legendFontsize);

% i
secondplotObject = findobj( iFigNum2, 'type', 'line' );
copyobj( secondplotObject, findobj( iFigNum1, 'type', 'axes' ) );
figure( iFigNum1 );
hold on
grid on
legend('\boldmath{$i$}, 2BP', '\boldmath{$i$}, J2', 'interpreter', 'latex', 'fontsize', legendFontsize);

% OM
secondplotObject = findobj( OMFigNum2, 'type', 'line' );
copyobj( secondplotObject, findobj( OMFigNum1, 'type', 'axes' ) );
figure( OMFigNum1 );
hold on
grid on
legend('\boldmath{$\Omega$}, 2BP', '\boldmath{$\Omega$}, J2', 'interpreter', 'latex', 'fontsize', legendFontsize)

% om
secondplotObject = findobj( omFigNum2, 'type', 'line' );
copyobj( secondplotObject, findobj( omFigNum1, 'type', 'axes' ) );
figure( omFigNum1 );
hold on
grid on
legend('\boldmath{$\omega$}, 2BP', '\boldmath{$\omega$}, J2', 'interpreter', 'latex', 'fontsize', legendFontsize);

% theta
secondplotObject = findobj( thetaFigNum2, 'type', 'line' );
copyobj( secondplotObject, findobj( thetaFigNum1, 'type', 'axes' ) );
figure( thetaFigNum1 );
hold on
grid on
legend('\boldmath{$\theta$}, 2BP', '\boldmath{$\theta$}, J2', 'interpreter', 'latex', 'fontsize', legendFontsize);


%% Close extra figures

figure(energyFigNum2);
close
figure(angmomentumFigNum2);
close
figure(eccentricityFigNum2);
close
figure(edothFigNum2);
close
figure(vrvthetaFigNum2);
close
figure(aFigNum2);
close
figure(iFigNum2);
close
figure(OMFigNum2);
close
figure(omFigNum2);
close
figure(thetaFigNum2);
close

end