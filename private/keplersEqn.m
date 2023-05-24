function [f, f_bar, k] = keplersEqn(mu, a, e, t, t0, f0)

M0 = f2M(e,f0);
M = M0 + sqrt(mu/a^3)*(t-t0);

k = floor(M/(2*pi));
M_bar = M - k*2*pi;

kepler_equation = @(E) E - e*sin(E) - M_bar;

E_bar_guess = M_bar + ( e*sin(M_bar) / ( 1 - sin(M_bar+e) + sin(M_bar) ) );
E_bar = fzero(kepler_equation, E_bar_guess);

f_bar = atan2( sqrt(1-e^2)*sin(E_bar), (cos(E_bar)-e) );
f_bar = wrapTo2Pi(f_bar);

f = f_bar + k*2*pi;

end