function [kep, EpochMJD] = TLE2kep(TLE, Mu_planet)

I_star = cell2mat(TLE(3,3));
OM_star = cell2mat(TLE(3,4));
ecc = cell2mat(TLE(3,5));
om_star = cell2mat(TLE(3,6)); 
M_star = cell2mat(TLE(3,7));
n_star = cell2mat(TLE(3,8));

I = deg2rad(I_star);
OM = deg2rad(OM_star);
om = deg2rad(om_star);

n = n_star*2*pi/86400;
T = 2*pi/n;
a = nthroot(Mu_planet/n/n, 3);

M = deg2rad(M_star);
eqn = @(E) M - E + ecc*sin(E);
eqn2 = @(E) E - ecc*sin(E);
E = fzero(eqn, M);
theta = acos((cos(E)-ecc)/(1-ecc*cos(E)));

kep = [a, ecc, I, OM, om, theta];

% Two-digit Epoch Years from 57-99 correspond to 1957-1999 and those 
% from 00-56 correspond to 2000-2056.

EpochYr_star = cell2mat(TLE(2,7));
if EpochYr_star > 56 && EpochYr_star < 100
    EpochYr = EpochYr_star + 1900;
else
    EpochYr = EpochYr_star + 2000;
end
EpochDy = cell2mat(TLE(2,8));

if EpochYr >= 2000
    EpochMJD2000 = (EpochYr - 2000)*365.25 + EpochDy - 0.5;
    EpochMJD = EpochMJD2000 + 51544.5;
else
    EpochMJD = (EpochYr - 1859)*365.25 + EpochDy - 44;
end

end
