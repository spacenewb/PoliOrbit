function E = f2E(e,f)

eqn = @(E) cos(f) - (cos(E)-e)/(1-e*cos(E)) - tan(f/2) + sqrt((1+e)/(1-e))*tan(E/2);
E = fzero(eqn,f);

end