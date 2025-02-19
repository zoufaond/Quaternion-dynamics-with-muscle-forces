function force = muscle_force(act, lm, fmax, lceopt)
    % simple muscle force using properties only based on muscle length

    f_gauss = 0.25;
    kpe = 5;
    epsm0 = 0.6;
    fpe = (exp(kpe * (lm / lceopt - 1) / epsm0)-1) / (exp(kpe) - 1);
    flce = (exp( - (lm / lceopt - 1)^2 / f_gauss));
    force = (flce * act +  fpe) * fmax;
end