function rhoVal = computeRho(phi,k)
    syms rho
    rhoFunc = @(theta) 2*rho^3*cos(theta)-k*rho^2+k-2;
    rhoVal = vpasolve(rhoFunc(phi) == 0, rho, [0,1]);
end