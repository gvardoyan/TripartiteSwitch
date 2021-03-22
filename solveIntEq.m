function phiCurrent = solveIntEq(k,M,epsilon)
    thetaVals = linspace(-pi,pi,M);
    
    % compute cotangents
    cotMat = zeros(M,M-1);
    for j = 1:M
        cotMat(j,:) = cot((thetaVals([1:j-1 j+1:end])-thetaVals(j))/2);
    end
    
    maxNorm = Inf;
    phiCurrent = thetaVals;
    % solve iteratively
    while maxNorm > epsilon
        % first, solve for rho's using current phi vals:
        rhoVals = zeros(1,M);
        for idx = 1:M
            rhoVals(idx) = computeRho(phiCurrent(idx),k);
        end
        % compute the new phi vals:
        phiNew = zeros(1,M);
        for j = 1:M
            S1 = log(rhoVals([1:j-1 j+1:end]))-log(rhoVals(j));
            S2 = cotMat(j,:);
            phiNew(j) = thetaVals(j) - (1/M)*sum(S1*S2');
        end
        % update the max norm
        maxNorm = max(abs(phiNew - phiCurrent));
        phiCurrent = phiNew;
    end
end