close all; clear all;

k=7;
% number of intervals to split [-pi,pi]
M = 15000*2+1;
theta0Idx = (M-1)/2+1;
epsilon=0.0001;
% solve the integral equation \phi(\theta):
phiVals = solveIntEq(k,M,epsilon)';

% compute rho again, using the final values of phi:
rhoVals = zeros(1,M)';
for idx = 1:M
    rhoVals(idx) = computeRho(phiVals(idx),k);
end

%%
% compute G(u) and g(u):
rhoCubed = rhoVals.^3;
denom = (1-rhoCubed.*exp(-1i*phiVals));
GVals = -(1-rhoCubed.*exp(1i*phiVals))./denom;
gVals = (1-rhoCubed.*cos(phiVals))./denom;
GVals(theta0Idx) = 1;
gVals(theta0Idx) = 0;

% compute log(w^{-1}G(w))
thetaVals = linspace(-pi,pi,M)';
wVals = exp(1i*thetaVals);
logwGw = log((1./wVals).*GVals);
%%
% compute H(u):
HVals = zeros(1,M)';
deltas = ones(1,M)*2*pi/M;
deltas = deltas';
for idx = 1:M
    if idx == 1
       % don't include the last point (singularity) 
       inds = [1:idx-1 idx+1:M-1];
    else
        inds = [1:idx-1 idx+1:M];
    end
    HVals(idx) = logwGw(idx) + (1/(2*pi*1i))*sum(deltas(inds)*1i.*wVals(inds).*(logwGw(inds)-logwGw(idx))./(wVals(inds)-wVals(idx)));
end

%% sanity check: should be 0.5
inds = [1:theta0Idx-1 theta0Idx+1:M];
san_check = (1/(2*pi*1i))*sum(deltas(inds)*1i.*wVals(inds)./(wVals(inds)-1));

%%
% compute the integrals I1, I2, I3, and I4
% ignore the value of theta=0
inds = [1:(M-1)/2 (M-1)/2+2:M];
I1 = exp((1/(2*pi*1i))*sum(deltas(inds)*1i.*wVals(inds).*logwGw(inds)./(wVals(inds)-1)));
I2 = exp((1/(2*pi*1i))*sum(deltas(inds)*1i.*wVals(inds).*logwGw(inds)./wVals(inds)));
I3 = (1/(2*pi*1i))*sum(deltas(inds)*1i.*wVals(inds).*exp(-HVals(inds)).*gVals(inds)./(wVals(inds)-1));
I4 = (1/(2*pi*1i))*sum(deltas(inds)*1i.*wVals(inds).*exp(-HVals(inds)).*gVals(inds)./wVals(inds));

%%
% compute the objectives
Delta = 2*I1*I2*(I3-I4)+2*I1*(1+I2)-I2;
F00 = (1/Delta)*(2/3)*((k-3)/(k-2))*I2;
F10 = (I1/Delta)*(2/3)*((k-3)/(k-2))*(I2*(I3-I4)+I2+1);
%2*F10-F00

%%
% compute alpha' and some other useful values
gamma0pr1 = 8*sqrt(2*(k-3))/(sqrt(6-k+sqrt((k+6)*(k-2)))*(3*sqrt(k-2)+sqrt(k+6)));
rhoSquared = rhoVals.^2;
rhoQuadred = rhoVals.^4;
% Eq. (60)
rhoprDenom = 3.*rhoVals.*cos(phiVals)-k;
rhoprVals = rhoSquared.*sin(phiVals)./rhoprDenom;
% Eqs. (55) and (56)
gprDenom = (1-rhoCubed.*exp(-1i*phiVals)).^2;
gprNumerPart = (rhoQuadred-rhoVals.*cos(phiVals)-3*rhoprVals.*sin(phiVals));
GprVals = -2i*rhoSquared.*gprNumerPart./gprDenom;
% Gpr(0) is special:
GprVals(theta0Idx) = 1i*k/(k-3);
gprVals = 1i*rhoSquared.*gprNumerPart./gprDenom;
% gpr(0) is also special:
gprVals(theta0Idx) = -1i*k/(2*(k-3));
% compute gamma_0'(w)
K=k-2;
omegaVals = rhoVals.*exp(1i*phiVals);
omegaSquared = omegaVals.^2;
omegaCubed = omegaVals.^3;
omegaQuadred = omegaVals.^4;
delta = (k-2-sqrt((k-2)^2+8*(k-2)))/4;
gamma01Pr = (omegaVals-1).*(K*omegaVals-2*omegaSquared+K)./omegaSquared;
gamma02Pr = (-K*delta^3*omegaSquared+3*K*delta^2*omegaCubed+2*delta^3*omegaCubed-4*delta^2*omegaQuadred-K*delta^3*omegaVals-K*delta^2*omegaSquared+2*delta^2*omegaCubed-2*K^2*omegaSquared-K*delta*omegaSquared+3*K*omegaCubed+2*K^2*delta-K*delta*omegaVals-K*omegaSquared)./(2*omegaSquared.*sqrt((K-omegaVals*delta^2).*(K-omegaVals)));
gamma0Pr = -2*(delta*gamma01Pr+gamma02Pr)/((1-delta)^2*(k-2-delta));
% Eq. (62)
phiprVals = 1i*exp(1i*(thetaVals-phiVals))./((rhoprVals+1i*rhoVals).*(gamma0Pr));
% Eq. (51): alpha'
alphaPrVals = -1i*exp(-1i*thetaVals).*(-1i+phiprVals.*GprVals./GVals);

%%
% compute beta'
% Eq. (47)
HprVals = zeros(1,M)';
deltas = ones(1,M)*2*pi/M;
deltas = deltas';
for idx = 1:M
    if idx == 1
       % don't include the last point (singularity) 
       inds = [1:idx-1 idx+1:M-1];
    elseif idx == M
        % don't inclide the first or the last point (singularities)
        inds = 2:M-1;
    else
        inds = [1:idx-1 idx+1:M];
    end
    HprVals(idx) = alphaPrVals(idx) + (1/(2*pi*1i))*sum(deltas(inds)*1i.*wVals(inds).*(logwGw(inds)-logwGw(idx)-alphaPrVals(idx).*(wVals(inds)-wVals(idx)))./(wVals(inds)-wVals(idx)).^2);
end
% beta'
betaPrVals = -1i*exp(-1i*thetaVals).*(phiprVals.*gprVals-1i*wVals.*gVals.*HprVals).*exp(-HVals);

%%
% compute I5 and I6
% phi'(0)
phi0 = 0; rho0 = 1;
phipr0 = 1/gamma0pr1;
% alpha'(1)
alphapr1 = phipr0*k/(k-3)-1;
% g'(1)
gpr1 = (1/2)*(k/(k-3))^2;
% H(0)
H0 = HVals(theta0Idx);
% beta'(1)
betapr1 = betaPrVals(theta0Idx);

inds = [1:(M-1)/2 (M-1)/2+2:M];
I5 = (1/(2*pi*1i))*sum(deltas(inds)*1i.*wVals(inds).*(alphaPrVals(inds)-alphapr1)./(wVals(inds)-1));
I6 = (1/(2*pi*1i))*sum(deltas(inds)*1i.*wVals(inds).*(betaPrVals(inds)-betapr1)./(wVals(inds)-1));
%%
% compute d/dx F(x,0)|x=1
dFx0dxat1 = gamma0pr1*(F10*(alphapr1+I5)+F00*I1*(betapr1+I6+1));

%% compute gamma''
gammaPrPr01 = 2*(K-omegaCubed)./omegaCubed;
gammaPrPr02denom = 4*omegaCubed.*((K-omegaVals*delta^2).*(K-omegaVals)).^(3/2);
omegaQuinted = omegaVals.^5;
omegaSexted = omegaVals.^6;
gammaPrPr02numer = K^2*delta^5*omegaCubed + 3*K^2*delta^4*omegaQuadred - 12*K*delta^4*omegaQuinted + 8*delta^4*omegaSexted + 3*K^2*delta^5*omegaSquared + K^2*delta^4*omegaCubed - 4*K*delta^5*omegaCubed - 4*K^3*delta^2*omegaCubed - 2*K^2*delta^3*omegaCubed + 18*K^2*delta^2*omegaQuadred - 12*K*delta^2*omegaQuinted - 12*K^3*delta^3*omegaVals + 18*K^2*delta^3*omegaSquared - 2*K^2*delta^2*omegaCubed - 4*K*delta^3*omegaCubed - 4*K^3*omegaCubed +K^2*delta*omegaCubed + 3*K^2*omegaQuadred + 8*K^4*delta - 12*K^3*delta*omegaVals +3*K^2*delta*omegaSquared + K^2*omegaCubed;
gammaPrPr02 =  -gammaPrPr02numer./gammaPrPr02denom;
gamma0PrPr = (-2/((1-delta)^2*(K-delta)))*(delta*gammaPrPr01 + gammaPrPr02);

%% compute rho''
rhoPrPr_denom = 3*rhoVals.*cos(phiVals) - k;
rhoPrPr_1 = 5*rhoVals.*(sin(phiVals).^2)./rhoPrPr_denom;
rhoPrPr_2 = 3*rhoSquared.*(sin(phiVals).^2).*cos(phiVals)./(rhoPrPr_denom.^2);
rhoPrPr = (rhoSquared./rhoPrPr_denom).*(cos(phiVals) + rhoPrPr_1 - rhoPrPr_2);

%% compute phi''(theta)
expThetaMinPhi = exp(1i*(thetaVals-phiVals));
phiPrPrPart1 = expThetaMinPhi./((gamma0Pr.*(rhoprVals+1i*rhoVals).^2).^2);
phiPrPrPart2 = wVals.*phiprVals.*(-1i*rhoSquared - 2*rhoprVals.*rhoVals + 1i*rhoprVals.^2).*gamma0PrPr + (-2*phiprVals.*rhoprVals + 1i*phiprVals.*rhoPrPr - 1i*phiprVals.*rhoVals + 1i*rhoVals + rhoprVals).*gamma0Pr;
phiPrPr = phiPrPrPart1.*phiPrPrPart2;

%% compute G'' and g''
expPhi = exp(-1i*phiVals);
expPhi2 = exp(-2i*phiVals);
GPrPr_denom = (1-rhoCubed.*expPhi).^3;
GPrPr_1 = rhoVals./GPrPr_denom;
GPrPr_2 = -4*(rhoVals.^8).*expPhi + (rhoVals.^5).*(3+expPhi2) + 3*rhoQuadred.*(2i*rhoprVals.*(expPhi2 - 3)-rhoPrPr.*(1-expPhi2)) + 12*rhoCubed.*(rhoprVals.^2).*(1-expPhi2) - 2i*rhoSquared.*sin(phiVals) + 6i*rhoVals.*(2*rhoprVals.*cos(phiVals) + rhoPrPr.*sin(phiVals)) + 12i*(rhoprVals.^2).*sin(phiVals);
GPrPr = GPrPr_1.*GPrPr_2;

gPrPr = -GPrPr/2;

% special values of G'' and g'' at 0:
GPrPr(theta0Idx) = -(k/(k-3))^2;
gPrPr(theta0Idx) = (k/(k-3))^2/2;

%% compute alpha''
alphaPrPr = -exp(-1i*thetaVals).*alphaPrVals - exp(-2i*thetaVals).*(phiPrPr.*GprVals./GVals + (phiprVals.^2).*GPrPr./GVals - (phiprVals.^2).*(GprVals./GVals).^2);
alphaPrPr(theta0Idx) = -(phipr0*k/(k-3) - 1)+1i*phiPrPr(theta0Idx)*k/(k-3);

% compute H''
HPrPrVals = zeros(1,M)';
for idx = 1:M
    if idx == 1
        % don't include the last point (singularity) 
        inds = [1:idx-1 idx+1:M-1];
    elseif idx == M
    % don't inclide the first or the last point (singularities)
    inds = 2:M-1;
    else
        inds = [1:idx-1 idx+1:M];
    end
    HPrPrVals(idx) = alphaPrPr(idx) + (1/(2*pi*1i))*sum(deltas(inds)*1i.*wVals(inds).*(2*(logwGw(inds)-logwGw(idx))-2*alphaPrVals(idx)*(wVals(inds)-wVals(idx))-alphaPrPr(idx)*(wVals(inds)-wVals(idx)).^2)./(wVals(inds)-wVals(idx)).^3);
end

%% compute beta''
betaPrPr = -exp(-2i*thetaVals).*(-1i*phiprVals.*gprVals + phiPrPr.*gprVals + (phiprVals.^2).*gPrPr - 1i*wVals.*phiprVals.*gprVals.*HprVals + exp(2i*thetaVals).*gVals.*HPrPrVals - wVals.*(1i*phiprVals.*gprVals + wVals.*gVals.*HprVals).*HprVals).*exp(-HVals);

% compute I7 and I8
inds = [1:(M-1)/2 (M-1)/2+2:M];
I7 = (1/(2*pi*1i))*sum(deltas(inds)*1i.*wVals(inds).*(alphaPrPr(inds)-alphaPrPr(theta0Idx))./(wVals(inds)-1));
I8 = (1/(2*pi*1i))*sum(deltas(inds)*1i.*wVals(inds).*(betaPrPr(inds)-betaPrPr(theta0Idx))./(wVals(inds)-1));
% compute d^2 F(x,0)/dx^2 at x=1
d2Fx0dx2at1 = (gamma0PrPr(theta0Idx)/gamma0pr1)*k/(6*(k-2)) + gamma0pr1*((k/(6*(k-2)))*(alphapr1 + I5) + F10*gamma0pr1*(alphaPrPr(theta0Idx)+I7) + F00*gamma0pr1*I1*((alphapr1+I5)*(betapr1+I6+1)+betaPrPr(theta0Idx)+I8));

% compute d^2 F(x,x)/dx^2 at x=1
EN = k/(k-3);
d2Fxxdx2at1 = (1/(k-3))*(2*k/3 - (k-6)*EN + 6*(k-2)*dFx0dxat1 + 3*(k-2)*d2Fx0dx2at1);

% compute the variance
varN = d2Fxxdx2at1 - 3*k/(k-3)^2;