%% ASSIGNMENT 1 EXERCISE 3 (GUIDANCE; CONTINUOUS GUIDANCE)

% Author: Massini Alessandro - 10990413

clearvars; close all; clc;
format long g
rng default

applyDefaultSettings();

tic
%% ------------------------ EXERCISE 3.1 --------------------------------%%
% --------------- Initial and Final state computation --------------------%

% PARAMETERS INITIALIZATION
hi = 800;                   %[km]
hf = 1000;                  %[km]
di = deg2rad(0.75);         %[degrees]
Re = 6378.1366;             %[km]
mu = 398600.435;            %[km^3/s^2]
rho0 = 750 + Re;            %{km]
k1 = 1e-05;                 %[-]
k2 = 1e-04;                 %[-]
mass = 1000;                %[kg]
T = 3;                      %[N]
Isp = 3120;                 %[s]
g0 = 9.81;                  %[m/s]

DU = 7178.1366;             %[km]
MU = mass;                    %[kg]

% VECTOR OF DISTANCE FROM THE EARTH
rho = linspace(hi+Re-100, hf+Re+100,1000);

% COMPUTATION OF DEBRIES SPATIAL DENSITY AS A FUNCTION OF RHO
q = k1 ./ (k2 + ((rho - rho0)./DU).^2);
[~,imax] = max(q);

% PLOT OF DEBRIES SPATIAL DENSITY FUNCTION
figure
plot(rho-Re,q, 'k')
hold on
xline(hi,'--', LineWidth=2, Color= [0 0.4470 0.7410])
xline(hf,'--', LineWidth=2, Color=[0.4660 0.6740 0.1880])
xline(rho(imax)-Re, 'r--', LineWidth=2)
xlabel('Altitude [km]')
ylabel('Debris spatial density[$\frac{1}{DU^3}$]')
legend('Debris spatial density','Initial Altitude', 'Final Altitude', 'Maximum Density')


% INITIAL STATE COMPUTATION
initialState = [Re+hi;0;0;0;sqrt(mu/(Re+hi));0];

% FINAL VELOCITY COMPONENTS COMPUTATION
finalVy = sqrt(mu/(Re+hf)) * cos(di);
finalVz = sqrt(mu/(Re+hf)) * sin(di);

% FINAL STATE COMPUTATION
finalState = [Re+hf;0;0;0;finalVy;finalVz];

%% ------------------------ EXERCISE 3.2 --------------------------------%%
% ----------------------- Adimensionalization ----------------------------%

% COMPUTE TU KNOWING THAT MU ADIM = 1
TU = sqrt(DU^3 / mu);
VU = DU/TU;

% SET MU = 1
mu = 1;

% ADIMENSIONALIZE INITIAL AND FINAL STATE
initialState = [initialState(1:3)./DU;initialState(4:6)./VU];
finalState   = [finalState(1:3)./DU;finalState(4:6)./VU];

% ADIMENSIONALIZE PARAMETERS
mass = mass/MU;
Isp = Isp/TU;
T = T * TU ^ 2 / (MU * DU * 10^3);
g0 = g0 * TU^2 / (DU * 10^3);
rho0 = rho0 / DU;

% CREATE PARAMETERS STRUCTURE
parameters = parametersDefinition(mu,T,Isp,g0,k1,k2,rho0);

%% -------------------- EXERCISE 3.3 and 3.4 --------------------------- %%
% ------------------ Write PMP and solve the Problem ---------------------%

% INITIALIZE INITIAL TIME AND FIRST GUESS OF FINAL TIME
t0 = 0;
initialFinalTime = 20 * pi;

% GENERATE RANDOM LAMBDA WITHIN THE APPROPRIATE RANGE
lambdaRange = [-250 250];
initialLambda = lambdaRange(1)+(lambdaRange(2)-lambdaRange(1))*rand(6,1);
lambdaM = 250 * rand(1,1);
initialLambda = [initialLambda;lambdaM];

% SET THE INITIAL GUESS
initialGuess = [initialLambda;initialFinalTime];

% TARGET [FINAL STATE; FINAL LAMBDA_M; FINAL HAMILTONIAN]
finalTarget = [finalState;0;0];

% ADD MASS TO INITIAL STATE
initialState = [initialState;mass];

% SET INITIAL EXITFLAG TO 0
exitFlag = 0;

% ITERATION UNTIL COVERGENCE IS REACHED BY FSOLVE
while exitFlag <= 0

    % DEFINE THE ZERO FINDING FUNCTION
    fun = @(initialGuess) zeroFindingFunction(t0, initialState, initialGuess, ...
        finalState, parameters);
 
    % SET THE OPTIONS
    options = optimoptions('fsolve', 'MaxFunctionEvaluations', 10000, ...
        'MaxIterations', 10000, 'Display','iter-detailed',...
    'Algorithm', 'levenberg-marquardt', 'FunctionTolerance', 1e-10, ...
    'OptimalityTolerance', 1e-10, 'StepTolerance',1e-10,'SpecifyObjectiveGradient',true);
    
    % SOLVER
    [solution, fval,exitFlag] = fsolve(fun,initialGuess,options);

    % IF CONVERGENCE IS NOT REACHED OR THE SOLUTION IS THE SUB-OPTIMAL
    if exitFlag <= 0 || solution(end) > 22*pi

        % NEW GENERATION OF LAMBDA
        lambdaRange = [-250 250];
        initialLambda = lambdaRange(1)+(lambdaRange(2)-lambdaRange(1))*rand(6,1);
        lambdaM = 250 * rand(1,1);
        initialLambda_new = [initialLambda;lambdaM];

        % UPDATE THE INITIAL GUESS
        initialGuess = [initialLambda_new;initialFinalTime];

        % UPDATE EXITFLAG FOR THE SUB-OPTIMAL SOLUTIONS
        exitFlag = 0;

        rng shuffle

    end
end

% CONCATENATION
state = [initialState;solution(1:end-1)];
tf = solution(end);

% COMPUTATION OF FINAL TIME IN MINUTES
finalTime = tf * TU / 60;

% PROPAGATION OF THE OPTIMAL SOLUTION
[~, ~, ~, propagate, tt] = propagatorEL(t0, state, tf, parameters);

% COMPUTE REMAINING MASS
finalMass = propagate(end,7) * MU;

% COMPUTE POSITION ERRORS IN KM
errPos = norm(propagate(end,1:3)' - finalState(1:3)) * DU;
fprintf('Position Error: %.10f [km]\n',errPos)

% COMPUTE VELOCITY ERRORS IN KM/S
errVel = norm(propagate(end,4:6)' - finalState(4:6)) * VU * 1000;
fprintf('Velocity Error: %.15f [km/s]\n',errVel)

% HAMILTONIAN INITIALIZATION
H = zeros(length(tt),1);

% COMPUTE HAMILTONIAN
for i = 1:length(propagate)
    
    H(i) = hamiltonian(propagate(i,:),parameters);

end

% PLOT HAMILTONIAN
figure
plot(tt,H, 'k')
xlabel('Time[TU]')
ylabel('Hamiltonian[-]')

% COMPUTE PRIMER VECTOR AT EACH TIME
alpha = - propagate(:,11:13)'./ vecnorm(propagate(:,11:13)');

% COMPUTE NORM OF PRIMER VECTOR AT EACH TIME
normalp = vecnorm(alpha);

% INITIALIZE PRIMER VECTOR IN NTW
alphaNTW = zeros(length(alpha), 3);
for i = 1:length(alpha)

    % COMPUTE PRIMER VECTOR IN NTW REFERENCE FRAME
    alphaNTW(i,:) = ECI2NTW(propagate(i,1:6),alpha(:,i));

end

% PLOT PRIMER VECTOR IN TIME IN NTW REFERENCE FRAME
figure
plot(tt, alphaNTW(:,1))
hold on
plot(tt, alphaNTW(:,2))
plot(tt, alphaNTW(:,3))
plot(tt, normalp)
xlabel('Time[TU]')
ylabel('$\alpha$[-]')
legend('$\alpha_N$','$\alpha_T$', '$\alpha_W$', 'norm($\alpha$)')

% PLOT ORBIT IN ECI
figure
plot3(propagate(2*156:end-2*156,1)*DU,propagate(2*156:end-2*156,2)*DU,propagate(2*156:end-2*156,3)*DU, 'k')
hold on 
plot3(propagate(1:2*156,1)*DU,propagate(1:2*156,2)*DU,propagate(1:2*156,3)*DU, 'r')
plot3(propagate(end-2*300:end,1)*DU,propagate(end-2*300:end,2)*DU,propagate(end-2*300:end,3)*DU, 'b')
plot3(propagate(1,1)*DU,propagate(1,2)*DU,propagate(1,3)*DU,'o','MarkerSize',10, 'MarkerEdgeColor','k','MarkerFaceColor','r');
plot3(propagate(end,1)*DU,propagate(end,2)*DU,propagate(end,3)*DU,'o','MarkerSize',12, 'MarkerFaceColor','b','MarkerEdgeColor','k')
%plot3(0,0,0,'Marker','o','MarkerSize',40, 'MarkerFaceColor','k','MarkerEdgeColor','k')
xlabel('x[km]')
ylabel('y[km]')
zlabel('z[km]')
legend('Trajectory', 'First Orbit', 'Last Orbit','Initial State', 'Final State')
title('@Earth J2000')

%% ------------------------ EXERCISE 3.5 --------------------------------%%
% --------------------- Numerical Continuation ---------------------------%

% DEFINE FINAL THRUST
Tend = 2.860 * TU ^ 2 / (MU * DU * 10^3);

% VECTOR OF DIFFERENT THRUST
TVec = linspace(T, Tend,10);

% INITIALIZE THE GUESS
guess = solution;

% OPTIONS FOR NUMERICAL CONTINUATION
optionsContinuation = optimoptions('fsolve', 'MaxFunctionEvaluations', 10000, ...
    'MaxIterations', 1000, 'Display','iter-detailed', 'Algorithm', 'trust-region-dogleg', ...
    'FunctionTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10, ...
    'SpecifyObjectiveGradient',true);

% ITERATE FOR EACH THRUST LEVEL
for i = 2 : length(TVec)

    % UPDATE THE VALUE OF THE THRUST
    parameters.T = TVec(i);

    % ZERO FINDING FUNCTION
    fun = @(guess) zeroFindingFunction(t0, initialState, guess, finalState, parameters);

    % SOLVER
    [solution(:,i), fval(:,i)] = fsolve(fun,guess,optionsContinuation);

    % UPDATE THE INITIAL GUESS
    guess = solution(:,i);

end

% SAVE THE FINAL STATE FOR T = 2.860
state = [initialState;solution(1:end-1, end)];
tfContinuation = solution(end, end);

% PROPAGATE
[~, ~, ~, propagateContinuation, ttContinuation] = propagatorEL(t0, state, tfContinuation, parameters);

%COMPUTE FINAL MASS
finalMassContinuation = propagateContinuation(end,7) * MU;

% INITIALIZE THE HAMILTONIAN
HContinuation = zeros(length(ttContinuation),1);
for i = 1:length(ttContinuation)
    
    % COMPUTE HAMILTONIAN FOR EACH TIME
    HContinuation(i) = hamiltonian(propagateContinuation(i,:),parameters);

end

% COMPUTE PRIMER VECTOR AT EACH TIME
alphaContinuation = - propagateContinuation(:,11:13)'./ vecnorm(propagateContinuation(:,11:13)');

% COMPUTE NORM OF ALPHA
normalpContinuation = vecnorm(alphaContinuation);

% INITIALIZE PRIMER VECTOR IN NTW
alphaNTWContinuation = zeros(length(alphaContinuation), 3);
for i = 1:length(alphaContinuation)

    % COMPUTE PRIMER VECTOR IN NTW REFERENCE FRAME
    alphaNTWContinuation(i,:) = ECI2NTW(propagateContinuation(i,1:6),alphaContinuation(:,i));

end

% PLOT PRIMER VECTOR IN TIME IN NTW REFERENCE FRAME
figure
plot(ttContinuation, alphaNTWContinuation(:,1))
hold on
plot(ttContinuation, alphaNTWContinuation(:,2))
plot(ttContinuation, alphaNTWContinuation(:,3))
plot(ttContinuation, normalpContinuation)
xlabel('Time[TU]')
ylabel('$\alpha$[-]')
legend('$\alpha_N$','$\alpha_T$', '$\alpha_W$', 'norm($\alpha$)')

% PLOT THE HAMILTONIAN
figure
plot(ttContinuation,HContinuation,'k')
xlabel('Time[TU]')
ylabel('Hamiltonian[-]')

% COMPUTE POSITION ERROR
errPos = norm(propagateContinuation(end,1:3)' - finalState(1:3)) * DU;
fprintf('Position Error: %.10f [km]\n',errPos)

% COMPUTE VELOCITY ERRORS IN KM/S
errVel = norm(propagateContinuation(end,4:6)' - finalState(4:6)) * VU * 1000;
fprintf('Velocity Error: %.15f [km/s]\n',errVel)
toc

%% ---------------------- Functions ------------------------------------ %%

% SETTINGS
function applyDefaultSettings()
    set(groot, 'defaultTextInterpreter', 'latex')
    set(groot,'defaultAxesXMinorGrid','on','defaultAxesXMinorGridMode','manual');
    set(groot,'defaultAxesYMinorGrid','on','defaultAxesYMinorGridMode','manual');
    set(groot, 'defaultLegendLocation', 'northeast')
    set(groot, 'defaultLegendInterpreter', 'latex')
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
    set(groot, 'defaultAxesFontWeight', 'bold')
    set(groot, 'defaultFigurePosition', [470, 360, 900, 530]-100)
    set(groot, 'defaultFigureColormap', turbo(256));
    set(groot, 'defaultAxesFontName', 'Palatino Linotype', 'defaultTextFontName', 'Palatino Linotype');
    set(groot, 'defaultSurfaceEdgeAlpha', 0.3);
    set(groot, 'defaultLineLineWidth', 2);
    set(groot, 'defaultFigureColor', [1; 1; 1]);
    set(groot, 'defaultAxesColor', 'none');
    set(groot,'DefaultAxesYLimitMethod',"padded")
    set(groot, 'defaultAxesFontSize',20);
end

% EULER LAGRANGE EQUATIONS WITH STM
function dxdt = eulerLagrangianRHS_STM(~, xx, parameters)

% eulerLagrangianRHS_STM.m - Returns the euler lagrange equations given the
%                            state and the parameters
%
%
% Prototype:  dxdt = eulerLagrangianRHS_STM(~, xx, parameters)
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs:     
%
%
% xx[210, 1]:  State of the system (rr,vv,m,llx,llv,llm,STM)
%
% t[nx1]:       Vector of the time grid
%
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                        mu   = Earth gravitational constant;
%                        T    = Thrust;
%                        g0   = Gravity Acceleration;
%                        Isp  = Specific Impulse;
%                        k1   = Debris Spatial Density Constant 1;
%                        k2   = Debris Spatial Density Constant 1;
%                        rho0 = Reference Radius;
%
% 
%
%
% Outputs:
%
%
% dxdt[210x1]:    Derivative of the state
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% Initialization
T = parameters.T;
mu = parameters.mu;
Isp = parameters.Isp;
g0 = parameters.g0;
k1 = parameters.k1;
k2 = parameters.k2;
rho0 = parameters.rho0;

% State
rr = xx(1:3);
vv = xx(4:6);
m  = xx(7);

% Costate
lambdaR = xx(8:10);
lambdaV = xx(11:13);

% STM
Phi = reshape(xx(15:end),14,14);

% Compute Norm
r = norm(rr);
lv = norm(lambdaV);

% Derivative of q with repect to variables
dqdx = [(2*k1*rr(1)*(rho0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2)))/((k2 + (rho0 - (rr(1)^2 + rr(2)^2 ...
        + rr(3)^2)^(1/2))^2)^2*(rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2));
        (2*k1*rr(2)*(rho0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2)))/((k2 + (rho0 - (rr(1)^2 + rr(2)^2 ...
        + rr(3)^2)^(1/2))^2)^2*(rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2));
        (2*k1*rr(3)*(rho0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2)))/((k2 + (rho0 - (rr(1)^2 + rr(2)^2 ...
        + rr(3)^2)^(1/2))^2)^2*(rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2))];
 
% Compute the Jacobian of the dynamics
J = jacobianOfDynamics(xx(1:14),parameters);

% Compute the derivative of the STM
Phidot = J*Phi;

% Initialize dxdt
dxdt = zeros(210,1);

% Compute the derivatives
dxdt(1:3)   = vv;

dxdt(4:6)   = - mu / r^3 * rr - T / m * lambdaV / lv;

dxdt(7)     = - T / (Isp * g0);

dxdt(8:10)  = - dqdx - 3 * mu / r^5 * (dot(rr,lambdaV)) * rr + mu / r^3 * lambdaV;

dxdt(11:13) = - lambdaR;

dxdt(14)    = - T / m^2 * lv;

dxdt(15:end) = Phidot(:);



end

% EULER LAGRANGE EQUATIONS
function dxdt = eulerLagrangianRHS(~, xx, parameters)

% eulerLagrangianRHS.m - Returns the euler lagrange equations given the
%                        state and the parameters
%
%
% Prototype:  dxdt = eulerLagrangianRHS_STM(~, xx, parameters)
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs:     
%
%
% xx[14, 1]:  State of the system (rr,vv,m,llx,llv,llm)
%
% t[nx1]:       Vector of the time grid
%
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                        mu [1]   = Earth gravitational constant;
%                        T [1]    = Thrust;
%                        g0 [1]   = Gravity Acceleration;
%                        Isp [1]  = Specific Impulse;
%                        k1  [1]  = Debris Spatial Density Constant 1;
%                        k2 [1]   = Debris Spatial Density Constant 1;
%                        rho0 [1] = Reference Radius;
%
% 
%
%
% Outputs:
%
%
% dxdt[14x1]:    Derivative of the state
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% Initialize
T = parameters.T;
mu = parameters.mu;
Isp = parameters.Isp;
g0 = parameters.g0;
k1 = parameters.k1;
k2 = parameters.k2;
rho0 = parameters.rho0;

% State
rr = xx(1:3);
vv = xx(4:6);
m  = xx(7);

% Costate
lambdaR = xx(8:10);
lambdaV = xx(11:13);

% Norm
r = norm(rr);
lv = norm(lambdaV);

% Compute derivative of q
dqdx = [(2*k1*rr(1)*(rho0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2)))/((k2 + (rho0 - (rr(1)^2 + rr(2)^2 ...
        + rr(3)^2)^(1/2))^2)^2*(rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2));
        (2*k1*rr(2)*(rho0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2)))/((k2 + (rho0 - (rr(1)^2 + rr(2)^2 ...
        + rr(3)^2)^(1/2))^2)^2*(rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2));
        (2*k1*rr(3)*(rho0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2)))/((k2 + (rho0 - (rr(1)^2 + rr(2)^2 ...
        + rr(3)^2)^(1/2))^2)^2*(rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2))];
 
% Compute derivative of state
dxdt = zeros(14,1);

dxdt(1:3)   = vv;

dxdt(4:6)   = - mu / r^3 * rr - T / m * lambdaV / lv;

dxdt(7)     = - T / (Isp * g0);

dxdt(8:10)  = - dqdx - 3 * mu / r^5 * (dot(rr,lambdaV)) * rr + mu / r^3 * lambdaV;

dxdt(11:13) = - lambdaR;

dxdt(14)    = - T / m^2 * lv;

end

% PROPAGATOR 
function [xf, tf, PHIf, xx, tt] = propagatorEL(t0, x0, tf, parameters)

% propagatorEL.m - Perform propagation of the state integrating the
%                  Euler-Lagrange Equations
%
%
% Prototype: [xf, tf, PHIf, xx, tt] = propagatorEL(t0, x0, tf, parameters)
%
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs:     
%
%
% t0[1]:    Initial time of propagation 
%
% x0[14x1]:  Initial state of the system (r,v,m,llr,llv,llm) 
%
% tf[1]:    Final time of propagation 
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                        mu [1]   = Earth gravitational constant;
%                        T [1]    = Thrust;
%                        g0 [1]   = Gravity Acceleration;
%                        Isp [1]  = Specific Impulse;
%                        k1  [1]  = Debris Spatial Density Constant 1;
%                        k2 [1]   = Debris Spatial Density Constant 1;
%                        rho0 [1] = Reference Radius;
%            
%
% Outputs:
%
%
% xf[14x1]:   Final state of the system (r,v,m,llr,llv,llm)
% 
% tf[1]:     Final time of propagation
%
% PHIf[14x14]: State Transition Matrix of the final state
%
% xx[nx210]:   Full propagation of the initial state and of the elements of the STM (x,y,vx,vy)
%
% tt[nx1]:   Vector of time
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025


 % Initialize
 Phi0 = eye(14);

 % Initialize
 x0Phi0 = [x0; Phi0(:)];

 % Integration
 options = odeset('reltol', 1e-12, 'abstol', 1e-12);
 [tt, xx] = ode113(@(t,x) eulerLagrangianRHS_STM(t,x,parameters), [t0 tf], x0Phi0, options);

 % Extraction
 xf = xx(end,1:14)';
 PHIf = reshape(xx(end,15:end),14,14);
 tf = tt(end);

end

% HAMOLTONIAN COMPUTATION
function H = hamiltonian(x0, parameters)

% hamiltonian.m -  Compute the Hamiltonian of the state
%
%
% Prototype: H = hamiltonian(x0, parameters)
%
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs:     
%
%
% x0[14x1]:  Initial state of the system (r,v,m,llr,llv,llm) 
%
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                        mu [1]   = Earth gravitational constant;
%                        T [1]    = Thrust;
%                        g0 [1]   = Gravity Acceleration;
%                        Isp [1]  = Specific Impulse;
%                        k1  [1]  = Debris Spatial Density Constant 1;
%                        k2 [1]   = Debris Spatial Density Constant 1;
%                        rho0 [1] = Reference Radius;
%            
%
% Outputs:
% 
%
% H[1] = Hamiltonian
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% Initialization
T = parameters.T;
mu = parameters.mu;
Isp = parameters.Isp;
g0 = parameters.g0;
k1 = parameters.k1;
k2 = parameters.k2;
rho0 = parameters.rho0;

% State
rr = x0(1:3);
vv = x0(4:6);
m  = x0(7);

% Costate
lambdaR = x0(8:10);
lambdaV = x0(11:13);
lambdaM = x0(14);

% Norm
r = norm(rr);

% q as a function of r
q = k1 / (k2 + (r - rho0)^2);

% Compute Hamiltonian
H = q + dot(lambdaR,vv) - mu / r^3 * dot(lambdaV,rr) - T / m * norm(lambdaV) ...
    - lambdaM * T / (Isp * g0);


end

% DEFINITION OF THE PARAMETERS STRUCT
function parameters = parametersDefinition (mu,T,Isp,g0,k1,k2,rho0)

% Inputs =               mu [1]   = Earth gravitational constant;
%                        T [1]    = Thrust;
%                        g0 [1]   = Gravity Acceleration;
%                        Isp [1]  = Specific Impulse;
%                        k1  [1]  = Debris Spatial Density Constant 1;
%                        k2 [1]   = Debris Spatial Density Constant 1;
%                        rho0 [1] = Reference Radius;
%
%
% Outputs:            parameters [1x1 struct]: Struct with fields


 parameters.mu = mu;                
 parameters.T = T;
 parameters.g0 = g0;
 parameters.Isp = Isp;
 parameters.k1 = k1;
 parameters.k2= k2;
 parameters.rho0 = rho0;

end

% OBJECTIVE FUNCTION DEFINITION
function [f, grad] = zeroFindingFunction(t0, initialState, initialGuess, finalState, parameters)

% zeroFindingFunction.m - Compute the residuals both with their gradients
%                         of the objective function of the zero finding
%                         problem
%
%
% Prototype: [f, grad] = zeroFindingFunction(t0, initialState, initialGuess, finalState, parameters)
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs:     
%
%
% t0[1]: Initial time of integration
%
%
% initialState[7, 1]:  State of the system (rr0,vv0,m0)
%
%
% initialGuess[8, 1]:  Guess of the system (llr, llv,llm, tf)
%
%
% finalState[6x1]: final target of the problem (rrf,vvf)
%
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                        mu   = Earth gravitational constant;
%                        T    = Thrust;
%                        g0   = Gravity Acceleration;
%                        Isp  = Specific Impulse;
%                        k1   = Debris Spatial Density Constant 1;
%                        k2   = Debris Spatial Density Constant 1;
%                        rho0 = Reference Radius;
%
% 
% Outputs:
%
%
% f[8x1]:    Residuals of the zero finding problem
%
%
% grad[8x8]: jacobian matrix of the residuals with respect to the variables
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025


% Extract the variables
tf = initialGuess(end);
initialGuess = initialGuess(1:end-1);

% Propagate the guess
[propagate, ~, Phif, ~, ~] = propagatorEL(t0, [initialState;initialGuess], tf, parameters);

% Compute the Hamiltonian
H = hamiltonian(propagate, parameters);

% Define the function
f = [propagate(1:6) - finalState(1:6);
     propagate(end);
     H];

% Jacobian of the function
grad =  jacobianOfFunction(tf, propagate, Phif, parameters);

end

% COMPUTE THE JACOBIAN OF THE RESIDUALS FUNCTION
function grad = jacobianOfFunction(t, state, Phi, parameters)

% jacobianOfFunction.m - Compute the Jacobian of the zero finding function
%                        with respect to the variables 
%
%
% Prototype: grad = jacobianOfFunction(t, state, Phi, parameters)
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs:     
%
%
% t[1]: Time of integration
%
%
% state[14, 1]:  State of the system (rr,vv,m, llr, llv, llm)
%
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                        mu   = Earth gravitational constant;
%                        T    = Thrust;
%                        g0   = Gravity Acceleration;
%                        Isp  = Specific Impulse;
%                        k1   = Debris Spatial Density Constant 1;
%                        k2   = Debris Spatial Density Constant 1;
%                        rho0 = Reference Radius;
%
% 
% Outputs:
%
%
% grad[8x8]: jacobian matrix of the residuals with respect to the variables
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% Extraction of the submatrices
STM_rr_lr = Phi(1:3,8:10);
STM_rr_lv = Phi(1:3,11:13);
STM_rr_lm = Phi(1:3,14);

STM_vv_lr = Phi(4:6,8:10);
STM_vv_lv = Phi(4:6,11:13);
STM_vv_lm = Phi(4:6,14);

STM_lm_ll = [Phi(14,8:10) Phi(14,11:13) Phi(14,14)];

% Extraction of the derivative of the state
dxdt =  eulerLagrangianRHS(t, state, parameters);

f = zeros(7,1);
f(1:6) = dxdt(1:6);
f(7) = dxdt(14);

dHdll =  dxdt(1:7)' * Phi(8:14,8:14) - dxdt(8:14)' * Phi(1:7,8:14);

% Assemble of the Jacobian matrix
grad = [STM_rr_lr STM_rr_lv STM_rr_lm f(1:3);
        STM_vv_lr STM_vv_lv STM_vv_lm f(4:6);
        STM_lm_ll f(7);
        dHdll     0  ];

end

% ROTATION OF ALPHA IN THE NTW REFERENCE FRAME
function alphaNTW = ECI2NTW(stateECI, alphaECI)

% ECI2NTW.m - Rotation from ECI to NTW reference frame
%
%
% Prototype: grad = jacobianOfFunction(t, state, Phi, parameters)
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs:     
%
%
% state[6, 1]:  State of the system (rr,vv)
%
%
% alphaECI[3x1]: Primer vector on ECI frame
%
% 
% Outputs:
%
%
% alphaNTW[3x1]: Primer vector on NTW frame
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025
    
    % Initialization
    rr = stateECI(1:3);
    vv = stateECI(4:6);

    % Computation of the director cosine
    T = vv/norm(vv);
    N = cross(rr,vv)/norm(cross(rr,vv));
    W = cross(N,T)/norm(cross(N,T));

    % Rotation Matrix
    rotMat = [N; T; W];

    % Rotation
    alphaNTW = rotMat * alphaECI;
    
end

% JACOBIAN OF THE DYNAMICS
function J = jacobianOfDynamics(state, parameters)

% jacobianOfDynamics.m -  Compute the jacobian of the dynamics for the
%                         variational approach
%
%
% Prototype: J = jacobianOfDynamics(state, parameters)
%
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs:     
%
%
% state[14x1]:  Initial state of the system (r,v,m,llr,llv,llm) 
%
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                        mu [1]   = Earth gravitational constant;
%                        T [1]    = Thrust;
%                        g0 [1]   = Gravity Acceleration;
%                        Isp [1]  = Specific Impulse;
%                        k1  [1]  = Debris Spatial Density Constant 1;
%                        k2 [1]   = Debris Spatial Density Constant 1;
%                        rho0 [1] = Reference Radius;
%            
%
% Outputs:
% 
%
% J[14x14] = Jacobian of the dynamics with respect to the state
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

    % Initialization
    x = state(1);
    y = state(2);
    z = state(3);
    m = state(7);
    lvx = state(11);
    lvy = state(12);
    lvz = state(13);

    T = parameters.T;
    mu = parameters.mu;
    k1 = parameters.k1;
    k2 = parameters.k2;
    rho0 = parameters.rho0;

    J = zeros(14,14);

    % Assemble the matrix, the components have been retrieved through syms
    % toolbox
    J(1,4) = 1;
    J(2,5) = 1;
    J(3,6) = 1;

    J(4,1) = -(mu*(- 2*x^2 + y^2 + z^2))/(x^2 + y^2 + z^2)^(5/2);
    J(4,2) = (3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2);
    J(4,3) = (3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2);
    J(4,7) = (T*lvx)/(m^2*(lvx^2 + lvy^2 + lvz^2)^(1/2));
    J(4,11) = -(T*(lvy^2 + lvz^2))/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2));
    J(4,12) = (T*lvx*lvy)/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2));
    J(4,13) = (T*lvx*lvz)/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2));

    J(5,1) = (3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2);
    J(5,2) = -(mu*(x^2 - 2*y^2 + z^2))/(x^2 + y^2 + z^2)^(5/2);
    J(5,3) = (3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2);
    J(5,7) = (T*lvy)/(m^2*(lvx^2 + lvy^2 + lvz^2)^(1/2));
    J(5,11) = (T*lvx*lvy)/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2));
    J(5,12) = -(T*(lvx^2 + lvz^2))/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2));
    J(5,13) = (T*lvy*lvz)/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2));

    J(6,1) = (3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2);
    J(6,2) = (3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2);
    J(6,3) = -(mu*(x^2 + y^2 - 2*z^2))/(x^2 + y^2 + z^2)^(5/2);
    J(6,7) = (T*lvz)/(m^2*(lvx^2 + lvy^2 + lvz^2)^(1/2));
    J(6,11) = (T*lvx*lvz)/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2));
    J(6,12) = (T*lvy*lvz)/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2));
    J(6,13) = -(T*(lvx^2 + lvy^2))/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2));

    J(8,1) = (15*mu*x^2*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) ...
        - (3*mu*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(5/2) - (2*k1*(rho0 ...
        - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - ...
        (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(1/2)) + (2*k1*x^2)/((k2...
        + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)) - ...
        (6*lvx*mu*x)/(x^2 + y^2 + z^2)^(5/2) + (2*k1*x^2*(rho0 - ...
        (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + ...
        z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*x^2*(rho0 - ...
        (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2));
    
    J(8,2) = (2*k1*x*y)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2))...
         - (3*lvx*mu*y)/(x^2 + y^2 + z^2)^(5/2) - (3*lvy*mu*x)/(x^2 + y^2 + z^2)^(5/2)...
         + (15*mu*x*y*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) + (2*k1*x*y*(rho0...
         - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 +...
         y^2 + z^2)^(3/2)) - (8*k1*x*y*(rho0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho0 -...
         (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2));
    
    J(8,3) = (2*k1*x*z)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2))...
        - (3*lvx*mu*z)/(x^2 + y^2 + z^2)^(5/2) - (3*lvz*mu*x)/(x^2 + y^2 + z^2)^(5/2)...
        + (15*mu*x*z*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) + (2*k1*x*z*(rho0...
        - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 +...
        y^2 + z^2)^(3/2)) - (8*k1*x*z*(rho0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho0...
        - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2));
    
    J(8,11) = (mu*(- 2*x^2 + y^2 + z^2))/(x^2 + y^2 + z^2)^(5/2);
    J(8,12) = -(3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2);
    J(8,13) = -(3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2);
    
    J(9,1) = (2*k1*x*y)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 +...
        y^2 + z^2)) - (3*lvx*mu*y)/(x^2 + y^2 + z^2)^(5/2) - (3*lvy*mu*x)/(x^2 +...
        y^2 + z^2)^(5/2) + (15*mu*x*y*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2)...
        + (2*k1*x*y*(rho0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 +...
        z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*x*y*(rho0 - (x^2 + y^2 +...
        z^2)^(1/2))^2)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2));
     
    J(9,2) = (15*mu*y^2*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2)...
        - (3*mu*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(5/2) - (2*k1*(rho0 -...
        (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 +...
        y^2 + z^2)^(1/2)) + (2*k1*y^2)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 +...
        y^2 + z^2)) - (6*lvy*mu*y)/(x^2 + y^2 + z^2)^(5/2) + (2*k1*y^2*(rho0 - (x^2 + y^2 +...
        z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2))...
        - (8*k1*y^2*(rho0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho0 - (x^2 + y^2 + ...
        z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2));
    
    J(9,3) = (2*k1*y*z)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 +...
        y^2 + z^2)) - (3*lvy*mu*z)/(x^2 + y^2 + z^2)^(5/2) - (3*lvz*mu*y)/(x^2 +...
        y^2 + z^2)^(5/2) + (15*mu*y*z*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) +...
        (2*k1*y*z*(rho0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + ...
        z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*y*z*(rho0 - (x^2 + y^2 + ...
        z^2)^(1/2))^2)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2));
     
    J(9,11) = -(3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2);
    J(9,12) = (mu*(x^2 - 2*y^2 + z^2))/(x^2 + y^2 + z^2)^(5/2);
    J(9,13) = -(3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2);
    
    J(10,1) = (2*k1*x*z)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 +...
        y^2 + z^2)) - (3*lvx*mu*z)/(x^2 + y^2 + z^2)^(5/2) - (3*lvz*mu*x)/(x^2 +...
        y^2 + z^2)^(5/2) + (15*mu*x*z*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 +...
        z^2)^(7/2) + (2*k1*x*z*(rho0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 -...
        (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*x*z*(rho0 -...
        (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2));
     
    J(10,2) = (2*k1*y*z)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 +...
        y^2 + z^2)) - (3*lvy*mu*z)/(x^2 + y^2 + z^2)^(5/2) - (3*lvz*mu*y)/(x^2 +...
        y^2 + z^2)^(5/2) + (15*mu*y*z*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 +...
        z^2)^(7/2) + (2*k1*y*z*(rho0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 -...
        (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*y*z*(rho0 -...
        (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2));
    
    J(10,3) = (15*mu*z^2*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2)...
        - (3*mu*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(5/2) - (2*k1*(rho0 -...
        (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 +...
        y^2 + z^2)^(1/2)) + (2*k1*z^2)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 +...
        y^2 + z^2)) - (6*lvz*mu*z)/(x^2 + y^2 + z^2)^(5/2) + (2*k1*z^2*(rho0 - (x^2 +...
        y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 +...
        z^2)^(3/2)) - (8*k1*z^2*(rho0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho0 -...
        (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2));
    
    J(10,11) = -(3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2);
    J(10,12) = -(3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2);
    J(10,13) = (mu*(x^2 + y^2 - 2*z^2))/(x^2 + y^2 + z^2)^(5/2);
    
    J(11,8) = -1;
    J(12,9) = -1;
    J(13,10) = -1;
    
    J(14,7) = (2*T*(lvx^2 + lvy^2 + lvz^2)^(1/2))/m^3;
    J(14,11) = -(T*lvx)/(m^2*(lvx^2 + lvy^2 + lvz^2)^(1/2));
    J(14,12) = -(T*lvy)/(m^2*(lvx^2 + lvy^2 + lvz^2)^(1/2));
    J(14,13) = -(T*lvz)/(m^2*(lvx^2 + lvy^2 + lvz^2)^(1/2));


end

