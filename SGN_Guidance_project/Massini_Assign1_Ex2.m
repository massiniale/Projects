%% ASSIGNMENT 1 EXERCISE 2 (GUIDANCE; IMPULSIVE GUIDANCE)

% Author: Massini Alessandro - 10990413

clearvars; close all; clc

cspice_kclear()
applyDefaultSettings();

format long

cspice_furnsh('ex02.tm')

% INITIALIZATION OF THE PARAMETERS
GM_E = cspice_bodvrd('EARTH', 'GM',1);  % gravitational constant
GM_M = cspice_bodvrd('MOON', 'GM',1);   % gravitational constant

R_E_vec = cspice_bodvrd('EARTH', 'RADII',3);  % Radius
R_E = R_E_vec(1);
 
R_M_vec = cspice_bodvrd('MOON', 'RADII',3);   % Radius
R_M = R_M_vec(1);

mu = GM_M / (GM_E + GM_M);              % scaled mass constant
m_s  = 3.28900541e05;                   % scaled sun mass
rho  = 3.88811143e02;                   % scaled sun-earth distance
om_s = -9.25195985e-1;                  % scaled sun angular velocity

om_em = 2.66186135e-06;                 % scaled moon amgular velocity
l_em = 3.84405e08;                      % earth-moon distance
 
DU = 3.84405000e05;                     % distance unit
VU = 1.02454018e03;                     % velocity unit
TU = 4.34256461;                        % time unit

hi = 167;                               % initial altitude
hf = 100;                               % final altitude

alpha = 0.2 * pi;
beta = 1.41;
delta = 4;
ti = 2;

r0 = (R_E + hi)/DU;
rf = (R_M + hf)/DU;
v0 = beta * sqrt((1 - mu)/r0);

parameters = parametersDefinition(mu,m_s,rho,om_s,om_em,l_em,...
        r0,rf,R_E,R_M,DU,VU,TU);


%% ----------------------EXERCISE 2.1----------------------------------- %%

% COMPUTATION OF INITIAL CARTESIAN STATE
x0 = r0 * cos(alpha) - mu;
y0 = r0 * sin(alpha);
vx0 = -(v0 - r0) * sin(alpha);
vy0 =  (v0 - r0) * cos(alpha);

initialState = [x0;y0;vx0;vy0];
tf = ti + delta;

% PROPAGATION IN THE ROTATING FRAME
[~,~, state, tt] = fourBodyPropagator(ti, initialState, tf, parameters);

% ROTATION OF THE TRAJECTORY IN THE EARTH-CENTERED INTERTIAL FRAME
stateECI = rotation2ECI(tt, state, parameters);

% PLOT OF THE INITIAL GUESS IN E-M ROTATING FRAME
theta = linspace(0, 2*pi, 1000);

xTargetOrbit = rf*cos(theta) + 1 - mu;
yTargetOrbit = rf*sin(theta);

xMoon = cos(theta);
yMoon = sin(theta);

xInitialOrbit = r0*cos(theta) - mu;
yInitialOrbit = r0*sin(theta);

figure
plot(state(:,1),state(:,2), 'k', 'LineWidth', 2)
hold on
plot(xTargetOrbit, yTargetOrbit, 'r--', 'LineWidth', 2)
plot(xInitialOrbit, yInitialOrbit, 'r--', 'LineWidth', 2)
plot(-mu,0,'b', 'Marker','o', 'MarkerSize',8,'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'b')
plot(1-mu,0, 'Marker','o', 'MarkerSize',6, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5])
grid on
axis equal
xlabel('x[DU]')
ylabel('y[DU]')
legend('Trajectory', '','','Earth', 'Moon','FontSize',18)
title('@EMB Earth-Moon Rotating Frame')

% PLOT OF THE INITIAL GUESS IN ECI FRAME
XXmoon = rotation2ECI(tf,[1-mu,0,0,0],parameters);

figure
plot(stateECI(:,1), stateECI(:,2),'k', 'LineWidth', 2)
hold on
plot(xMoon,yMoon,'k--', 'LineWidth', 2)
plot(xInitialOrbit + mu,yInitialOrbit,'r--', 'LineWidth', 2)
plot(0, 0, 'b', 'Marker','o', 'MarkerSize',8,'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'b')
plot(XXmoon(1), XXmoon(2), 'Marker','o', 'MarkerSize',6, 'MarkerEdgeColor', ...
    'k', 'MarkerFaceColor', [0.5, 0.5, 0.5])
grid on
axis equal
xlabel('x[DU]')
ylabel('y[DU]')
legend('Trajectory', 'Moon Orbit','','Earth', 'Moon','FontSize',18)
ylim([-1.6 1.2])
title('@Earth Earth-Centered Inerital Frame')


%% ---------------------EXERCISE 2.2------------------------------------ %%
%-----------------------SIMPLE SHOOTING-----------------------------------%
%-----------------------WITHOUT GRADIENTS---------------------------------%

% CONSTRUCTION OF THE INITIAL GUESS
initialGuess = [initialState; ti; tf];

% OBJECTIVE FUNCTION
objFun = @(initialGuess) costFunction(initialGuess, parameters);

% NON LINEAR CONSTRAINTS
c =  @(initialGuess) nonLinearConstraints(initialGuess, parameters);

% FMINCON OPTIONS
options = optimoptions('fmincon', 'Algorithm', 'active-set',...
    'Display','iter-detailed','ConstraintTolerance', 1e-10, 'FiniteDifferenceStepSize',1e-12);

% COMPUTATION OF THE OPTIMAL INITIAL STATE 
[optGuess,deltaV] = fmincon(objFun,initialGuess,[],[],[],[],[],[], c, options);


% PROPAGATION OF THE OPTIMAL GUESS
[~,~, optState, ttOpt] = fourBodyPropagator(optGuess(5),optGuess(1:4), optGuess(6), parameters);

% ROTATION IN ECI
optStateEci = rotation2ECI(ttOpt, optState, parameters);
XXMoonSS = rotation2ECI(ttOpt(end),[1-mu,0,0,0],parameters);
xTargetOrbitSS = rf*cos(theta) + XXMoonSS(1);
yTargetOrbitSS = rf*sin(theta) + XXMoonSS(2);

% PLOT IN ROTATING FRAME
figure
plot(optState(:,1),optState(:,2), 'k', 'LineWidth', 2)
hold on
plot(xTargetOrbit, yTargetOrbit, 'r--', 'LineWidth', 2)
plot(xInitialOrbit, yInitialOrbit, 'r--', 'LineWidth', 2)
plot(-mu,0,'b', 'Marker','o', 'MarkerSize',8,'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'b')
plot(1-mu,0, 'Marker','o', 'MarkerSize',6, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5])
grid on
axis equal
xlabel('x[DU]')
ylabel('y[DU]')
legend('Trajectory', '','','Earth', 'Moon','FontSize',18)
xlim([-1.3 1.45])
ylim([-0.1 1.7])
title('@EMB Earth-Moon Rotating Frame')

% PLOT IN ECI FRAME
figure
plot(optStateEci(:,1), optStateEci(:,2),'k', 'LineWidth',2)
hold on
plot(xMoon,yMoon,'k--', 'LineWidth', 2)
plot(xInitialOrbit + mu,yInitialOrbit,'r--', 'LineWidth', 2)
plot(xTargetOrbitSS, yTargetOrbitSS, 'r--', 'LineWidth', 2)
plot(0,0,'b', 'Marker','o', 'MarkerSize',8,'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'b')
plot(XXMoonSS(1), XXMoonSS(2), 'Marker','o', 'MarkerSize',6, 'MarkerEdgeColor', ...
    'k', 'MarkerFaceColor', [0.5, 0.5, 0.5])
grid on
axis equal
xlabel('x[DU]')
ylabel('y[DU]')
legend('Trajectory', 'Moon Orbit','','','Earth', 'Moon','FontSize',18)
xlim([-1.7 1.7])
ylim([-1.2 1.2])
title('@Earth Earth-Centered Inertial Frame')

%--------------------------WITH GRADIENTS---------------------------------%

gradientSelector = 'true';

% OBJECTIVE FUNCTION
objFunGrad = @(initialGuess) costFunction(initialGuess,parameters, gradientSelector);

% NON LINEAR CONSTRAINTS
cGrad =  @(initialGuess) nonLinearConstraints(initialGuess,parameters,gradientSelector);

% OPTION FOR THE 'WITH GRADIENTS' CASE
optionsGrad = optimoptions('fmincon', 'Algorithm', 'active-set',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
    'Display','iter-detailed', 'ConstraintTolerance', 1e-10);

% COMPUTATION OF THE OPTIMAL INITIAL STATE 
[optGuessGrad,deltaVGrad] = fmincon(objFunGrad,initialGuess,[],[],[],[],[], ...
    [], cGrad, optionsGrad);

% PROPAGATION OF THE OPTIMAL GUESS
[~,~, optStateGrad, ttOptGrad] = fourBodyPropagator(optGuessGrad(5), ...
    optGuessGrad(1:4), optGuessGrad(6), parameters);

% ROTATION IN ECI
optStateEciGrad = rotation2ECI(ttOptGrad, optStateGrad, parameters);
XXMoonSSGrad = rotation2ECI(ttOptGrad(end),[1-mu,0,0,0],parameters);
xTargetOrbitSSGrad = rf*cos(theta) + XXMoonSSGrad(1);
yTargetOrbitSSGrad = rf*sin(theta) + XXMoonSSGrad(2);

% PLOT IN ROTATING FRAME
figure
plot(optStateGrad(:,1),optStateGrad(:,2), 'k', 'LineWidth', 2)
hold on
plot(xTargetOrbit, yTargetOrbit, 'r--', 'LineWidth', 2)
plot(xInitialOrbit, yInitialOrbit, 'r--', 'LineWidth', 2)
plot(-mu,0,'b', 'Marker','o', 'MarkerSize',8,'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'b')
plot(1-mu,0, 'Marker','o', 'MarkerSize',6, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5])
grid on
axis equal
xlabel('x[DU]')
ylabel('y[DU]')
legend('Trajectory', '','','Earth', 'Moon','FontSize',18)
xlim([-1.3 1.45])
ylim([-0.1 1.7])
title('@EMB Earth-Moon Rotating Frame')

% PLOT IN ECI FRAME
figure
plot(optStateEciGrad(:,1), optStateEciGrad(:,2),'k', 'LineWidth',2)
hold on
plot(xMoon,yMoon,'k--', 'LineWidth', 2)
plot(xInitialOrbit + mu,yInitialOrbit,'r--', 'LineWidth', 2)
plot(xTargetOrbitSSGrad, yTargetOrbitSSGrad, 'r--', 'LineWidth', 2)
plot(0,0,'b', 'Marker','o', 'MarkerSize',8,'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'b')
plot(XXMoonSSGrad(1), XXMoonSSGrad(2), 'Marker','o', 'MarkerSize',6, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.5, 0.5, 0.5])
grid on
axis equal
xlabel('x[DU]')
ylabel('y[DU]')
legend('Trajectory', 'Moon Orbit','','','Earth', 'Moon','FontSize',18)
xlim([-1.7 1.7])
ylim([-1.2 1.2])
title('@EMB Earth-Centered Inertial Frame')

% COMPUTE ERRORS
[~,ceq] = c(optGuess);

% Position errors in [m]:
errorsPosition = ceq(3) * DU * 1000;

% Velocity errors in [m/s]: 
errorsVelocity = ceq(4) * VU;

errors = [errorsPosition; errorsVelocity];

[c,ceqGrad] = cGrad(optGuessGrad);

% Position errors in [m]:
errorsPositionGrad = ceqGrad(3) * DU * 1000; 

% Velocity errors in [m/s]: 
errorsVelocityGrad = ceqGrad(4) * VU; 

errorsGrad = [errorsPositionGrad; errorsVelocityGrad];


%% ---------------------EXERCISE 2.3------------------------------------ %%
%-----------------------MULTIPLE SHOOTING---------------------------------%

% INITIALIZATION
N = 4;
t = zeros(N,1);
ti = 2;
tf = 6;
t(1) = ti;
initialStateMS = zeros(N,4);
initialStateMS(1,:) = initialState;

% DEFINITION OF NODES AND VARIABLES
for i = 2:N
    t(i) = ti + (tf - ti) * (i - 1)/(N - 1);
    [initialStateMS(i,:),~,~,~]= fourBodyPropagator(t(i-1),initialStateMS(i-1,:), ...
        t(i), parameters);
end

initialStateMS = initialStateMS';
initialGuessMS = [initialStateMS(:);ti;tf];

% OBJECTIVE FUNCTION
objFunMS =  @(initialGuessMS) costFunction_MS(initialGuessMS,parameters, N);

% CONSTRAINTS
constraints = @(initialGuessMS) constraints_MS(initialGuessMS, parameters, N);

% COMPUTATION OF THE OPTIMAL INITIAL STATE 
optionsMS = optimoptions('fmincon', 'Algorithm', 'active-set',...
    'SpecifyObjectiveGradient',true, 'SpecifyConstraintGradient',true,...
    'Display','iter-detailed', 'ConstraintTolerance', 1e-7,'MaxFunctionEvaluations',5000, ...
    'MaxIterations',500);
[optGuessMS,deltaVMS] = fmincon(objFunMS,initialGuessMS,[],[],[],[],[],[], constraints, optionsMS);

% PROPAGATION AND PLOT OF THE SOLUTION IN ROTATING FRAME
tt = zeros(N,1);
tt(1) = optGuessMS(end-1);
stateMS = cell(N, 2);

figure
for i = 1:N-1
    tt(i+1) = optGuessMS(end-1) + (optGuessMS(end) - optGuessMS(end-1)) * (i)/(N - 1);
    [~,~,propagatedMS,ttMS] = fourBodyPropagator(tt(i), optGuessMS(4*i-3:4*i), ...
        tt(i+1), parameters);
    stateMS{i,1} = propagatedMS;
    stateMS{i,2} = ttMS;
    plot(propagatedMS(:,1),propagatedMS(:,2), 'k')
    hold on
    plot(optGuessMS(4*i-3),optGuessMS(4*i-2),  'Marker','o', 'MarkerSize',8, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
end
plot(optGuessMS(end-5),optGuessMS(end-4),'Marker','o', 'MarkerSize',8, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
plot(-mu,0,'b', 'Marker','o', 'MarkerSize',8,'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'b')
plot(1-mu,0, 'Marker','o', 'MarkerSize',6, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5])
grid on
axis equal
xlabel('x[DU]')
ylabel('y[DU]')
legend('Trajectory','Nodes','', '','','','','Earth', 'Moon','FontSize',18)
xlim([-1.3 1.45])
ylim([-0.1 1.7])
title('@EMB Earth-Moon Rotating Frame')

% ROTATION AND PLOT IN ECI
XXMoonMS = rotation2ECI(ttMS(end),[1-mu,0,0,0],parameters);
xTargetOrbitMS = rf*cos(theta) + XXMoonMS(1);
yTargetOrbitMS = rf*sin(theta) + XXMoonMS(2);

figure
hold on
plot(xMoon,yMoon,'k--', 'LineWidth', 2)
plot(xInitialOrbit + mu,yInitialOrbit,'r--', 'LineWidth', 2)
plot(xTargetOrbitMS, yTargetOrbitMS, 'r--', 'LineWidth', 2)
grid on
axis equal
xlabel('x[DU]')
ylabel('y[DU]')


for i = 1:N-1
    optStateEciMS = rotation2ECI(stateMS{i,2}, stateMS{i,1}, parameters);
    plot(optStateEciMS(:,1),optStateEciMS(:,2), 'k')
    plot(optStateEciMS(1,1),optStateEciMS(1,2),'Marker','o', 'MarkerSize',8, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
end
plot(optStateEciMS(end,1),optStateEciMS(end,2),'Marker','o', 'MarkerSize',8, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
plot(0,0,'b', 'Marker','o', 'MarkerSize',8,'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'b')
plot(XXMoonMS(1), XXMoonMS(2), 'Marker','o', 'MarkerSize',6, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.5, 0.5, 0.5])
grid on
axis equal
legend('Moon Orbit','','','Trajectory','','','Nodes','','','', 'Earth', 'Moon','FontSize',18)
xlim([-1.7 1.72])
ylim([-1.2 1.22])
title('@ECI Earth-Centered Inertial Frame')

% COMPUTE ERRORS
[~, ceqMS] = constraints(optGuessMS);

% Position errors in [m]:
errorsPositionEarthMS = ceqMS(13) * DU * 1000;  
% Velocity errors in [m/s]: 
errorsVelocityEarthMS = ceqMS(14) * VU; 

% Position errors in [m]:
errorsPositionMoonMS = ceqMS(15) * DU * 1000;  
% Velocity errors in [m/s]: 
errorsVelocityMoonMS = ceqMS(16) * VU; 

% Position errors in [m]:
errorsXnode2MS = ceqMS(1) * DU * 1000;
errorsYnode2MS = ceqMS(2) * DU * 1000;

% Velocity errors in [m/s]: 
errorsVXnode2MS = ceqMS(3) * VU;
errorsVYnode2MS = ceqMS(4) * VU;

% Position errors in [m]:
errorsXnode3MS = ceqMS(5) * DU * 1000;
errorsYnode3MS = ceqMS(6) * DU * 1000;

% Velocity errors in [m/s]: 
errorsVXnode3MS = ceqMS(7) * VU;
errorsVYnode3MS = ceqMS(8) * VU;

% Position errors in [m]:
errorsXnode4MS = ceqMS(9) * DU * 1000;
errorsYnode4MS = ceqMS(10) * DU * 1000;

% Velocity errors in [m/s]: 
errorsVXnode4MS = ceqMS(11) * VU;
errorsVYnode4MS = ceqMS(12) * VU;

errorsMS = [sqrt(errorsXnode2MS^2 + errorsYnode2MS^2);
            sqrt(errorsVXnode2MS^2 + errorsVYnode2MS^2);
            sqrt(errorsXnode3MS^2 + errorsYnode3MS^2);
            sqrt(errorsVXnode3MS^2 + errorsVYnode3MS^2);
            sqrt(errorsXnode4MS^2 + errorsYnode4MS^2);
            sqrt(errorsVXnode4MS^2 + errorsVYnode4MS^2);
            errorsPositionEarthMS;
            errorsVelocityEarthMS
            errorsPositionMoonMS;  
            errorsVelocityMoonMS; ];




%% ---------------------EXERCISE 2.4------------------------------------ %%

% DATA INITIALIZATION
tformat = 'YYYY-MON-DD-HR:MN:SC.####::TDB';
t_i_tdb = 'September 28 00:00:00.000 TDB 2024';

labels = {'SUN';
          'EARTH';
          'MOON';
          'MERCURY';
          'VENUS';
          'MARS BARYCENTER';
          'JUPITER BARYCENTER';
          'SATURN BARYCENTER';
          'URANUS BARYCENTER';
          'NEPTUNE BARYCENTER';
          'PLUTO BARYCENTER'};

bodies = nbody_init(labels);
frame.theta = 'J2000';
frame.nbody = 'J2000';
center = 'EARTH';

moonRevPeriod = 2 * pi / om_em;

et_l = cspice_str2et(t_i_tdb);
et_u = et_l + moonRevPeriod;

% DEFINITION OF THE TARGET ANGLE NEEDED TO IDENTIFY THE INITIAL EPOCH
theta_target = om_s * optGuessGrad(5); 
if theta_target <= 0
    theta_target = theta_target + 2 * pi;
end

et = linspace(et_l, et_u, 100);

% FIND THE INITIAL AND FINAL EPOCH 
fun =  @(et) thetaFinder(theta_target, et, bodies,frame, center);
et_guess = (et_l + et_u) / 2;

et_i = fzero(fun,[et_guess, et_u]);
et_f = et_i + (optGuessGrad(6)-optGuessGrad(5)) * TU * 24 * 3600;

initialEpoch = cspice_et2utc(et_i,'C', 3);
disp(['initial Epoch: ', initialEpoch])
finalEpoch = cspice_et2utc(et_f,'C', 3);
disp(['initial Epoch: ', finalEpoch])

% PROPAGATION WITH NBODY PROPAGATOR OF THE INITIAL STATE IN J200 EC FRAME
initialStateECIscaled = [optStateEciGrad(1,1);optStateEciGrad(1,2);0;...
    optStateEciGrad(1,3);optStateEciGrad(1,4);0];
initialStateECI = [initialStateECIscaled(1:3) * DU; initialStateECIscaled(4:6) * VU*10^-3];


options = odeset('reltol', 1e-12, 'abstol', [1e-6*ones(1,3),1e-9*ones(1,3)]);
[tt, statenBody] = ode113(@(t,x) nbodyShiftRHS(t,x,bodies,frame,center), ...
    [et_i et_f], initialStateECI, options);

etProp = linspace(et_i,et_f,1000);

% PLOT OF THE TWO TRAJECTORIES
figure
hold on
plot(statenBody(:,1),statenBody(:,2),'LineWidth', 2)
plot(optStateEciGrad(:,1) * DU, optStateEciGrad(:,2) * DU, 'k')
plot(xMoon * DU, yMoon * DU,'k--', 'LineWidth',1)
plot(0,0,'b', 'Marker','o', 'MarkerSize',6,'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b')
axis equal
grid on
set(gca, 'FontSize',18)
xlim([-5e05, 7e05 ])
ylim([-5e05, 5e05 ])
xticks([-5e5 -2.5e5 0 2.5e5 5e5]); % Specifica i tick sugli assi X
yticks([-5e5 -2.5e5 0 2.5e5 5e5]); % Specifica i tick sugli assi Y
legend('NBody','PBRFBP','Earth')
xlabel('x[km]')
ylabel('y[km]')
title('@ECI Earth-Centered Inertial Frame')


%%
% CLEAR KERNEL POOL
cspice_kclear()

%% ---------------------FUNCTIONS----------------------------------------%%
 
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

% PARAMETERS DEFINITION 
function parameters = parametersDefinition (mu,m_s,rho,om_s,om_em,l_em,...
        r0,rf,R_E,R_M,DU,VU,TU)

% parametersDefinition.m - it assigns the parameter to the respective field
%                          of the structure 'parameters' in output, see the
%                          main text for the meaning of each parameter

 parameters.mu = mu;
 parameters.m_s = m_s;
 parameters.rho = rho;
 parameters.om_s = om_s;
 parameters.om_em = om_em;
 parameters.l_em = l_em;
 parameters.VU = VU;
 parameters.DU = DU;
 parameters.TU = TU;
 parameters.r0 = r0;
 parameters.rf = rf;
 parameters.R_E = R_E/DU;
 parameters.R_M = R_M/DU;

end

% PROPAGATOR FOR THE PBRFBP, WITHOUT THE COMPUTATION OF STM
function [xf,tf, xx, tt]  = fourBodyPropagator(t0,x0,tf,parameters)

% fourBodyPropagator.m - Perform propagation of the state integrating the
%                        equations of motion in the Planar Bicircular
%                        restricted four body problem
%
%
% Prototype:       [xf,tf,xx,tt] = fourBodyPropagator(t0,x0,tf,parameters)
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
% x0[4x1]:  Initial state of the system (x,y,vx,vy) 
%
% tf[1]:    Final time of propagation 
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                         mu[1]:    Mass Ratio 
%                         m_s[1]:   Scaled Mass of the Sun 
%                         rho[1]:   Scaled Sun-(Earth+Moon) distance
%                         om_s[1]:  Scaled Angular velocity of the Sun
%            
%
% Outputs:
%
%
% xf[4x1]:   Final state of the system (x,y,vx,vy)
% 
% tf[1]:     Final time of propagation 
%
% xx[nx4]:   Full propagation of the initial state (x,y,vx,vy)
%
% tt[nx1]:   Vector of time
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025
                          

% Perform integration
options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[tt, xx] = ode78(@(t,x) PBRFBP(t,x,parameters), [t0 tf], x0, options);


% Extract final state vector and time
xf = xx(end,1:4)';
tf = tt(end);

end

% EQUATION OF MOTIONS FOR THE PBRFBP, WITHOUT STM COMPUTATION
function [dxdt] = PBRFBP(t,xx, parameters)

% PBRFBP.m - Returns the derivative of the state in the Planar Bicircular
%            Restricted Four Body Problem, in the rotating frame components
%
%
% Prototype:       [dxdt] = PBRFBP(t,xx, parameters)
%
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs: 

% t[1]:     Current time of integration
% xx[4x1]:  Current state of the system (x,y,vx,vy) 
% parameters[1x1 struct]: Parameters necessary for the computation
%                         
%                         mu[1]:    Mass Ratio 
%                         m_s[1]:   Scaled Mass of the Sun 
%                         rho[1]:   Scaled Sun-(Earth+Moon) distance
%                         om_s[1]:  Scaled Angular velocity of the Sun
%            
%
% Outputs:
%
%
% dxdt[4x1]:  Derivative state of the system
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025
  

% Initialization 
x  = xx(1);
y  = xx(2);
vx = xx(3);
vy = xx(4);   

mu = parameters.mu;
m_s = parameters.m_s;
rho = parameters.rho;
om_s = parameters.om_s;    

% Compute derivative of the potential
dUdx = x - (m_s*cos(om_s*t))/rho^2 - (mu*(mu + x - 1))/((mu + x - 1)^2 ...
       + y^2)^(3/2) + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2))...
       - (m_s*(2*x - 2*rho*cos(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y ...
       - rho*sin(om_s*t))^2)^(3/2));

dUdy = y - (m_s*sin(om_s*t))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) ...
       - (m_s*(2*y - 2*rho*sin(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - ...
       rho*sin(om_s*t))^2)^(3/2)) + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2);

% Assemble the derivative of the state
dxdt = zeros(4,1);
dxdt(1:2) = xx(3:4);
dxdt(3)   = dUdx + 2*vy;
dxdt(4)   = dUdy - 2*vx;

end

% EQUATION OF MOTIONS FOR THE PBRFBP, WITH STM COMPUTATION
function [dxdt] = PBRFBP_STM(t,xx, parameters)

% PBRFBP.m - Returns the derivative of the state in the Planar Bicircular
%            Restricted Four Body Problem, in the rotating frame components
%
%
% Prototype:       [dxdt] = PBRFBP(t,xx, parameters)
%
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs: 

% t[1]:     Current time of integration
% xx[20x1]:  Current state of the system (x,y,vx,vy) 
% parameters[1x1 struct]: Parameters necessary for the computation
%                         
%                         mu[1]:    Mass Ratio 
%                         m_s[1]:   Scaled Mass of the Sun 
%                         rho[1]:   Scaled Sun-(Earth+Moon) distance
%                         om_s[1]:  Scaled Angular velocity of the Sun
%            
%
% Outputs:
%
%
% dxdt[20x1]:  Derivative state of the system
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025
    
    % Extract variables
    x  = xx(1);
    y  = xx(2);
    vx = xx(3);
    vy = xx(4);

    mu = parameters.mu;
    m_s = parameters.m_s;
    rho = parameters.rho;
    om_s = parameters.om_s;

    % Put PHI in matrix form
    Phi = reshape(xx(5:end),4,4);

    % Compute derivative of the potential
    dUdx = x - (m_s*cos(om_s*t))/rho^2 - (mu*(mu + x - 1))/((mu + x - 1)^2 ...
        + y^2)^(3/2) + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2))...
        - (m_s*(2*x - 2*rho*cos(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y ...
        - rho*sin(om_s*t))^2)^(3/2));

    dUdy = y - (m_s*sin(om_s*t))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) ...
        - (m_s*(2*y - 2*rho*sin(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - ...
        rho*sin(om_s*t))^2)^(3/2)) + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2);
 
    % Assemble the matrix A(t)=dfdx 4x4 matrix
    dfdx = [0, 0, 1, 0;
            0, 0, 0, 1;
            (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - m_s/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2) - (3*(mu + x)^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*m_s*(x - rho*cos(om_s*t))^2)/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2) + (3*mu*(mu + x - 1)^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1, (3*m_s*(2*x - 2*rho*cos(om_s*t))*(2*y - 2*rho*sin(om_s*t)))/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)), 0, 2;
            (3*m_s*(2*x - 2*rho*cos(om_s*t))*(2*y - 2*rho*sin(om_s*t)))/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)), (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - m_s/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*m_s*(y - rho*sin(om_s*t))^2)/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1, -2, 0];
 
    % Compute the derivative of the STM
    Phidot = dfdx*Phi;

    % Assemble right-hand side
    dxdt = zeros(20,1);

    dxdt(1:2) = xx(3:4);
    dxdt(3)   = dUdx + 2*vy;
    dxdt(4)   = dUdy - 2*vx;
    dxdt(5:end) = Phidot(:);
    
end

% PROPAGATOR FOR THE PBRFBP, WITH THE COMPUTATION OF STM
function [xf,tf, PHIf, xx, tt]  = fourBodyPropagator_STM(t0,x0,tf,parameters)

% fourBodyPropagator_STM.m - Perform propagation of the state integrating the
%                        equations of motion in the Planar Bicircular
%                        restricted four body problem
%
%
% Prototype: [xf,tf,PHIf, xx,tt] = fourBodyPropagator_STM(t0,x0,tf,parameters)
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
% x0[4x1]:  Initial state of the system (x,y,vx,vy) 
%
% tf[1]:    Final time of propagation 
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                         mu[1]:    Mass Ratio 
%                         m_s[1]:   Scaled Mass of the Sun 
%                         rho[1]:   Scaled Sun-(Earth+Moon) distance
%                         om_s[1]:  Scaled Angular velocity of the Sun
%            
%
% Outputs:
%
%
% xf[4x1]:   Final state of the system (x,y,vx,vy)
% 
% tf[1]:     Final time of propagation
%
% PHIf[4x4]: State Transition Matrix of the final state
%
% xx[nx20]:   Full propagation of the initial state and of the elements of the STM (x,y,vx,vy)
%
% tt[nx1]:   Vector of time
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

    % Initialize State Transition Matrix at t0
    Phi0 = eye(4);

    % Append to initial conditions the conditions for the STM
    x0Phi0 = [x0; Phi0(:)];
    
    % Perform integration
    options_STM = odeset('reltol', 1e-11, 'abstol', 1e-11);
    [tt, xx] = ode78(@(t,x) PBRFBP_STM(t,x,parameters), [t0 tf], x0Phi0, options_STM);

    % Extract state vector and State Transition Matrix
    xf = xx(end,1:4)';
    PHIf = reshape(xx(end,5:end),4,4);
    tf = tt(end);

end

% ROTATION OF THE TRAJECTORY IN ECI FRAME
function XX = rotation2ECI(t, xx, parameters)

% rotation2ECI.m - Perform the rotation of the trajectory from the rotating
%                  Earth-Moon frame to the Earth-centered inertial frame
%
%
% Prototype:       XX = rotation2ECI(t, xx, parameters)
%
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs:     
%
%
% t[nx1]:    time 
%
% xx[nx4]:  Trajectory in Earth-Moon rotating frame (x,y,vx,vy) 
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                         mu[1]:     Mass Ratio 
%            
%
% Outputs:
%
%
% XX[N,4]:   Trajectory in Earth-centered intertial frame (X,Y,VX,VY)
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025


% Initialization
x   = xx(:,1);
y   = xx(:,2);
vx  = xx(:,3);
vy  = xx(:,4);

mu = parameters.mu;

% Rotation in ECI
XX(:,1) = (x + mu) .* cos(t) - y .* sin(t);
XX(:,2) = (x + mu) .* sin(t) + y .* cos(t);
XX(:,3) = (vx - y) .* cos(t) - (vy + x + mu) .* sin(t);
XX(:,4) = (vx - y) .* sin(t) + (vy + x + mu) .* cos(t);

end

% COST FUNCTION FOR SIMPLE SHOOTING
function [totalDeltaV, grad] = costFunction(variables, parameters, selector)

% costFunction.m - Perform the computation of the cost function (DeltaV)
%                  for simple shooting optimization
%
%
% Prototype:     [totalDeltaV, grad] = costFunction(variables, parameters, selector)
%
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs:     
%
%
% variables[6x1]:  variables of the problem, respectively initialState and 
%                  initial and Final time 
%
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                         mu[1]:     Mass Ratio 
%                         r0[1]:     norm of position of departure orbit
%                         rf[1]:     norm of position of arrival orbit
%                         m_s[1]:    Scaled Mass of the Sun 
%                         rho[1]:    Scaled Sun-(Earth+Moon) distance
%                         om_s[1]:   Scaled Angular velocity of the Sun
%
% selector[string]: True if Gradients are needed, False or ~ if not
%            
%
% Outputs:
%
%
% totalDeltaV[1]:   DeltaV of the Orbit
% 
% grad[6x1]:        Gradient of the Cost Function with respect to the variables
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025


% Switch for cases with or without gradients
if nargin < 3
    selector = 'false';
end

% Initialization
xi  = variables(1);
yi  = variables(2);
vxi = variables(3);
vyi = variables(4);
ti  = variables(5);
tf  = variables(6);

mu = parameters.mu;
ri = parameters.r0;
rf = parameters.rf;
 
% Propagation and Gradients Computation
if strcmpi(selector,'false')
    [finalState,~,~, ~]  = fourBodyPropagator(ti,[xi;yi;vxi;vyi],tf,parameters);
else 
    if strcmpi(selector,'true')
        [finalState,~,PHI,~, ~]  = fourBodyPropagator_STM(ti,[xi;yi;vxi;vyi],tf,parameters);
        grad = gradientCostFunction(variables,finalState,PHI,parameters);
    end
end

xf =  finalState(1);
yf =  finalState(2);
vxf = finalState(3);
vyf = finalState(4);

% Assembly of the Cost Function
deltaV_i = sqrt((vxi - yi)^2 + (vyi + xi + mu)^2) - sqrt((1 - mu) / ri);
deltaV_f = sqrt((vxf - yf)^2 + (vyf + xf + mu - 1)^2) - sqrt(mu / rf);

totalDeltaV = deltaV_f + deltaV_i ;

end

% NON LINEAR CONSTRAINTS FOR SIMPLE SHOOTING
function [c,ceq, dc, dceq] = nonLinearConstraints(variables, parameters,selector)

% nonLinearConstraints.m - Perform the computation of the non linear
%                          equality and inequality constraints for simple 
%                          shooting optimization
%
%
% Prototype: [c,ceq, dc, dceq] = nonLinearConstraints(variables, parameters,selector)
%
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs:     
%
%
% variables[6x1]:  variables of the problem, respectively initialState and 
%                  initial and Final time 
%
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                         mu[1]:     Mass Ratio 
%                         r0[1]:     norm of position of departure orbit
%                         rf[1]:     norm of position of arrival orbit
%                         m_s[1]:    Scaled Mass of the Sun 
%                         rho[1]:    Scaled Sun-(Earth+Moon) distance
%                         om_s[1]:   Scaled Angular velocity of the Sun
%
% selector[string]: True if Gradients are needed, False or ~ if not
%            
%
% Outputs:
%
%
% c[1]:         Equality Constraint
%
% ceq[4x1]:     Inequality Constraints
% 
% dc[6x1]:      Gradient of the Eq Constraint with respect to the variables
%
% dceq[6x4]:    Jacobian of the Ineq Constraint with respect to the variables
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025


% Switch for cases without gradients
if nargin < 3
    selector = 'false';
end

% Initialization 
xi  = variables(1);
yi  = variables(2);
vxi = variables(3);
vyi = variables(4);
ti  = variables(5);
tf  = variables(6);

mu = parameters.mu;
ri = parameters.r0;
rf = parameters.rf;

% Propagation and Jacobian Computation
if strcmpi(selector,'false')
    [finalState,tf, ~, ~]  = fourBodyPropagator(ti,[xi;yi;vxi;vyi],tf,parameters);
else
    if strcmpi(selector,'true')
        [finalState,tf,PHI,~, ~]  = fourBodyPropagator_STM(ti,[xi;yi;vxi;vyi],tf,parameters);
        [dc,dceq] = jacobianConstraints(variables,finalState,PHI,parameters);
    end
end

xf =  finalState(1);
yf =  finalState(2);
vxf = finalState(3);
vyf = finalState(4);

% Assembly of the inequality and equality constraints
c = ti - tf;

ceq1 = (xi + mu)^2 + yi^2 - ri^2;
ceq2 = (xi + mu) * (vxi - yi) + yi * (vyi + xi + mu);
ceq3 = (xf + mu - 1)^2 + yf^2 - rf^2;
ceq4 = (xf + mu - 1) * (vxf - yf) + yf * (vyf + xf + mu - 1);

ceq = [ceq1; ceq2; ceq3; ceq4];

end

% GRADIENT OF THE OBJECTIVE FUNCTION FOR SIMPLE SHOOTING
function grad = gradientCostFunction (variables, propagatedVariables, PHI, parameters)

% gradientCostFunction.m - Compute the gradients of the cost function with
%                          respect to the variables of the problem for the 
%                          simple shooting optimization
%
%
% Prototype:  grad = gradientCostFunction (variables, propagatedVariables, PHI, parameters)
%
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs:     
%
%
% variables[6x1]:  Variables of the problem, respectively initialState and 
%                  initial and Final time 
%
% propagatedVariables[4x1]: Fimal State of the trajectory
%
% PHI[4x4]:        State Transition Matrix of the final state with respect
%                  to the initial one
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                         mu[1]:     Mass Ratio 
%                         r0[1]:     norm of position of departure orbit
%                         rf[1]:     norm of position of arrival orbit
%                         m_s[1]:    Scaled Mass of the Sun 
%                         rho[1]:    Scaled Sun-(Earth+Moon) distance
%                         om_s[1]:   Scaled Angular velocity of the Sun           
%
%
% Outputs:
%
%
% grad[6x1]:   Gradients of the objective Function with respect to the
%              variables 
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025


    % Initialization
    xi = variables(1);
    yi = variables(2);
    vxi = variables(3);
    vyi = variables(4);
    ti = variables(5);
    tf = variables(6);

    xf = propagatedVariables(1);
    yf = propagatedVariables(2);
    vxf = propagatedVariables(3);
    vyf = propagatedVariables(4);

    mu = parameters.mu;

    % Computation of derivative of the cost function with respect to the initial state
    dDVidxi = 1/sqrt((vxi-yi)^2+(vyi+xi+mu)^2) * [vyi+xi+mu; yi-vxi; vxi-yi; vyi+xi+mu];
    dDVfdxi = PHI' * 1/sqrt((vxf-yf)^2+(vyf+xf+mu-1)^2) * [vyf+xf+mu-1; yf-vxf; vxf-yf; vyf+xf+mu-1];
    dJdxi = dDVidxi + dDVfdxi;

    dDvfdxf = (1/sqrt((vxf-yf)^2+(vyf+xf+mu-1)^2)*[vyf+xf+mu-1; yf-vxf; vxf-yf; vyf+xf+mu-1]);
    
    % Computation of derivative of cost function with respect to the initial time
    dfdti = PBRFBP(ti,[xi;yi;vxi;vyi],parameters);
    dJdti = - dDvfdxf' * PHI * dfdti;

    % Computation of derivative of cost function with respect to the final time
    dfdtf = PBRFBP(tf,[xf;yf;vxf;vyf],parameters);
    dJdtf = dDvfdxf' * dfdtf;
    
    % Assembly of Gradient
    grad = [dJdxi;dJdti;dJdtf];
end

% JACOBIAN OF THE CONSTRAINTS FOR SIMPLE SHOOTING
function [dC,dCeq] = jacobianConstraints(variables, propagatedVariables, PHI, parameters)

% jacobianConstraints.m -  Perform the computation of the jacobian of the 
%                          non linear equality and inequality constraints 
%                          for simple shooting optimization
%
%
% Prototype: [dC,dCeq] = jacobianConstraints(variables, propagatedVariables, PHI, parameters)
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs:     
%
%
% variables[6x1]:  Variables of the problem, respectively initialState and 
%                  initial and Final time 
%
% propagatedVariables[4x1]: Fimal State of the trajectory
%
% PHI[4x4]:        State Transition Matrix of the final state with respect
%                  to the initial one
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                         mu[1]:     Mass Ratio 
%                         r0[1]:     norm of position of departure orbit
%                         rf[1]:     norm of position of arrival orbit
%                         m_s[1]:    Scaled Mass of the Sun 
%                         rho[1]:    Scaled Sun-(Earth+Moon) distance
%                         om_s[1]:   Scaled Angular velocity of the Sun           
%
%
% Outputs:
%
%
% c[1]:         Equality Constraint
%
% ceq[4x1]:     Inequality Constraints
% 
% dc[6x1]:      Gradient of the Eq Constraint with respect to the variables
%
% dceq[6x4]:    Jacobian of the Ineq Constraint with respect to the variables
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

    % Initialization
    xi = variables(1);
    yi = variables(2);
    vxi = variables(3);
    vyi = variables(4);
    ti = variables(5);
    tf = variables(6);

    xf = propagatedVariables(1);
    yf = propagatedVariables(2);
    vxf = propagatedVariables(3);
    vyf = propagatedVariables(4);

    mu = parameters.mu;
    
    % Computation of derivatives of the Ineq. constraints on the departure orbit
    % with respect to state and times
    dC12dxi = [2 * (xi + mu), 2 * yi, 0, 0;
              vxi, vyi, xi + mu, yi];
    dC12dti = [0;0];
    dC12dtf = [0;0];

    dC34dxf = [2 * (xf + mu - 1), 2 * yf, 0, 0;
              vxf, vyf, xf + mu - 1, yf];
    dC34dxi =dC34dxf * PHI; 

    % Computation of derivatives of the Ineq. constraints on the arrival orbit
    % with respect to state and times
    dfdti = PBRFBP(ti,[xi;yi;vxi;vyi],parameters);
    dC34dti = - dC34dxf * PHI * dfdti;

    dfdtf = PBRFBP(tf,[xf;yf;vxf;vyf],parameters);
    dC34dtf = dC34dxf * dfdtf;

    % Computation of derivatives of the Eq. constraint with respect to variables
    dC = [0;0;0;0;1;-1];

    % Assembly of the Jacobian
    dCeq = [dC12dxi, dC12dti, dC12dtf;
            dC34dxi, dC34dti, dC34dtf];
    dCeq = dCeq';

end

% COST FUNCTION FOR MULTIPLE SHOOTING
function [totalDeltaV,grad] = costFunction_MS(variablesVector, parameters, N)

% costFunction_MS.m - Perform the computation of the cost function (DeltaV)
%                     for multiple shooting optimization
%
%
% Prototype:   [totalDeltaV, grad] = costFunction_MS(variables, parameters, selector)
%
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs:     
%
%
% variables[4N+2, 1]:  variables of the problem, respectively state of the 
%                      system at the nodes and initial and Final time 
%
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                         mu[1]:     Mass Ratio 
%                         r0[1]:     norm of position of departure orbit
%                         rf[1]:     norm of position of arrival orbit
%                         m_s[1]:    Scaled Mass of the Sun 
%                         rho[1]:    Scaled Sun-(Earth+Moon) distance
%                         om_s[1]:   Scaled Angular velocity of the Sun
%
% N[1]:                   Number of nodes
%            
%
% Outputs:
%
%
% totalDeltaV[1]:    DeltaV of the Orbit
% 
% grad[4N+2, 1]:     Gradient of the Cost Function with respect to the variables
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% Initialization
xi  = variablesVector(1);
yi  = variablesVector(2);
vxi = variablesVector(3);
vyi = variablesVector(4);

xf  = variablesVector(end-5);
yf  = variablesVector(end-4);
vxf = variablesVector(end-3);
vyf = variablesVector(end-2);

mu = parameters.mu;
ri = parameters.r0;
rf = parameters.rf;
 
% Computation of the cost function
deltaV_i = sqrt((vxi - yi)^2 + (vyi + xi + mu)^2) - sqrt((1 - mu) / ri);
deltaV_f = sqrt((vxf - yf)^2 + (vyf + xf + mu - 1)^2) - sqrt(mu / rf);

% Assembly of Cost Function and Gradient
totalDeltaV = deltaV_f + deltaV_i ;
grad = gradientCostFunction_MS(variablesVector,parameters, N);

end

% GRADIENT OF THE OBJECTIVE FUNCTION FOR MULTIPLE SHOOTING
function grad = gradientCostFunction_MS(variables, parameters, N)

% gradientCostFunction_MS.m - Compute the gradients of the cost function with
%                             respect to the variables of the problem for the 
%                             multiple shooting optimization
%
%
% Prototype:  grad = gradientCostFunction_MS(variables, parameters, N)
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs:     
%
%
% variables[4N+2, 1]:  Variables of the problem, respectively state of the
%                      system at the nodes and initial and final time
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                         mu[1]:     Mass Ratio 
%                         r0[1]:     norm of position of departure orbit
%                         rf[1]:     norm of position of arrival orbit
%                         m_s[1]:    Scaled Mass of the Sun 
%                         rho[1]:    Scaled Sun-(Earth+Moon) distance
%                         om_s[1]:   Scaled Angular velocity of the Sun           
%
% N[1]:                   Number of Nodes
%
%
% Outputs:
%
%
% grad[4N+2, 1]:   Gradients of the objective Function with respect to the
%                  variables 
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

    % Initialization
    xi = variables(1);
    yi = variables(2);
    vxi = variables(3);
    vyi = variables(4);
    xN  = variables(end-5);
    yN  = variables(end-4);
    vxN = variables(end-3);
    vyN = variables(end-2);

    mu = parameters.mu;

    % Computation of derivative wrt the variables
    dDVdxi = 1/sqrt((vxi-yi)^2+(vyi+xi+mu)^2) * [vyi+xi+mu; yi-vxi; vxi-yi; vyi+xi+mu];
    dDVdxN = 1/sqrt((vxN-yN)^2+(vyN+xN+mu-1)^2)*[vyN+xN+mu-1; yN-vxN; vxN-yN; vyN+xN+mu-1];
    
    % Assembly of the gradient
    grad = zeros(4*N+2,1);
    grad(1:4) = dDVdxi;
    grad(end-5:end-2) = dDVdxN;
end

% NON LINEAR CONSTRAINTS FOR MULTIPLE SHOOTING
function [c, ceq, dc, dceq] = constraints_MS(variables, parameters, N)

% constraints_MS.m - Compute the Inequality and the Equality, Non Linear 
%                    Constraints of the problem for the multiple shooting 
%                    optimization and their analytical derivatives
%
%
% Prototype:  [c, ceq, dc, dceq] = constraints_MS(variables, parameters, N)
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs:     
%
%
% variables[4N+2, 1]:  Variables of the problem, respectively state of the
%                      system at the nodes and initial and final time
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                         mu[1]:     Mass Ratio 
%                         r0[1]:     norm of position of departure orbit
%                         rf[1]:     norm of position of arrival orbit
%                         m_s[1]:    Scaled Mass of the Sun 
%                         rho[1]:    Scaled Sun-(Earth+Moon) distance
%                         om_s[1]:   Scaled Angular velocity of the Sun           
%
% N[1]:                   Number of Nodes
%
%
% Outputs:
%
%
% c[2N+1,1]:         Inequality Constraints
% 
% dc[4N+2, 2N+1]:    Jacobian of the Inequality Constraints
%
% ceq[4N,1]:         Equality Constraints
%
% dceq[4N+2, 4N]:    Jacobian of the Equality Constraints
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% Initialization
t1 = variables(end-1);
tN = variables(end);

xi  = variables(1);
yi  = variables(2);
vxi = variables(3);
vyi = variables(4);

xf  = variables(end-5);
yf  = variables(end-4);
vxf = variables(end-3);
vyf = variables(end-2);

mu = parameters.mu;
ri = parameters.r0;
rf = parameters.rf;
Re = parameters.R_E;
Rm = parameters.R_M;

% Reshape variables for integration
variablesRSP = reshape(variables(1:end-2),4,N).';
t = zeros(N,1);

for i = 1:N

    t(i) = t1 + (tN - t1) * (i - 1)/(N - 1);

end


ceq = zeros(4*N,1);
c = zeros(2*N + 1,1);

% Subdivision of dCeq in submatrices, for simplicity reasons
topLeft = zeros(3*N,4*N);
topRight = zeros(3*N,2);
bottomLeft = zeros(4,4*N);
bottomRight = zeros(4,2);
dc = zeros(2*N+1,4*N+2);


% Definition of the first two Inequality constraints
c(1:2) = [ Re^2 - (variables(1) + mu)^2 - variables(2)^2;
           Rm^2 - (variables(1) + mu - 1)^2 - variables(2)^2];


for i = 1:N-1
    
    % Propagation
    [propagate,~,PHI,~,~] = fourBodyPropagator_STM(t(i),variables(4*i-3:4*i),t(i+1),parameters);
    
    % Computation of the Defects
    ceq(4*i-3:4*i) = propagate - variables(4*i+1:4*i+4);

    % Inequality Constraints
    c(2*i+1:2*i+2) = [ Re^2 - (variables(4*i + 1) + mu)^2 - variables(4*i + 2)^2;
                      Rm^2 - (variables(4*i + 1) + mu - 1)^2 - variables(4*i + 2)^2];
    
    % Derivatives of the state
    dfidti = PBRFBP(t(i), variablesRSP(i,:), parameters);
    df2dt2 = PBRFBP(t(i+1),propagate, parameters);

    % Computation of the elements of the submatrices
    Q1i = -(N - i)/(N - 1) * PHI * dfidti + (N - i - 1)/(N - 1) * df2dt2;
    QNi = -(i - 1)/(N - 1) * PHI * dfidti + i/(N - 1) * df2dt2;

    % Computation of the submatrices of dceq
    topLeft(4*i-3:4*i,4*i-3:4*i) = PHI;
    topLeft(4*i-3:4*i,4*i+1:4*i+4) = -eye(4);

    topRight(4*i-3:4*i,1) = Q1i;
    topRight(4*i-3:4*i,2) = QNi;


end

% Computation of the submatrices
bottomLeft(1:2,1:4) = [2 * (xi + mu), 2 * yi, 0, 0;
             vxi, vyi, xi + mu, yi];
bottomLeft(3:4,end-3:end) = [2 * (xf + mu - 1), 2 * yf, 0, 0;
             vxf, vyf, xf + mu - 1, yf];

% Assembly of the Jacobian of the equality constraints
dceq = [topLeft, topRight;
        bottomLeft, bottomRight];
dceq = dceq';

for j = 1:N

    xj = variablesRSP(j,1);
    yj = variablesRSP(j,2);

    % Computation of the elements of the jacobian of the Inequality
    % constraints
    Sj = [-2*(xj + mu) -2*yj 0 0;
          -2*(xj + mu -1) -2*yj 0 0];
    dc(2*j-1:2*j,4*j-3:4*j) = Sj; 

end

dc(end,end-1:end) = [1 -1];
dc = dc';

% Assembly
c(end) = t1 - tN;

ceq(end - 3) = (variables(1) + mu)^2 + variables(2)^2 - ri^2;

ceq(end - 2) = (variables(1) + mu) * (variables(3) - variables(2)) ...
    + variables(2) * (variables(4) + variables(1) + mu);

ceq(end - 1) = (variables(3*N + 1) + mu - 1)^2 + variables(3*N + 2)^2 - rf^2;

ceq(end) = (variables(3*N + 1) + mu - 1) * (variables(3*N + 3) - variables(3*N + 2)) ...
    + variables(3*N + 2) * (variables(3*N + 4) + variables(3*N + 1) + mu - 1);

end

% INITIALIZATION OF THE NBODY PROPAGATOR
function [bodies] = nbody_init(labels)

%NBODY_INIT Initialize planetary data for n-body propagation
%   Given a set of labels of planets and/or barycentres, returns a
%   cell array populated with structures containing the body label and the
%   associated gravitational constant.
%
%
% Author
%   Name: ALESSANDRO 
%   Surname: MORSELLI
%   Research group: DART
%   Department: DAER
%   University: Politecnico di Milano 
%   Creation: 26/09/2021
%   Contact: alessandro.morselli@polimi.it
%   Copyright: (c) 2021 A. Morselli, Politecnico di Milano. 
%                  All rights reserved.
%
%
% Notes:
%   This material was prepared to support the course 'Satellite Guidance
%   and Navigation', AY 2021/2022.
%
%
% Inputs:
%   labels : [1,n] cell-array with object labels
%
% Outputs:
%   bodies : [1,n] cell-array with struct elements containing the following
%                  fields
%                  |
%                  |--bodies{i}.name -> body label
%                  |--bodies{i}.GM   -> gravitational constant [km**3/s**2]
%
%
% Prerequisites:
%   - MICE (Matlab SPICE)
%   - Populated kernel pool (PCK kernels)
%

% Initialize output
bodies = cell(size(labels));

% Loop over labels
for i = 1:length(labels)
    % Store body label
    bodies{i}.name = labels{i};
    % Store body gravitational constant
    bodies{i}.GM   = cspice_bodvrd(labels{i}, 'GM', 1);
end

end

% FIND THE ANGLE THETA AT EACH GIVEN EPOCH
function fun = thetaFinder(theta_target, et_vec, bodies, frame, center)

% thetaFinder.m - At each given epoch, returns the value of the difference
%                 between theta (angle between Sun and Moon in J2000) and
%                 the target angle
%
%
% Prototype:  fun = thetaFinder(theta_target, et_vec, bodies, frame, center)
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs:     
%
%
% theta_target[1]:     Target Angle [rad]
%
%
% et_vec[nx1]:         Vector containing epoch, could be any length
%
%
% bodies[1,11]:    cell-array with struct elements containing the following
%                  fields
%                  |
%                  |--bodies{i}.name -> body label
%                  |--bodies{i}.GM   -> gravitational constant [km**3/s**2]
%
% frame[struct]:   field frame.theta has to be 'J2000'
%
%
% center[string]: 'EARTH'
%
%
% Outputs:
%
%
% fun[1,length(et_Vec)]: Difference between the angle at the epoch and the
%                         target one
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% Initialization
sunPosEMB = zeros(3, length(et_vec));
moonPosECI = zeros(3, length(et_vec));

% Computation of the position of Sun and Moon in J200 with respect to the
% EMB
for i = 1:length(et_vec)

    sunPosEMB(:,i) = cspice_spkpos(bodies{1}.name, et_vec(i), frame.theta, 'NONE', 'EMB');
    moonPosECI(:,i) = cspice_spkpos(bodies{3}.name, et_vec(i), frame.theta, 'NONE', center);

end

% Extraction of the bidimensional components
sunPosEMB_xy = sunPosEMB(1:2, :);
moonPosECI_xy = moonPosECI(1:2,:); 

% Computation of the angles 
theta_sun = atan2(sunPosEMB_xy(2, :), sunPosEMB_xy(1, :));
theta_EMB = atan2(moonPosECI_xy(2, :), moonPosECI_xy(1, :));

% Wrap to 2pi
theta = theta_sun - theta_EMB;
theta = wrapTo2Pi(theta);

% Computation of the difference with the target
fun = theta - theta_target;

end

% NBODY PROPAGATOR
function [dxdt] = nbodyShiftRHS(t, x, bodies, frame, center)
%NBODYSHIFTRHS Evaluates the right-hand-side of a N-body propagator
%   Evaluates the right-hand-side of a newtonian N-body propagator.
%   The integration centre is the can be decided by the user changing the 
%   body contained into the struct center.
%
% Inputs:
%   t      : [1,1] ephemeris time (ET SPICE), seconds past J2000 (TDB)
%   x      : [6,1] cartesian state vector wrt desired object
%   bodies : [1,n] cell-array created with function nbody_init
%
% Outputs:
%   dxdt   : [6,1] RHS, newtonian gravitational acceleration only
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

if not( strcmpi(frame.nbody, 'ECLIPJ2000') || strcmpi(frame.nbody, 'J2000') )
    msg = 'Invalid integration reference frame, select either J2000 or ECLIPJ2000';
    error(msg);
end

dxdt = zeros(6,1);
dxdt(1:3) = x(4:6);

% Extract the object position from state x

rr_center_obj = x(1:3);

% Extract GM of central body GM0

for i = 1:length(bodies)
    value = strcmpi(center,bodies{i}.name);
    if value == 1
        center_index = i;
        GM0 = bodies{i}.GM;
    end
end


% Compute contribution to acceleration of central body
dist2 = dot(rr_center_obj, rr_center_obj);
dist = sqrt(dist2);

aa_grav_center =  - GM0 * rr_center_obj /(dist*dist2);

% Loop over all bodies (except GM0)
    
for i = 1:length(bodies)

    if i ~= center_index

        % Retrieve position and velocity of i-th celestial body wrt desired
        % center in inertial frame:
        rv_center_body = cspice_spkezr(bodies{i}.name, t, frame.nbody, 'NONE', center);
        
        % Extract object position wrt. i-th celestial body
        rho = rv_center_body(1:3);

        rr_body_obj = rr_center_obj - rv_center_body(1:3);

        % Compute non-inertial terms as in the slides (d, rho, q, f):
        d = rr_body_obj ;

        q = dot(rr_center_obj,(rr_center_obj-2*rho)) / dot(rho, rho);

        f = q * (3 + 3*q + q^2) / (1 + (1 + q) ^ (3/2));

        dist2_d = dot(d,d);
        dist_d = sqrt(dist2_d);
        
        % Compute their contribution to acceleration:
        aa_grav_i = -bodies{i}.GM / (dist2_d * dist_d) * (rr_center_obj + rho * f) ;

        dxdt(4:6) = dxdt(4:6) + aa_grav_i;
    end
    
end

% Sum up acceleration to right-hand-side
dxdt(4:6) = dxdt(4:6) + aa_grav_center;

end











