%% ASSIGNMENT 2 EXERCISE 1 ( NAVIGATION; UNCERTAINTIES PROPAGATION )

% Author: Massini Alessandro - 10990413

clc; clearvars; close all

applyDefaultSettings();
rng default

format long g

cspice_kclear()

cspice_furnsh('assignment02.tm')

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

alpha = 1;
beta = 2;

parameters = parametersDefinition(mu,m_s,rho,om_s,om_em,l_em,R_E,R_M,DU,VU,TU);

% INITIALIZATION OF MEAN, COVARIANCE AND TIME GRID
initialState = [-0.011965533749906; -0.017025663128129; 10.718855256727338; 0.116502348513671];

ti = 1.282800225339865;
tf = 9.595124551366348;
tvec = linspace(ti,tf,5);

P0 = [1.041e-15 6.026e-17 5.647e-16 4.577e-15;
      6.026e-17 4.287e-18 4.312e-17 1.855e-16;
      5.647e-16 4.312e-17 4.432e-16 1.455e-15;
      4.577e-15 1.855e-16 1.455e-15 2.822e-14];

%% ----------------------EXERCISE 1.1------------------------------------%%
% ----------UNCERTAINTY PROPAGATION WITH LINCOV AND UT -------------------%

% LINCOV
linCov = LinCovMethod(initialState, P0, tvec, parameters);

% UT
UT = UTmethod(initialState, P0, tvec, parameters, alpha, beta);

% COMPUTATION OF ELLIPSES COORDINATES
ellipsCoorLin = errorEllipse(linCov.meanState(end,1:2), linCov.P(1:2,1:2,end));
ellipsCoorUT  = errorEllipse(UT.meanState(1:2,end)', UT.P(1:2,1:2,end));

% PLOT OF MEAN STATE AND COVARIANCE ELLIPSE FOR UT AND LINCOV
figure()
hold on
plot(ellipsCoorLin(:,1),ellipsCoorLin(:,2),'r--', 'Color', [0.4660 0.6740 0.1880], 'LineWidth',2)
plot(ellipsCoorUT(:,1),ellipsCoorUT(:,2),'k--','LineWidth',2)
plot(UT.meanState(1,1,end), UT.meanState(2,1,end), 'kx', MarkerSize=15, LineWidth=2.5)
plot(linCov.meanState(5,1), linCov.meanState(5,2), 'x', Color= [0.4660 0.6740 0.1880], MarkerSize=15, LineWidth=2.5)
xlabel('x[-]')
ylabel('y[-]')
legend('Covariance Ellipse 3$\sigma$(LinCov)', 'Covariance Ellipse 3$\sigma$(UT)',...
    'Mean Position (LinCov)', 'Mean Position (UT)', fontsize= 16)

%% ----------------------EXERCISE 1.2------------------------------------%%
% ----------UNCERTAINTY PROPAGATION WITH MONTECARLO -------------------%

% MONTECARLO
population = 1000;
montecarlo = montecarloMethod(initialState, P0, tvec, parameters, population);
ellipsCoorMonte = errorEllipse(montecarlo.meanState(1:2,end)', montecarlo.P(1:2,1:2,end));
%%
% PLOT OF THE TIME EVOLUTION OF THE MAXIMUM STANDARD DEVIATION FOR LINCOV,
% UT AND MONTECARLO METHODS FOR POSITION
figure()

% PLOT OF THE TIME EVOLUTION OF THE MAXIMUM STANDARD DEVIATION FOR POSITION
subplot(2,1,1)
hold on
plot(tvec, linCov.tSdrr, '-^','Color', [0.4660 0.6740 0.1880], 'MarkerSize', 8, 'LineWidth', 1.5)
plot(tvec, UT.tSdrr, ':ko', 'MarkerSize', 8, 'LineWidth', 1.5) 
plot(tvec, montecarlo.tSdrr, '--s','Color', [0.8500 0.3250 0.0980], 'MarkerSize', 8, 'LineWidth', 1.5) 
xticklabels([])
ylabel('3$\sqrt{max\lambda_i}$ $P_{rr}$[-]  ')
legend('LinCov', 'UT', 'MC', 'Location', 'best')
grid on

% PLOT OF THE TIME EVOLUTION OF THE MAXIMUM STANDARD DEVIATION FOR VELOCITY
subplot(2,1,2)
hold on
plot(tvec, linCov.tSdvv, '-^', 'Color', [0.4660 0.6740 0.1880], 'MarkerSize', 8, 'LineWidth', 1.5) % Rosso, triangolo
plot(tvec, UT.tSdvv, ':ko', 'MarkerSize', 8, 'LineWidth', 1.5) % Verde, cerchio
plot(tvec, montecarlo.tSdvv, '--s', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 8, 'LineWidth', 1.5) % Nero, quadrato
xlabel('TU [-]')
ylabel('3$\sqrt{max\lambda_i}$ $P_{vv}$[-]  ')
legend('LinCov', 'MC', 'UT', 'Location', 'best')
grid on

% PLOT OF MEAN STATE AND COVARIANCE ELLIPSE FOR UT, LINCOV AND MONTECARLO
% WITH MONTECARLO SAMPLES
figure()
hold on
plot(montecarlo.samples(1,:,end),montecarlo.samples(2,:,end), 'b.')
plot(ellipsCoorMonte(:,1),ellipsCoorMonte(:,2), '--',Color=[0.8500 0.3250 0.0980], LineWidth=2)
plot(ellipsCoorLin(:,1),ellipsCoorLin(:,2),'--', 'Color', [0.4660 0.6740 0.1880], 'LineWidth',2)
plot(ellipsCoorUT(:,1),ellipsCoorUT(:,2),'k--','LineWidth',2)
plot(linCov.meanState(5,1), linCov.meanState(5,2), 'x', Color= [0.4660 0.6740 0.1880], MarkerSize=15, LineWidth=2.5)
plot(UT.meanState(1,1,end), UT.meanState(2,1,end), 'kx', MarkerSize=15, LineWidth=2.5)
plot(montecarlo.meanState(1,1,end), montecarlo.meanState(2,1,end),'x', MarkerSize=15,Color=[0.8500 0.3250 0.0980], LineWidth=2.5)
xlabel('x[-]')
ylabel('y[-]')
legend('MC samples', 'Covariance Ellipse 3$\sigma$(MC)','Covariance Ellipse 3$\sigma$(LinCov)', 'Covariance Ellipse 3$\sigma$(UT)',...
    'Mean Position (LinCov)', 'Mean Position (UT)', 'Mean Position (MC)')
title('Mean and Covariance of LinCov, UT and Montecarlo')

% QQPLOT OF THE COMPONENTS OF THE MONTECARLO SAMPLES
figure()
tiledlayout(2, 2);
stateLabels = {'x', 'y', '$v_x$', '$v_y$'};
for i = 1:4
    nexttile;
    qqplot(montecarlo.samples(i,:,end))
    title('');
    legend('Normal Quantile', '', ['MC ',stateLabels{i}, ' component'],'Location', 'best', 'FontSize', 15,'Box', 'off')
    ylabel('Input Quantiles', FontSize=18)
    xlabel('Standard Normal Quantiles', FontSize=18)
    grid on;
end

%%
% CLEAR KERNEL POOL
cspice_kclear()

%% -----------------------FUNCTIONS----------------------------------------

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
function parameters = parametersDefinition (mu,m_s,rho,om_s,om_em,l_em,R_E,R_M,DU,VU,TU)

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
 parameters.R_E = R_E/DU;
 parameters.R_M = R_M/DU;


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
    options_STM = odeset('reltol', 1e-12, 'abstol', 1e-12);
    [tt, xx] = ode78(@(t,x) PBRFBP_STM(t,x,parameters), [t0 tf], x0Phi0, options_STM);


    % Extract state vector and State Transition Matrix
    xf = xx(end,1:4)';
    PHIf = reshape(xx(end,5:end),4,4);
    tf = tt(end);

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

% UNCERTAINTY PROPAGATION WITH THE LINCOV METHOD
function linCov = LinCovMethod(initialState, initialCovariance, tvec, parameters)

% LinCovMethod.m - Returns a structure containing the mean and the
%                  covariance of the state at each time instant of the time grid
%                  and the maximum eigenvalue for the final position and
%                  velocity, obtained with LinCov method
%
% Prototype: linCov = LinCovMethod(initialState, initialCovariance, tvec, parameters)
%
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs: 
%
%
% initialState[4x1]:     Initial state of the system
% 
% initialCovarianc[4x4]: Initial covariance of the state
%
% tvec[5x1]:             Time grid
%
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
% linCov[1x1 struct]:   
%
%                       meanState[5x4]: Mean state of the system
%                       P[4x4x5]: Covariance of the system
%                       tSdrr[5x1]: three times the max eigevalue of the
%                                 position submatrix 
%                       tSdvv[5x1]: three times the max eigevalue of the
%                                 velocity submatrix 
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025
  

linCov.meanState(1,:) = initialState;
linCov.P(:,:,1) = initialCovariance;

% PROPAGATION AND COMUT.PATION OF MEAN AND COVARIANCE WITH LINCOV METHOD
for i = 2:length(tvec)
    [linCov.meanState(i,:), ~, PHIf, ~, ~]  = fourBodyPropagator_STM(tvec(1),initialState,tvec(i),parameters);
    linCov.P(:,:,i) = PHIf * initialCovariance * PHIf';
end

% MAX EIGENVALUES OF THE LINCOV FOR POSTITION AND VELOCITY 
[maxEigrrLin,maxEigvvLin] = findMaxEigenvalue(linCov.P);
linCov.tSdrr = 3 * sqrt(maxEigrrLin);
linCov.tSdvv = 3 * sqrt(maxEigvvLin);
end

% UNCERTAINTY PROPAGATION WITH THE UT METHOD
function UT = UTmethod(initialState, initialCovariance, tvec, parameters, alpha, beta)

% UTmethod.m -  Returns a structure containing the mean and the
%               covariance of the state at each time instant of the time grid
%               and the maximum eigenvalue for the final position and
%               velocity, obtained with UT method
%
% Prototype: UT = UTmethod(initialState, initialCovariance, tvec, parameters, alpha, beta)
%
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs: 
%
%
% initialState[4x1]:     Initial state of the system
% 
% initialCovarianc[4x4]: Initial covariance of the state
%
% tvec[5x1]:             Time grid
%
% parameters[1x1 struct]: Parameters necessary for the computation
%                         
%                         mu[1]:    Mass Ratio 
%                         m_s[1]:   Scaled Mass of the Sun 
%                         rho[1]:   Scaled Sun-(Earth+Moon) distance
%                         om_s[1]:  Scaled Angular velocity of the Sun
%
% alpha[1], beta[1]:      Parameters needed for the UT propagation
%            
%
% Outputs:
%
%
% UT[1x1 struct]:   
%
%                       meanState[4x1x5]: Mean state of the system
%                       P[4x4x5]: Covariance of the system
%                       tSdrr[5x1]: three times the max eigevalue of the
%                                 position submatrix 
%                       tSdvv[5x1]: three times the max eigevalue of the
%                                 velocity submatrix 
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% Initialization
n = length(initialState);
lambda = alpha^2 * n - n;
sqrtP = sqrtm((lambda + n) * initialCovariance);

% SIGMA POINTS GENERATION 
sigmaPoints = generateSigmaPoints(initialState, sqrtP);

% INITIALIZATION FOR THE PROPAGATION
sigmaProp = zeros(length(sigmaPoints),4,length(tvec));
sigmaProp(:,:,1) = sigmaPoints';

% PROPAGATION
for i = 1:length(sigmaPoints)
    for j = 2:length(tvec)
        [sigmaProp(i,:,j),~, ~, ~]  = fourBodyPropagator(tvec(1),sigmaPoints(:,i),tvec(j),parameters);
    end
end
sigmaProp = permute(sigmaProp,[2,1,3]);

% WEIGHT DEFINITION
w = 1/(2*(n+lambda)) * ones(8,1);
wm = [lambda/(n+lambda); w];
wc = [lambda/(n+lambda)+(1-alpha^2+beta); w];

% COMPUTATION OF MEAN AND COVARIANCE WITH SAMPLE MEAN AND SAMPLE COVARIANCE
UT.meanState = sum(sigmaProp .*wm', 2);
UT.P(:,:,:) = pagemtimes((sigmaProp - UT.meanState) .* wc','none',(sigmaProp - UT.meanState),'transpose');

% MAX EIGENVALUES OF THE UT FOR POSITION AND VELOCITY 
[maxEigrrUT,maxEigvvUT] = findMaxEigenvalue(UT.P);
UT.tSdrr = 3 * sqrt(maxEigrrUT);
UT.tSdvv = 3 * sqrt(maxEigvvUT);
end

% UNCERTAINTY PROPAGATION WITH THE MONTECARLO METHOD
function montecarlo = montecarloMethod(initialState, initialCovariance, tvec, parameters, population)

% montecarloMethod.m -  Returns a structure containing the mean and the
%                       covariance of the state at each time instant of the time grid
%                       and the maximum eigenvalue for the final position and
%                       velocity, obtained with UT method
%
% Prototype: montecarlo = montecarloMethod(initialState, initialCovariance, tvec, parameters, population)
%
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs: 
%
%
% initialState[4x1]:     Initial state of the system
% 
% initialCovarianc[4x4]: Initial covariance of the state
%
% tvec[5x1]:             Time grid
%
% parameters[1x1 struct]: Parameters necessary for the computation
%                         
%                         mu[1]:    Mass Ratio 
%                         m_s[1]:   Scaled Mass of the Sun 
%                         rho[1]:   Scaled Sun-(Earth+Moon) distance
%                         om_s[1]:  Scaled Angular velocity of the Sun
%
% population[1]:          Samples of the montecarlo
%            
%
% Outputs:
%
%
% montecarlo[1x1 struct]:   
%
%                       meanState[4x1x5]: Mean state of the system
%                       P[4x4x5]: Covariance of the system
%                       tSdrr[5x1]: three times the max eigevalue of the
%                                 position submatrix 
%                       tSdvv[5x1]: three times the max eigevalue of the
%                                 velocity submatrix
%                       samples[4x1000x5]: state of each sample
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% INITIALIZATION OF THE SAMPLES
montePoints = mvnrnd(initialState,initialCovariance,population);
montePoints = montePoints';

% INITIALIZATION FOR PROPAGATION
sampleMonte = zeros(population,4,5);
sampleMonte(:,:,1) = montePoints';

t0 = tvec(1);

% PROPAGATION
for j = 2:length(tvec)
    tj = tvec(j);
   for i = 1:population
        [sampleMonte(i,:,j),~, ~, ~]  = fourBodyPropagator(t0,montePoints(:,i),tj,parameters);
   end
end
sampleMonte = permute(sampleMonte,[2,1,3]);

% COMPUTATION OF MEAN AND COVARIANCE WITH SAMPLE MEAN AND SAMPLE COVARIANCE
montecarlo.meanState = mean(sampleMonte,2);
montecarlo.P = pagemtimes(sampleMonte-montecarlo.meanState, 'none', ...
    sampleMonte-montecarlo.meanState, 'transpose') / (population- 1);

% MAX EIGENVALUES OF THE MONTECARLO FOR POSITION AND VELOCITY 
[maxEigrrMonte,maxEigvvMonte] = findMaxEigenvalue(montecarlo.P);
montecarlo.tSdrr = 3 * sqrt(maxEigrrMonte);
montecarlo.tSdvv = 3 * sqrt(maxEigvvMonte);

montecarlo.samples = sampleMonte;

end

% CREATION OF COVARIANCE ELLIPSE
function ellipseCoordinates = errorEllipse(meanState, P)

% errorEllipse.m -  Return the covariance ellipse for a given mean state
%                   and covariance, the ellipse is computed with a trust
%                   region of 99.73%
%
% Prototype: ellipseCoordinates = errorEllipse(meanstate, P)
%
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs: 
%
%
% meanState[4x1]:        Mean state of the system
% 
% P[4x4]:                Covariance of the state
%            
%
% Outputs:
%
%
% ellipseCoordinates[100x2]: Cartesian Coordinates of the ellipse
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025


% Extraction of eigenvectors and eigenvalue
[eigvec, eigval] = eig(P);

% Generation of the unitary circumference
theta = linspace(0,2*pi,100);
x = cos(theta);
y = sin(theta);

% Projection, rotation and traslation of the components into the
% covariance ellipsem space
ellipseCoordinates = 3 * eigvec * sqrt(eigval) * [x;y] + meanState';
ellipseCoordinates = ellipseCoordinates';

end

% SIGMA POINTS GENERATION
function sigmaPoints = generateSigmaPoints(mean, M)

% generateSigmaPoints.m - Generates the sigma points starting from a given
%                         mean and square root of the covariance matrix
%
%
% Prototype: sigmaPoints = generateSigmaPoints(mean, M)
%
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs: 
%
%
% mean[4x1]:        Mean state of the system
% 
% M[4x4]:           Square root of the covariance of the state
%            
%
% Outputs:
%
%
% sigmaPoints[4x9]: Sigma Points state 
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025


% Initialization
sigmaPoints = zeros(4,2*length(mean));

% Generation of sigma points
for i = 1:length(mean)
    sigmaPoints(:,i) = mean + M(:,i);
    sigmaPoints(:,length(mean) + i) = mean - M(:,i);
end

sigmaPoints = [mean sigmaPoints];

end

% MAXIMUM EIGENVALUE
function [maxEigrr,maxEigvv] = findMaxEigenvalue(P)

% findMaxEigenvalue.m -   Starting from the covariance matrix of the state
%                         it computes the maximum eigenvalue of position and 
%                         velocity submatrices
%
% Prototype: [maxEigrr,maxEigvv] = findMaxEigenvalue(P)
%
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs: 
%
%
% P[4x4]:   Covariance of the state
%            
%
% Outputs:
%
%
% maxEigrr[5x1]: Max Eigenvalue for Position submatrix
% maxEigvv[5x1]: Max Eigenvalue for Velocity submatrix
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025


[~,~,k] = size(P);

% Division of the covariance matrix in Sub Matrices
Prr = P(1:2,1:2,:);
Pvv = P(3:4,3:4,:);

% Initialization
maxEigrr = zeros(k,1);
maxEigvv = zeros(k,1);

for i = 1:k

    % Extraction of the eigenvalues
    [~,eigVrr]  = eig(Prr(:,:,i));
    [~,eigVvv]  = eig(Pvv(:,:,i));

    % Find the maximum
    maxEigrr(i) = max(max(eigVrr));
    maxEigvv(i) = max(max(eigVvv));
end

end
