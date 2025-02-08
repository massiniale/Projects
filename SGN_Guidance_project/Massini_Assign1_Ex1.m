%%  ASSIGNMENT 1 EXERCISE 1 ( GUIDANCE; PERIODIC ORBIT )

% Author: Massini Alessandro - 10990413

clc; clearvars; close all;
tic

format long g
applyDefaultSettings()

%% --------------------- EXERCISE 1.1 -----------------------------------%%
% ------------------- Find Lagrangian Points ------------------------------

% INITIALIZATION
mu = 0.012150;

% DEFINE THE POTENTIAL
U = @(x,y) 0.5*(x.^2 + y.^2) + (1-mu)./sqrt((x+mu).^2 + y.^2) + ...
mu./sqrt((x-(1-mu)).^2 + y.^2) + 0.5*mu*(1 - mu);

% DEFINE DERIVATIVES WITH RESPECT TO X AND Y 
dUdx = @(x, y) x - (1-mu)*(x+mu)./((x+mu).^2 + y.^2).^(3/2) - ...
    mu*(x-(1-mu))./((x-(1-mu)).^2 + y.^2).^(3/2);

dUdy = @(x, y) y - (1-mu)*y./((x+mu).^2 + y.^2).^(3/2) -...
    mu*y./((x-(1-mu)).^2 + y.^2).^(3/2);

% FUNCTION TO FIND COLLINEAR
fun = @(pos)[dUdx(pos(1),pos(2)); dUdy(pos(1),pos(2))];

% DEFINE THE INITIAL GUESSES
L_guess(:,1) = [0.5, 0];            
L_guess(:,2) = [1.5, 0];            
L_guess(:,3) = [-0.5, 0];           
L_guess(:,4) = [0.5, sqrt(3/2)];    
L_guess(:,5) = [0.5, -sqrt(3/2)];   

% SOLVER OPTION TO HAVE THE 10th DECIMAL ACCURACY
options = optimoptions('fsolve','FunctionTolerance',1e-10,'OptimalityTolerance',1e-10);

% FIND ALL THE LAGRANGIAN POINTS
for i = 1:length(L_guess)
    lagrangian.coordinates(i,:)  = fsolve(fun,L_guess(:,i),options);
    lagrangian.jacobiConstant(i) = 2*U(lagrangian.coordinates(i,1),lagrangian.coordinates(i,2));
end

% SET TO ZERO Y-COMPONENT OF THE COLLINEAR
lagrangian.coordinates(1:3,2) = 0;

% DEFINE XVEC TO PLOT dUdx
xvec = linspace(-2,2,1000);


% TO HAVE AN UNDERSTANDABLE PLOT
dUdx = xvec - (1-mu)*(xvec+mu)./abs(xvec+mu).^3 - mu*(xvec+mu-1)./abs(xvec+mu-1).^3;
dUdx(abs(dUdx) > 5) = NaN;

% CIRCLES FOR THE TRIANGULAR POINTS
theta = linspace(0, 2*pi, 100); 
circle1_x = -mu +  cos(theta);
circle1_y =  sin(theta);

circle2_x = 1- mu +  cos(theta);
circle2_y =  sin(theta);


% PLOT OF THE LAGRANGIAN POINTS 
figure
plot(xvec, dUdx,'k--', LineWidth=1)
hold on
plot(circle1_x, circle1_y, 'k--', 'LineWidth', 1); 
plot(circle2_x, circle2_y, 'k--', 'LineWidth', 1); 
for i = 1:length(lagrangian.coordinates)
    plot(lagrangian.coordinates(i,1), lagrangian.coordinates(i,2), 'Marker', ...
        'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize',10)
end
plot(-mu,0,'b', 'Marker','o', 'MarkerSize',18,'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'b')
plot(1-mu,0, 'Marker','o', 'MarkerSize',12, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5])
axis equal
grid on 
xlabel('x[DU]')
ylabel('y[DU]')
axis([-1.5 2.3 -1.3 1.3])

% ADD TEXT
text(lagrangian.coordinates(1,1) - 0.22, lagrangian.coordinates(1,2), ...
    'L_1', 'FontSize', 23, 'Color', 'k', 'Interpreter', 'tex');
text(lagrangian.coordinates(2,1) + 0.1, lagrangian.coordinates(2,2), ...
    'L_2', 'FontSize', 23, 'Color', 'k', 'Interpreter', 'tex');
text(lagrangian.coordinates(3,1) - 0.23, lagrangian.coordinates(3,2), ...
    'L_3', 'FontSize', 23, 'Color', 'k', 'Interpreter', 'tex');
text(lagrangian.coordinates(4,1) + 0.07, lagrangian.coordinates(4,2) + 0.15, ...
    'L_4', 'FontSize', 23, 'Color', 'k', 'Interpreter', 'tex');
text(lagrangian.coordinates(5,1) + 0.07, lagrangian.coordinates(5,2) - 0.18, ...
    'L_5', 'FontSize', 23, 'Color', 'k', 'Interpreter', 'tex');
text(-mu - 0.1, -0.19, 'Earth', 'FontSize', 20, 'Color', 'k', 'Interpreter', 'tex');
text(1 - mu -0.1, -0.19, 'Moon', 'FontSize', 20, 'Color', 'k', 'Interpreter', 'tex');


%% ------------------------ EXERCISE 1.2 --------------------------------%%
% --------------------- Find Periodic Halo Orbit ------------------------- %

% STATE INITIALIZATION
x0  = 1.068792441776;
y0  = 0;
z0  = 0.071093328515;
vx0 = 0;
vy0 = 0.319422926485;
vz0 = 0;

% CONCATENATION
initialStateWrong = [x0;y0;z0;vx0;vy0;vz0];

% PROPAGATION TIME SET A HIGH VALUE TO SHOW NON PERIODICITY
tf = 5;

% PROPAGATION FORWARD WITH STOPPING CRITERIA TO EXTRACT THE EVENT TIME
[~,~,te,trajectoryForw]  = tdPropagator(0,initialStateWrong,tf,mu,true);

% PROPAGATION BACKWARD WITHOUT STOPPING CRITERIA TO SHOW NON PERIODICITY
[~,~,~,trajectoryBack]  = tdPropagator(0,initialStateWrong,-tf,mu,false);

% PLOT OF THE NON-CORRECTED AND NON-PERIODIC ORBIT
figure
hold on
plot3(trajectoryForw(:,1),trajectoryForw(:,2),trajectoryForw(:,3),'LineWidth',2)
plot3(trajectoryBack(:,1),trajectoryBack(:,2),trajectoryBack(:,3),'LineWidth',2)
xlabel('x[DU]')
ylabel('y[DU]')
zlabel('z[DU]')
legend('Forward Propagation', 'Backward Propagation')
grid on
title('Non Periodic Orbit @EMB Earth-Moon rotating frame')


% INITIALIZATION OF THE ERRORS TO A HIGH VALUE
err_y   = 1;  
err_vxf = 1;    
err_vzf = 1;
err_J   = 1;

% MAXIMUM ITERATIONS AND TOLERANCE
Nmax    = 100;   
iter    = 0;    
tol     = 1e-12; 
options_findPeriodic = [Nmax iter tol];

% UPDATED COMPONENTS INITIALLY SETTED TO THE INITIAL GUESS
x0New  = x0;
z0New  = z0;
vy0New = vy0;
tfNew  = te;

% TARGET VALUE TO BE REACH AT THE END OF THE CORRECTION
y_ref   = 0;
vxf_ref = 0;
vzf_ref = 0;
J_ref   = 3.09;
target = [y_ref; vxf_ref; vzf_ref; J_ref];

% INITIALIZE VECTOR
initialGuess = [x0New;y0;z0New;vx0;vy0New;vz0];
errors = [err_y;err_vxf;err_vzf;err_J];

% FIND THE PERIODIC HALO ORBIT
[initialStateRight, tfNew, J, iter] = findPeriodic(initialGuess,tfNew, errors, target, ...
    mu, options_findPeriodic);

% PROPAGATE WITH NO STOPPING CRITERIA TO SHOW THE PERIODICITY 
[~,~,te, trajectoryPeriodic, tt_F] = tdPropagator(0,initialStateRight, tf,mu,false);

% PLOT OF THE PERIODIC ORBIT
figure
plot3(trajectoryPeriodic(:,1),trajectoryPeriodic(:,2),trajectoryPeriodic(:,3),' k')
hold on
plot3(lagrangian.coordinates(2,1), lagrangian.coordinates(2,2), 0, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 12, 'HandleVisibility', 'off');
plot3(1-mu, 0, 0, 'o', 'MarkerFaceColor', [0.8, 0.8, 0.8], 'MarkerEdgeColor', 'k', 'MarkerSize', 17, 'HandleVisibility', 'off');
text(1 - mu -0.02, -0.07, 'Moon', 'FontSize', 16, 'Color', 'k', 'Interpreter', 'tex');
text(lagrangian.coordinates(2,1) + 0.02,lagrangian.coordinates(2,2), 'L_2', 'FontSize', 16, 'Color', 'k', 'Interpreter', 'tex');
xlabel('x[DU]')
ylabel('y[DU]')
zlabel('z[DU]')
title('Halo Orbit @EMB Earth-Moon Rotating Frame')

%% ------------------------ EXERCISE 1.3 --------------------------------%%
% --------------------- Numerical Continuation ---------------------------%

% INITIALIZE JACOBI VECTOR WITH 0.05 STEP
J_vec = linspace(3.09,3.04,21);

% INITIALIZE VECTOR OF INITIAL CONDITIONS
initialStateMat = zeros(length(initialStateRight),length(J_vec));
initialStateMat(:,1) = initialStateRight;
finalTimeMat = zeros(1,length(J_vec));
finalTimeMat(1) = tfNew;

% VECTOR OF JACOBI CONSTANTS
C = zeros(length(J_vec),1);
C(1) = J;

% VECTOR OF ITERATIONS
iterMat = zeros(length(J_vec),1);
iterMat(1) = iter;

% EXPLOIT NUMERICAL CONTINUATION TO FIND THE TARGET VALUE
for i = 2:length(J_vec)
    
    % UPDATE FOR THE NUMERICAL CONTINUATION
    target(4) = J_vec(i);

    % FIND PERIODIC HALO ORBITS
    [initialStateMat(:,i),finalTimeMat(i),C(i),iterMat(i)] = findPeriodic( ...
        initialStateMat(:,i-1),finalTimeMat(i-1),errors,target, mu, options_findPeriodic);
end

% INITIALIZE CELLS
propagateCellForw = cell(length(J_vec), 1);
propagateCellBack = cell(length(J_vec), 1);

% PROPAGATE ALL THE SOLUTIONS
for i = 1:length(J_vec)

    % PROPAGATE FORWARD AND BACKWARD
    [~,~,te, trajectoryPeriodicF, tt_F] = tdPropagator(0,initialStateMat(:,i), finalTimeMat(i),mu,true);
    [~,~,~, trajectoryPeriodicB, tt_B] = tdPropagator(0,initialStateMat(:,i), -finalTimeMat(i),mu,true);
    
    % SAVE THE RESULT INTO THE CELLS
    propagateCellForw{i} =  trajectoryPeriodicF(:,1:6);
    propagateCellBack{i} =  trajectoryPeriodicB(:,1:6);
end 

% CONCATENATION
propagateForw = vertcat(propagateCellForw{:});
propagateBack = vertcat(propagateCellBack{:});


% PLOT THE DIFFERENT HALO ORBITS
figure
clr = linspace(C(1),C(end), max(length(propagateForw)));
scatter3(propagateForw(:,1),propagateForw(:,2),propagateForw(:,3),5,clr)
hold on
scatter3(propagateBack(:,1),propagateBack(:,2),propagateForw(:,3),5,clr)
cb = colorbar;
grid on 
xlabel('X')
ylabel('Y')
zlabel('Z')
plot3(lagrangian.coordinates(2,1), lagrangian.coordinates(2,2), 0, 'o', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 12);
plot3(1-mu, 0, 0, 'o', 'MarkerFaceColor', [0.8, 0.8, 0.8], 'MarkerEdgeColor', ...
    'k', 'MarkerSize', 17, 'HandleVisibility', 'off');
text(1 - mu -0.02, -0.07, 'Moon', 'FontSize', 16, 'Color', 'k', 'Interpreter', 'tex');
text(lagrangian.coordinates(2,1) + 0.02,lagrangian.coordinates(2,2), 'L_2', ...
    'FontSize', 16, 'Color', 'k', 'Interpreter', 'tex');
xlabel('x[DU]')
ylabel('y[DU]')
zlabel('z[DU]')
title(cb, 'Jacobi Constant')
cb.Position = [0.85, 0.18, 0.03, 0.7];
title('Halo Orbit families @EMB Earth-Moon Rotating Frame')

%% ---------------------- FUNCTIONS -------------------------------------%%


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
    set(groot, 'defaultAxesFontSize',19);
end

% PROPAGATOR 3D CIRCULAR RESTRICTED THREE BODY PROBLEM
function [xf,PHIf,tf, xx, tt]  = tdPropagator(t0,x0,tf,mu,varargin)

% tdPropagator.m -       Perform propagation of the state integrating the
%                        equations of motion in the  3D Circular Restricted
%                        three body problem
%
%
% Prototype: [xf,PHIf,tf, xx, tt]  = tdPropagator(t0,x0,tf,mu,varargin)
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
% x0[6x1]:  Initial state of the system (x,y,z,vx,vy,vz) 
%
% tf[1]:    Final time of propagation 
%
% mu[1]:    Mass Ratio                         
%                        
% varargin['char']: can be true or false, if true the integration stopped 
%                   when the event is reached         
%
% Outputs:
%
%
% xf[6x1]:   Final state of the system (x,y,z,vx,vy,vz)
% 
% tf[1]:     Final time of propagation
%
% PHIf[6x6]: State Transition Matrix of the final state
%
% xx[nx42]:  Full propagation of the initial state and of the elements of the STM 
%
% tt[nx1]:   Vector of time
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

    if nargin>4
        evtFlag=varargin{1};
    else
        evtFlag=true;
    end

    tof = tf - t0;

    % Initialize State Transition Matrix at t0
    Phi0 = eye(6);

    % Append to initial conditions the conditions for the STM
    x0Phi0 = [x0; Phi0(:)];
    
    % Perform integration
    options_STM = odeset('reltol', 1e-12, 'abstol', 1e-12,'Events',@(x,y) xz_plane_crossing(x,y,evtFlag));
    [tt, xx] = ode78(@(t,x) CR3BP_STM(t,x,mu), [0 tof], x0Phi0, options_STM);

    % Extract state vector and State Transition Matrix
    xf = xx(end,1:6)';
    PHIf = reshape(xx(end,7:end),6,6);
    tf = tt(end);

end

% COMPUTES THE DERIVATIVE OF THE STATES AND OF THE STM 
function [dxdt] = CR3BP_STM(~,xx, mu)

% CR3BP_STM.m - Returns the derivative of the state in the Circular Restricted
%               three body problem, in the rotating frame components
%
%
% Prototype:       [dxdt] = CR3BP_STM(~,xx, mu)
%
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs: 

% t[1]:         Current time of integration
%
% xx[42x1]:     Current state of the system (x,y,z,vx,vy,vz) 
%
% mu[1]:        Mass Ratio 
%                       
%                       
% Outputs:
%
%
% dxdt[42x1]:  Derivative state of the system
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025
    
    % Extract variables
    x  = xx(1);
    y  = xx(2);
    z  = xx(3);
    vx = xx(4);
    vy = xx(5);

    % Put PHI in matrix form
    Phi = reshape(xx(7:end),6,6);

    % Compute distances from bodies 1 and 2
    r1 = sqrt((x + mu)^2 + y^2 + z^2);
    r2 = sqrt((x + mu - 1)^2 + y^2 + z^2);

    % Compute derivative of the potential
    dUdx = x - (1-mu)/r1^3*(mu+x) + mu/r2^3*(1-mu-x);
    dUdy = y - (1-mu)/r1^3*y - mu/r2^3*y;
    dUdz = - (1-mu)/r1^3*z - mu/r2^3*z;

    % Assemble the matrix A(t)=dfdx 4x4 matrix
    A = [0,0,0,1,0,0;
         0,0,0,0,1,0;
         0,0,0,0,0,1;
         (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (3*mu*(2*mu + 2*x - 2)*(mu + x - 1))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(2*mu + 2*x)*(mu + x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)) + 1, (3*mu*y*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2), (3*mu*z*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*z*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2),  0, 2, 0;
         (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)), (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) + 1, (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2), -2, 0, 0;
         (3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)), (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2), (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*z^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*z^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2), 0, 0, 0];

    % Compute the derivative of the STM
    Phidot = A*Phi;

    % Assemble right-hand side
    dxdt = zeros(42,1);

    dxdt(1:3) = xx(4:6);
    dxdt(4)   = dUdx + 2*vy;
    dxdt(5)   = dUdy - 2*vx;
    dxdt(6)   = dUdz;
    dxdt(7:end) = Phidot(:);
    
end

% FLAG THAT CHECKS FOR XZ PLANE CROSSING OF THE ORBIT 
function [value, isterminal, direction] = xz_plane_crossing(~,y,isTerminal)
    value = y(2);
    isterminal = isTerminal;
    direction = 0;
end

% IT CORRECTS INITIAL STATE IN ORDER TO FIND A PERIODIC ORBIT WITH A 
% TARGETED JACOBI CONSTANT  
function [xx0_new, tf_new, J, iter] = findPeriodic(initialGuess, tf_initial, err, target , mu, options)

% findPeriodic.m - Computes the new initial state of the system in order to
%                  find a periodic halo orbit.
%
%
% Prototype: [xx0_new, tf_new, J, iter] = findPeriodic(initialGuess, tf_initial, err, target , mu, options)
%
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs: 
%
%
% initialGuess[6x1]: Initial State that returns a non-periodic orbit (x,y,z,vx,vy,vz) 
%
% tf_initial[1]:   Initial Time at which the orbit cross the xz plane
%
% err[4x1]: Vector containing the errors between the actual values after the
%           propagation and the targeted ones, initially they are set to 1
%
% target[4x1]: Targeted value in the final values of y, vx, vz and J 
%  
% mu[1]: Mass ratio of the system
%
% options[3x1]: It contains the max number of iteration for the while, the
%               current iteration number and the tolerance
%
%                       
% Outputs:
%
%
% xx0_New[6x1]:  Vector of the initial state that returns a periodic orbit
%
% tf_new[1]: Time at which the orbit crosses the xz plane
%
% J[1]: Final value of the orbit Jacobi Constant
%
% iter[1]: Number of iterations needed to find the solution
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% Initialization
x0_new = initialGuess(1);
y0 = initialGuess(2);
z0_new = initialGuess(3);
vx0 = initialGuess(4);
vy0_new = initialGuess(5);
vz0 = initialGuess(6);

xx0_new = [x0_new;y0;z0_new;vx0;vy0_new;vz0];

tf_new = tf_initial;

y_targ = target(1);
vx_targ = target(2);
vz_targ = target(3);
J_targ = target(4);

Nmax = options(1);
iter = options(2);
tol = options(3);

err_y = err(1);
err_vxf = err(2);
err_vzf = err(3);
err_J = err(4);

% Doesn't stop untile all the errors are below the tolerance or the
% iterations exceed the maximum number
while (abs(err_vxf)>tol || abs(err_vzf)>tol || abs(err_y)>tol || ...
        abs(err_J)>tol) && iter < Nmax

    % Propagator
    [xf,PHI,te,~]  = tdPropagator(0,xx0_new,tf_new,mu, false);

    % Jacobi Constant computation with current values
    J = jacobiConstant(xx0_new,mu);

    % Gradient of the Jacobian with respect to the state components
    gradJ = gradientJ(xx0_new,mu);

    % Extraction of the derivatives of the dynamic at the event time
    xxphi = [xf;PHI(:)];
    dxdt = CR3BP_STM(te,xxphi,mu);

    % Compute the errors
    err_y   =  xf(2) - y_targ;
    err_vxf =  xf(4) - vx_targ;
    err_vzf =  xf(6) - vz_targ;
    err_J   =  J - J_targ;

    err = [err_y;err_vxf;err_vzf;err_J];

    % Build the system matrix
    A = systemMatrixBuild(PHI,dxdt,gradJ);

    % Retrieve the differential correction
    corrections = A\err;

    % Update the corrected components 
    x0_new = xx0_new(1) - corrections(1);
    z0_new = xx0_new(3) - corrections(2);
    vy0_new = xx0_new(5) - corrections(3);
    tf_new = tf_new - corrections(4);

    % Assemble the new initial state
    xx0_new = [x0_new;y0;z0_new;vx0;vy0_new;vz0];

    % Update iteration count
    iter = iter+1;

end
end

% COMPUTE THE GRADIENT OF THE JACOBI CONSTANT
function gradJ = gradientJ(xx,mu)

% gradientJ.m - Compute gradient of J with respect to the state components
%
%
% Prototype: gradJ = gradientJ(xx,mu)
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs: 
%
%
% xx[6x1]: State components (x,y,z,vx,vy,vz) 
%  
% mu[1]: Mass ratio of the system
%
%                       
% Outputs:
%
%
% gradJ[6x1]:  Gradient of the Jacobi Constant with respect to the state
%              components
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% Initialization
x   = xx(1);
y   = xx(2);
z   = xx(3);
vx  = xx(4);
vy  = xx(5);
vz  = xx(6);

% Computation of the Gradient
gradJ = [2*x + ((2*mu + 2*x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) ...
    - (mu*(2*mu + 2*x - 2))/((mu + x - 1)^2 + y^2 + z^2)^(3/2);
         2*y - (2*mu*y)/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (2*y*(mu ...
         - 1))/((mu + x)^2 + y^2 + z^2)^(3/2);
        (2*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) - (2*mu*z)/((mu + x ...
        - 1)^2 + y^2 + z^2)^(3/2);
        -2*vx;
        -2*vy;
        -2*vz];
end

% ASSEMBLE THE DIFFERENTIAL CORRECTION SYSTEM MATRIX
function A = systemMatrixBuild(PHI, dxdt, gradJ)

% systemMatrixBuild.m - Assemble the matrix of the differential correction
%                       system
%
%
% Prototype: A = systemMatrixBuild(PHI, dxdt, gradJ)
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs: 
%
%
% PHI[6x6]: State Transition Matrix of the system
%  
% dxdt[6x1]: Derivative of the state
%
% gradJ[6x1]: Gradient of the Jacobi Constant
%
%                       
% Outputs:
%
%
% gradJ[6x1]:  Gradient of the Jacobi Constant with respect to the state
%              components
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025



    A_left_up = [PHI(2,1) PHI(2,3) PHI(2,5);
         PHI(4,1) PHI(4,3) PHI(4,5);
         PHI(6,1) PHI(6,3) PHI(6,5)];
     
    A_right_up = [dxdt(2); dxdt(4); dxdt(6)];

    A_left_down = [gradJ(1) gradJ(3) gradJ(5)];

    A_right_down = 0;

    A = [A_left_up A_right_up;
         A_left_down A_right_down];

end

% JACOBI CONSTANT COMPUTATION
 function J = jacobiConstant(xx, mu)

% jacobiConstant.m - Given the state it computes the jacobi constant
%
%
% Prototype:J = jacobiConstant(xx, mu)
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs: 
%
%
%  
% xx [6x1]: State of the System
%
% mu[1]: Mass Ratio
%
%                       
% Outputs:
%
%
% J[1]:  Jacobi constant of the state
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% Compute Potential
U = 0.5 * (xx(1)^2 + xx(2)^2) + (1 - mu)/sqrt((xx(1) + mu)^2 + xx(2)^2 ...
    + xx(3)^2) + mu/sqrt((xx(1) + mu - 1)^2 + xx(2)^2 + xx(3)^2) + ...
    0.5 * mu * (1 - mu) ;

% Compute Jacobi Constant
J = 2 * U - dot(xx(4:6),xx(4:6)) ;

end
