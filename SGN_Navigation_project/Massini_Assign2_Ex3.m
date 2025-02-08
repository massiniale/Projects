%% ASSIGNMENT 2 EXERCISE 3 (NAVIGATION; UNSCENTED KALMAN FILTER)

% Author: Massini Alessandro - 10990413

clc; clearvars; close all


format long g
applyDefaultSettings()
rng default

cspice_kclear()
cspice_furnsh('assignment02.tm')

% INITIALIZATION
parameters.mu = cspice_bodvrd('MOON', 'GM', 1);  % [km^3/s^2]
Rm = cspice_bodvrd('MOON','RADII', 3);
parameters.Rm = Rm(1);                           % [km]

t0 = '2024-11-18T16:30:00.000'; % UTC
tf = '2024-11-18T20:30:00.000'; % UTC
timeStep = 30;                  % [s]

% NOISE STANDARD DEVIATION
sigmaNoise = 0.1;               % [km]

alpha = 0.01;
beta  = 2;

% INITIAL COVARIANCE MATRIX
P0 =  diag([10, 1, 1, 0.001, 0.001, 0.001, 0.00001, 0.00001]); % [km^2, km^2/s^2, rad^2]

initialEpoch = cspice_str2et(t0);
finalEpoch   = cspice_str2et(tf);

orbiter.initialState = [4307.844185282820;-1317.980749248651; 2109.210101634011;   % [km]
                       -0.110997301537882; -0.509392750828585; 0.815198807994189]; % [km/s]

lander.latitude  = 78;      % [deg]
lander.longitude = 15;      % [deg]
lander.altitude  = 0;       % [km]

% COMPUTATION OF TIMESPAN STARTING FROM THE TIME STEP
npoints = round((finalEpoch - initialEpoch)/timeStep)+1;
etvec = linspace(initialEpoch, finalEpoch, npoints);

%% --------------------------EXERCISE 3.1--------------------------------%%
% --------------------Check the visibility window-------------------------%

% COMPUTATION OF STATE OF ORBITER AND LANDER, AND RANGE, AZ, EL
[orbiter, lander] = landerPointing(orbiter, lander, etvec, parameters);

% SET TICK IN UTC
time_labels = cell(length(etvec), 1); 
for i = 1:length(etvec)
    utc_full = cspice_et2utc(etvec(i), 'C', 0); 
    time_labels{i} = utc_full(12:end);
end
num_ticks = 5;
tick_indices = round(linspace(1, length(etvec), num_ticks)); 
tick_values = etvec(tick_indices); 
tick_labels = time_labels(tick_indices); 


% PLOT OF ELEVATION OVER TIME AND CHECK THAT IT IS ALWAYS GREATER THAN THE
% MINIMUM ONE
figure
plot(etvec,lander.el,'k')
yline(0,'r--', 'LineWidth',1.5)
xticks(tick_values)
xticklabels(tick_labels)
ylabel('Elevation [deg]')
xlabel('18-NOV-2024')
xlim([etvec(1), etvec(end)])
legend('Elevation', 'Minimum Elevation Threshold', 'Location','best')
title('Elevation of Orbiter over the Lander')

visibility = lander.el >= 0;

if any(~visibility)
    error('condition not met')
else
    fprintf('Relative Visibility between Lander and Orbiter respected\n')
end

%% --------------------------EXERCISE 3.2--------------------------------%%
% -----------------------Simulate Measurements----------------------------%

% SIMULATION OF THE MEASUREMENTS (POSITION COMPONENTS AND RANGE) ADDING NOISE
[measurements, landerReal] = measurmentsSimulation(orbiter, etvec, sigmaNoise);

% PLOT OF THE ORBIT AROUND THE MOON
figure
[X, Y, Z] = sphere;
hSurface = surf(X*(Rm(1)), Y*(Rm(1)), Z*(Rm(1)));
set(hSurface,'FaceColor',[0.5 0.5 0.5])
hold on; grid on; axis equal;
plot3(orbiter.trajectory(:,1), orbiter.trajectory(:,2), orbiter.trajectory(:,3), ...
    'k', 'LineWidth',2.5)
plot3(orbiter.trajectory(1,1), orbiter.trajectory(1,2), orbiter.trajectory(1,3), ...
    'rs', 'MarkerSize',10)
plot3(orbiter.trajectory(end,1), orbiter.trajectory(end,2), orbiter.trajectory(end,3), ...
    'gs', 'MarkerSize',10)
xlabel('X[km]')
ylabel('Y[km]')
zlabel('Z[km]')
legend('', 'Orbiter Trajectory in MCI', 'Initial position $@t_0$', 'Final Position $@t_f$')
title('Orbiter Trajectory around Moon @Moon J2000')

% PLOT OF THE MEASUREMENTS
figure
plot(etvec,orbiter.trajectory(:,1))
hold on
plot(etvec,orbiter.trajectory(:,2))
plot(etvec,orbiter.trajectory(:,3))
plot(etvec, landerReal.range)
xticks(tick_values)
xticklabels(tick_labels)
ylabel('Measurements [km]')
xlabel('18-NOV-2024')
xlim([etvec(1), etvec(end)])
legend('X Component', 'Y Component', 'Z Component', 'Range', 'Location','best')
title('Measurements over time')

% PLOT OF THE DIFFERENCE BETWEEN THE RANGE MEASUREMENTS AND THE NOMINAL ONE
% AMONG WITH THE STANDARD DEVIATION
figure
plot(etvec,orbiter.trajectory(:,1) - measurements(:,1), 'LineWidth',1.3)
hold on
plot(etvec,orbiter.trajectory(:,2) - measurements(:,2),'LineWidth',1.3)
plot(etvec,orbiter.trajectory(:,3) - measurements(:,3), 'LineWidth',1.3)
plot(etvec, landerReal.range' - measurements(:,4), 'LineWidth',1.3)
plot(etvec, 3 * sigmaNoise * ones(length(etvec),1), 'r--', LineWidth=2)
plot(etvec, -3 * sigmaNoise * ones(length(etvec),1), 'r--', LineWidth=2)
xticks(tick_values)
xticklabels(tick_labels)
ylabel('Noise [km]')
xlabel('18-NOV-2024')
xlim([etvec(1), etvec(end)])
legend('X Component', 'Y Component', 'Z Component', 'Range','$3\sigma$ Boundary', 'Location','best')


% VERIFICATION THAT THE NOISE IS IN THE RIGHT BOUNDARIES
inRangeSigma = (landerReal.range' - measurements(:,4) <= sigmaNoise) & ...
    (landerReal.range' - measurements(:,4) >= -sigmaNoise);
inRangeTwoSigma = (landerReal.range' - measurements(:,4) <= 2 * sigmaNoise) & ...
    (landerReal.range' - measurements(:,4) >= -2 * sigmaNoise);
inRangeThreeSigma = (landerReal.range' - measurements(:,4) <= 3 * sigmaNoise) & ...
    (landerReal.range' - measurements(:,4) >= -3 * sigmaNoise);

inRangePerc = sum(inRangeSigma) / length(inRangeSigma) * 100;
inRangePerc2 = sum(inRangeTwoSigma) / length(inRangeTwoSigma) * 100;
inRangePerc3 = sum(inRangeThreeSigma) / length(inRangeThreeSigma) * 100;

fprintf('Percentage of points in one sigma: %.2f%%\n',inRangePerc)
fprintf('Percentage of points in two sigma: %.2f%%\n',inRangePerc2)
fprintf('Percentage of points in three sigma: %.2f%%\n',inRangePerc3)

%% --------------------------EXERCISE 3.3--------------------------------%%
% -------------Estimate the lunar orbiter absolute state------------------%

% INITIALIZATION OF THE FILTER
initialState = (mvnrnd(orbiter.initialState, P0(1:6,1:6)))';

% INITIALIZATION OF A POSTERIORI ESTIMATE AND COVARIANCE
estimatePost   = initialState;
covariancePost = P0(1:6,1:6);

% INITIALIZATION OF THE VECTOR (CELL) OF ALL THE A POSTERIORI ESTIMATE (COVARIANCE)
estimate   = zeros(length(estimatePost), length(etvec));
covariance = cell(1, length(etvec)); 

estimate(:,1) = estimatePost;
covariance{1} = covariancePost;

% UNSCENTED KALMAN FILTER
for i = 2:length(etvec)

    % GENERATION OF SIGMA POINTS
    [sigma] = generateSigmaPoints(estimatePost, covariancePost, alpha, beta);

    % TIME VECTOR BETWEEN TWO DIFFERENT MEASUREMENTS FOR PROPAGATION
    tvec = [etvec(i-1) etvec(i)];

    % UNSCENTED KALMAN FILTER
    [estimatePost, covariancePost] = ukf(tvec, sigma, measurements(i,1:3), sigmaNoise, parameters);

    % FORCE SIMMETRY
    covariancePost = (covariancePost + covariancePost') / 2 ;

    % SAVE DATA INTO THE ARRAYS
    estimate(:, i) = estimatePost;
    covariance{i}  = covariancePost;      
end

% COMPUTATION OF ERROR AND STANDARD DEVIATION 
errorPosition = zeros(length(estimate),1);
errorVelocity = zeros(length(estimate),1);
stdPosition   = zeros(length(estimate),1);
stdVelocity   = zeros(length(estimate),1);

for i = 1:length(estimate)
   errorPosition(i) = sqrt((orbiter.trajectory(i,1) - estimate(1,i)')^2 + (orbiter.trajectory(i,2)- estimate(2,i)')^2 + (orbiter.trajectory(i,3) - estimate(3,i)')^2);
   stdPosition(i) = 3 * sqrt(trace(covariance{i}(1:3,1:3)));
   errorVelocity(i) = sqrt((orbiter.trajectory(i,4) - estimate(4,i)')^2 + (orbiter.trajectory(i,5)- estimate(5,i)')^2 + (orbiter.trajectory(i,6) - estimate(6,i)')^2);
   stdVelocity(i) = 3 * sqrt(trace(covariance{i}(4:6,4:6)));
end

% PLOT OF ERROR AND STANDARD DEVIATION FOR POSITION
figure
semilogy(etvec, errorPosition,'k')
hold on
semilogy(etvec, stdPosition, 'r--', LineWidth=2)
xticks(tick_values)
xticklabels(tick_labels)
ylabel('Position Error [km]')
xlabel('18-NOV-2024')
xlim([etvec(1), etvec(end)])
grid on
legend('Error', '$3\sigma$  Bound','Location','best','FontSize', 36)
title('Position error and $3\sigma$')

% PLOT OF ERROR AND STANDARD DEVIATION FOR VELOCITY

figure
semilogy(etvec, errorVelocity,'k')
hold on
semilogy(etvec, stdVelocity,'r--',LineWidth=2)
xticks(tick_values)
xticklabels(tick_labels)
ylabel('Velocity Error [km/s]')
xlabel('18-NOV-2024')
xlim([etvec(1), etvec(end)])
grid on
legend('Error', '$3\sigma$  Bound', 'Location','best','FontSize', 36)
title('Velocity error and $3\sigma$')



%% --------------------------EXERCISE 3.4--------------------------------%%
% ---------------Estimate the lunar lander coordinates--------------------%

% INITIALIZATION OF THE FILTER
orbiter.initialStateFull = [orbiter.initialState; ...
    lander.latitude * cspice_rpd(); lander.longitude * cspice_rpd()];

initialStateFull = (mvnrnd(orbiter.initialStateFull, P0))';

% INITIALIZATION OF A POSTERIORI ESTIMATE AND COVARIANCE
estimatePostFull   = initialStateFull;
covariancePostFull = P0;

% INITIALIZATION OF THE VECTOR (CELL) OF ALL THE A POSTERIORI ESTIMATES (COVARIANCES)
estimateFull   = zeros(length(estimatePostFull), length(etvec));
covarianceFull = cell(1, length(etvec)); 

estimateFull(:,1) = estimatePostFull;
covarianceFull{1} = covariancePostFull;

% UNSCENTED KALMAN FILTER
for i = 2:length(etvec)

    % SIGMA POINTS GENERATION 
    [sigma] = generateSigmaPoints(estimatePostFull, covariancePostFull, alpha, beta);

    % TIME VECTOR BETWEEN TWO DIFFERENT MEASUREMENTS FOR PROPAGATION
    tvec = [etvec(i-1) etvec(i)];

    % UNSCETED KALMAN FILTER
    [estimatePostFull, covariancePostFull] = ukf(tvec, sigma, measurements(i,:), sigmaNoise, parameters);
    
    % FORCE SIMMETRY
    covariancePostFull = (covariancePostFull + covariancePostFull') / 2;

    % SAVE RESULTS INTO THE ARRAYS
    estimateFull(:, i) = estimatePostFull;
    covarianceFull{i}  = covariancePostFull;      
end

% COMPUTATION OF ERRORS AND STANDARD DEVIATIONS
errorPositionFull = zeros(length(estimateFull),1);
errorVelocityFull = zeros(length(estimateFull),1);
stdPositionFull   = zeros(length(estimateFull),1);
stdVelocityFull   = zeros(length(estimateFull),1);
errorLatFull      = zeros(length(estimateFull),1);
stdLatFull        = zeros(length(estimateFull),1);
errorLongFull     = zeros(length(estimateFull),1);
stdLongFull       = zeros(length(estimateFull),1);

for i = 1:length(estimateFull)

   % ERROR AND STANDARD DEVIATION FOR POSITION
   errorPositionFull(i) = sqrt((orbiter.trajectory(i,1) - estimateFull(1,i)')^2 ...
       + (orbiter.trajectory(i,2)- estimateFull(2,i)')^2 + (orbiter.trajectory(i,3) - estimateFull(3,i)')^2);
   stdPositionFull(i)   = 3 * sqrt(trace(covarianceFull{i}(1:3,1:3)));

   % ERROR AND STANDARD DEVIATION FOR VELOCITY
   errorVelocityFull(i) = sqrt((orbiter.trajectory(i,4) - estimateFull(4,i)')^2 ...
       + (orbiter.trajectory(i,5)- estimateFull(5,i)')^2 + (orbiter.trajectory(i,6) - estimateFull(6,i)')^2);
   stdVelocityFull(i)   = 3 * sqrt(trace(covarianceFull{i}(4:6,4:6)));

   % ERROR AND STANDARD DEVIATION FOR LATITUDE
   errorLatFull(i) = abs(landerReal.latitude(i) - estimateFull(7,i) * cspice_dpr');
   stdLatFull(i)   = 3 * sqrt((covarianceFull{i}(7,7))) * cspice_dpr;

   % ERROR AND STANDARD DEVIATION FOR LONGITUDE
   errorLongFull(i) = abs(landerReal.longitude(i) - estimateFull(8,i) * cspice_dpr');
   stdLongFull(i)   = 3 * sqrt((covarianceFull{i}(8,8))) * cspice_dpr;
end

% PLOT ERRORS OF THE FILTER
figure
semilogy(etvec, errorPositionFull,'k')
hold on
semilogy(etvec, errorPosition, Color= [0.9290 0.6940 0.1250])
semilogy(etvec, stdPositionFull, 'r--', LineWidth=2)
semilogy(etvec, stdPosition, 'b--', LineWidth=2)
xticks(tick_values)
xticklabels(tick_labels)
ylabel('Position Error [km]')
xlabel('18-NOV-2024')
xlim([etvec(1), etvec(end)])
grid on
legend('Error with Lander Measurements', 'Error without Lander Measurements', ...
    '$3\sigma$ Bound with lander Measurements',['$3\sigma$ Bound without lander' ...
    ' Measurements'],'Location','best')
title('Position error and $3\sigma$')

% PLOT OF ERROR AND STANDARD DEVIATION FOR VELOCITY

figure
semilogy(etvec, errorVelocityFull,'k')
hold on
semilogy(etvec, errorVelocity, Color= [0.9290 0.6940 0.1250])
semilogy(etvec, stdVelocityFull, 'r--', LineWidth=2)
semilogy(etvec, stdVelocity, 'b--', LineWidth=2)
xticks(tick_values)
xticklabels(tick_labels)
ylabel('Velocity Error [km/s]')
xlabel('18-NOV-2024')
xlim([etvec(1), etvec(end)])
grid on
legend('Error with Lander Measurements', 'Error without Lander Measurements', ...
    '$3\sigma$ Bound with lander Measurements',['$3\sigma$ Bound without lander' ...
    ' Measurements'],'Location','best')
title('Velocity error and $3\sigma$')


% PLOT OF ERROR AND STANDARD DEVIATION LATITUDE
figure
plot(etvec, errorLatFull,'k')
hold on
plot(etvec, stdLatFull,'r--')
xticks(tick_values)
xticklabels(tick_labels)
ylabel('Latitude Error [deg]')
xlabel('18-NOV-2024')
xlim([etvec(1), etvec(end)])
grid on
legend('Error', '$3\sigma$ Bound', 'Location','best')
title('Latitude error and $3\sigma$')

% PLOT OF ERROR AND STANDARD DEVIATION LONGITUDE
figure
plot(etvec, errorLongFull,'k')
hold on
plot(etvec, stdLongFull,'r--')
xticks(tick_values)
xticklabels(tick_labels)
ylabel('Longitude Error [deg]')
xlabel('18-NOV-2024')
xlim([etvec(1), etvec(end)])
grid on
legend('Error', '$3\sigma$  Bound', 'Location','best')
title('Longitude error and $3\sigma$')

%%

% CLEAR  KERNEL POOL 
cspice_kclear()

%% ------------------------ FUNCTIONS -----------------------------------%%

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

% KEPLERIAN PROPAGATOR
function [xf, tf, xx, tt] = keplerianPropagator(x0, tvec, parameters)

% keplerianPropagator.m - Perform propagation of the state integrating the
%                         equations of motion in the two body problem
%
% Prototype: [xf,tf, xx,tt] = keplerianPropagator(x0, tvec, parameters)
%
%
% Inputs:     
%
%
% tvec[nx1]:    time vector of integration [s]
%
% x0[6x1]:  Initial state of the system (x,y,z,vx,vy,vz) [km, km/s]
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                         mu[1]:    Moon Gravitational Constant [km^3/s^2]
%                         Rm[1]:    Moon Radius [km]
%                         
%
% Outputs:
%
%
% xf[6x1]:   Final state of the system (x,y,z,vx,vy,vz) [km, km/s]
% 
% tf[1]:     Final time of propagation [s]
%
% xx[nx6]:   Full propagation of the initial state and of the elements of
%            the STM (x,y,z,vx,vy,z) [km, km/s]
%
% tt[nx1]:   Vector of time [s]
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% Perform integration
options= odeset('reltol', 1e-12, 'abstol', 1e-12);
[tt, xx] = ode113(@(t,x) keplerianRHS(t,x,parameters), tvec, x0, options);

% Extract final state and time 
xf = xx(end,1:6)';
tf = tt(end);

end

% TWO BODY PROBLEM RIGHT HAND SIDE
function dy = keplerianRHS(~, y, parameters)

% keplerianRHS.m - Returns the derivative state of the system in 
%                  case of an unperturbed two-body problem 
%
% PROTOTYPE:    dy = keplerianRHS(t, y, parameters)
%
%
% INPUTS:
%
% t[1]:       time (can be omitted as the system is autonomous)        [s]
%
% y[6x1]:     Cartesian state of the system     (rx,ry,rz,vx,vy,vz)    [km,km/s]
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                         mu[1]:    Moon Gravitational Constant [km^3/s^2]
%                         Rm[1]:    Moon Radius [km]
%                        
%
% OUTPUTS:
%
% dy[6x1]:                 Vector containing the derivative state of the           
%                          system in cartesian components  [km/s, km/s^2]

% Initialize position and velocity 
r = y(1:3) ;
v = y(4:6) ;

% Compute norm of position
rnorm = norm(r);

% Extract parameters
mu = parameters.mu;

% Create Right Hand Side
dy = [v
        (-mu/rnorm^3)*r];

end

% COMPUTE ORBITER AND LANDER STATE, AND RELATIVE AZ, EL AND RANGE
function [orbiter, lander] = landerPointing(orbiter, lander, etvec, parameters)

% landerPointing.m - Given the initial state of the orbiter in MCI and
%                    coordinates of the lander (LAT, LONG, ALT), it returns
%                    the propagated state of orbiter and lander in MCI,
%                    both with relative range azimuth and elevation
%
% Prototype: [orbiter, lander] = landerPointing(orbiter, lander, etvec, parameters)
%
%
% Inputs:     
%
%
% orbiter[1x1 structs]:    Contains properties of the orbiter:
%                          - initialState[6x1]: Initial state of the system in MCI
%                               (x,y,z,vx,vy,vz) [km, km/s]
%
% lander[1x1 structs]:    Contains properties of the lander:
%                          - latitude [1x1]: Latitude  [deg]
%                          - longitude[1x1]: Longitude [deg]
%                          - altitude [1x1]: Altitude  [km]
%
% etvec[1xn]: It contains the time vector in ephemeris time [s]
%
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                         mu[1]:    Moon Gravitational Constant [km^3/s^2]
%                         Rm[1]:    Moon Radius [km]
%                         
%
% Outputs:
%
%
% orbiter[1x1 structs]:    Updated properties of the orbiter:
%                          - initialState[6x1]: Initial state of the system in MCI
%                               (x,y,z,vx,vy,vz) [km, km/s]
%
%                          - finalState[6x1]: final state of the system in MCI
%                               (x,y,z,vx,vy,vz) [km, km/s]
%
%                         - trajectory[nx6]: trajectory in MCI
%                               (x,y,z,vx,vy,vz) [km, km/s]
%
%
% lander[1x1 structs]:    Contains properties of the lander:
%                          - latitude [1x1]: Latitude  [deg]
%                          - longitude[1x1]: Longitude [deg]
%                          - altitude [1x1]: Altitude  [km]
%                          - range    [1xn]: Relative Range [km]
%                          - az       [1xn]: Relative Azimuth [deg]
%                          - el       [1xn]: Relative Elevation [deg]
%                          - stateMci [6xn]: Lander State in Mci [km, km/s]
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% Initialization
Rm = parameters.Rm;

initialStateOrbiter = orbiter.initialState;

latitude  = lander.latitude * cspice_rpd();
longitude = lander.longitude * cspice_rpd();
altitude  = lander.altitude * cspice_rpd();
radius    = altitude + Rm;

% Propagation of the orbiter state
[orbiter.finalState, ~, orbiter.trajectory, ~] = keplerianPropagator...
    (initialStateOrbiter, etvec, parameters);

% Lander state in Moon Centered Moon Fixed reference frame
landerPosIAUmoon = cspice_latrec(radius, longitude, latitude);
landerStateIAUmoon = [landerPosIAUmoon; zeros(3,1)];

% Rotation Matrix from MCMF to MCI
rotMatIau2Mci = cspice_sxform('IAU_MOON', 'J2000', etvec);

% Rotation in MCI
landerStateMci = zeros(6, length(etvec));
for i = 1:length(etvec)
    landerStateMci(:,i) = rotMatIau2Mci(:,:,i)* landerStateIAUmoon;
end

% Extraction of the lander state
lander.stateMci = landerStateMci;

% Computation of relative state
relativeStateMci = orbiter.trajectory' - landerStateMci;

% Rotation from MCI to TOPOCENTRIC Reference Frame
rotMatMci2Topo = cspice_sxform('J2000', 'MOONLANDER_TOPO', etvec);

% Rotation in TOPOCENTRIC
relativeStateTopo = zeros(6,length(etvec));
for i = 1:length(etvec)
    relativeStateTopo(:,i) = rotMatMci2Topo(:,:,i)*relativeStateMci(:,i);
end

% Computation of Range, Azimuth and Elevation
rll = cspice_xfmsta(relativeStateTopo,'RECTANGULAR','LATITUDINAL','MOON');

% Extraction
lander.range   = rll(1,:);   
lander.az      = rll(2,:) * cspice_dpr();   
lander.el      = rll(3,:) * cspice_dpr(); 

end

% SIMULATE MEASUREMENTS WITH NOISE
function [measurements, landerReal] = measurmentsSimulation(orbiter, etvec, sigmaNoise)

% measurementsSimulation.m - Given the initial state of the orbiter in MCI and
%                            the noise standard deviation, computes
%                            measurements of orbiter position and range. It
%                            returns also the range, azimuth, elevation,
%                            longitude and latitude starting from the
%                            Lander Kernel
%
%
% Prototype: [measurements, landerReal] = measurmentsSimulation(orbiter, etvec, sigmaNoise)
%
%
% Inputs:     
%
%
% orbiter[1x1 struct]:    Contains properties of the orbiter:
%                          - initialState[6x1]: Initial state of the system in MCI
%                               (x,y,z,vx,vy,vz) [km, km/s]
%
%                          - finalState[6x1]: final state of the system in MCI
%                               (x,y,z,vx,vy,vz) [km, km/s]
%
%                         - trajectory[nx6]: trajectory in MCI
%                               (x,y,z,vx,vy,vz) [km, km/s]
%
%
% etvec[1xn]: It contains the time vector in ephemeris time [s]
%
%
% sigmaNoise[1x1]: Standard deviation of the measurements nois [km]
%                         
%
% Outputs:
%
%
% measurements[nx4]: Matrix which columns are measurements of, respectively, 
%       [km]         x-pos, y-pos, z-pos, of the Orbiter in MCI and Range
%                    between Orbiter and Lander             
%
%
% landerReal[1x1 structs]:    Contains properties of the lander:
%                               - latitude [1x1]: Latitude  [deg]
%                               - longitude[1x1]: Longitude [deg]
%                               - altitude [1x1]: Altitude  [km]
%                               - range    [1xn]: Relative Range [km]
%                               - az       [1xn]: Relative Azimuth [deg]
%                               - el       [1xn]: Relative Elevation [deg]
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% Creation of the Noice Covariance Matrix
covarianceNoise = diag(sigmaNoise^2 * ones(4,1));

% Initialization
trajectory = orbiter.trajectory;

% State of lander in MCI
stateLanderReal = cspice_spkezr('MOONLANDER', etvec, 'J2000', 'NONE', 'MOON');

% Computation of relative state
relativeStateMci = trajectory' - stateLanderReal;

% Rotation matrix from MCI to TOPO
rotMatMci2Topo = cspice_sxform('J2000', 'MOONLANDER_TOPO', etvec);

% Rotation in TOPO frame
relativeStateTopo = zeros(6,length(etvec));
for i = 1:length(etvec)
    relativeStateTopo(:,i) = rotMatMci2Topo(:,:,i)*relativeStateMci(:,i);
end

% Computation of Range, Azimuth and Elevation
rll = cspice_xfmsta(relativeStateTopo,'RECTANGULAR','LATITUDINAL','MOON');

landerReal.range   = rll(1,:); 
landerReal.az      = rll(2,:) * cspice_dpr(); 
landerReal.el      = rll(3,:) * cspice_dpr(); 

% Simulation of measurements
measurements = [trajectory(:,1:3) landerReal.range'];
measurements = mvnrnd(measurements,covarianceNoise);

% State of lander in MCMF
stateLanderIau =   cspice_spkezr('MOONLANDER', etvec, 'IAU_MOON', 'NONE', 'MOON');

% Longitude and Latitude
[~, long, lat] = cspice_reclat(stateLanderIau(1:3,:));

landerReal.longitude = long * cspice_dpr();
landerReal.latitude  = lat * cspice_dpr();

end

% GENERATION OF SIGMA POINTS
function [sigma] = generateSigmaPoints(mean, covariance, alpha, beta)

% generateSigmaPoints.m - Generates the sigma points starting from a given
%                         mean and covariance matrix
%
%
% Prototype: sigmaPoints = generateSigmaPoints(mean, covariance, alpha, beta)
%
%
%
% Reminder:  All the quantities, but when indicated, are scaled and so
%            adimensional
%
% Inputs: 
%
%
% mean[nx1]:           Mean state of the system  [km, km/s], if 8x8 [rad]
% 
% covariance[nxn]:     Covariance of the state  [km^2,km^2/s^2], if 8x8 [rad^2]
%            
%
% Outputs:
%
%
% sigma[1x1 struct]: Structure containing the sigma points and the weights:
%
%                        - sigmaPoints[nx(2n+1)]: Sigma Points state 
%                                       [km, km/s], if 8x8 [rad]
%                        - wc, wm[(2n+1)x1]: weights for mean and covariance
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% Initialization
n = length(mean);

% Computation of square root of Covariance matrix
lambda = alpha^2 * n - n;
M = sqrtm((lambda + n) * covariance);

% Generation of sigma points
sigmaPoints = zeros(n,2*n);
for i = 1:n
    sigmaPoints(:,i) = mean + M(:,i);
    sigmaPoints(:,n + i) = mean - M(:,i);
end

% Concatenation with the first sigma point
sigma.sigmaPoints = [mean sigmaPoints];

% Computation of weights
w = 1/(2*(n+lambda)) * ones(2*n,1);

sigma.wm = [lambda/(n+lambda); w];
sigma.wc = [lambda/(n+lambda)+(1-alpha^2+beta); w];

end

% UNSCENTED KALMAN FILTER
function [estimatePost, covariancePost] = ukf(tvec, sigma, measurements, sigmaNoise, parameters)

% ukf.m - Sequential Filter, starting from the sigma points, given in
%         input, it compute the a priori estimate, the kalman gain and it
%         updates estimate and covariance.
%
%
% Prototype: [estimatePost, covariancePost] = ukf(tvec, sigma, measurements, sigmaNoise, parameters)
%
%
% Inputs:     
%
%
% etvec[1x2]: it contains time at step i-1 and at step i [s]
%
%
% sigma[1x1 struct]: Structure containing the sigma points and the weights:
%
%                        - sigmaPoints[nx(2n+1)]: Sigma Points state 
%                                       [km, km/s], if 8x8 [rad]
%                        - wc, wm[(2n+1)x1]: weights for mean and covariance
%
%
% measurements[1x3] or [1x4]: Vector which columns are measurements of, 
%                             respectively,  x-pos, y-pos, z-pos, of the 
%       [km]                  Orbiter in MCI and Range between Orbiter and
%                             Lander. If it is [1x3] the filter computes
%                             just the estimate of the orbiter, if it is
%                             [1x4] it computes even the estimate of the
%                             latitude and longitude of the lander
%       
%
% sigmaNoise[1x1]: Standard deviation of the measurements nois [km]
%                         
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                         mu[1]:    Moon Gravitational Constant [km^3/s^2]
%                         Rm[1]:    Moon Radius [km]
%
%
% Outputs:
%
%
% estimatePost[1x6] or [1x8]: It contains the a posteriori estimate of the 
%  [km, km/s, rad]            state of the lander (if there are 3 
%                             measurements) and of the latitude and 
%                             longitude (if there are 4 measurements)          
%
%
% covariancePost[6x6] or [8x8]: It contains the a posteriori covariance
%  [km^2, km^2/s^2, rad^2]      matrix of the state of the lander 
%                               (if there are 3 measurements) and of the 
%                               latitude and longitude (if there are 4 measurements)  
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025
% Computation of Noise Covariance
covarianceNoise = diag(sigmaNoise^2 * ones(4,1));

% Check the number of sigma points, since they depend on the dimension of
% the state, it uses it as a discriminant to check if there are present
% longitude and latitude
if length(sigma.sigmaPoints) < 14

    % Extract sigma points
    sigmaPoints = sigma.sigmaPoints;
    wm          = sigma.wm;
    wc          = sigma.wc;

    % Selector for the computation of the range
    selector    = 0;
else

    % Extract sigma points
    sigmaPoints = sigma.sigmaPoints(1:6,:);

    % Extract longitude and latitude
    latitude    = sigma.sigmaPoints(7,:);
    longitude   = sigma.sigmaPoints(8,:);
    radius      = parameters.Rm;

    % Extract weights
    wm          = sigma.wm;
    wc          = sigma.wc;

    % Selector for computation of the range
    selector    = 1;
    
end

% Initialization
sigmaProp = zeros(size(sigmaPoints,2),  size(sigmaPoints,1));
range     = zeros(1, size(sigmaPoints,2));

% Propagation 
for i = 1:length(sigmaPoints)
    [sigmaProp(i,:),~, ~, ~]  = keplerianPropagator(sigmaPoints(:,i), tvec, parameters);

    % Range computation 
    if selector 
        range(i) = rangeComputation(latitude(i), longitude(i), radius, sigmaProp(i,:), tvec);
    end
end

% Concatenation in case of Range computation
if selector
    sigmaProp = [sigmaProp'; latitude; longitude];
else
    sigmaProp = sigmaProp';
end

% Computation of a priori estimate
estimatePriori = sigmaProp * wm;

% Computation of a priori covariance
covariancePriori = zeros(length(estimatePriori), length(estimatePriori));
for i = 1:size(sigmaProp, 2)
    covariancePriori = covariancePriori + wc(i) * (sigmaProp(:, i) - ...
        estimatePriori) * (sigmaProp(:, i) - estimatePriori)';
end

% Computation of the sigma points into the measurements space
if selector
    sigmaMeas = [sigmaProp(1:length(measurements)-1,:); range];
else
    sigmaMeas = sigmaProp(1:length(measurements),:);
end

% Computation of the estimate of the predicted measurements
predictedMeas = sigmaMeas * wm;

% Computation of the covariance of the predicted measurements
covarianceMeasurements = covarianceNoise(1:length(measurements),1:length(measurements));
for i = 1:size(sigmaMeas, 2)
    covarianceMeasurements = covarianceMeasurements + wc(i) * ...
        (sigmaMeas(:, i) - predictedMeas) * (sigmaMeas(:, i) - predictedMeas)';
end

% Computation of the cross-covariance matrix
crossCovariance = zeros(length(estimatePriori), length(predictedMeas));
for i = 1:size(sigmaProp, 2)
    crossCovariance = crossCovariance + wc(i) * (sigmaProp(:, i) - ...
        estimatePriori) * (sigmaMeas(:, i) - predictedMeas)';
end

% Computation of the kalman gain
kalmanGain = crossCovariance / covarianceMeasurements;

% Update of estimate and covariance, computation of a posteriori estimate
% and covariance
estimatePost = estimatePriori + kalmanGain * (measurements' - predictedMeas);
covariancePost = covariancePriori - kalmanGain * covarianceMeasurements * kalmanGain';

end

% COMPUTE RELATIVE RANGE BETWEEN ORBITER AND LANDER
function range = rangeComputation(latitude, longitude, radius, orbiterState, tvec)

% rangeComputation.m - Given latitude, longitude and orbiter state, it
%                      computes the relative range
%
%
% Prototype: range = rangeComputation(latitude, longitude, radius, orbiterState, tvec)
%
%
% Inputs:     
%
%
% latitude[1]: Lander Latitude [rad]
%
%
% longitude[1]: Lander longitude [rad]
%
%
% radius[1]: Distance of the lander from the center of the moon [km]
%
%
% orbiterState[6x1]: State of the orbiter [km, km/s]
%                         
%
% tvec[1]: Current time of evaluation [s]
%
%
% Outputs:
%
%
% range[1]: Relative range between orbiter and lander [km]
%
% 
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% Computation of lander state in Moon Centered Moon Fixed RF
landerPosIAUmoon = cspice_latrec(radius, longitude, latitude);
landerStateIAUmoon = [landerPosIAUmoon; zeros(3,1)];

% Rotation matrix from MCMF to MCI
rotMatIau2Mci = cspice_sxform('IAU_MOON', 'J2000', tvec(end));

% Computation of lander state in MCI
landerStateMci = rotMatIau2Mci * landerStateIAUmoon;

% Computation of relative state
relativeStateMci = orbiterState' - landerStateMci;

% Rotation matrix from MCI to TOPOCENTRIC
rotMatMci2Topo = cspice_sxform('J2000', 'MOONLANDER_TOPO', tvec(end));

% Rotation into TOPOCENTRIC RF
relativeStateTopo = rotMatMci2Topo*relativeStateMci;

% Range computation
rll = cspice_xfmsta(relativeStateTopo,'RECTANGULAR','LATITUDINAL','MOON');
range   = rll(1);

end
