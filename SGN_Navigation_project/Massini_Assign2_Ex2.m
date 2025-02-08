%% ASSIGNMENT 2 EXERCISE 2 ( NAVIGATION; BATCH FILTER )

% Author: Massini Alessandro - 10990413

clc; clearvars; close all 

format long g
applyDefaultSettings()

cspice_kclear()

cspice_furnsh('assignment02.tm')
addpath('sgp4');

typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)
arcsec2rad = pi / (180*3600);
sigmaRange  = 0.01;
sigmaAngles = 0.125;

% INITIALIZATION
stations = stationInitialization();

t0 = '2024-11-18T20:30:00.000';  % UTC
initialEpoch = cspice_str2et(t0);

tf = '2024-11-18T22:15:00.000';  % UTC
finalEpoch = cspice_str2et(tf);

spacecraftID = 36036;

SMOS = read_3LE(spacecraftID, 'tle\36036.3le', whichconst);

parameters.mu = SMOS.mu;
Re = cspice_bodvrd('EARTH', 'RADII', 3);
parameters.Re = Re(1);

[year,mon,day,hr,min,sec] = invjday(SMOS.jdsatepoch, SMOS.jdsatepochf);

SMOS_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
referenceEpoch = cspice_str2et(SMOS_epoch_str);

% CORRECTION FOR NUTATION, FROM EOP
dPsi = -0.114752 * arcsec2rad;     % [rad]
dEpsilon = -0.007529 * arcsec2rad; % [rad]
corrections.dPsi = dPsi;
corrections.dEpsilon = dEpsilon;

% COMPUTATION OF THE REFERENCE STATE AT REFERENCE EPOCH IN ECI
[stateRefEci, elements] = refStateEciComputation(SMOS, referenceEpoch, corrections);

% PROPAGATION FROM REFERENCE EPOCH TO INITIAL EPOCH
etvec0 = linspace(referenceEpoch,initialEpoch,100);
[initialState,~,~,~] = keplerianPropagator(stateRefEci,etvec0,parameters);

%% ------------------------EXERCISE 2.1 ---------------------------------%%
% -------------------COMPUTE VISIBILITY WINDOWS --------------------------%

% COMPUTE VISIBILITY WINDOW USING KEPLERIAN MOTION
stations = myAntennaPointing(initialEpoch, initialState, finalEpoch, stations, parameters);

% ASSIGN TO THE PROPER STATION THE RESULTS
kourou   = stations{1};
troll    = stations{2};
svalbard = stations{3};

visibilityK = kourou.el >= kourou.minimumElevation;
visibilityT = troll.el >= troll.minimumElevation;
visibilityS = svalbard.el >= svalbard.minimumElevation;

% COMPUTE VISIBILITY WINDOW
etvecVisibleK = kourou.etvec(visibilityK);
kourou.timeWindow = [cspice_et2utc(etvecVisibleK(1), 'C', 3); cspice_et2utc(etvecVisibleK(end), 'C', 3)];

etvecVisibleT = troll.etvec(visibilityT);
troll.timeWindowT = [cspice_et2utc(etvecVisibleT(1), 'C', 3); cspice_et2utc(etvecVisibleT(end), 'C', 3)];

etvecVisibleS = svalbard.etvec(visibilityS);
svalbard.timeWindowS = [cspice_et2utc(etvecVisibleS(1), 'C', 3); cspice_et2utc(etvecVisibleS(end), 'C', 3)];


% SET THE TICK IN UTC TIME
time_labels = cell(length(kourou.etvec), 1); 
for i = 1:length(kourou.etvec)
    utc_full = cspice_et2utc(kourou.etvec(i), 'C', 0); 
    time_labels{i} = utc_full(12:end);
end

num_ticks = 5;
tick_indices = round(linspace(1, length(kourou.etvec), num_ticks)); 
tick_values = kourou.etvec(tick_indices); 
tick_labels = time_labels(tick_indices); 

% PLOT OF AZIMUTH AND ELEVATION AS A FUNCTION OF TIME
figure()
subplot(1,2,1)
yyaxis right
hold on
plot(kourou.etvec, kourou.el)
plot(kourou.etvec(visibilityK), kourou.el(visibilityK),'Color',  [1, 0.84, 0], 'LineStyle','-','LineWidth',3)
xticks(tick_values)
xticklabels(tick_labels)

ylabel('Elevation [deg]')
ax = gca;
ax.YColor = 'k';

yyaxis left
plot(kourou.etvec, kourou.az)
plot(kourou.etvec(visibilityK), kourou.az(visibilityK),'Color', [0.5 0.5 0.5], 'LineStyle','-','LineWidth',3)
ylabel('Azimuth [deg]')
xlabel('18-NOV-2024')
ax.YColor = 'k';
legend( 'Azimuth','Visible Azimuth','Elevation', 'Visible Elevation','Location','southeast')

% PLOT OF VISIBLE ELEVATION AND AZIMUTH
subplot(1,2,2)
plot(kourou.az(visibilityK), kourou.el(visibilityK), '+', 'MarkerSize',10)
xlabel('Azimuth [deg]')
ylabel('Elevation [deg]')
axis([-180,180,0, 90])
legend('SMOS Visibility Window')

% PLOT OF AZIMUTH AND ELEVATION AS FUNCTION OF TIME
figure()
subplot(1,2,1)
yyaxis right
hold on
plot(troll.etvec, troll.el)
plot(troll.etvec(visibilityT), troll.el(visibilityT),'Color',  [1, 0.84, 0], 'LineStyle','-','LineWidth',3)
xticks(tick_values)
xticklabels(tick_labels)

ylabel('Elevation [deg]')
ax = gca;
ax.YColor = 'k';

yyaxis left
plot(troll.etvec, troll.az)
plot(troll.etvec(visibilityT), troll.az(visibilityT),'Color', [0.5 0.5 0.5], 'LineStyle','-','LineWidth',3)
ylabel('Azimuth [deg]')
xlabel('18-NOV-2024')
ax.YColor = 'k';

legend( 'Azimuth','Visible Azimuth','Elevation', 'Visible Elevation','Location','southeast')

% SET THE THICK IN UTC TIME
time_labels = cell(length(troll.etvec), 1); 
for i = 1:length(troll.etvec)
    utc_full = cspice_et2utc(troll.etvec(i), 'C', 0); 
    time_labels{i} = utc_full(12:end);
end

% PLOT OF VISIBLE ELEVATION AND AZIMUTH
subplot(1,2,2)
plot(troll.az(visibilityT), troll.el(visibilityT), '+', 'MarkerSize',10)
xlabel('Azimuth [deg]')
ylabel('Elevation [deg]')
axis([-180,180,0, 90])
legend('SMOS Visibility Window')

% PLOT OF AZIMUTH AND ELEVATION AS FUNCTION OF TIME
figure()
subplot(1,2,1)
yyaxis right
hold on
plot(svalbard.etvec, svalbard.el)
plot(svalbard.etvec(visibilityS), svalbard.el(visibilityS),'Color',  [1, 0.84, 0], 'LineStyle','-','LineWidth',3)
xticks(tick_values)
xticklabels(tick_labels)

ylabel('Elevation [deg]')
ax = gca;
ax.YColor = 'k';

yyaxis left
plot(svalbard.etvec, svalbard.az)
plot(svalbard.etvec(visibilityS), svalbard.az(visibilityS),'Color', [0.5 0.5 0.5], 'LineStyle','-','LineWidth',3)
ylabel('Azimuth [deg]')
xlabel('18-NOV-2024')
ax.YColor = 'k';

legend( 'Azimuth','Visible Azimuth','Elevation', 'Visible Elevation','Location','southeast')

% SET THE TICK IN UTC TIME
time_labels = cell(length(svalbard.etvec), 1); 
for i = 1:length(svalbard.etvec)
    utc_full = cspice_et2utc(svalbard.etvec(i), 'C', 0); 
    time_labels{i} = utc_full(12:end);
end

% PLOT OF VISIBLE ELEVATION AND AZIMUTH
subplot(1,2,2)
plot(svalbard.az(visibilityS), svalbard.el(visibilityS), '+', 'MarkerSize',10)
xlabel('Azimuth [deg]')
ylabel('Elevation [deg]')
axis([-180,180,0, 90])
legend('SMOS Visibility Window')

%% ------------------------EXERCISE 2.2 ---------------------------------%%
% ------------------SIMULATE MEASUREMENTS --------------------------------%

% COMPUTE THE EXPECTED MEASUREMENTS THROUGH SGP4
[stationsSGP4] = myAntennaPointingSGP4(referenceEpoch, initialEpoch, SMOS, finalEpoch, stations, corrections);

% ASSIGN TO THE PROPER STATION THE RESULTS
kourouSGP4   = stationsSGP4{1};
trollSGP4    = stationsSGP4{2};
svalbardSGP4 = stationsSGP4{3};

% IDEAL MEASUREMENTS ( WITHOUTH NOISE ) 
idealMeasK = [kourouSGP4.range; kourouSGP4.az; kourouSGP4.el]';
idealMeasT = [trollSGP4.range; trollSGP4.az; trollSGP4.el]';
idealMeasS = [svalbardSGP4.range; svalbardSGP4.az; svalbardSGP4.el]';

% NOISE DIAGONAL MATRIX
covarianceNoise = diag([sigmaRange^2;sigmaAngles^2;sigmaAngles^2]);

% ADD NOISE TO HAVE THE REAL MEASUREMENTS
kourouSGP4.realMeas   = mvnrnd(idealMeasK, covarianceNoise);
trollSGP4.realMeas    = mvnrnd(idealMeasT, covarianceNoise);
svalbardSGP4.realMeas = mvnrnd(idealMeasS, covarianceNoise);

% COMPUTE AND SAVE THE VISIBILITY REGION
kourouSGP4.visibilityMeas = kourouSGP4.realMeas(:,3) >= kourouSGP4.minimumElevation;
trollSGP4.visibilityMeas = trollSGP4.realMeas(:,3) >= trollSGP4.minimumElevation;
svalbardSGP4.visibilityMeas = svalbardSGP4.realMeas(:,3) >= svalbardSGP4.minimumElevation;

% SAVE ONLY THE MEASUREMENTS RESPECTING THE VISIBILITY CONDITION
kourouSGP4.realMeas   = kourouSGP4.realMeas(kourouSGP4.visibilityMeas,:);
trollSGP4.realMeas    = trollSGP4.realMeas(trollSGP4.visibilityMeas,:);
svalbardSGP4.realMeas = svalbardSGP4.realMeas(svalbardSGP4.visibilityMeas,:);

% SAVE ONLY THE VISIBLE TIMES
kourouSGP4.etvec = kourouSGP4.etvec(kourouSGP4.visibilityMeas);
trollSGP4.etvec = trollSGP4.etvec(trollSGP4.visibilityMeas);
svalbardSGP4.etvec = svalbardSGP4.etvec(svalbardSGP4.visibilityMeas);

% CREATE A NEW CELLS ARRAY WITH THE NEW INFORMATIONS
stationsMeas{1} = kourouSGP4;
stationsMeas{2} = trollSGP4;
stationsMeas{3} = svalbardSGP4;

% PLOT OF DIFFERENCE BETWEEN SIMULATED MEASUREMENTS AND KEPLERIAN PREDICTION
figure()
hold on
plot(kourou.az(visibilityK), kourou.el(visibilityK), '*','MarkerSize', 15)
plot(kourouSGP4.realMeas(:,2), kourouSGP4.realMeas(:,3),'x', 'MarkerSize', 15)
plot(troll.az(visibilityT), troll.el(visibilityT), '*', 'MarkerSize', 15)
plot(trollSGP4.realMeas(:,2), trollSGP4.realMeas(:,3),'x','MarkerSize', 15)
plot(svalbard.az(visibilityS), svalbard.el(visibilityS), '*', 'MarkerSize', 15)
plot(svalbardSGP4.realMeas(:,2), svalbardSGP4.realMeas(:,3),'x', 'MarkerSize', 15)
axis([-180,180,0, 90])
xlabel('Azimuth [deg]')
ylabel('Elevation [deg]', 'Color', 'k')
legend([kourou.name, ' Kepler'], [kourou.name, ' Real'], ...
    [troll.name, ' Kepler'], [troll.name, ' Real'], [svalbard.name, ' Kepler'], ...
    [svalbard.name, ' Real'])
title('Difference between simulated measurements and keplerian')

%% -----------------------EXERCISE 2.3-----------------------------------%%
% -----------------SOLVE NAVIGATION PROBLEM ------------------------------%

% INITIALIZATION
initialGuess = initialState;
weightMatrix = diag(1./[sigmaRange; sigmaAngles; sigmaAngles]);

%-------------------------------------------------------------------------%
% a) Navigation Problem with One station measurements


% COST FUNCTION AS FUNCTION OF THE VARIABLES
fun = @(initialGuess) costFunction(initialEpoch, initialGuess, finalEpoch,...
    stationsMeas{1}, weightMatrix, parameters);

% OPTIONS AND SOLVER
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', ...
    'Display', 'iter', 'StepTolerance',1e-10);
[solution,resnorm,residual,~,~,~,jacobian] = lsqnonlin(fun, initialGuess, [], [], options);

% COMPUTATION OF THE COVARIANCE MATRIX
jacobian = full(jacobian);
P = resnorm/(length(residual)-length(initialGuess)).*inv(jacobian.'*jacobian);

% COMPUTE STANDARD DEVIATION FOR POSITION AND VELOCITY
navSolutionK.pos= sqrt(trace(P(1:3,1:3)));
navSolutionK.vel = sqrt(trace(P(4:6,4:6)));

% LINEAR MAPPING FROM CARTESIAN STATE TO KEPLERIAN ELEMENTS AND ESTRACTION
% OF COVARIANCE MATRIX FOR KEP AND ST.DEV FOR SMAXIS AND INCLINATION
elementsSol = car2kep(solution, parameters);
[Pelements, navSolutionK] = linearMappingCovariance(solution, P, parameters, navSolutionK);

%-------------------------------------------------------------------------%
% 2.3b Navigation Problem with All stations measurements

% COST FUNCTION AS A FUNCTION OF TH VARIABLES
funAll = @(initialGuess) costFunction(initialEpoch, initialGuess, finalEpoch,...
    stationsMeas, weightMatrix, parameters);

% SOLVER
[solutionAll,resnormAll,residualAll,~,~,~,jacobianAll] = lsqnonlin(funAll, ...
    initialGuess, [], [], options);

% COMPUTATION OF THE COVARIANCE MATRIX
jacobianAll = full(jacobianAll);
PAll = resnormAll/(length(residualAll)-length(initialGuess)).*inv(jacobianAll.'*jacobianAll);

% COMPUTE STANDARD DEVIATION FOR POSITION AND VELOCITY
navSolutionAll.pos = sqrt(trace(PAll(1:3,1:3)));
navSolutionAll.vel = sqrt(trace(PAll(4:6,4:6)));

% LINEAR MAPPING FROM CARTESIAN STATE TO KEPLERIAN ELEMENTS AND ESTRACTION
% OF COVARIANCE MATRIX FOR KEP AND ST.DEV FOR SMAXIS AND INCLINATION
elementsSolAll = car2kep(solutionAll, parameters);
[PelementsAll, navSolutionAll] = linearMappingCovariance(solutionAll, PAll, ...
    parameters, navSolutionAll);

% ------------------------------------------------------------------------%
% 2.3c Navigation Problem with j2 and All stations measurements

% ADD J2 TO THE PARAMETERS STRUCT
parametersj2  = parameters;
parametersj2.j2 = 0.0010826269;

% COMPUTE NEW INITIAL GUESS THROUGH THE KEPLERIAN SOLVER WITH J2
[initialGuessj2,~,~,~] = keplerianPropagator(stateRefEci,etvec0,parametersj2);

% COST FUNCTION
funj2  = @(initialGuessj2) costFunction(initialEpoch, initialGuessj2, finalEpoch,...
    stationsMeas, weightMatrix, parametersj2);

% SOLVER
[solutionj2,resnormj2,residualj2,~,~,~,jacobianj2] = lsqnonlin(funj2, ...
    initialGuessj2, [], [], options);

% COMPUTATION OF THE COVARIANCE MATRIX
jacobianj2 = full(jacobianj2);
Pj2 = resnormj2/(length(residualj2)-length(initialGuessj2)).*inv(jacobianj2.'*jacobianj2);

% COMPUTE STANDARD DEVIATION FOR POSITION AND VELOCITY
navSolutionj2.pos = sqrt(trace(Pj2(1:3,1:3)));
navSolutionj2.vel = sqrt(trace(Pj2(4:6,4:6)));

% LINEAR MAPPING FROM CARTESIAN STATE TO KEPLERIAN ELEMENTS AND ESTRACTION
% OF COVARIANCE MATRIX FOR KEP AND ST.DEV FOR SMAXIS AND INCLINATION
[elementsj2] = car2kep(solutionj2, parametersj2);
[Pelementsj2, navSolutionj2] = linearMappingCovariance(solutionj2, Pj2, ...
    parametersj2, navSolutionj2);

%% ------------------------EXERCISE 2.4 ---------------------------------%%
% -----------------------TRADE OFF ANALYSIS ------------------------------%

% INITIALIZE THE CELL ARRAYS NEEDED
kourouTroll{1}    = kourouSGP4;
kourouTroll{2}    = trollSGP4;

kourouSvalbard{1} = kourouSGP4;
kourouSvalbard{2} = svalbardSGP4;

trollSvalbard{1}  = trollSGP4;
trollSvalbard{2}  = svalbardSGP4;

%-------------------------------------------------------------------------%
% Case 1: Kourou-Troll

% COST FUNCTION
funKT  = @(initialGuessj2) costFunction(initialEpoch, initialGuessj2, finalEpoch,...
    kourouTroll, weightMatrix, parametersj2);

% SOLVER
[solutionKT,resnormKT,residualKT,~,~,~,jacobianKT] = lsqnonlin(funKT, ...
    initialGuessj2, [], [], options);

% INITIALIZE THE RESULT STRUCTURE
navSolutionKT = struct();

% COMPUTE COVARIANCE MATRIX
jacobianKT = full(jacobianKT);
PKT = resnormKT/(length(residualKT)-length(initialGuessj2)).*inv(jacobianKT.'*jacobianKT);

% LINEAR MAPPING FROM CARTESIAN STATE TO KEPLERIAN ELEMENTS AND ESTRACTION
% OF COVARIANCE MATRIX FOR KEP AND ST.DEV FOR SMAXIS AND INCLINATION
[PelementsKT, navSolutionKT] = linearMappingCovariance(solutionKT, PKT, ...
    parametersj2, navSolutionKT);

% COMPUTE TOTAL COST IN DOLLARS
totalCostKT = kourouTroll{1}.costPerPass + kourouTroll{2}.costPerPass;

%-------------------------------------------------------------------------%
% Case 2: Kourou-Svalbard

% COST FUNCTION
funKS  = @(initialGuessj2) costFunction(initialEpoch, initialGuessj2, finalEpoch,...
    kourouSvalbard, weightMatrix, parametersj2);

% SOLVER
[solutionKS,resnormKS,residualKS,~,~,~,jacobianKS] = lsqnonlin(funKS, ...
    initialGuessj2, [], [], options);

% INITIALIZE THE RESULT STRUCTURE
navSolutionKS = struct();

% COMPUTE COVARIANCE MATRIX
jacobianKS = full(jacobianKS);
PKS = resnormKS/(length(residualKS)-length(initialGuessj2)).*inv(jacobianKS.'*jacobianKS);

% LINEAR MAPPING FROM CARTESIAN STATE TO KEPLERIAN ELEMENTS AND ESTRACTION
% OF COVARIANCE MATRIX FOR KEP AND ST.DEV FOR SMAXIS AND INCLINATION
[PelementsKS, navSolutionKS] = linearMappingCovariance(solutionKS, PKS, ...
    parametersj2, navSolutionKS);

% COMPUTE TOTAL COST IN DOLLARS
totalCostKS = kourouSvalbard{1}.costPerPass + kourouSvalbard{2}.costPerPass;

%-------------------------------------------------------------------------%
% Case 3: Troll-Svalbard

% COST FUNCTION
funTS  = @(initialGuessj2) costFunction(initialEpoch, initialGuessj2, finalEpoch,...
    trollSvalbard, weightMatrix, parametersj2);

% SOLVER
[solutionTS,resnormTS,residualTS,~,~,~,jacobianTS] = lsqnonlin(funTS, ...
    initialGuessj2, [], [], options);

% INITIALIZE RESULTS STRUCTURE
navSolutionTS = struct();

% COMPUTE COVARIANCE MATRIX
jacobianTS = full(jacobianTS);
PTS = resnormTS/(length(residualTS)-length(initialGuessj2)).*inv(jacobianTS.'*jacobianTS);

% LINEAR MAPPING FROM CARTESIAN STATE TO KEPLERIAN ELEMENTS AND ESTRACTION
% OF COVARIANCE MATRIX FOR KEP AND ST.DEV FOR SMAXIS AND INCLINATION
[PelementsTS, navSolutionTS] = linearMappingCovariance(solutionTS, PTS, ...
    parametersj2, navSolutionTS);

% COMPUTE TOTAL COST IN DOLLARS
totalCostTS = trollSvalbard{1}.costPerPass + trollSvalbard{2}.costPerPass;

%% ---------------------- EXERCISE 2.5 ----------------------------------%%
% --------------------- Long-Term Analysis -------------------------------%

% DEFINE LONG TERM ANALYSIS WINDOW
referenceEpochMinus = referenceEpoch - 5 * 3600;
referenceEpochPlus  = referenceEpoch + 5 * 3600;

% REMOVE THE SAVED ETVEC ON THE STATIONS
for i = 1:numel(stationsSGP4)
    stationsSGP4{i} = rmfield(stationsSGP4{i}, "etvec");
end

% CREATE TWO NEW CELLS PER FORWARD AND BACKWARD
stationsSGP4Plus  = stationsSGP4;
stationsSGP4Minus = stationsSGP4;

% COMPUTE ELEVATION
stationsSGP4Plus = myAntennaPointingSGP4(referenceEpoch, referenceEpoch, ...
    SMOS, referenceEpochPlus, stationsSGP4Plus, corrections);
stationsSGP4Minus = myAntennaPointingSGP4(referenceEpoch, referenceEpoch, ...
    SMOS, referenceEpochMinus, stationsSGP4Minus, corrections);

% PLOT ELEVATION ON TIME AND REFERENCE EPOCH
figure
plot(stationsSGP4Plus{1}.etvec/cspice_spd(), stationsSGP4Plus{1}.el,'k')
hold on
plot(stationsSGP4Minus{1}.etvec/cspice_spd(), stationsSGP4Minus{1}.el, 'k')
xline(referenceEpoch/cspice_spd(), 'r--', 'LineWidth', 1.5)
yline(kourou.minimumElevation, 'b--', 'LineWidth',1.5)
legend('Elevation over Kourou', '', 'Reference Time', 'Minimum Elevation', ...
    'Location','best', 'Fontsize', 18)
xlabel('Epoch [MJD2000]')
ylabel('Elevation [deg]')
title('SMOS elevation over Kourou')

figure
plot(stationsSGP4Plus{2}.etvec/cspice_spd(), stationsSGP4Plus{2}.el,'k')
hold on
plot(stationsSGP4Minus{2}.etvec/cspice_spd(), stationsSGP4Minus{2}.el, 'k')
xline(referenceEpoch/cspice_spd(), 'r--', 'LineWidth', 1.5)
yline(troll.minimumElevation, 'b--', 'LineWidth',1.5)
legend('Elevation over Troll', '', 'Reference Time', 'Minimum Elevation', ...
    'Location','best', 'Fontsize', 18)
xlabel('Epoch [MJD2000]')
ylabel('Elevation [deg]')
title('SMOS elevation over Troll')

figure
plot(stationsSGP4Plus{3}.etvec/cspice_spd(), stationsSGP4Plus{3}.el,'k')
hold on
plot(stationsSGP4Minus{3}.etvec/cspice_spd(), stationsSGP4Minus{3}.el, 'k')
xline(referenceEpoch/cspice_spd(), 'r--', 'LineWidth', 1.5)
yline(svalbard.minimumElevation, 'b--', 'LineWidth',1.5)
legend('Elevation over Svalbard', '', 'Reference Time', 'Minimum Elevation', ...
    'Location','best', 'Fontsize', 18)
xlabel('Epoch [MJD2000]')
ylabel('Elevation [deg]')
title('SMOS elevation over Svalbard')


%%

% CLEAR  KERNEL POOL 
cspice_kclear()

%% --------------------------FUNCTIONS-----------------------------------%%

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
    set(groot, 'defaultAxesFontSize',18);
end

% STATIONS INITIALIZATION
function station = stationInitialization()

% stationInitialization.m - Return a cell array, every cell represent a 
%                           station which is represented by a struct with fields
%
%
% Prototype: station = stationInitialization()
%
%
% Inputs:     
%
%
% Outputs:
%
%
% station[1,n] cell array with struct elements containing:
%                   name             -> station name
%                   minimumFrequency -> station frequency of acquisition
%                   minimumElevation -> minimum visible elevation
%                   costPerPass      -> cost in dollars for each pass 
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025


station{1}.name = 'KOUROU';
station{1}.minimumFrequency = 60;      % [s]
station{1}.minimumElevation = 6;       % [deg]
station{1}.costPerPass = 30000;        % [$]

station{2}.name = 'TROLL';
station{2}.minimumFrequency = 30;      % [s]
station{2}.minimumElevation = 0;       % [deg] 
station{2}.costPerPass = 35000;        % [$]

station{3}.name = 'SVALBARD';          
station{3}.minimumFrequency = 60;      % [s]
station{3}.minimumElevation = 8;       % [deg]
station{3}.costPerPass = 35000;        % [$]

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
%                         mu[1]:    Earth Gravitational Constant [km^3/s^2]
%                         Re[1]:    Earth Radius [km]
%                         j2[1]:    Earth j2 constant
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

% Check if the j2 is present
if ~isfield(parameters, "j2")

    % Perform integration in the unperturbed case
    options= odeset('reltol', 1e-12, 'abstol', 1e-12);
    [tt, xx] = ode113(@(t,x) keplerianRHS(t,x,parameters), tvec, x0, options);
else
    % Perform the integration in the j2 perturbed case
    options= odeset('reltol', 1e-12, 'abstol', 1e-12);
    [tt, xx] = ode113(@(t,x) keplerianRHSPerturbed(t,x,parameters), tvec, x0, options);
end

 % Extract final state and time 
 xf = xx(end,1:6)';
 tf = tt(end);

end

% KEPLERIAN RIGHT HAND SIDE FOR UNPERTURBED
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
%                         mu[1]:    Earth Gravitational Constant [km^3/s^2]
%                         Re[1]:    Earth Radius [km]
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

% KEPLERIAN RIGHT HAND SIDE FOR J2 PERTURBED 
function dy = keplerianRHSPerturbed(t, y, parameters)

% keplerianRHS.m - Returns the derivative state of the system in 
%                          case of j2 perturbed two-body problem
%
% PROTOTYPE:    dy = keplerianRHSPerturbed(t, y, parameters)
%
%
% INPUTS:
%
% t[1]:       time    [s]
%
% y[6x1]:     Cartesian state of the system     (rx,ry,rz,vx,vy,vz)    [km,km/s]
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                         mu[1]:    Earth Gravitational Constant [km^3/s^2]
%                         Re[1]:    Earth Radius [km]
%                         j2[1]:    Earth j2 constant
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
R_E = parameters.Re;
j2 = parameters.j2;

% Rotation from J2000 to ITRF93
rotm_eci2ecef = cspice_pxform('J2000', 'ITRF93', t);
    
% Position in ecef
recef = rotm_eci2ecef * r;
   
% Computation of J2 acceleration components 
aj2ecef(1) = (3/2*(j2*mu*R_E^2)/rnorm^4)*(recef(1)/rnorm*(5*recef(3)^2/rnorm^2-1));
aj2ecef(2) = (3/2*(j2*mu*R_E^2)/rnorm^4)*(recef(2)/rnorm*(5*recef(3)^2/rnorm^2-1));
aj2ecef(3) = (3/2*(j2*mu*R_E^2)/rnorm^4)*(recef(3)/rnorm*(5*recef(3)^2/rnorm^2-3));

% Transposition and rotation in eci
aj2ecef = aj2ecef';
aj2 = rotm_eci2ecef' * aj2ecef;

% Create Right Hand Side
dy = [v
        (-mu/rnorm^3)*r + aj2];

end

% REFERENCE STATE COMPUTATION
function [stateRefEci, elements] = refStateEciComputation(satrec, referenceEpoch, corrections)

% refStateEciComputation.m - Given the satrec structure of the satellite,
%                            the reference epoch and the needed corrections
%                            for nutation, returns the reference state of
%                            the satellite (at reference epoch) and its
%                            osculating elements
%
% PROTOTYPE: [stateRefEci, elements] = refStateEciComputation(satrec, referenceEpoch, corrections)
%
%
% INPUTS:
%
% satrec[1x1 struct]:  satrec structure of the satellite
%
% referenceEpoch:  Epoch in which TLE is given [s]
%
% corrections[1x1 struct]: Corrections needed for nutation
%                         
%                           ddPsi[1]: [rad]
%                           ddEps[1]: [rad]
%
% OUTPUTS:
%
% stateRefEci[6x1]:        Vector containing the state of the           
%                          system in cartesian components  [km, km/s]
%
% elements[1x1 struct]:    Osculating element of the state
%                           
%                           - semi major axis       [km]
%                           - eccentricity          [-]
%                           - inclination           [deg]
%                           - raan                  [deg]
%                           - argument of perigee   [deg]
%                           - mean anomaly          [deg]
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% Initialization
[satrec,posTeme,velTeme] = sgp4(satrec, 0.0);

% Computation of osculating elements
oscElements = cspice_oscelt([posTeme;velTeme], referenceEpoch, satrec.mu );

smAxis = oscElements(1) / (1 - oscElements(2));
eccentricity = oscElements(2);
inclination = oscElements(3) * cspice_dpr();
raan = oscElements(4) * cspice_dpr();
perigeeArg = oscElements(5) * cspice_dpr();
meanAnomaly = oscElements(6) * cspice_dpr();

elements.smAxis = smAxis;
elements.eccentricity = eccentricity;
elements.inclination = inclination;
elements.raan = raan;
elements.perigeeArg = perigeeArg;
elements.meanAnomaly = meanAnomaly;


% Correction for Nutation, from EOP
dPsi = corrections.dPsi;
dEpsilon = corrections.dEpsilon;

% Correction for Precession 
ttt = cspice_unitim(referenceEpoch, 'ET', 'TDT')/cspice_jyear()/100;

% Computation of state in Eci
[posRefEci, velRefEci, ~] = teme2eci(posTeme, velTeme, [0;0;0], ttt, dPsi, dEpsilon);

% Concatenation
stateRefEci = [posRefEci;velRefEci];

end

% COMPUTATION OF RANGE, AZIMUTH AND ELEVATION WITH KEPLERIAN MOTION
function [station] = myAntennaPointing(t0, initialState, tf, station, parameters)

% myAntennaPointing.m - Given a station (as a struct) or multiple stations
%                       as a cell array containing the stations structs,
%                       returns the same structs with new fields range,
%                       azimuth and elevation of a given satellite in a
%                       given time interval. Propagation with Keplerian motion
%
% PROTOTYPE: [station] = myAntennaPointing(t0, initialState, tf, station, parameters)
%
%
% INPUTS:
%
% t0[1]:  Initial Epoch [s]
%
% initialState:  Vector containing the state of system in cartesian 
%                components at the initial epoch [km, km/s]
%
% tf[1]:  Final Epoch [s]
%
% station: can be either a struct with fields or a [1xn] cell array containing
%          n-struct with fields: 
%
%                      name             -> station name
%                      minimumFrequency -> station frequency of acquisition
%                      minimumElevation -> minimum visible elevation
%                      costPerPass      -> cost in dollars for each pass 
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                       mu[1]:    Earth Gravitational Constant [km^3/s^2]
%                       Re[1]:    Earth Radius [km]
%                       j2[1]:    Earth j2 constant
%
% OUTPUTS:
%
% station: can be either a struct with fields or a [1xn] cell array containing
%          n-struct with updated fields with respect to the input one: 
%
%                      Updated Fields:
%
%                      range     -> satellite range [km]
%                      azimuth   -> satellite azimuth [deg]
%                      elevation -> satellite elevation [deg]
%                      etvec     -> time vector [s]
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% If station is a struct, it computes just one station, if not it iterates
% on the cell array computing range, azimuth and elevation for all the
% stations
if isstruct(station)

    % Initialization
    stationName = station.name;
    topoFrame = [stationName, '_TOPO'];
    minimumFrequency = station.minimumFrequency;

    initialEpoch = t0;
    finalEpoch = tf;

        % Check if the visibility window time vector is already present, if
        % and in case utilize it propagating only at the visibile time
        if ~isfield(station, "etvec")
             % Creation of vector of time based on minimum frequency acquisition of
             % the station
            npoints = round((finalEpoch - initialEpoch)/minimumFrequency)+1;
            etvec = linspace(initialEpoch, finalEpoch, npoints);
        else
            % Concatenation with the initial epoch
            etvec = [initialEpoch station.etvec];
        end

    % Propagation
    [~,~,propagate] = keplerianPropagator(initialState, etvec, parameters);

    % Extraction of the station position in ECI
    stateStationEci = cspice_spkezr(stationName, etvec, 'J2000', 'NONE', 'EARTH');

    % Station to Satellite vector computation
    stateStationSatEci = propagate' - stateStationEci;

    % Rotation matrix from J2000 to TOPO
    rotMatEci2Topo = cspice_sxform('J2000', topoFrame, etvec);

    % Rotation in TOPOCENTRIC frame
    stateStationSatTopo = zeros(6,length(etvec));
    for i = 1:length(stateStationSatEci)
        stateStationSatTopo(:,i) = rotMatEci2Topo(:,:,i)*stateStationSatEci(:,i);
    end

    % Computation of Range, Azimuth and Elevation
    rll = cspice_xfmsta(stateStationSatTopo,'RECTANGULAR','LATITUDINAL','EARTH');

    % Update of the input Structure
    station.range   = rll(1,:);   
    station.az = rll(2,:) * cspice_dpr();   
    station.el = rll(3,:) * cspice_dpr(); 

    % Save the time vector only if there isn't already one
    if ~isfield(station, "etvec")
        station.etvec = etvec;
    end

else
    % Repeating of all the steps iterating on the different input stations

    for i = 1:numel(station)
        stationName = station{i}.name;
        topoFrame = [stationName, '_TOPO'];
        minimumFrequency = station{i}.minimumFrequency;

        initialEpoch = t0; 
        finalEpoch = tf;
        
        if ~isfield(station{i}, "etvec")
            npoints = round((finalEpoch - initialEpoch)/minimumFrequency)+1;
            etvec = linspace(initialEpoch, finalEpoch, npoints);
        else
                etvec = [initialEpoch station{i}.etvec];
        end

        [~,~,propagate] = keplerianPropagator(initialState, etvec, parameters);

        stateStationEci = cspice_spkezr(stationName, etvec, 'J2000', 'NONE', 'EARTH');

        stateStationSatEci = propagate' - stateStationEci;

        rotMatEci2Topo = cspice_sxform('J2000', topoFrame, etvec);

        stateStationSatTopo = zeros(6,length(etvec));

        for j = 1:length(stateStationSatEci)
            stateStationSatTopo(:,j) = rotMatEci2Topo(:,:,j)*stateStationSatEci(:,j);
        end

        rll = cspice_xfmsta(stateStationSatTopo,'RECTANGULAR','LATITUDINAL','EARTH');

        station{i}.range   = rll(1,:);   
        station{i}.az = rll(2,:) * cspice_dpr();   
        station{i}.el = rll(3,:) * cspice_dpr(); 

        % Save the time vector only if there isn't already one
        if ~isfield(station{i}, "etvec")
            station{i}.etvec = etvec;
        end
    end
end
end

% COMPUTATION OF RANGE, AZIMUTH AND ELEVATION WITH SGP4 PROPAGATION
function [station] = myAntennaPointingSGP4(tref, t0, satrec, tf, station, corrections)

% myAntennaPointingSGP4.m - Given a station (as a struct) or multiple stations
%                           as a cell array containing the stations structs,
%                           returns the same structs with new fields range,
%                           azimuth and elevation of a given satellite in a
%                           given time interval. Propagation with SGP4
%
% PROTOTYPE: [station] = myAntennaPointingSGP4(tref, t0, satrec, tf, station, corrections)
%
%
% INPUTS:
%
% tref[1]: Reference Epoch on which the TLE is defined [s]
%
% t0[1]:  Initial Epoch [s]
%
% satrec[1x1 struct]:  satrec structure of the satellite
%
% tf[1]:  Final Epoch [s]
%
% station: can be either a struct with fields or a [1xn] cell array containing
%          n-struct with fields: 
%
%                      name             -> station name
%                      minimumFrequency -> station frequency of acquisition
%                      minimumElevation -> minimum visible elevation
%                      costPerPass      -> cost in dollars for each pass 
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                       mu[1]:    Earth Gravitational Constant [km^3/s^2]
%                       Re[1]:    Earth Radius [km]
%                       j2[1]:    Earth j2 constant
%
% OUTPUTS:
%
% station: can be either a struct with fields or a [1xn] cell array containing
%          n-struct with updated fields with respect to the input one: 
%
%                      Updated Fields:
%
%                      range     -> satellite range [km]
%                      azimuth   -> satellite azimuth [deg]
%                      elevation -> satellite elevation [deg]
%                      etvec     -> time vector [s]
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% If station is a struct, it computes just one station, if not it iterates
% on the cell array computing range, azimuth and elevation for all the
% stations
if isstruct(station)

    % Initialization
    stationName     = station.name;
    topoFrame       = [stationName, '_TOPO'];
    minimumFreq     = station.minimumFrequency;
    
    referenceEpoch  = tref;
    initialEpoch    = t0;
    finalEpoch      = tf;
    dPsi            = corrections.dPsi;
    dEpsilon        = corrections.dEpsilon;

    % Creation of time vector based on the minimum frequency acquisition of
    % the station
    npoints = round((finalEpoch - initialEpoch)/minimumFreq)+1;
    etvec = linspace(initialEpoch, finalEpoch, npoints);

    % Initialization of the position and velocity vectors
    posEci = zeros(3, npoints);
    velEci = zeros(3,npoints);

    for i = 1: npoints

        % Compute time elapsed from reference epoch in minutes
        tsince = (etvec(i) - referenceEpoch)/60;

        % Compute position and velocity at actual time
        [~,posTeme,velTeme] = sgp4(satrec, tsince);

        % Correction for precession
        ttt = cspice_unitim(etvec(i), 'ET', 'TDT')/cspice_jyear()/100;

        % Transformation from TEME to ECI
        [posEci(:,i), velEci(:,i)] = teme2eci(posTeme, velTeme, [0.0;0.0;0.0],  ttt, dPsi, dEpsilon);
    end

    % Concatenation
    stateEci = [posEci;velEci];

    % Computation of Station state in J2000
    stateStationEci = cspice_spkezr(stationName, etvec, 'J2000', 'NONE', 'EARTH');

    % Computation of Station to Satellite vector
    stateStationSatEci = stateEci - stateStationEci;

    % Rotation matrix from J2000 to TOPO
    rotMatEci2Topo = cspice_sxform('J2000', topoFrame, etvec);

    % Rotation of the vector in TOPOCENTRIC reference frame
    stateStationSatTopo = zeros(6,length(etvec));

    for i = 1:length(stateStationSatEci)
        stateStationSatTopo(:,i) = rotMatEci2Topo(:,:,i)*stateStationSatEci(:,i);
    end

    % COmputation of Range, Azimuth and Elevation
    rll = cspice_xfmsta(stateStationSatTopo,'RECTANGULAR','LATITUDINAL','EARTH');

    % Update of the input struct
    station.range   = rll(1,:);   
    station.az = rll(2,:) * cspice_dpr();   
    station.el = rll(3,:) * cspice_dpr(); 
    station.etvec = etvec;

else
    % Repeat the steps iterating on the different stations saved in one
    % element of the cell array
     for i = 1:numel(station)
        stationName = station{i}.name;
        topoFrame = [stationName, '_TOPO'];
        minimumFreq = station{i}.minimumFrequency;

        referenceEpoch  = tref;
        initialEpoch    = t0;
        finalEpoch      = tf;
        dPsi            = corrections.dPsi;
        dEpsilon        = corrections.dEpsilon;

        
        npoints = abs(round((finalEpoch - initialEpoch)/minimumFreq))+1;
        etvec = linspace(initialEpoch, finalEpoch, npoints);

        posEci = zeros(3, npoints);
        velEci = zeros(3,npoints);

        for j = 1: npoints
            tsince = (etvec(j) - referenceEpoch)/60;
            [~,posTeme,velTeme] = sgp4(satrec, tsince);
            ttt = cspice_unitim(etvec(j), 'ET', 'TDT')/cspice_jyear()/100;
            [posEci(:,j), velEci(:,j)] = teme2eci(posTeme, velTeme, [0.0;0.0;0.0],  ttt, dPsi, dEpsilon);
        end

        stateEci = [posEci;velEci];
    
        stateStationEci = cspice_spkezr(stationName, etvec, 'J2000', 'NONE', 'EARTH');

        stateStationSatEci = stateEci - stateStationEci;

        rotMatEci2Topo = cspice_sxform('J2000', topoFrame, etvec);

        stateStationSatTopo = zeros(6,length(etvec));

        for j = 1:length(stateStationSatEci)
            stateStationSatTopo(:,j) = rotMatEci2Topo(:,:,j)*stateStationSatEci(:,j);
        end

        rll = cspice_xfmsta(stateStationSatTopo,'RECTANGULAR','LATITUDINAL','EARTH');

        station{i}.range   = rll(1,:);   
        station{i}.az = rll(2,:) * cspice_dpr();   
        station{i}.el = rll(3,:) * cspice_dpr(); 
        station{i}.etvec = etvec;
     end
end

end

% RESIDUALS OF THE MEASUREMENTS
function residuals = costFunction(t0, vars, tf, station, weightMatrix, parameters)

% costFunction.m - Computes the residual between the measurements computed 
%                  through the variables of the problem and a set of real
%                  measurements given as input.  
%
% PROTOTYPE: residuals = costFunction(t0, vars, tf, station, weightMatrix, parameters)
%
%
% INPUTS:
%
%
% t0[1]:  Initial Epoch [s]
%
% vars[6x1]:  variables of the problem, state of the system at the initial
%             epoch [km, km/s]
%
% tf[1]:  Final Epoch [s]
%
% station: can be either a struct with fields or a [1xn] cell array containing
%          n-struct with fields: 
%
%                      name             -> station name
%                      minimumFrequency -> station frequency of acquisition
%                      minimumElevation -> minimum visible elevation
%                      costPerPass      -> cost in dollars for each pass 
%                      range            -> satellite range [km]
%                      azimuth          -> satellite azimuth [deg]
%                      elevation        -> satellite elevation [deg]
%                      etvec            -> time vector [s]
%                      realMeas         -> noised Range, Azimuth,
%                                          Elevation, computed in the
%                                          visibility window
%                      visibility       -> Logical array definining the
%                                          visibility window
%
% weightMatrix[3x3]: Weight Matrix for the measurements [1/km, 1/deg, 1/deg]
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                       mu[1]:    Earth Gravitational Constant [km^3/s^2]
%                       Re[1]:    Earth Radius [km]
%                       j2[1]:    Earth j2 constant
%
% OUTPUTS:
%
% 
% residuals[nx3]: difference between the measurements computed starting
%                 from the variables and the real measurements in input,
%                 the number of rows depends on how many measurements are
%                 available. The columns are respectively the errors on
%                 Range, Azimuth and Elevation. 
%                   
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% If station is a struct, it computes just one station, if not it iterates
% on the cell array computing residuals for all the measurements and then
% vertcat them

if isstruct(station)
   
   % Initialization
   realMeas = station.realMeas;
   diff = zeros(length(realMeas), 3);

   % Computation of Range, Azimuth and Elevation from variables
   stationPredicted = myAntennaPointing(t0, vars, tf, station, parameters);

   % Selection of the measurements in the visibility window
   predictedMeas = [stationPredicted.range(2:end); ...
       stationPredicted.az(2:end); stationPredicted.el(2:end)]';

   % Compute differences
   diff(:,1)   = predictedMeas(:,1) - realMeas(:,1);
   diff(:,2:3) = angdiff(predictedMeas(:,2:3) * cspice_rpd, realMeas(:,2:3) * cspice_rpd) * cspice_dpr;
   
   % Compute residuals
   residuals = (weightMatrix * diff')'; 

else
    
    % Repeat the steps iterating on the different stations saved in one
    % element of the cell array

    % Initialize cell array for the residuals
    residualsCell = cell(1,numel(station));

    % Compute predicted measurements
    stationPredicted = myAntennaPointing(t0, vars, tf, station, parameters);

    for i = 1:numel(station)
        realMeas = station{i}.realMeas;
        diff = zeros(length(realMeas), 3);

        predictedMeas = [stationPredicted{i}.range(2:end); ...
            stationPredicted{i}.az(2:end); stationPredicted{i}.el(2:end)]';
        
        diff(:,1)   = predictedMeas(:,1) - realMeas(:,1);
        diff(:,2:3) = angdiff(predictedMeas(:,2:3) * cspice_rpd, realMeas(:,2:3) * cspice_rpd) * cspice_dpr;
        residualsCell{i} = (weightMatrix * diff')'; 
      
    end

    % Vertcat of the elements of the cell array
    residuals = vertcat(residualsCell{:}); 
 
end

end

% TRANSFORMATION FROM CARTESIAN STATE TO KEPLERIAN ELEMENTS
function [elements, elementsVec] = car2kep(state, parameters)

% car2kep.m - Compute the Keplerian Elements of the state starting from the
%             cartesian state
%
% PROTOTYPE: [elements, elementsVec] = car2kep(state, parameters)
%
%
% INPUTS:
%
%
% state[6x1]:              Vector containing the state of the           
%                          system in cartesian components  [km, km/s]
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                         mu[1]:    Earth Gravitational Constant [km^3/s^2]
%                         Re[1]:    Earth Radius [km]
%                         j2[1]:    Earth j2 constant
%
% OUTPUTS:
%
%
% elements[1x1 struct]:    Osculating element of the state
%                           
%                           - semi major axis       [km]
%                           - eccentricity          [-]
%                           - inclination           [deg]
%                           - raan                  [deg]
%                           - argument of perigee   [deg]
%                           - true anomaly          [deg]
%
% elementsVec[6x1]: It contains the same quantities of the strucure but in
%                   a vectorial form [km,-,deg,deg,deg,deg]
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% Initialization
position = state(1:3);
velocity = state(4:6);

mu = parameters.mu;
rnorm = norm(position);

% Computation of angular momentum vector and norm
hvec = cross(position, velocity);
h = norm(hvec);

% Computation of eccentricity vector and norm
evec = cross(velocity, hvec) / mu - position / rnorm;
e = norm(evec);

% Computation of normal vector and norm
nvec = cross([0;0;1], hvec);
n = norm(nvec);

% Computation of true anomaly
if dot(position, velocity) >= 0
   trueAn = acos(dot(evec, position) / (e * rnorm));
else
    trueAn = 2*pi - acos(dot(evec, position) / (e * rnorm));
end

% Computation of inclination
i = acos(hvec(3) / h);

% Computation of right ascension of the ascending node
if nvec(2) >= 0
    raan = acos(nvec(1) / n);
else
    raan = 2 * pi - acos(nvec(1) / n);
end

% Computation of argument of perigee
if evec(3) >= 0
    perigeeArg = acos(dot(nvec,evec) / (n * e));
else
    perigeeArg = 2 * pi - acos(dot(nvec,evec) / (n * e));
end

% Computation of Semi Major Axis
smAxis = 1 / (2 / rnorm - dot(velocity,velocity) / mu);

% Creation of output structure and vector
elements.smAxis       = smAxis;
elements.eccentricity = e;
elements.inclination  = i * cspice_dpr;
elements.raan         = raan * cspice_dpr;
elements.perigeeArg   = perigeeArg * cspice_dpr;
elements.trueAn       = trueAn * cspice_dpr;

elementsVec = [smAxis; e; i * cspice_dpr; raan * cspice_dpr;
                perigeeArg * cspice_dpr; trueAn * cspice_dpr];

end

% LINEAR MAPPING FOR KEPLERIAN ELEMENTS COVARIANCE MATRIC
function [Pelements, navSolution] = linearMappingCovariance(nominalState, P, parameters, navSolution)

% linearMappingCovariance.m - Given the covariance matrix of the state and
%                             the nominal state, it computes the linear
%                             mapping given in output the covariance matrix
%                             of the elements and the standard deviation
%                             of the inclination and semi major axis
%
% PROTOTYPE: [Pelements, navSolution] = linearMappingCovariance(nominalState, P, parameters, navSolution)
%
%
% INPUTS:
%
%
% nominalState[6x1]:       Vector containing the state of the           
%                          system in cartesian components  [km, km/s]
%
% P[6x6]:    Covariance Matrix of the state [km^2, km^2/s^2]
%
% parameters[1x1 struct]: Parameters necessary for the integration
%                         
%                         mu[1]:    Earth Gravitational Constant [km^3/s^2]
%                         Re[1]:    Earth Radius [km]
%                         j2[1]:    Earth j2 constant
%
% navSolution[1x1 struct]: struct to update, can be empty
%
% OUTPUTS:
%
%
% Pelements[6x6]:   Covariance Matrix of the Keplerian Elements
%                   [km^2,-,deg^2]
%
% navSolution[1x1 struct]: Updated struct with fields:
%
%                     stdDeviationSmAxis      -> standard deviation of the  
%                                                semi major axis [km]     
%                     stdDeviationInclination -> standard deviation of
%                                                the inclination [deg]
%
%
% Contributors:
%
% Alessandro Massini, Spacecraft Guidance and Navigation course 2024/2025

% Initialization of the Elements in nominal Case
[~, nominalElements] = car2kep(nominalState, parameters);

J = zeros(length(nominalElements), length(nominalState));
for i = 1:length(nominalState)
    
    % Initialize every iteration the state as the nominal state
    state = nominalState;

    % Perturbation for the Forward Finite Differences
    epsPert = sqrt(eps) * max(1,abs(nominalState(i)));

    % Perturb only the i-th component
    state(i) = state(i) + epsPert;

    % Find the variation of the element due to the perturbation
    [~, elementsPert] = car2kep(state, parameters);

    % Compute the i-th column of the Jacobian
    J(:,i) = (nominalElements - elementsPert) / epsPert; 

end

% Compute the covariance matrix
Pelements = J * P * J';

% Extract the needed elements and compute standard deviation
navSolution.stdDeviationSmAxis = sqrt(Pelements(1,1));
navSolution.stdDeviationInclination = sqrt(Pelements(3,3));

end
