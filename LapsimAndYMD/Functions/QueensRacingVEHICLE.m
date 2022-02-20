function [veh] = QueensRacingVEHICLE(filename)
%% Code to generate a vehicle model for the Queen's Racing LapSim
%
% Initially written by Maurice Nayman
%
% Vehicle model generation code for use in QueensRacingLAP
% January 2022

%% Clearing Memory
tic;

clc
close all force
diary('off')
fclose('all') ;


%% Declaration of global variables
%Variables imported from files
global tireData
global M
global cgx
global cgh
global L
global FT
global RT
global Cl
global Cd
global CoP
global A
global rho
global FRCH
global RRCH
global kF
global kR
global rG
global fTP
global fIA
global rTP
global rIA

%variables derived from imported variables
global a
global b
global cgMomentArm

%% Generate car model from excel file and load in tire data
carInfo = readInfo(filename,'Info') ;
engineInfo = readTorqueCurve(filename, 'Torque Curve');

name = table2array(carInfo(1,2)) ;
type = table2array(carInfo(2,2)) ;
% index
i = 3 ;
% mass
M = str2double(table2array(carInfo(i,2))) ; i = i+1 ; % [kg]
cgx = str2double(table2array(carInfo(i,2)))/100 ; i = i+1 ; % [-]
cgh = str2double(table2array(carInfo(i,2)))/1000 ; i = i+1 ; % [m]
% wheelbase and tracks
L = str2double(table2array(carInfo(i,2)))/1000 ; i = i+1 ; % [m]
FT = str2double(table2array(carInfo(i,2)))/1000 ; i = i+1 ; % [m]
RT = str2double(table2array(carInfo(i,2)))/1000 ; i = i+1 ; % [m]
% steering
rackRatio = str2double(table2array(carInfo(i,2)))/1000 ; i = i+1 ; % [-]
% aerodynamics
Cl = str2double(table2array(carInfo(i,2))) ; i = i+1 ; % [-]
Cd = str2double(table2array(carInfo(i,2))) * -1 ; i = i+1 ; % [-]
CoP = str2double(table2array(carInfo(i,2)))/100 ; i = i+1 ; % [-]
A = str2double(table2array(carInfo(i,2))) ; i = i+1 ; % [m2]
rho = str2double(table2array(carInfo(i,2))) ; i = i+1 ; % [kg/m3]
% kinematics and stiffnesses
FRCH = str2double(table2array(carInfo(i,2)))/1000 ; i = i+1 ; % [m]
RRCH = str2double(table2array(carInfo(i,2)))/1000 ; i = i+1 ; % [m]
kF = str2double(table2array(carInfo(i,2))) * 180/pi ; i = i+1 ; % [Nm/rad]
kR = str2double(table2array(carInfo(i,2))) * 180/pi ; i = i+1 ; % [Nm/rad]
% tyre information
rG = str2double(table2array(carInfo(i,2))) ; i = i+1 ; % [-]
fTP = str2double(table2array(carInfo(i,2))) ; i = i+1 ; % [psi]
fIA = str2double(table2array(carInfo(i,2))) ; i = i+1 ; % [deg]
rTP = str2double(table2array(carInfo(i,2))) ; i = i+1 ; % [psi]
rIA = str2double(table2array(carInfo(i,2))) ; i = i+1 ; % [deg]
tyreRadius = str2double(table2array(carInfo(i,2)))/1000 ; i = i+1 ; % [m]
Cr = str2double(table2array(carInfo(i,2))) ; i = i+1 ; % [-]
vWeight = str2double(table2array(carInfo(i,2))) ; i = i+1 ; % [m/s]
tireFile = (table2array(carInfo(i,2))) ; i = i+1 ; % [-]
% engine
relativePower = str2double(table2array(carInfo(i,2))) ; i = i+1 ; % [-]
nuTherm = str2double(table2array(carInfo(i,2))) ; i = i+1 ; %[-]
fLHV = str2double(table2array(carInfo(i,2))) ; i = i+1 ; %[J/kg]
% transmission
driveType = table2array(carInfo(i,2)) ; i = i+1 ; %[-]
shiftTime = str2double(table2array(carInfo(i,2))) ; i = i+1 ; %[s]
nuPG = str2double(table2array(carInfo(i,2))) ; i = i+1 ; %[-]
nuFG = str2double(table2array(carInfo(i,2))) ; i = i+1 ; %[-]
nuGB = str2double(table2array(carInfo(i,2))) ; i = i+1 ; %[-]
PGR = str2double(table2array(carInfo(i,2))) ; i = i+1 ; %[-]
FGR = str2double(table2array(carInfo(i,2))) ; i = i+1 ; %[-]
GBGR = str2double(table2array(carInfo(i:end,2))) ;
nog = length(GBGR) ;

%Derived variables
a = (1 - cgx) * L;
b = cgx * L;
cgMomentArm = cgh - ((RRCH - FRCH)/L * a + FRCH);

tireData = load(tireFile);

%% HUD

[folder_status,folder_msg] = mkdir('QueensRacingVEHICLE Vehicles') ;
vehname = "QueensRacingVEHICLE Vehicles/"+name ;
delete(vehname+".log") ;
diary(vehname+".log") ;
disp([...
    ' _______                               ______                __           ';...
    ' __  __ \____  ___________________ ___ __  __ \_____________/_/___________';...
    ' _  / / / / / / /  ___/ ___/ __  / ___/_  /_/ / / __  / ___/ / __  / __  /';...
    ' / / / / / / / / ___/ ___// / / /_ /_  /   __/ / / / / /  / / / / / /_/ / ';...
    '/ /_/ / / /_/ / /__  /__ / / / /__/ / / /\ \  / /_/ / /__/ / / / /___  /  ';...
    '\___\_\/_____/\___/\___//_/ /_//___/ /_/  \_\/___/_/____/_/_/ /_/___/ /   ';...
    '                                                                /____/    '...
    ]) ;
disp('====================================================================================')
disp(filename)
disp('File read successfully')
disp('====================================================================================')
disp("Name: "+name)
disp("Type: "+type)
disp("Date: "+datestr(now,'dd/mm/yyyy'))
disp("Time: "+datestr(now,'HH:MM:SS'))
disp('====================================================================================')
disp('Vehicle generation started.')

%% Steering Model
C = pi/90 * [kF,kF+kR;kF*a,kF*a+kR*b] ; % steering model matrix, converting from rad back to degrees
% HUD
disp('Steering model generated successfully.')

%% Driveline Model

% fetching engine curves
enSpeedCurve = table2array(engineInfo(:,1)) ; % [rpm]
enTorqueCurve = table2array(engineInfo(:,2)) ; % [N*m]
enPowerCurve = enTorqueCurve.*enSpeedCurve*2*pi/60 ; % [W]
% memory preallocation
% wheel speed per gear for every engine speed value
wheelSpeedGear = zeros(length(enSpeedCurve),nog) ;
% vehicle speed per gear for every engine speed value
vehicleSpeedGear = zeros(length(enSpeedCurve),nog) ;
% wheel torque per gear for every engine speed value
wheelTorqueGear = zeros(length(enTorqueCurve),nog) ;
% calculating values for each gear and engine speed
for i=1:nog
    wheelSpeedGear(:,i) = enSpeedCurve/PGR/GBGR(i)/FGR ;
    vehicleSpeedGear(:,i) = wheelSpeedGear(:,i)*2*pi/60*tyreRadius ;
    wheelTorqueGear(:,i) = enTorqueCurve*PGR*GBGR(i)*FGR*nuPG*nuGB*nuFG ;
end
% minimum and maximum vehicle speeds
v_min = min(vehicleSpeedGear,[],'all') ;
v_max = max(vehicleSpeedGear,[],'all') ;
% new speed vector for fine meshing
dv = 0.5/3.6 ;
vehicleSpeed = linspace(v_min,v_max,(v_max-v_min)/dv)' ;
% memory preallocation
% gear
gear = zeros(length(vehicleSpeed),1) ;
% engine tractive force
fxEngine = zeros(length(vehicleSpeed),1) ;
% engine tractive force per gear
fx = zeros(length(vehicleSpeed),nog) ;
% optimising gear selection and calculating tractive force
for i=1:length(vehicleSpeed)
    % going through the gears
    for j=1:nog
        fx(i,j) = interp1(vehicleSpeedGear(:,j),wheelTorqueGear(:,j)/tyreRadius,vehicleSpeed(i),'linear',0) ;
    end
    % getting maximum tractive force and gear
    [fxEngine(i),gear(i)] = max(fx(i,:)) ;
end
% adding values for 0 speed to vectors for interpolation purposes at low speeds
vehicleSpeed = [0;vehicleSpeed] ;
gear = [gear(1);gear] ;
fxEngine = [fxEngine(1);fxEngine] ;
% final vectors
% engine speed
engineSpeed = FGR*GBGR(gear)*PGR.*vehicleSpeed/tyreRadius*60/2/pi ;
% wheel torque
wheelTorque = fxEngine*tyreRadius ;
% engine torque
engineTorque = wheelTorque/FGR./GBGR(gear)/PGR/nuPG/nuGB/nuFG ;
% engine power
enginePower = engineTorque.*engineSpeed*2*pi/60 ;
% HUD
disp('Driveline model generated successfully.')

%% Shifting Points and Rev Drops

% finding gear changes
gearChange = diff(gear) ; % gear change will appear as 1
% getting speed right before and after gear change
gearChange = logical([gearChange;0]+[0;gearChange]) ;
% getting engine speed at gear change
engineSpeedGearChange = engineSpeed(gearChange) ;
% getting shift points
shiftPoints = engineSpeedGearChange(1:2:length(engineSpeedGearChange)) ;
% getting arrive points
arrivePoints = engineSpeedGearChange(2:2:length(engineSpeedGearChange)) ;
% calculating revdrops
revDrops = shiftPoints-arrivePoints ;
% creating shifting table
rownames = cell(nog-1,1) ;
for i=1:nog-1
    rownames(i) = {[num2str(i,'%d'),'-',num2str(i+1,'%d')]} ;
end
shifting = table(shiftPoints,arrivePoints,revDrops,'RowNames',rownames) ;
% HUD
disp('Shift points calculated successfully.')

%% Tire coefficient generation

steerMax = 14; % deg, hardcoded
bodySlipMax = 14; % deg, hardcoded
stepSize = 1; % Points, hardcoded
dv = 2; % Delta speed between each point used to calculate mu
vMaxMu = ceil(120/3.6); % Maximum desired speed by rules. No need to evaluate mu past it

vMuCalc = linspace(0, vMaxMu, vMaxMu/dv + 1); % Speed range used to calculate averaged friction coefficient

%HUD
disp('Starting friction coefficient derivation...')
startYMD = tic;
mu = calculateSteadyMu(vMuCalc, vWeight, steerMax, bodySlipMax, stepSize);
toc(startYMD)

disp('Friction coefficient calculated successfully...')

%% Force model

% gravitational constant
g = 9.81 ;
% drive and aero factors
switch driveType
    case 'RWD'
        factorDrive = (1-cgx) ; % weight distribution
        factorAero = (1-CoP) ; % aero distribution
        drivenWheels = 2 ; % number of driven wheels
    case 'FWD'
        factorDrive = cgx ;
        factorAero = CoP ;
        drivenWheels = 2 ;
    otherwise % AWD
        factorDrive = 1 ;
        factorAero = 1 ;
        drivenWheels = 4 ;
end
% Z axis
fzMass = -M*g ;
fzAero = 1/2*rho*Cl*A*vehicleSpeed.^2 ;
fzTotal = fzMass+fzAero ;
fzTyre = (factorDrive*fzMass+factorAero*fzAero)/drivenWheels ;
% x axis
fxAero = 1/2*rho*Cd*A*vehicleSpeed.^2 ;
fxRoll = Cr*abs(fzTotal) ;
fxTyre = drivenWheels * mu .* abs(fzTyre) ;
% HUD
disp('Forces calculated successfully.')

%% GGV Map

% tyre coefficients
W = M*g;
% speed map vector
dv = 2 ;
v = (0:dv:v_max)' ;
if v(end)~=v_max
    v = [v;v_max] ;
end
% friction ellipse points
N = 45 ;
% map preallocation
GGV = zeros(length(v),2*N-1,3) ;
for i=1:length(v)
    % aero forces
    AeroDf = 1/2*rho*Cl*A*v(i)^2 ;
    AeroDr = 1/2*rho*Cd*A*v(i)^2 ;
    % rolling resistance
    RollDr = Cr*abs(-AeroDf+W) ;
    % normal load on driven wheels
    Wd = (factorDrive*W+(-AeroDf * factorAero))/drivenWheels ;
    % drag acceleration
    axDrag = (AeroDr+RollDr)/M ;
    % maximum lat acc available from tyres
    ayMax = 1/M*(mu*(W-AeroDf)) ;
    % max long acc available from tyres
    axTyreMaxAcc = 1/M*mu*Wd*drivenWheels ;
    % max long acc available from tyres
    axTyreMaxDec = -1/M*mu*(W-AeroDf) ;
    % getting power limit from engine
    axPowerLimit = 1/M*(interp1(vehicleSpeed,relativePower*fxEngine,v(i))) ;
    axPowerLimit = axPowerLimit*ones(N,1) ;
    % lat acc vector
    ay = ayMax*cosd(linspace(0,180,N))' ;
    % long acc vector
    axTyreAcc = axTyreMaxAcc*sqrt(1-(ay/ayMax).^2) ; % friction ellipse
    axAcc = min(axTyreAcc,axPowerLimit)+axDrag ; % limiting by engine power
    axDec = axTyreMaxDec*sqrt(1-(ay/ayMax).^2)+axDrag ; % friction ellipse
    % saving GGV map
    GGV(i,:,1) = [axAcc',axDec(2:end)'] ;
    GGV(i,:,2) = [ay',flipud(ay(2:end))'] ;
    GGV(i,:,3) = v(i)*ones(1,2*N-1) ;
end
% HUD
disp('GGV map generated successfully.')

%% Saving vehicle

% saving and reloading into useful form for lapsim
save(vehname+".mat")
veh = load(vehname+".mat");

%% Plot

% figure
set(0,'units','pixels') ;
SS = get(0,'screensize') ;
H = 900-90 ;
W = 900 ;
Xpos = floor((SS(3)-W)/2) ;
Ypos = floor((SS(4)-H)/2) ;
f = figure('Name','Vehicle Model','Position',[Xpos,Ypos,W,H]) ;
sgtitle(name)

% rows and columns
rows = 4 ;
cols = 2 ;

% engine curves
subplot(rows,cols,1)
hold on
title('Engine Curve')
xlabel('Engine Speed [rpm]')
yyaxis left
plot(enSpeedCurve,relativePower*enTorqueCurve)
ylabel('Engine Torque [Nm]')
grid on
xlim([enSpeedCurve(1),enSpeedCurve(end)])
yyaxis right
plot(enSpeedCurve,relativePower*enPowerCurve/745.7)
ylabel('Engine Power [Hp]')

% gearing
subplot(rows,cols,3)
hold on
title('Gearing')
xlabel('Speed [m/s]')
yyaxis left
plot(vehicleSpeed,engineSpeed)
ylabel('Engine Speed [rpm]')
grid on
xlim([vehicleSpeed(1),vehicleSpeed(end)])
yyaxis right
plot(vehicleSpeed,gear)
ylabel('Gear [-]')
ylim([gear(1)-1,gear(end)+1])

% traction model
subplot(rows,cols,[5,7])
hold on
title('Traction Model')
plot(vehicleSpeed,relativePower*fxEngine,'k','LineWidth',4)
plot(vehicleSpeed,min([relativePower*fxEngine';fxTyre']),'r','LineWidth',2)
plot(vehicleSpeed,-fxAero)
plot(vehicleSpeed,-fxRoll)
plot(vehicleSpeed,fxTyre)
for i=1:nog
    plot(vehicleSpeed(2:end),fx(:,i),'k--')
end
grid on
xlabel('Speed [m/s]')
ylabel('Force [N]')
xlim([vehicleSpeed(1),vehicleSpeed(end)])
legend({'Engine tractive force','Final tractive force','Aero drag','Rolling resistance','Max tyre tractive force','Engine tractive force per gear'},'Location','southoutside')

% ggv map
subplot(rows,cols,[2,4,6,8])
hold on
title('GGV Map')
surf(GGV(:,:,2),GGV(:,:,1),GGV(:,:,3))
grid on
xlabel('Lat acc [m/s^2]')
ylabel('Long acc [m/s^2]')
zlabel('Speed [m/s]')
view(105,5)
set(gca,'DataAspectRatio',[1 1 0.8])

% saving figure
savefig(vehname+".fig")
% HUD
disp('Plots created and saved.')

%% HUD

% HUD
disp('Vehicle generated successfully.')
% diary
diary('off') ;

end