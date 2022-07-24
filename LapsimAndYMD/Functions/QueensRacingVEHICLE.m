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


%% Read in car model from excel
veh = readCarFromExcel(filename);

%% HUD

[folder_status,folder_msg] = mkdir('QueensRacingVEHICLE Vehicles') ;
vehname = "QueensRacingVEHICLE Vehicles/"+veh.name ;
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
disp("Name: "+veh.name)
disp("Type: "+veh.type)
disp("Date: "+datestr(now,'dd/mm/yyyy'))
disp("Time: "+datestr(now,'HH:MM:SS'))
disp('====================================================================================')
disp('Vehicle generation started.')

%% Steering Model
veh.C = pi/90 * [veh.kF,veh.kF+veh.kR;veh.kF*veh.a,veh.kF*veh.a+veh.kR*veh.b] ; % steering model matrix, converting from rad back to degrees
% HUD
disp('Steering model generated successfully.')

%% Driveline Model

% fetching engine curves
veh.enSpeedCurve = table2array(veh.engineInfo(:,1)) ; % [rpm]
veh.enTorqueCurve = table2array(veh.engineInfo(:,2)) ; % [N*m]
veh.enPowerCurve = veh.enTorqueCurve.*veh.enSpeedCurve*2*pi/60 ; % [W]
% memory preallocation
% wheel speed per gear for every engine speed value
veh.wheelSpeedGear = zeros(length(veh.enSpeedCurve),veh.nog) ;
% vehicle speed per gear for every engine speed value
veh.vehicleSpeedGear = zeros(length(veh.enSpeedCurve),veh.nog) ;
% wheel torque per gear for every engine speed value
veh.wheelTorqueGear = zeros(length(veh.enTorqueCurve),veh.nog) ;
% calculating values for each gear and engine speed
for i=1:veh.nog
    veh.wheelSpeedGear(:,i) = veh.enSpeedCurve/veh.PGR/veh.GBGR(i)/veh.FGR ;
    veh.vehicleSpeedGear(:,i) = veh.wheelSpeedGear(:,i)*2*pi/60*veh.tyreRadius ;
    veh.wheelTorqueGear(:,i) = veh.enTorqueCurve*veh.PGR*veh.GBGR(i)*veh.FGR*veh.nuPG*veh.nuGB*veh.nuFG ;
end
% minimum and maximum vehicle speeds
veh.v_min = min(veh.vehicleSpeedGear,[],'all') ;
veh.v_max = max(veh.vehicleSpeedGear,[],'all') ;
% new speed vector for fine meshing
veh.dv = 0.5/3.6 ;
veh.vehicleSpeed = linspace(veh.v_min,veh.v_max,(veh.v_max-veh.v_min)/veh.dv)' ;
% memory preallocation
% gear
veh.gear = zeros(length(veh.vehicleSpeed),1) ;
% engine tractive force
veh.fxEngine = zeros(length(veh.vehicleSpeed),1) ;
% engine tractive force per gear
veh.fx = zeros(length(veh.vehicleSpeed),veh.nog) ;
% optimising gear selection and calculating tractive force
for i=1:length(veh.vehicleSpeed)
    % going through the gears
    for j=1:veh.nog
        veh.fx(i,j) = interp1(veh.vehicleSpeedGear(:,j),veh.wheelTorqueGear(:,j)/veh.tyreRadius,veh.vehicleSpeed(i),'linear',0) ;
    end
    % getting maximum tractive force and gear
    [veh.fxEngine(i),veh.gear(i)] = max(veh.fx(i,:)) ;
end
% adding values for 0 speed to vectors for interpolation purposes at low speeds
veh.vehicleSpeed = [0;veh.vehicleSpeed] ;
veh.gear = [veh.gear(1);veh.gear] ;
veh.fxEngine = [veh.fxEngine(1);veh.fxEngine] ;
% final vectors
% engine speed
veh.engineSpeed = veh.FGR*veh.GBGR(veh.gear)*veh.PGR.*veh.vehicleSpeed/veh.tyreRadius*60/2/pi ;
% wheel torque
veh.wheelTorque = veh.fxEngine*veh.tyreRadius ;
% engine torque
veh.engineTorque = veh.wheelTorque/veh.FGR./veh.GBGR(veh.gear)/veh.PGR/veh.nuPG/veh.nuGB/veh.nuFG ;
% engine power
veh.enginePower = veh.engineTorque.*veh.engineSpeed*2*pi/60 ;
% HUD
disp('Driveline model generated successfully.')

%% Shifting Points and Rev Drops

% finding gear changes
veh.gearChange = diff(veh.gear) ; % gear change will appear as 1
% getting speed right before and after gear change
veh.gearChange = logical([veh.gearChange;0]+[0;veh.gearChange]) ;
% getting engine speed at gear change
veh.engineSpeedGearChange = veh.engineSpeed(veh.gearChange) ;
% getting shift points
veh.shiftPoints = veh.engineSpeedGearChange(1:2:length(veh.engineSpeedGearChange)) ;
% getting arrive points
veh.arrivePoints = veh.engineSpeedGearChange(2:2:length(veh.engineSpeedGearChange)) ;
% calculating revdrops
veh.revDrops = veh.shiftPoints-veh.arrivePoints ;
% creating shifting table
veh.rownames = cell(veh.nog-1,1) ;
for i=1:veh.nog-1
    veh.rownames(i) = {[num2str(i,'%d'),'-',num2str(i+1,'%d')]} ;
end
veh.shifting = table(veh.shiftPoints,veh.arrivePoints,veh.revDrops) ;
% HUD
disp('Shift points calculated successfully.')

%% Tire coefficient generation

veh.steerMax = 20; % deg, hardcoded
veh.bodySlipMax = 15; % deg, hardcoded
veh.stepSize = 1; % Points, hardcoded
veh.dv = 2; % Delta speed between each point used to calculate mu
veh.vMaxMu = ceil(120/3.6); % Maximum desired speed by rules. No need to evaluate mu past it

veh.vMuCalc = linspace(0, veh.vMaxMu, veh.vMaxMu/veh.dv + 1); % Speed range used to calculate averaged friction coefficient

%HUD
disp('Starting friction coefficient derivation...')
startYMD = tic;
veh.mu = calculateSteadyMu(veh);
toc(startYMD)

disp('Friction coefficient calculated successfully...')

%% Force model

% gravitational constant
veh.g = 9.81 ;
% drive and aero factors
switch veh.driveType
    case 'RWD'
        veh.factorDrive = (1-veh.cgx) ; % weight distribution
        veh.factorAero = (1-veh.CoP) ; % aero distribution
        veh.drivenWheels = 2 ; % number of driven wheels
    case 'FWD'
        veh.factorDrive = veh.cgx ;
        veh.factorAero = veh.CoP ;
        veh.drivenWheels = 2 ;
    otherwise % AWD
        veh.factorDrive = 1 ;
        veh.factorAero = 1 ;
        veh.drivenWheels = 4 ;
end
% Z axis
veh.fzMass = -veh.M*veh.g ;
veh.fzAero = 1/2*veh.rho*veh.Cl*veh.A*veh.vehicleSpeed.^2 ;
veh.fzTotal = veh.fzMass+veh.fzAero ;
veh.fzTyre = (veh.factorDrive*veh.fzMass+veh.factorAero*veh.fzAero)/veh.drivenWheels ;
% x axis
veh.fxAero = 1/2*veh.rho*veh.Cd*veh.A*veh.vehicleSpeed.^2 ;
veh.fxRoll = veh.Cr*abs(veh.fzTotal) ;
veh.fxTyre = veh.drivenWheels * veh.mu .* abs(veh.fzTyre) ;
% HUD
disp('Forces calculated successfully.')

%% GGV Map

% tyre coefficients
veh.W = veh.M*veh.g;
% speed map vector
veh.dv = 2 ;
veh.v = (0:veh.dv:veh.v_max)' ;
if veh.v(end)~=veh.v_max
    veh.v = [veh.v;veh.v_max] ;
end
% friction ellipse points
veh.N = 45 ;
% map preallocation
veh.GGV = zeros(length(veh.v),2*veh.N-1,3) ;
for i=1:length(veh.v)
    % aero forces
    veh.AeroDf = 1/2*veh.rho*veh.Cl*veh.A*veh.v(i)^2 ;
    veh.AeroDr = 1/2*veh.rho*veh.Cd*veh.A*veh.v(i)^2 ;
    % rolling resistance
    veh.RollDr = veh.Cr*abs(-veh.AeroDf+veh.W) ;
    % normal load on driven wheels
    veh.Wd = (veh.factorDrive*veh.W+(-veh.AeroDf * veh.factorAero))/veh.drivenWheels ;
    % drag acceleration
    veh.axDrag = (veh.AeroDr+veh.RollDr)/veh.M ;
    % maximum lat acc available from tyres
    veh.ayMax = 1/veh.M*(veh.mu*(veh.W-veh.AeroDf)) ;
    % max long acc available from tyres
    veh.axTyreMaxAcc = 1/veh.M*veh.mu*veh.Wd*veh.drivenWheels ;
    % max long acc available from tyres
    veh.axTyreMaxDec = -1/veh.M*veh.mu*(veh.W-veh.AeroDf) ;
    % getting power limit from engine
    veh.axPowerLimit = 1/veh.M*(interp1(veh.vehicleSpeed,veh.relativePower*veh.fxEngine,veh.v(i))) ;
    veh.axPowerLimit = veh.axPowerLimit*ones(veh.N,1) ;
    % lat acc vector
    veh.ay = veh.ayMax*cosd(linspace(0,180,veh.N))' ;
    % long acc vector
    veh.axTyreAcc = veh.axTyreMaxAcc*sqrt(1-(veh.ay/veh.ayMax).^2) ; % friction ellipse
    veh.axAcc = min(veh.axTyreAcc,veh.axPowerLimit)+veh.axDrag ; % limiting by engine power
    veh.axDec = veh.axTyreMaxDec*sqrt(1-(veh.ay/veh.ayMax).^2)+veh.axDrag ; % friction ellipse
    % saving GGV map
    veh.GGV(i,:,1) = [veh.axAcc',veh.axDec(2:end)'] ;
    veh.GGV(i,:,2) = [veh.ay',flipud(veh.ay(2:end))'] ;
    veh.GGV(i,:,3) = veh.v(i)*ones(1,2*veh.N-1) ;
end
% HUD
disp('GGV map generated successfully.')

%% Saving vehicle

% saving vehicle model for later use
save(vehname+".mat")

%% Plot

% figure
set(0,'units','pixels') ;
SS = get(0,'screensize') ;
H = 900-90 ;
W = 900 ;
Xpos = floor((SS(3)-W)/2) ;
Ypos = floor((SS(4)-H)/2) ;
f = figure('Name','Vehicle Model','Position',[Xpos,Ypos,W,H]) ;
sgtitle(veh.name)

% rows and columns
rows = 4 ;
cols = 2 ;

% engine curves
subplot(rows,cols,1)
hold on
title('Engine Curve')
xlabel('Engine Speed [rpm]')
yyaxis left
plot(veh.enSpeedCurve,veh.relativePower*veh.enTorqueCurve)
ylabel('Engine Torque [Nm]')
grid on
xlim([veh.enSpeedCurve(1),veh.enSpeedCurve(end)])
yyaxis right
plot(veh.enSpeedCurve,veh.relativePower*veh.enPowerCurve/745.7)
ylabel('Engine Power [Hp]')

% gearing
subplot(rows,cols,3)
hold on
title('Gearing')
xlabel('Speed [m/s]')
yyaxis left
plot(veh.vehicleSpeed,veh.engineSpeed)
ylabel('Engine Speed [rpm]')
grid on
xlim([veh.vehicleSpeed(1),veh.vehicleSpeed(end)])
yyaxis right
plot(veh.vehicleSpeed,veh.gear)
ylabel('Gear [-]')
ylim([veh.gear(1)-1,veh.gear(end)+1])

% traction model
subplot(rows,cols,[5,7])
hold on
title('Traction Model')
plot(veh.vehicleSpeed,veh.relativePower*veh.fxEngine,'k','LineWidth',4)
plot(veh.vehicleSpeed,min([veh.relativePower*veh.fxEngine';veh.fxTyre']),'r','LineWidth',2)
plot(veh.vehicleSpeed,-veh.fxAero)
plot(veh.vehicleSpeed,-veh.fxRoll)
plot(veh.vehicleSpeed,veh.fxTyre)
for i=1:veh.nog
    plot(veh.vehicleSpeed(2:end),veh.fx(:,i),'k--')
end
grid on
xlabel('Speed [m/s]')
ylabel('Force [N]')
xlim([veh.vehicleSpeed(1),veh.vehicleSpeed(end)])
legend({'Engine tractive force','Final tractive force','Aero drag','Rolling resistance','Max tyre tractive force','Engine tractive force per gear'},'Location','southoutside')

% ggv map
subplot(rows,cols,[2,4,6,8])
hold on
title('GGV Map')
surf(veh.GGV(:,:,2),veh.GGV(:,:,1),veh.GGV(:,:,3))
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