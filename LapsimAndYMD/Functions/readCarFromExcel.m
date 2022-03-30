function [veh] = readCarFromExcel(filename)
%% Generate car model from excel file and load in tire data
veh.carInfo = readInfo(filename,'Info') ;
veh.engineInfo = readTorqueCurve(filename, 'Torque Curve');

veh.name = table2array(veh.carInfo(1,2)) ;
veh.type = table2array(veh.carInfo(2,2)) ;
% index
i = 3 ;
% mass
veh.M = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; % [kg]
veh.cgx = str2double(table2array(veh.carInfo(i,2)))/100 ; i = i+1 ; % [-]
veh.cgh = str2double(table2array(veh.carInfo(i,2)))/1000 ; i = i+1 ; % [m]
% wheelbase and tracks
veh.L = str2double(table2array(veh.carInfo(i,2)))/1000 ; i = i+1 ; % [m]
veh.FT = str2double(table2array(veh.carInfo(i,2)))/1000 ; i = i+1 ; % [m]
veh.RT = str2double(table2array(veh.carInfo(i,2)))/1000 ; i = i+1 ; % [m]
% steering
veh.rackRatio = str2double(table2array(veh.carInfo(i,2)))/1000 ; i = i+1 ; % [-]
% aerodynamics
veh.Cl = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; % [-]
veh.Cd = str2double(table2array(veh.carInfo(i,2))) * -1 ; i = i+1 ; % [-]
veh.CoP = str2double(table2array(veh.carInfo(i,2)))/100 ; i = i+1 ; % [-]
veh.A = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; % [m2]
veh.rho = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; % [kg/m3]
% kinematics and stiffnesses
veh.FRCH = str2double(table2array(veh.carInfo(i,2)))/1000 ; i = i+1 ; % [m]
veh.RRCH = str2double(table2array(veh.carInfo(i,2)))/1000 ; i = i+1 ; % [m]
veh.kF = str2double(table2array(veh.carInfo(i,2))) * 180/pi ; i = i+1 ; % [Nm/rad]
veh.kR = str2double(table2array(veh.carInfo(i,2))) * 180/pi ; i = i+1 ; % [Nm/rad]
% tyre information
veh.rG = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; % [-]
veh.fTP = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; % [psi]
veh.fIA = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; % [deg]
veh.rTP = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; % [psi]
veh.rIA = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; % [deg]
veh.tyreRadius = str2double(table2array(veh.carInfo(i,2)))/1000 ; i = i+1 ; % [m]
veh.Cr = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; % [-]
veh.vWeight = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; % [m/s]
veh.tireFile = (table2array(veh.carInfo(i,2))) ; i = i+1 ; % [-]
% engine
veh.relativePower = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; % [-]
veh.nuTherm = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; %[-]
veh.fLHV = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; %[J/kg]
% transmission
veh.driveType = table2array(veh.carInfo(i,2)) ; i = i+1 ; %[-]
veh.shiftTime = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; %[s]
veh.nuPG = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; %[-]
veh.nuFG = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; %[-]
veh.nuGB = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; %[-]
veh.PGR = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; %[-]
veh.FGR = str2double(table2array(veh.carInfo(i,2))) ; i = i+1 ; %[-]
veh.GBGR = str2double(table2array(veh.carInfo(i:end,2))) ;
veh.nog = length(veh.GBGR) ;

%Derived variables
veh.a = (1 - veh.cgx) * veh.L;
veh.b = veh.cgx * veh.L;
veh.cgMomentArm = veh.cgh - ((veh.RRCH - veh.FRCH)/veh.L * veh.a + veh.FRCH);

veh.tireData = load(veh.tireFile);
end