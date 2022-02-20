function [YMD] = QueensRacingYMD(filename, vx, steerMax, bodySlipMax, steerResolution, bodySlipResolution)
%% Code to generate a vehicle model for the Queen's Racing LapSim
%
% Initially written by Maurice Nayman
%
% YMD generation code for sensitivity analyses
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

%Derived variables
a = (1 - cgx) * L;
b = cgx * L;
cgMomentArm = cgh - ((RRCH - FRCH)/L * a + FRCH);

tireData = load(tireFile);

%% HUD

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
disp('YMD generation started.')

%% Generating Yaw Moment Diagram
[YMD] = extractYMD(vx, steerMax, bodySlipMax, steerResolution, bodySlipResolution);

toc
disp('YMD generated successfully.')

end