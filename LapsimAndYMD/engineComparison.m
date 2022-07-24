%ADD PATH
%Add the template and function folders to the path

addpath('Q:\Formula\VehicleDynamics_Modelling\LapsimAndYMD\Templates')
addpath('Q:\Formula\VehicleDynamics_Modelling\LapsimAndYMD\Functions')

%FILES
%Input vehicle and track filenames
%Each file has motor weight and torque curves changed. 
vehicleFile = 'FSAE_QueensRacingVEHICLE_AMKDT7-75.xlsx';
% vehicleFile = 'FSAE_QueensRacingVEHICLE_EMRAX208.xlsx';
% vehicleFile = 'FSAE_QueensRacingVEHICLE_EMRAX228.xlsx';
% vehicleFile = 'FSAE_QueensRacingVEHICLE_EMRAX268.xlsx';
% vehicleFile = 'FSAE_QueensRacingVEHICLE_YASAP400.xlsx';

trackFile = 'FSAE_Michigan2019_AutoX.xlsx';

%VEHICLE
%Run the QueensRacingVehicle file to generate the vehicle data
vehicle = QueensRacingVEHICLE(vehicleFile);

%TRACK
%Run the QueensRacingTRACK to generate the track data
track = QueensRacingTRACK(trackFile, 'shape data');

%LAP
%Run QueensRacingLAP to run a given vehicle through the track
sim = QueensRacingLAP(vehicle, track);

%MISC
%Useful information we want to extract
%wheel torque, engine torque, engine power, engine speed, gear

