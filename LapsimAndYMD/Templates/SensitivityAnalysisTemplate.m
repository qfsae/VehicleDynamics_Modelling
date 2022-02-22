%% Script to generate a sensitivity analysis for use in QueensRacingLAP

%% Clearing memory

clear
clc
close all force
diary('off')
fclose('all') ;

addpath 'pathToLapSimFunctions'
%% Initialization of simInput variable

%Excel workbook that will be modified in the sensitivity analysis and read
%into QueensRacingLAP
carFile = 'excelWorkbook.xlsx';
trackFile = 'pathToTrack' ;

% This assumes the track is already generated. Could be replaced with a
% function call of QueensRacingTRACK
tr = load(trackFile);

%Definition of parameters that you wish to be changed, along with their
%cell in the input sheet
CL = [-3.0 -3.5 -4.0 -4.5 -5.0];
CLCell = 'C11';

CoP = [40 42 44 46 48 50 52 54 56 58 60];
CoPCell = 'C13';

laptimes = zeros(length(CL),length(CoP));

%% Loop to generate and execute LapSim models
% This first model acts as the baseline solution that other laptimes will
% be compared to
baseCL = 0.081;
baseCD = 0.427;
baseCoP = 174.1;
name = 'FSAE-V1';
xlswrite(carFile, cellstr(name), 'Info', 'C2');
xlswrite(carFile, baseCL, 'Info', CLCell);
xlswrite(carFile, baseCD, 'Info', 'C12');
xlswrite(carFile, baseCoP, 'Info', CoPCell);

[veh] = QueensRacingVEHICLE(carFile);
[sim] = QueensRacingLAP(veh, tr);

baseLap = sim.laptime.data;

% Nested loops to update the Excel workbook with the desired parameters on
% each loop
for i = 1:length(CL)
    xlswrite(carFile, CL(i), 'Info', CLCell);
    CD = CL(i)/-2.5;
    xlswrite(carFile, CD, 'Info', 'C12');
    
    for j = 1:length(CoP)
        name = strcat('FSAE-V1_CL',num2str(CL(i)),'_CoP', num2str(CoP(j)));
        xlswrite(carFile, cellstr(name), 'Info', 'C2');
        xlswrite(carFile, CoP(j), 'Info', CoPCell);
        
        [veh] = QueensRacingVEHICLE(carFile);
        %Generate .mat file variable that will get used in QueensRacingLAP
        [sim] = QueensRacingLAP(veh, tr);
        
        %Extract laptimes as a percentage
        laptimes(i,j) = sim.laptime.data/baseLap * 100;
    end
end

% Plotting routines.. The figure positioning is set to get the figure
% nicely sized.
set(0,'units','pixels') ;
SS = get(0,'screensize') ;
H = 900-90 ;
W = 1080 ;
Xpos = floor((SS(3)-W)/2) ;
Ypos = floor((SS(4)-H)/2) ;
figure('Name','CL and CoP Sensitivity','Position',[Xpos,Ypos,W,H]) ;

% Line chart plotting
for i =1:length(CL)
    hold on
    plot(CoP,laptimes(i,:), '-*', 'DisplayName', strcat("CL = ", num2str(-CL(i))))
end

legend('Location', 'northeastoutside')
ylabel('Laptime Relative to Baseline [%]')
xlabel('Centre of Pressure [% Front]')
grid on

set(0,'units','pixels') ;
SS = get(0,'screensize') ;
H = 900-90 ;
W = 1080 ;
Xpos = floor((SS(3)-W)/2) ;
Ypos = floor((SS(4)-H)/2) ;
figure('Name','CL and CoP Sensitivity','Position',[Xpos,Ypos,W,H]) ;

% Surface plotting. Make sure when using meshgrid that each axis
% corresponds to the correct column
[X, Y] = meshgrid(CoP, -CL);
surf(X, Y, laptimes)
xlabel('Centre of Pressure [% Front]')
ylabel('CL [-]')
cb = colorbar;
caxis([min(laptimes,[], 'all'), max(laptimes,[], 'all')])
ylabel(cb, 'Laptime Relative to Baseline [%]', 'FontSize',12)
