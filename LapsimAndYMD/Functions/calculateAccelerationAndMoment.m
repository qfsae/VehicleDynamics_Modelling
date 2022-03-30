function [ayCalc, mz, mu] = calculateAccelerationAndMoment(ay, vx, SA, BS, veh)
    wz = ay ./ vx;
    vy = vx .* sin(BS);
    
    v2 = vx .^2;
    
    staticFront = veh.M * 9.81 * veh.cgx / 2 - 0.25 * veh.rho * veh.Cl * veh.A .* v2 * veh.CoP; %Static weight including aero contributions
    staticRear = veh.M * 9.81 * (1 - veh.cgx) / 2 - 0.25 * veh.rho * veh.Cl * veh.A .* v2 * (1 - veh.CoP); %Static weight including aero contributions

    delWF = ay .* veh.M ./ veh.FT .* ((veh.cgMomentArm * veh.kF)/(veh.kF + veh.kR) + veh.b/veh.L * veh.FRCH);
    delWR = ay .* veh.M ./ veh.FT .* ((veh.cgMomentArm * veh.kR)/(veh.kF + veh.kR) + veh.a/veh.L * veh.RRCH);
    
    weightFR = staticFront + delWF;
    weightFL = staticFront - delWF;
    weightRR = staticRear + delWR; 
    weightRL = staticRear - delWR; 
    
    alphaFR = (vy + wz .* veh.a)./(vx + wz .* veh.FT/2) - SA;
    alphaFL = (vy + wz .* veh.a)./(vx - wz .* veh.FT/2) - SA;
    alphaRR = (vy - wz .* veh.b)./(vx + wz .* veh.RT/2);
    alphaRL = (vy - wz .* veh.b)./(vx - wz .* veh.RT/2);
    
    [FY_FR, muFR] = lateralForce(-weightFR, veh.fIA, veh.fTP, alphaFR, veh);
    FY_FR = -FY_FR;
    muFR = -muFR;
    
    [FY_FL, muFL] = lateralForce(-weightFL, veh.fIA, veh.fTP, -alphaFL, veh);
        
    [FY_RR, muRR] = lateralForce(-weightRR, veh.rIA, veh.rTP, alphaRR, veh);
    FY_RR = -FY_RR;
    muRR = -muRR;
    
    [FY_RL, muRL] = lateralForce(-weightRL, veh.rIA, veh.rTP, -alphaRL, veh);
    
    ayCalc = (FY_FR + FY_FL + FY_RR + FY_RL) / veh.M;
    mz = (FY_FR + FY_FL)*veh.a*cos(SA) - (FY_RR + FY_RL)*veh.b;
    mu = (muFR + muFL + muRR + muRL) ./ 4;
end