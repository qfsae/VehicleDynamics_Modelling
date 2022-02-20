function [ayCalc, mz, mu] = calculateAccelerationAndMoment(ay, vx, SA, BS)
    global M
    global FT
    global RT
    global cgMomentArm
    global cgx
    global L
    global a
    global b
    global kF
    global kR
    global FRCH
    global RRCH
    global Cl
    global CoP
    global A
    global rho
    global fTP
    global fIA
    global rTP
    global rIA
    
    wz = ay ./ vx;
    vy = vx .* BS';
    
    v2 = vx .^2 + vy .^2;
    
    staticFront = M * 9.81 * cgx / 2 - 0.25 * rho * Cl * A .* v2 * CoP; %Static weight including aero contributions
    staticRear = M * 9.81 * (1 - cgx) / 2 - 0.25 * rho * Cl * A .* v2 * (1 - CoP); %Static weight including aero contributions

    delWF = ay .* M ./ FT .* ((cgMomentArm * kF)/(kF + kR) + b/L * FRCH);
    delWR = ay .* M ./ FT .* ((cgMomentArm * kR)/(kF + kR) + a/L * RRCH);
    
    weightFR = staticFront + delWF;
    weightFL = staticFront - delWF;
    weightRR = staticRear + delWR; 
    weightRL = staticRear - delWR; 
    
    alphaFR = (vy + wz .* a)./(vx + wz .* FT/2) - SA;
    alphaFL = (vy + wz .* a)./(vx - wz .* FT/2) - SA;
    alphaRR = (vy - wz .* b)./(vx + wz .* RT/2);
    alphaRL = (vy - wz .* b)./(vx - wz .* RT/2);
    
    [FY_FR, muFR] = lateralForce(-weightFR, fIA, fTP, alphaFR);
    FY_FR = -FY_FR;
    muFR = -muFR;
    
    [FY_FL, muFL] = lateralForce(-weightFL, fIA, fTP, -alphaFL);
        
    [FY_RR, muRR] = lateralForce(-weightRR, rIA, rTP, alphaRR);
    FY_RR = -FY_RR;
    muRR = -muRR;
    
    [FY_RL, muRL] = lateralForce(-weightRL, rIA, rTP, -alphaRL);
    
    ayCalc = (FY_FR + FY_FL + FY_RR + FY_RL) / M;
    mz = (FY_FR + FY_FL)*a*cos(SA) - (FY_RR + FY_RL)*b;
    mu = (muFR + muFL + muRR + muRL) ./ 4;
end