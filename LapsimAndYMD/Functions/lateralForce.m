function [Fy, mu] = lateralForce(Fz,Ia,P,slipAngle, veh)
    %ensure no positive values values
    Fz = min(0,Fz);
    
    B = zeros(1,length(Fz));
    C = zeros(1,length(Fz));
    D = zeros(1,length(Fz));
    Ep = zeros(1,length(Fz));
    En = zeros(1,length(Fz));
    x0 = zeros(1,length(Fz));
    y0 = zeros(1,length(Fz));
    
    %do some weird stuff
    [x,y,z] = meshgrid(veh.tireData.IA_domain, veh.tireData.FZ_domain, veh.tireData.P_domain);
    
    %interpolate to get coefficients for magic formula
    B = interp3(x,y,z,veh.tireData.processed_data(:,:,:,1),Ia,Fz,P,'makima')';
    C = interp3(x,y,z,veh.tireData.processed_data(:,:,:,2),Ia,Fz,P,'makima')';
    D = interp3(x,y,z,veh.tireData.processed_data(:,:,:,3),Ia,Fz,P,'makima')';
    Ep = interp3(x,y,z,veh.tireData.processed_data(:,:,:,4),Ia,Fz,P,'makima')';
    En = interp3(x,y,z,veh.tireData.processed_data(:,:,:,5),Ia,Fz,P,'makima')';
    x0 = interp3(x,y,z,veh.tireData.processed_data(:,:,:,6),Ia,Fz,P,'makima')';
    y0 = interp3(x,y,z,veh.tireData.processed_data(:,:,:,7),Ia,Fz,P,'makima')';
    
    slipAngle = slipAngle - x0;

    if slipAngle > 0
        mu = D.*sin(C.*atan(B.*slipAngle-Ep.*(B.*slipAngle-atan(B.*slipAngle))));
    else
        mu = D.*sin(C.*atan(B.*slipAngle-En.*(B.*slipAngle-atan(B.*slipAngle)))); 
    end

    mu = mu + y0;
    
    %solve for lateral load
    Fy = veh.rG.*Fz.*mu;
    mu = veh.rG.*mu;
   
end