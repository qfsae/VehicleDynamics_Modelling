function mu = calculateSteadyMu(veh)
    %Initialization of mu and weights vectors, as well as initial angle vectors
    muVals = zeros(1,length(veh.vMuCalc));
    weights = veh.vWeight./max(abs(veh.vMuCalc - veh.vWeight), 1);
    
    bodySlipMin = -veh.bodySlipMax;
    steerMin = -veh.steerMax;
    
    %Modification of speed vector to find if any speeds are duplicated
    veh.vMuCalc = max(5, veh.vMuCalc);
    
    %Angle vector generation. Body slip vector is not changed during the code execution
    bsRes = ceil((veh.bodySlipMax - bodySlipMin)/veh.stepSize) + 1;
    bodySlip = linspace(bodySlipMin, veh.bodySlipMax, bsRes);
    saRes = ceil(veh.steerMax/veh.stepSize) + 1;
    steerAngle = linspace(0, veh.steerMax, saRes);
        
    %Initial mu calculation for speed vector
    [muVals(1), jMax] = steadyMu(veh.vMuCalc(1), steerAngle, bodySlip, veh);
    
    %Modification of the steer angle vector based on difference from bounds
    difference = max(min(abs(steerAngle(jMax) - veh.steerMax),abs(steerAngle(jMax))), 2 * veh.stepSize);
    steerLimMin = max(steerAngle(jMax) - difference, 0);
    steerLimMax = min(steerAngle(jMax) + difference, veh.steerMax);
        
    i = 2;
    while i <= length(veh.vMuCalc)
        if veh.vMuCalc(i) == veh.vMuCalc(i - 1)
            muVals(i) = muVals(i - 1);
        else
            %Generation of new angle limits for each speed and rerunning of steady-state mu calculation
            saRes = ceil((steerLimMax - steerLimMin)/veh.stepSize) + 1;
            steerAngle = linspace(steerLimMin, steerLimMax, saRes);
        
            [muVals(i), jMax] = steadyMu(veh.vMuCalc(i), steerAngle, bodySlip, veh);
            
            %Check to see if there were no sign changes in any columns, and to rerun the iteration if necessary.
            if jMax == 0
                i = i - 1;
                counter = counter + 1;
                difference = max(veh.steerMax/2, counter * veh.stepSize);
                steerLimMin = max(veh.steerMax/2 - difference, 0);
                steerLimMax = min(veh.steerMax/2 + difference, veh.steerMax);
            else
                counter = 2;                
                difference = max(min(abs(steerAngle(jMax) - steerLimMax),abs(steerAngle(jMax) - steerLimMin)), 2 * veh.stepSize);
                steerLimMin = max(steerAngle(jMax) - difference, 0);
                steerLimMax = min(steerAngle(jMax) + difference, veh.steerMax);
            end
        end
        i = i + 1;
    end
    %Weighted average of friction coefficients
    mu = sum(muVals.*weights)/sum(weights);

function [muSteady, jMax] = steadyMu(v, steer, bodySlip, veh)
    muSteady = 0;
    
    %Extraction of yaw moment diagram for the specified steer and body slip vectors
    [YMD] = extractYMD(v, steer, bodySlip, 0, 0, veh);
        
    jMax = 0;
    
    %Iteration over each column, will return jMax = 0 if there are no columns with sign changes. 
    for j = 1:saRes
        f = @(x) pchip(YMD.ayVals(:,j), YMD.mzVals(:,j), x);

        signChanges = find(sign(YMD.mzVals(1:end-1, j)) ~= sign(YMD.mzVals(2:end, j)));
                
        if isempty(signChanges)
                continue
        end
            
        diffs = diff(YMD.mzVals(:,j)');
        [~, idx] = min(diffs(signChanges));
        idx = signChanges(idx);
        ayInt = fzero(f, YMD.ayVals(idx, j));
        muInt = pchip(YMD.ayVals(:,j), YMD.muVals(:, j), ayInt);
    
        muSteady = max(abs(muInt), muSteady);
        if muSteady == abs(muInt)
            jMax = j;
        end
    end
end % End statement for steadyMu

end % End statement for calculateSteadyMu