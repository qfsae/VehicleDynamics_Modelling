function mu = calculateSteadyMu(v, vWeight, steerMax, bodySlipMax, stepSize)
    %Initialization of mu and weights vectors, as well as initial angle vectors
    muVals = zeros(1,length(v));
    weights = vWeight./max(abs(v - vWeight), 1);
    
    bodySlipMin = -bodySlipMax;
    steerMin = -steerMax;
    
    %Modification of speed vector to find if any speeds are duplicated
    v = max(5, v);
    
    %Angle vector generation. Body slip vector is not changed during the code execution
    bsRes = ceil((bodySlipMax - bodySlipMin)/stepSize) + 1;
    bodySlip = linspace(bodySlipMin, bodySlipMax, bsRes);
    saRes = ceil(steerMax/stepSize) + 1;
    steerAngle = linspace(0, steerMax, saRes);
        
    %Initial mu calculation for speed vector
    [muVals(1), jMax] = steadyMu(v(1), steerAngle, bodySlip);
    
    %Modification of the steer angle vector based on difference from bounds
    difference = max(min(abs(steerAngle(jMax) - steerMax),abs(steerAngle(jMax))), 2 * stepSize);
    steerLimMin = max(steerAngle(jMax) - difference, 0);
    steerLimMax = min(steerAngle(jMax) + difference, steerMax);
        
    i = 2;
    while i <= length(v)
        if v(i) == v(i - 1)
            muVals(i) = muVals(i - 1);
        else
            %Generation of new angle limits for each speed and rerunning of steady-state mu calculation
            saRes = ceil((steerLimMax - steerLimMin)/stepSize) + 1;
            steerAngle = linspace(steerLimMin, steerLimMax, saRes);
        
            [muVals(i), jMax] = steadyMu(v(i), steerAngle, bodySlip);
            
            %Check to see if there were no sign changes in any columns, and to rerun the iteration if necessary.
            if jMax == 0
                i = i - 1;
                counter = counter + 1;
                difference = max(min(abs(steerAngle(floor(j/2)) - steerLimMax),abs(steerAngle(floor(j/2)) - steerLimMin)), counter * stepSize);
                steerLimMin = max(steerAngle(floor(j/2)) - difference, 0);
                steerLimMax = min(steerAngle(floor(j/2)) + difference, steerMax);
            else
                counter = 2;                
                difference = max(min(abs(steerAngle(jMax) - steerLimMax),abs(steerAngle(jMax) - steerLimMin)), 2 * stepSize);
                steerLimMin = max(steerAngle(jMax) - difference, 0);
                steerLimMax = min(steerAngle(jMax) + difference, steerMax);
            end
        end
        i = i + 1;
    end
    %Weighted average of friction coefficients
    mu = sum(muVals.*weights)/sum(weights);

function [muSteady, jMax] = steadyMu(v, steer, bodySlip)
    muSteady = 0;
    
    %Extraction of yaw moment diagram for the specified steer and body slip vectors
    [YMD] = extractYMD(v, steer, bodySlip, 0, 0);
        
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