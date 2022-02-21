function [ay, mz, mu] = calculateLateralAccelInterp(vx, steerAngle, bodySlip)
    aMax = 40;
    aTMP = linspace(-aMax, aMax, 400);
    aOut = calculateAccelerationAndMoment(aTMP, vx, steerAngle, bodySlip) - aTMP;
    
    signChanges = find(sign(aOut(1:end-1)) ~= sign(aOut(2:end)));
    
    if isempty(signChanges) == 1
        %This will give an offset to ensure the zero is found correctly in
        %case the domain doesn't cross zero
        [~, idx] = min(abs(aOut));
        ay = aOut(idx) + aTMP(idx);
        [~, mz, mu] = calculateAccelerationAndMoment(ay, vx, steerAngle, bodySlip);
    else
        [~, idx] = max(abs(aTMP(signChanges)));
    
        f = @(x) pchip(aTMP, aOut, x);
        [xZero, ay] = fzero(f, aTMP(signChanges(idx)));
    
        ay = xZero + ay;
        [~, mz, mu] = calculateAccelerationAndMoment(ay, vx, steerAngle, bodySlip);
    end
end