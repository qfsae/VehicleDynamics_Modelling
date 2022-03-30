function [YMD] = extractYMD(vx, steer, bodySlip, SAresolution, BSresolution, veh)
    %Function that returns the maximum lateral acceleration that results in
    %zero yaw moment. Body slip and steer angle ranges are from -Max to
    %Max, with each having the same resolution 'N', to create a YMD
    %resolution of N x N.
    
    %Ensure vx is at least 5 m/s so a non-zero value is returned.
    vx = max([5, vx]);
    
    bodySlip = bodySlip * pi / 180;
    steer = steer * pi / 180;
    
    if length(bodySlip) == 1
        BSDom = linspace(-bodySlip, bodySlip, BSresolution);
    else
        BSDom = bodySlip;
    end
    
    if length(steer) == 1
        SADom = linspace(-steer, steer, SAresolution);
    else
        SADom = steer;
    end

    ayVals = zeros(length(BSDom), length(SADom));
    mzVals = zeros(length(BSDom), length(SADom));
    muVals = zeros(length(BSDom), length(SADom));

    for i = 1:length(BSDom)
        for j = 1:length(SADom)
            [ayVals(i,j), mzVals(i,j), muVals(i,j)] = calculateLateralAccelInterp(vx, SADom(j), BSDom(i), veh);
        end
    end
    
    [aMax, idx] = max(ayVals, [], 'all', 'linear');
    mzMax = mzVals(idx);
    muMax = max(abs(muVals(idx)));
    
    YMD.ayVals = ayVals;
    YMD.mzVals = mzVals;
    YMD.muVals = muVals;
    YMD.aMax = aMax;
    YMD.mzMax = mzMax;
    YMD.muMax = muMax;
end