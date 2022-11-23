
function cellColorMaps = RH_ColorMaps
%based on a script by Sylvia Schroeder

vecRed = [1 0 .5];
vecBlue = [0 .5 1];
vecBlack = [1 1 1].*0.5;
vecGrad = linspace(0,1,100)';
vecReds = vecRed.*flip(vecGrad) + [1 1 1].*vecGrad;
vecBlacks = vecBlack.*flip(vecGrad) + [1 1 1].*vecGrad;
vecON = [vecBlacks; flip(vecReds(1:end-1,:),1)];
vecBlues = vecBlue.*flip(vecGrad) + [1 1 1].*vecGrad;
vecOFF = [vecBlacks; flip(vecBlues(1:end-1,:),1)];
cellColorMaps = {vecON, vecOFF};