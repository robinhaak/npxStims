
function [vecSDF,vecT] = getSpikeDensityFunction(vecSpikes,vecTrialStarts,dblDur,dblTimeStep,dblSigma)
%GETSPIKEDENSITYFUNCTION
%2023, Robin Haak

vecT = 0:dblTimeStep:dblDur; 

%sort spikes
intTrials = numel(vecTrialStarts);
matSDF = zeros(intTrials,length(vecT));
for intTrial=1:intTrials
    %get spikes
    vecTheseSpikes = vecSpikes((vecSpikes >= vecTrialStarts(intTrial)) & vecSpikes < (vecTrialStarts(intTrial)+ dblDur));
    vecTheseSpikes = vecTheseSpikes - vecTrialStarts(intTrial);

    matGauss = zeros(numel(vecTheseSpikes),length(vecT));
    for intSpike = 1:numel(vecTheseSpikes)
        %center gaussian at spike time
        dblMu = vecTheseSpikes(intSpike);

        %compute gaussian
        p1 = -0.5 * ((vecT - dblMu) / dblSigma) .^2;
        p2 = (dblSigma * sqrt(2 * pi));
        matGauss(intSpike,:) = exp(p1 ./ p2);

    end
        
        %sum over all distributions to get spike density function
        matSDF(intTrial,:) = sum(matGauss,1);
end

%average over time
vecSDF = mean(matSDF);





