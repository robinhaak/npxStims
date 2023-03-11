function vecCorCumSpike = corCumFun( vecSpikes )
%coCumFun computes detrended cumulative spike count 
%
% vecCorCumSpike = corCumFun( vecSpikes )
%
% 2023, Alexander Heimel

vecCorCumSpike = (1:length(vecSpikes))' - (vecSpikes - vecSpikes(1)) * length(vecSpikes) / (vecSpikes(end)-vecSpikes(1));
