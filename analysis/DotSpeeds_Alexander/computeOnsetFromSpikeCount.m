function dblOnsetTime = computeOnsetFromSpikeCount( vecSpikeTimes, verbose)
%computeOnsetFromSpikeCount gives onset based on minimum of corrected cumultive spike count
%
% dblOnsetTime = computeOnsetFromSpikeCount( vecSpikeTimes, [verbose=false])
%
% 2023, Alexander Heimel

disp('DEPRECATED: use COMPUTE_ONSET_FROM_SPIKECOUNT instead.');

if nargin<2 || isempty(verbose)
    verbose = false;
end

dblOnsetTime = NaN;

if isempty(vecSpikeTimes)
    return
end

vecCorCurSpike = corCumFun(vecSpikeTimes);

[~,indMax] = max(vecCorCurSpike);

[dblMin,indMin] = min(vecCorCurSpike(1:indMax));
dblOnsetTime = vecSpikeTimes(indMin);

if verbose
    figure
    hold on
    plot(vecSpikeTimes,vecCorCurSpike,'-');
    plot(dblOnsetTime,dblMin,'o');
end
