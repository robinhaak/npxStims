function dblOnsetTime = computeOnsetFromSpikeCount( vecSpikeTimes, verbose)
if nargin<2 || isempty(verbose)
    verbose = false;
end

dblOnsetTime = NaN;

if isempty(vecSpikeTimes)
    return
end

vecCorCurSpike = (1:length(vecSpikeTimes))' - (vecSpikeTimes - vecSpikeTimes(1)) * length(vecSpikeTimes) / (vecSpikeTimes(end)-vecSpikeTimes(1));

[dblMax,indMax] = max(vecCorCurSpike);


[dblMin,indMin] = min(vecCorCurSpike(1:indMax));
dblOnsetTime = vecSpikeTimes(indMin);

if verbose
    figure
    hold on
    plot(vecSpikeTimes,vecCorCurSpike,'-');
    plot(dblOnsetTime,dblMin,'o');
end
