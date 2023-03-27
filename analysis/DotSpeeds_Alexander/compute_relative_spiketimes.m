function relspiketimes = compute_relative_spiketimes(spiketimes,eventtimes,maxduration)
%compute_relative_spiketimes Gives spiketimes relative to last trial start to duration
%
%  RELSPIKETIMES = compute_relative_spiketimes(SPIKETIMES,EVENTTIMES,[MAXDURATION])
%
%      SPIKETIMES is a vector with all spike times.
%      EVENTTIMES is a vector with all event start times.
%      MAXDURATION is the time to collect spikes after each event. If not
%      given, then MAXDURATION will be the minimum inter-event interval.
%
% 2021-2023, Alexander Heimel

if nargin<3 || isempty(maxduration)
    eventtimes = sort(eventtimes);
    maxduration = min(diff(eventtimes));
end

n_trialstarts = length(eventtimes);
spikespertrial = cell(n_trialstarts,1);
for i = 1:n_trialstarts
    start = eventtimes(i);
    stop = start + maxduration;
    spikespertrial{i} = spiketimes(spiketimes < stop & spiketimes > start) - start;
end
relspiketimes = vertcat(spikespertrial{:});
relspiketimes = sort(relspiketimes);
