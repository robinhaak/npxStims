function [onset_time,bootstrapped_error] = compute_onset_from_spikecount( spiketimes, eventtimes, maxduration, verbose)
%compute_onset_from_spikecount Gives onset based on minimum of corrected cumultive spike count
%
% [ONSET_TIME, STANDARD_ERROR] = compute_onset_from_spikecount(SPIKETIMES, [EVENTTIMES=0], [MAXDURATION], [VERBOSE=false])
%
%     SPIKETIMES is a vector with spike times relative to stimulus onset
%     EVENTTIMES is the time at start of the spike window. if EVENTTIMES is a
%     vector than the spiketimes are assumed to be from multiple events, 
%     each starting at an element of EVENTTIMES.
%     For the onset, the spiketimes will be then be computed relative to
%     the last event time before a spike.
%     MAXDURATION is maximum period to look at after event onset. If not 
%     given, then MAXDURATION is minimum inter-event interval.
%     if VERBOSE is true, a figure with the detrended count and onset is
%     shown.
%     STANDARD_ERROR is the standard deviation of the onset times computed
%     by bootstrapping the events. If not required, then do not ask for it, 
%     as then no bootstrapping needs to be done.
%
% 2023, Alexander Heimel

logmsg('DEPRECATED: USE COMPUTE_ONSET_LATENCY INSTEAD');

if nargin<2 || isempty(eventtimes)
    eventtimes = 0;
end
if nargin<3 || isempty(maxduration)
    eventtimes = sort(eventtimes);
    maxduration = min(diff(unique(eventtimes)));
end
if nargin<4 || isempty(verbose)
    verbose = false;
end

onset_time = NaN;

if isempty(spiketimes)
    return
end

if length(eventtimes)>1 
    relspiketimes = compute_relative_spiketimes(spiketimes,eventtimes,maxduration);
else
    relspiketimes = sort(spiketimes); % just for safety
end

% add spike time before first spike
% either at first spikes minus max isi, or at t_start, whichever
% comes first
isi = diff(relspiketimes);
max_isi = max(isi);
relspiketimes = [relspiketimes(1)-max_isi;relspiketimes(:)];
if relspiketimes(1)>eventtimes
    relspiketimes(1) = eventtimes;
end
    
detrended_count = detrend_count(relspiketimes);

[~,ind_max] = max(detrended_count(2:end));
ind_max = ind_max + 1;

[min_val,ind_min] = min(detrended_count(1:ind_max));
onset_time = relspiketimes(ind_min);

if nargout>1 % compute_standard_error 
    num_bootstraps = 30;
    bootstrapped_onsettimes = zeros(num_bootstraps,1);
    for i_bootstrap = 1:num_bootstraps
        bootstrap = randi(length(eventtimes),size(eventtimes));
        bootstrapped_onsettimes(i_bootstrap) = compute_onset_from_spikecount( spiketimes, eventtimes(bootstrap), maxduration, false);
    end
    bootstrapped_error = std(bootstrapped_onsettimes);
else
    bootstrapped_error = NaN;
end

if verbose
    figure
    hold on
    plot(relspiketimes,detrended_count,'-');
    plot(onset_time,min_val,'o');

    logmsg([ 'Onset_time = ' num2str(onset_time,'%.3f') ' +- ' num2str(bootstrapped_error,'%.3f')]);
end
