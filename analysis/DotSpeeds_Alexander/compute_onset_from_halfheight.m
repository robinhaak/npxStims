function [onset_time, bootstrapped_error] = compute_onset_from_halfheight( spiketimes, eventtimes, binwidth, prestim_duration, verbose)
%compute_onset_from_halfheight Returns the time that the rate or response reaches half its maximum
%
%  [ONSET_TIME, BOOTSTRAPPED_ERROR] = compute_onset_from_halfheight( SPIKETIMES, [EVENTTIMES=0], [BINWIDTH=0.1], [PRESTIM_DURATION=0], [VERBOSE=false])
%
%     SPIKETIMES is a vector with spike times relative to stimulus onset
%     EVENTTIMES is the time at start of the spike window. if EVENTTIMES is a
%     vector than the spiketimes are assumed to be from multiple events, 
%     each starting at an element of EVENTTIMES.
%     For the onset, the spiketimes will be then be computed relative to
%     the last event time before a spike.
%     BINWIDTH, in seconds, is size of spike bins to compute rate.
%     PRESTIM_DURATION is time after start of event that stimulus has not
%     started. if PRESTIM_DURATION is 0, then no period of spontaneous firing
%     is assumed, and not spontaneous firing rate is subtracted for the
%     rate.
%     STANDARD_ERROR is the standard deviation of the onset times computed
%     by bootstrapping the events. If not required, then do not ask for it, 
%     as then no bootstrapping needs to be done.
%
% 2023, Alexander Heimel

logmsg('DEPRECATED: USE COMPUTE_ONSET_LATENCY INSTEAD');


if nargin<2 || isempty(eventtimes)
    eventtimes = 0;
end
if nargin<3 || isempty(binwidth)
    binwidth = 0.1; % s
end
if nargin<4 || isempty(prestim_duration)
    prestim_duration = 0; % s
end
if nargin<5 || isempty(verbose)
    verbose = false;
end

onset_time = NaN;

if isempty(spiketimes)
    return
end

num_events = length(eventtimes);
if num_events>1 
    eventtimes = sort(eventtimes);
    maxduration = min(diff(unique(eventtimes)));
    relspiketimes = compute_relative_spiketimes(spiketimes,eventtimes,maxduration);
end

[psth_count,bin_edges] = histcounts(relspiketimes,'BinWidth',binwidth);
psth_rate = psth_count / binwidth / num_events;

bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2;

if prestim_duration>0
    rate_spont = length(find(relspiketimes<prestim_duration)) / prestim_duration / num_events;
else
    rate_spont = 0;
end

[max_rate,ind_peak] = max(psth_rate);
halfmax =  rate_spont + (max_rate-rate_spont)/2;
threshold = halfmax;

ind_onset = find(psth_rate>threshold,1);
onset_time = bin_centers(ind_onset);

if nargout>1 % compute_standard_error 
    num_bootstraps = 30;
    bootstrapped_onsettimes = zeros(num_bootstraps,1);
    for i_bootstrap = 1:num_bootstraps
        bootstrap = randi(length(eventtimes),size(eventtimes));
        bootstrapped_onsettimes(i_bootstrap) = compute_onset_from_halfheight( spiketimes, eventtimes(bootstrap), binwidth, prestim_duration, false);
    end
    bootstrapped_error = std(bootstrapped_onsettimes);
else
    bootstrapped_error = NaN;
end

if verbose 
   figure
   hold on
   bar(bin_centers,psth_rate);
   plot(xlim,rate_spont * [1 1],'--k');
   plot(onset_time*[1 1],ylim,'-');
   xlabel('Time (s)')
   ylabel('Rate (sp/s)')

   logmsg([ 'Onset_time = ' num2str(onset_time,'%.3f') ' +- ' num2str(bootstrapped_error,'%.3f')]);
end
