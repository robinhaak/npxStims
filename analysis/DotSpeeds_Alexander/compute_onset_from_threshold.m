function [onset_time, bootstrapped_error] = compute_onset_from_threshold( spiketimes, eventtimes, prestim_duration, binwidth, threshold_sigma, verbose)
%compute_onset_by_change_from_spontaneous Returns the time that the rate or response reaches half its maximum
%
%  [ONSET_TIME, BOOTSTRAPPED_ERROR] = compute_onset_from_threshold( SPIKETIMES, [EVENTTIMES=0], PRESTIM_DURATION, [BINWIDTH=0.1], [THRESHOLD_SIGMA=3], [VERBOSE=false])
%
%     SPIKETIMES is a vector with spike times relative to stimulus onset
%     EVENTTIMES is the time at start of the spike window. if EVENTTIMES is a
%     vector than the spiketimes are assumed to be from multiple events, 
%     each starting at an element of EVENTTIMES.
%     For the onset, the spiketimes will be then be computed relative to
%     the last event time before a spike.
%     PRESTIM_DURATION is time after start of event that stimulus has not
%     started. if PRESTIM_DURATION is 0, then no period of spontaneous firing
%     is assumed, and not spontaneous firing rate is subtracted for the
%     rate.
%     BINWIDTH, in seconds, is size of spike bins to compute rate.
%     THRESHOLD_SIGMA is the threshold in number of standard deviation of
%     rate in spontaneous period. 
%     STANDARD_ERROR is the standard deviation of the onset times computed
%     by bootstrapping the events. If not required, then do not ask for it, 
%     as then no bootstrapping needs to be done.
%
% 2023, Alexander Heimel

logmsg('DEPRECATED: USE COMPUTE_ONSET_LATENCY INSTEAD');


if nargin<2 || isempty(eventtimes)
    eventtimes = 0;
end
if nargin<3 || isempty(prestim_duration)
    logmsg('It is necessary to provide a prestimulus duration period to compute the spontaneous rate.');
    return
end
if nargin<4 || isempty(binwidth)
    binwidth = 0.1; % s
end
if nargin<5 || isempty(threshold_sigma)
    threshold_sigma = 3; % x std of rate during spontaneous period
end
if nargin<6 || isempty(verbose)
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

bin_edges = 0:binwidth:(max(relspiketimes)+binwidth);

if isempty(bin_edges)
    keyboard
end

[psth_count,bin_edges] = histcounts(relspiketimes,bin_edges);
psth_rate = psth_count / binwidth / num_events;

bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2;

ind_spont = find(bin_edges<=prestim_duration);
ind_spont(end) = []; 

rate_spont = mean(psth_rate(ind_spont));
rate_spont_std = std(psth_rate(ind_spont));
threshold = rate_spont + threshold_sigma * rate_spont_std;

ind_onset = find(psth_rate>threshold,1);
onset_time = bin_centers(ind_onset);

if nargout>1 % compute_standard_error 
    num_bootstraps = 30;
    bootstrapped_onsettimes = zeros(num_bootstraps,1);
    for i_bootstrap = 1:num_bootstraps
        bootstrap = randi(length(eventtimes),size(eventtimes));
        bootstrapped_onsettimes(i_bootstrap) = compute_onset_from_threshold( spiketimes, eventtimes(bootstrap), prestim_duration, binwidth, threshold_sigma, false);
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

