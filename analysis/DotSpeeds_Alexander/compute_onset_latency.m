function [onset_time, bootstrapped_error] = compute_onset_latency( spiketimes, eventtimes, prestim_duration, method, params, verbose)
%compute_onset_by_change_from_spontaneous Returns the time that the rate or response reaches half its maximum
%
%  [ONSET_TIME, BOOTSTRAPPED_ERROR] = compute_onset_from_threshold( SPIKETIMES, [EVENTTIMES=0], 
%               PRESTIM_DURATION, [METHOD='cusum_minimum'], [PARAMS], [VERBOSE=false])
%
%     SPIKETIMES is a vector with spike times relative to stimulus onset
%     EVENTTIMES is the time at start of the spike window. if EVENTTIMES is a
%         vector than the spiketimes are assumed to be from multiple events, 
%         each starting at an element of EVENTTIMES.
%         For the onset, the spiketimes will be then be computed relative to
%         the last event time before a spike.
%     PRESTIM_DURATION is time after start of event that stimulus has not
%         started. if PRESTIM_DURATION is 0, then no period of spontaneous firing
%         is assumed, and not spontaneous firing rate is subtracted for the
%         rate.
%     METHOD can be 'cusum_threshold', 'cusum_minimum', 'rate_threshold',
%         'rate_halfheight' 
%     PARAMS is a struct with method specific fields
%         For rate methods:
%         PARAMETERS.binwidth, in seconds, is size of bins to compute rate.
%         For threshold methods:
%         PARAMETERS.threshold_sigma is the threshold in number of standard deviation 
%         in spontaneous period. 
%
%     ONSET_TIME is the onset latency 
%     STANDARD_ERROR is the standard deviation of the onset times computed
%     by bootstrapping the events. If not required, then do not ask for it, 
%     as then no bootstrapping needs to be done.
%
% 2023, Alexander Heimel

if nargin<2 || isempty(eventtimes)
    eventtimes = 0;
end
if nargin<3 || isempty(prestim_duration)
    switch method
        case {'cusum_threshold','rate_threshold'}
            logmsg('It is necessary to provide a prestimulus duration period to compute the spontaneous rate.');
            return
    end
    prestim_duration = 0;
end
if nargin<4 || isempty(method)
    method = 'cusum_minimum';
end
if nargin<5 || isempty(params)
    switch method
        case 'cusum_minimum'
        case 'cusum_threshold'
            params.threshold_sigma = 3; % x std of rate during spontaneous period
        case 'rate_threshold'
            params.binwidth = 0.1; % s
            params.threshold_sigma = 3; % x std of rate during spontaneous period
        case 'rate_halfheight'
            params.binwidth = 0.1; % s
    end
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

switch method
    case 'rate_threshold' 
        bin_edges = 0:params.binwidth:(max(relspiketimes)+params.binwidth);
        [psth_count,bin_edges] = histcounts(relspiketimes,bin_edges);
        psth_rate = psth_count / params.binwidth / num_events;
        bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2;
        ind_spont = find(bin_edges<=prestim_duration);
        ind_spont(end) = [];
        rate_spont = mean(psth_rate(ind_spont));
        rate_spont_std = std(psth_rate(ind_spont));
        threshold = rate_spont + params.threshold_sigma * rate_spont_std;
        ind_onset = find(psth_rate>threshold,1);
        onset_time = bin_centers(ind_onset);
    case 'rate_halfheight'
        bin_edges = 0:params.binwidth:(max(relspiketimes)+params.binwidth);
        [psth_count,bin_edges] = histcounts(relspiketimes,bin_edges);
        psth_rate = psth_count / params.binwidth / num_events;
        bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2;
        ind_spont = find(bin_edges<=prestim_duration);
        ind_spont(end) = [];
        rate_spont = mean(psth_rate(ind_spont));
        max_rate = max(psth_rate);
        halfmax =  rate_spont + (max_rate-rate_spont)/2;
        threshold = halfmax;
        ind_onset = find(psth_rate>threshold,1);
        onset_time = bin_centers(ind_onset);
    case 'cusum_minimum'
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
end
        
        



if nargout>1 % compute_standard_error 
    num_bootstraps = 30;
    bootstrapped_onsettimes = zeros(num_bootstraps,1);
    for i_bootstrap = 1:num_bootstraps
        bootstrap = randi(length(eventtimes),size(eventtimes));
        bootstrapped_onsettimes(i_bootstrap) = compute_onset_latency( spiketimes, eventtimes(bootstrap), prestim_duration, method, params, false);
    end
    bootstrapped_error = std(bootstrapped_onsettimes);
else
    bootstrapped_error = NaN;
end

if verbose 
   figure
   hold on
   switch method
       case {'rate_threshold','rate_halfheight'}
           bar(bin_centers,psth_rate);
           plot(xlim,rate_spont * [1 1],'--k');
           ylabel('Rate (sp/s)')
       case {'cusum_threshold','cusum_minimum'}
           plot(relspiketimes,detrended_count,'-');
           plot(onset_time,min_val,'o');
           ylabel('Cusum (sp)')
   end
   
   plot(onset_time*[1 1],ylim,'-');
   xlabel('Time (s)')

   logmsg([ 'Onset_time = ' num2str(onset_time,'%.3f') ' +- ' num2str(bootstrapped_error,'%.3f')]);
end


