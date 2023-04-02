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
%         'rate_halfmax', 'poisson'
%     PARAMS is a struct with method specific fields
%         For rate methods:
%         PARAMETERS.bin_width, in seconds, is size of bins to compute rate.
%         For threshold methods:
%         PARAMETERS.threshold_sigma is the threshold in number of standard deviation
%         in spontaneous period.
%
%     ONSET_TIME is the onset latency
%     STANDARD_ERROR is the standard deviation of the onset times computed
%     by bootstrapping the events. If not required, then do not ask for it,
%     as then no bootstrapping needs to be done.
%
%
%   References
%   'poisson': Legendy and Salcman, J Neurophys 1985
%
% 2023, Alexander Heimel

if nargin<2 || isempty(eventtimes)
    eventtimes = 0;
end
if nargin<3 || isempty(prestim_duration)
    switch method
        case {'cusum_ouellette','rate_threshold'}
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
        case 'cusum_ouellette'
            params.bin_width = 0.111; % s,  Ouellete: 0.001 s
            params.threshold_sigma = 3; % Ouellete: 1 s.e.m.
        case 'rate_threshold'
            params.bin_width = 0.1; % s
            params.threshold_sigma = 3; % x std of rate during spontaneous period
        case 'rate_halfmax'
            params.bin_width = 0.1; % s
        case 'berenyi'
            params.bin_width = 0.005; % s, Berenyi: 0.005 s
            params.window_widths = [0.150  0.300]; % s, Berenyi: [0.150 0.200 0.250 0.300] s
            params.stepsize = params.bin_width; % s, Berenyi: stepsize = bin_width
            params.n_offset_sods = [15  30]; % Berenyi:  [15 20 25 30]
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

if isfield(params,'bin_width') % i.e. psth-based
    bin_edges = 0:params.bin_width:(max(relspiketimes)+params.bin_width);
    [psth_count,bin_edges] = histcounts(relspiketimes,bin_edges);
    psth_rate = psth_count / params.bin_width / num_events;
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2;
    ind_spont = find(bin_edges<=prestim_duration);
    ind_spont(end) = [];
    rate_spont = mean(psth_rate(ind_spont));
    rate_spont_std = std(psth_rate(ind_spont));
else
    % add spike time before first spike
    % either at first spikes minus max isi, or at t_start, whichever
    % comes first
    isi = diff(relspiketimes);
    max_isi = max(isi);
    relspiketimes = [relspiketimes(1)-max_isi;relspiketimes(:)];
    if relspiketimes(1)>eventtimes
        relspiketimes(1) = eventtimes;
    end
    isi = [relspiketimes(2) - relspiketimes(1) ; isi];
    detrended_count = detrend_count(relspiketimes);
    
    rate_spont = sum(relspiketimes < prestim_duration) / prestim_duration / num_events;
end

switch method
    case 'rate_threshold'
        threshold = rate_spont + params.threshold_sigma * rate_spont_std;
        ind_onset = find(psth_rate>threshold,1);
        onset_time = bin_centers(ind_onset);
    case 'rate_halfmax'
        threshold =  rate_spont + (max(psth_rate)-rate_spont)/2;
        ind_onset = find(psth_rate>threshold,1);
        onset_time = bin_centers(ind_onset);
    case 'cusum_ouellette'
        cusum = cumsum(psth_rate - rate_spont);
        threshold = params.threshold_sigma * rate_spont_std;  % Ouellette says to use s.e.m. instead
        ind_onset = find(cusum>threshold,1);
        onset_time = bin_centers(ind_onset);
    case 'cusum_minimum'
        [~,ind_max] = max(detrended_count(2:end));
        ind_max = ind_max + 1;
        [min_val,ind_min] = min(detrended_count(1:ind_max));
        onset_time = relspiketimes(ind_min);
    case 'berenyi'
        onset_times = NaN(1,length(params.window_widths)*length(params.n_offset_sods));
        i_onset_times = 1;
        for i_window_width = 1:length(params.window_widths)
            window_width = params.window_widths(i_window_width);
            
            num_bins_in_window = ceil(window_width/params.bin_width);
            rate_window = movmean(psth_rate,num_bins_in_window,'Endpoints','fill');
            [~,ind_max] = max(rate_window);
            if mod(num_bins_in_window,2)==0
                num_bins_before = num_bins_in_window/2;
                num_bins_after = num_bins_in_window/2-1;
            else
                num_bins_before = (num_bins_in_window-1)/2;
                num_bins_after = (num_bins_in_window-1)/2;
            end
            significance = NaN(1,length(psth_rate));
            for i_sample = (num_bins_before+1) : length(psth_rate)-num_bins_after
                [~,significance(i_sample)] = ttest2(...
                    psth_rate( i_sample-num_bins_before : i_sample+num_bins_after),...
                    psth_rate( ind_max-num_bins_before : ind_max+num_bins_after));
            end % i_sample
            for i_n_offset_sod = 1:length(params.n_offset_sods)
                sod = NaN(1,length(psth_rate));
                n_offset_sod = params.n_offset_sods(i_n_offset_sod);
                for i_sample = (n_offset_sod+1) : length(psth_rate)-n_offset_sod
                    sod(i_sample) = abs(significance(i_sample-n_offset_sod)-significance(i_sample)) - ...
                        abs(significance(i_sample+n_offset_sod)-significance(i_sample));
                end % i_sample
                [~,ind_onset] = min(sod);
                onset_times(i_onset_times) = bin_centers(ind_onset);
                i_onset_times = i_onset_times + 1;
            end % i_n_offset_sods
        end % i_window_widths
        onset_time = median(onset_times);
    case 'poisson'
        num_spikes = length(relspiketimes);
        
        % First find offset of response period
        surprise_index = NaN(num_spikes,1);
        ind_first_spike = 1; 
        for i_spike = (ind_first_spike+1):num_spikes
            interval = relspiketimes(i_spike) - relspiketimes(ind_first_spike);
            p = poisscdf(i_spike - ind_first_spike, rate_spont * num_events * interval,'upper');
            surprise_index( i_spike ) = -log(p);
        end
        [~,ind_offset] = max(surprise_index);
        
        % Find onset 
        surprise_index = NaN(num_spikes,1);
        ind_first_spike = 1; 
        for i_spike = (ind_first_spike+1):(ind_offset-1)
            interval = relspiketimes(ind_offset) - relspiketimes(i_spike);
            lambda = rate_spont * num_events * interval;
            k = ind_offset - i_spike;
            p = poisscdf(k, lambda,'upper');
            if p>0
                surprise_index( i_spike ) = -log(p);
            else
                % use poisson tail approximation
                surprise_index( i_spike ) =  (-k+lambda) - k .*log(lambda./k) - log(k+1) + ...
                    log(sqrt(2*pi*k)) + log(k+1-lambda);
            end
        end
        [~,ind_onset] = max(surprise_index);
        onset_time = relspiketimes(ind_onset);
        
        %         figure
        %         hold on
        %         plot(relspiketimes,surprise_index/max(surprise_index),'r');
        %         plot(relspiketimes,linspace(0,1,num_spikes),'k');
        %         plot(relspiketimes(ind_offset)*[1 1],ylim,'r-');
        %         plot(relspiketimes(ind_onset)*[1 1],ylim,'g-');
        
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
    if isfield(params,'bin_width')
        bar(bin_centers,psth_rate);
        plot(xlim,rate_spont * [1 1],'--k');
        ylabel('Rate (sp/s)')
    else
        plot(relspiketimes,detrended_count,'-');
        plot(onset_time*[1 1],ylim,'-k');
        ylabel('Cusum (sp)')
    end
    
    plot(onset_time*[1 1],ylim,'-');
    xlabel('Time (s)')

    logmsg([ 'Onset_time = ' num2str(onset_time,'%.3f') ' +- ' num2str(bootstrapped_error,'%.3f')]);
end


