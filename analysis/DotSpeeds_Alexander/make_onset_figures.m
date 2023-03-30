%MAKE_ONSET_FIGURES makes figures to show cast onset determination
%
% 2023, Alexander Heimel, Robin Haak

%%
trial_start_time = 0;
pretrial_duration = 5; % s
prestim_duration = 2; % s
stim_duration = 1; % s
poststim_duration = 0; % s
posttrial_duration = 5; % s
num_stimuli = 10;
event_duration = prestim_duration + stim_duration + poststim_duration;
event_times = trial_start_time + pretrial_duration + 0:event_duration:(num_stimuli-1)*event_duration;
trial_end_time = event_times(end) + event_duration + posttrial_duration;

peak_response = 100; % sp/s
response_fun = @(t) + peak_response*thresholdlinear(  -(t-1/2-prestim_duration).^2 + 1/4  );
rate_spont = 1; % sp/s
rate_fun = @(t) rate_spont + sum(response_fun( t - event_times'));

spike_times = generate_spiketimes_from_ratefunction(rate_fun, trial_start_time, trial_end_time,true);

[p_zeta,zeta] = zetatest(spike_times,event_times);
logmsg(['p_zeta = ' num2str(p_zeta) ]);

%% From detrended spikecount
[onset_time, bootstrapped_error] = compute_onset_latency( spike_times, event_times,[],'cusum_minimum',[],true);

%% From halfheight
[onset_time, bootstrapped_error] = compute_onset_latency( spike_times, event_times,prestim_duration,'rate_halfheight',[],true);

%% From threshold
[onset_time, bootstrapped_error] = compute_onset_latency( spike_times, event_times, prestim_duration, 'rate_threshold', [], true);

