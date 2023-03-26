%MAKE_ONSET_FIGURES makes figures to show cast onset determination
%
% 2023, Alexander Heimel, Robin Haak

%%
trial_start_time = 0;
pretrial_duration = 5; % s
prestimulus_duration = 2; % s
stimulus_duration = 1; % s
poststimulus_duration = 0; % s
posttrial_duration = 5; % s
num_stimuli = 10;
event_duration = prestimulus_duration + stimulus_duration + poststimulus_duration;
event_times = trial_start_time + pretrial_duration + 0:event_duration:(num_stimuli-1)*event_duration;
trial_end_time = event_times(end) + event_duration + posttrial_duration;

response_fun = @(t) + 50*thresholdlinear(  -(t-1/2-prestimulus_duration).^2 + 1/4  );
rate_spont = 1; % sp/s
rate_fun = @(t) rate_spont + sum(response_fun( t - event_times'));

spiketimes = generate_spiketimes_from_ratefunction(rate_fun, trial_start_time, trial_end_time,true);

[p_zeta,zeta] = zetatest(spiketimes,event_times);
onset_time = compute_onset_from_spikecount( zeta.vecSpikeT, 0, true);
logmsg(['p_zeta = ' num2str(p_zeta) ', onset_time = ' num2str(onset_time) ]);

