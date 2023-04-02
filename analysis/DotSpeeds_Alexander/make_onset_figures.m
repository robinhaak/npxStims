%MAKE_ONSET_FIGURES makes figures to show cast onset determination
%
% 2023, Alexander Heimel, Robin Haak

%%
trial_start_time = 0;
pretrial_duration = 5; % s
prestim_duration = 2; % s
stim_duration = 1; % s
poststim_duration = 1; % s
posttrial_duration = 5; % s
num_stimuli = 10; % 10
peak_response = 50; % sp/s
rate_spont = 1; % sp/s
response_fun = @(t) + peak_response*thresholdlinear(  -(t-1/2-prestim_duration).^2 + 1/4  );

event_duration = prestim_duration + stim_duration + poststim_duration;
event_times = trial_start_time + pretrial_duration + 0:event_duration:(num_stimuli-1)*event_duration;
trial_end_time = event_times(end) + event_duration + posttrial_duration;

rate_fun = @(t) rate_spont + sum(response_fun( t - event_times'));
%%
    spike_times = generate_spiketimes_from_ratefunction(rate_fun, trial_start_time, trial_end_time,false);
    compute_onset_latency( spike_times, event_times,prestim_duration,'poisson',[],true)
return
%%

CUSUM_MINIMUM = 1;
CUSUM_OUELLETTE = 2;
RATE_THRESHOLD = 3;
RATE_HALFMAX = 4;
POISSON = 5;
BERENYI = 6;

methods{CUSUM_MINIMUM} = 'cusum_minimum';
methods{CUSUM_OUELLETTE} = 'cusum_ouellette';
methods{RATE_THRESHOLD} = 'rate_threshold';
methods{RATE_HALFMAX} = 'rate_halfmax';
methods{POISSON} = 'poisson';
methods{BERENYI} = 'berenyi';

%[p_zeta,zeta] = zetatest(spike_times,event_times);
%logmsg(['p_zeta = ' num2str(p_zeta) ]);

num_samples = 100;

onsets = NaN(length(methods),num_samples);
onset_errors = NaN(length(methods),num_samples);

for i_sample = 1:num_samples
    logmsg(['Computing sample ' num2str(i_sample)]);
    spike_times = generate_spiketimes_from_ratefunction(rate_fun, trial_start_time, trial_end_time,false);
    for i_method = [CUSUM_MINIMUM POISSON RATE_THRESHOLD]
            logmsg(['Computing method ' methods{i_method}]);

        [onsets(i_method,i_sample), onset_errors(i_method,i_sample)] = ...
            compute_onset_latency( spike_times, event_times,prestim_duration,methods{i_method},[],false);
    end % i_method
end % i_sample

onset_means = mean(onsets,2);
onset_std = std(onsets,[],2);

for i_method = 1:length(methods)
    if isnan(onset_means(i_method))
        continue
    end
    
    logmsg([ methods{i_method} ': Onset = ' num2str(onset_means(i_method),3) ' +- ' num2str(onset_std(i_method),3) ' s (mean+-std)']);
end % i_method
