function spike_times = generate_spiketimes_from_constant_rate(rate, t_start, t_end)
%GENERATE_SPIKETIMES_FROM_CONSTANT_RATE
%
% generate example spike times for a Poisson neuron with the given firing rate
% in the interval [t_start, t_end]
%
% 2023, Alexander Heimel

% mean inter-spike interval (ISI) in seconds
mean_isi = 1 / rate;

% generate ISIs from an exponential distribution
num_isis = ceil((t_end ) * rate * 2); % generate twice as many ISIs as needed
isis = exprnd(mean_isi, [num_isis, 1]);

% cumulatively sum ISIs to obtain spike times
spike_times = cumsum(isis);

% remove spike times that fall outside the interval [t_start, t_end]
spike_times(spike_times < t_start | spike_times > t_end) = [];

end