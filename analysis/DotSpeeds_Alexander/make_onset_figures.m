%MAKE_ONSET_FIGURES makes figures to show cast onset determination
%
% 2023, Alexander Heimel, Robin Haak

%%


rate_fun = @(t) 100 + 100*sin(2*pi*t); % Hz

rate_fun = @(t) 1000 * (1 + 1*sign(t-3) + (1 + sign(t-5)) + (-1 - sign(t-7))+ (-1 - sign(t-9)) );
start_time = 3; % s
end_time = 10; %

num_samples = 1000;
num_spikes = zeros(num_samples,1);
for i_sample = 1:num_samples
    spiketimes = generate_spiketimes_from_ratefunction(rate_fun, start_time, end_time);
    num_spikes(i_sample) = length(spiketimes);
end % i_sample

disp([ num2str(mean(num_spikes)) ' +- ' num2str(sem(num_spikes))]);

spiketimes = generate_spiketimes_from_ratefunction(rate_fun, start_time, end_time,true);

 
 