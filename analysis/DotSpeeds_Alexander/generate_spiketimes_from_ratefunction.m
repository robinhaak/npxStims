function spike_times = generate_spiketimes_from_ratefunction(rate_fun, t_start, t_end, verbose)
%GENERATE_SPIKETIMES_FROM_RATEFUNCTION generates spike times from a poisson process with variable rate
%
% SPIKE_TIMES = generate_spiketimes_from_ratefunction(RATE_FUN, T_START, T_END, [VERBOSE=false])
%     RATE_FUN is a handle to a function of time that gives the rate
%     T_START is start of interval to generate spike times for
%     T_END is end of interval to generate spike times for
%     if VERBOSE is true, a figure with the detrended count is shown
%
% 2023, Alexander Heimel 

if nargin<4 || isempty(verbose)
    verbose = false;
end

num_steps = 10000;
dt = (t_end - t_start) / num_steps;
t = t_start:dt:t_end;

cdf = dt*cumtrapz(rate_fun(t));
mean_num_spikes = cdf(end);
num_spikes = poissrnd(mean_num_spikes);

cdf = cdf/mean_num_spikes;

% make sure that cdf(t) can be inverted
[cdf,ind]= unique(cdf);
t = t(ind);

u = rand(num_spikes, 1);
spike_times = interp1(cdf, t, u, 'linear', 'extrap');
spike_times = sort(spike_times);

if verbose
  figure;
  subplot(2,1,1);
  hold on
  plot(t, cdf*mean_num_spikes);
  plot(spike_times,1:num_spikes,'-b');
  xlabel('Time');
  ylabel('Cum. count');
  
  subplot(2,1,2);
  hold on
  binwidth = 0.1; % s
  histogram(spike_times,'BinWidth',binwidth,'Normalization','countdensity');
  t = linspace(t_start,t_end,100);
  plot( t,rate_fun(t),'k-');
  xlabel('Time');
  ylabel('Rate');
  
end
  
return
