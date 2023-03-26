function spike_times = generate_spiketimes_from_ratefunction(rate_fun, t_start, t_end, verbose)
%GENERATE_SPIKETIMES_FROM_RATEFUNCTION generates spike times from a poisson process with variable rate
%
% spike_times = generate_spiketimes_from_ratefunction(rate_fun, t_start, t_end, [verbose=false])
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
  n_bars = 20;
  histogram(spike_times,n_bars,'Normalization','countdensity');
  t = linspace(t_start,t_end,100);
  plot( t,rate_fun(t),'k-');
  xlabel('Time');
  ylabel('Rate');
  
end
  
return
