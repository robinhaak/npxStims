function onsettime = compute_onset_from_spikecount( spiketimes, t_start, verbose)
%compute_onset_from_spikecount gives onset based on minimum of corrected cumultive spike count
%
% ONSETTIME = compute_onset_from_spikecount(SPIKETIMES, [T_START=0], [VERBOSE=false])
%
%     SPIKETIMES is a vector with spike times
%     T_START is the time at start of the spike window
%     if VERBOSE
%
% 2023, Alexander Heimel

if nargin<2 || isempty(t_start)
    t_start = 0;
end
if nargin<3 || isempty(verbose)
    verbose = false;
end

onsettime = NaN;

if isempty(spiketimes)
    return
end

spiketimes = sort(spiketimes); % just for safety

% add spike time before first spike
% either at first spikes minus max isi, or at t_start, whichever
% comes first
isi = diff(spiketimes);
max_isi = max(isi);
spiketimes = [spiketimes(1)-max_isi;spiketimes(:)];
if spiketimes(1)>t_start
    spiketimes(1) = t_start;
end
    
detrended_count = detrend_count(spiketimes);

[~,ind_max] = max(detrended_count(2:end));
ind_max = ind_max + 1;

[min_val,ind_min] = min(detrended_count(1:ind_max));
onsettime = spiketimes(ind_min);

if verbose
    figure
    hold on
    plot(spiketimes,detrended_count,'-');
    plot(onsettime,min_val,'o');
end
