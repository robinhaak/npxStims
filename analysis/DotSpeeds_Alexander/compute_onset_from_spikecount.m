function onsettime = compute_onset_from_spikecount( spiketimes, verbose)
%computeOnsetFromSpikeCount gives onset based on minimum of corrected cumultive spike count
%
% onsettime = compute_onset_from_spikecount( spiketimes, [verbose=false])
%
% 2023, Alexander Heimel

if nargin<2 || isempty(verbose)
    verbose = false;
end

onsettime = NaN;

if isempty(spiketimes)
    return
end

detrended_count = detrend_count(spiketimes);

[~,ind_max] = max(detrended_count);

[min_val,ind_min] = min(detrended_count(1:ind_max));
onsettime = spiketimes(ind_min);

if verbose
    figure
    hold on
    plot(spiketimes,detrended_count,'-');
    plot(onsettime,min_val,'o');
end
