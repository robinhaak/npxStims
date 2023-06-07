function detrended_count = detrend_count( spiketimes )
%DETREND_COUNT computes detrended cumulative spike count 
%
% detrended_count = detrend_count( spiketimes )
%
% 2023, Alexander Heimel

detrended_count = (1:length(spiketimes))' - ...
    (spiketimes - spiketimes(1)) * length(spiketimes) / (spiketimes(end)-spiketimes(1));

