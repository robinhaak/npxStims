function vecStimPos_pix = timeToStimulusCenterPosition(vecTime,sStimuli,ind)
%timeToStimulusCenterPosition computes stimulus position at the spike time
%
%  vecStimPos_pix = spikeTimeToStimulusCenterPosition(vecSpikeTimes,sStimuli,ind)
%
%    ind should be integer or a vector with the same number of elements of
%    vecTime
%
% 2023, Alexander Heimel

if length(ind) == 1
    vecStimPos_pix = ...
        sStimuli.vecStimStartX_pix(ind) + cos(sStimuli.vecDirection(ind)/180*pi) * vecTime * sStimuli.vecSpeed_pix(ind);
elseif length(ind) == length(vecTime)
    vecStimPos_pix = ...
        sStimuli.vecStimStartX_pix(ind) + cos(sStimuli.vecDirection(ind)/180*pi) .* vecTime .* sStimuli.vecSpeed_pix(ind);
else
    logmsg('vecTime and ind have different number of elements and ind is not a single integer.')
end
