function dblRateSpontaneous = computeRateSpontaneous( vecSpikeTimes, vecStimOnTime,vecStimOffTime, sParams)
%computeRateSpontaneous computes spontaneous rate
%
% dblRateSpontaneous = computeRateSpontaneous( vecSpikeTimes, vecStimOnTime,vecStimOffTime, sParams)
%
% 2023, Alexander Heimel

intCount = 0;
dblPeriod = 0;
for i=1:length(vecStimOffTime)-1
    intCount = intCount + ...
        length(find(vecSpikeTimes>(vecStimOffTime(i)+sParams.separationFromPrevStimOff) & ...
        vecSpikeTimes<vecStimOnTime(i+1)));
    dblPeriod = dblPeriod + vecStimOnTime(i+1) - (vecStimOffTime(i)+sParams.separationFromPrevStimOff);
end
dblRateSpontaneous = intCount / dblPeriod;
if dblPeriod<5
    disp('Less than 5s to compute spontaneous rate.')
end
if intCount<10
    disp('Less than 10 spikes to compute spontaneous rate.')
end
end