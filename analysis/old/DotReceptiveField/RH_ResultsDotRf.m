
function RH_ResultsDotRF(record)

global measures %#ok<GVMIS> 
evalin('base','global measures');

set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesFontSize', 10);

% indRightW = record.sStimuli.sAllDots.stimID(record.sStimuli.sAllDots.vecDirection==0); %0= rightwards, 180= leftwards
% indLeftW = record.sStimuli.sAllDots.stimID(record.sStimuli.sAllDots.vecDirection==180);


%% azimuth-trajectories
indDownW  = record.sStimuli.sAllDots.stimID(record.sStimuli.sAllDots.vecDirection==90); %90= downwards, 270= upwards 
indUpW = record.sStimuli.sAllDots.stimID(record.sStimuli.sAllDots.vecDirection==270);

%calculate x-coordinate of each trajectory
indDownUpW = [indDownW indUpW];
vecX = zeros(1,length(indDownUpW));
for intTraj = 1:length(indDownUpW)
    vecX(intTraj) = record.sStimuli.sAllDots.vecBoundingRect{1,indDownUpW(intTraj)}(1,1) + ...
        record.sStimuli.sAllDots.vecBoundingRect{1,indDownUpW(intTraj)}(3,1)-record.sStimuli.sAllDots.vecBoundingRect{1,indDownUpW(intTraj)}(1,1);
end
%%
for intNeuron = 1:65
%mean rate
vecMeanRate = measures(intNeuron).dblZetaP(indDownUpW);
vecMeanMean = (vecMeanRate(1:16)+vecMeanRate(17:32))/2;
figure;
plot(vecX(1:16),vecMeanMean)
end






%% azimuth-trajectories
indRightW = record.sStimuli.sAllDots.stimID(record.sStimuli.sAllDots.vecDirection==0); %0= rightwards, 180= leftwards
indLeftW = record.sStimuli.sAllDots.stimID(record.sStimuli.sAllDots.vecDirection==180);

%calculate x-coordinate of each trajectory
indRightLeftW = [indRightW indLeftW];
vecX = zeros(1,length(indRightLeftW));
for intTraj = 1:length(indRightLeftW)
    vecX(intTraj) = record.sStimuli.sAllDots.vecBoundingRect{1,indRightLeftW(intTraj)}(2,1) + ...
        record.sStimuli.sAllDots.vecBoundingRect{1,indRightLeftW(intTraj)}(4,1)-record.sStimuli.sAllDots.vecBoundingRect{1,indRightLeftW(intTraj)}(2,1);
end

%%
figure;hold on
for intNeuron = 1:3
%mean rate
vecMeanRate = measures(intNeuron).dblZetaP(indRightLeftW);
vecMeanMean = (vecMeanRate(1:9)+vecMeanRate(10:18))/2
plot(vecX(1:9),vecMeanMean)
end



