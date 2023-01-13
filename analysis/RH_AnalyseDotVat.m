
%% set parameters
sParams = struct;
sParams.strStimulus = 'dot_variations';
sParams.strArea = 'superior colliculus optic';
sParams.boolUseOnlyGoodUnits = false;
sParams.dblMaxAbsNonStationarity = 0.25; %indicates how unit's spikes are distbuted wihtin recording (dbl, [0 1])
sParams.dblMaxViolations1ms = 0.25; %<1 means that unit shows a refractory period (the closer to zero, the better)
sParams.dblSecsFromPrevStimOff = 0.1; %s, time to stay clear of off-response for calculation of spontaneous rate

%% load data
strSessionPath = fullfile(record.path,record.project,'Data_analysis',record.dataset,record.subject,record.sessionid);
fprintf('[%s] Loading data from %s...\n',getTime,strSessionPath);
if ~exist(strSessionPath,'dir')
    msg = [strSessionPath ' does not exist']; error(msg);
end

%load acquipix (sAP) output file
strAP = [record.subject '_' record.date '_AP.mat'];
load(fullfile(strSessionPath,strAP),'sAP');

%% get units in area of interest
vecClustInArea = [];
cellClustAreas = {sAP.sCluster(:).Area};
for intClust = 1:numel(sAP.sCluster)
    strClustArea = cellClustAreas{intClust};
    if ~isempty(strClustArea) && contains(strClustArea,sParams.strArea,'IgnoreCase',true)
        vecClustInArea = [vecClustInArea,intClust]; %#ok<AGROW>
    end
end
sAggNeuron = sAP.sCluster(vecClustInArea);
fprintf('[%s] Found %d clusters in %s\n',getTime,numel(sAggNeuron),sParams.strArea);

%% plot some quality metrics, select 'good' units
figure; hold on
scatter([sAggNeuron.NonStationarity],[sAggNeuron.Violations1ms],'k');
scatter([sAggNeuron([sAggNeuron.KilosortGood]).NonStationarity],[sAggNeuron([sAggNeuron.KilosortGood]).Violations1ms],'k','filled'); %filled dots are marked as 'good'by kilosort
xline([-sParams.dblMaxAbsNonStationarity sParams.dblMaxAbsNonStationarity],'r--');
yline(sParams.dblMaxViolations1ms,'b--');
%set(gca,'YScale','log'); %N.B. zero-values disappear
title([record.subject ' - ' record.date]);
xlabel('Non-stationarity');
ylabel('Violations 1ms');
%fixfig;

%select clusters that pass quality control
if sParams.boolUseOnlyGoodUnits
    vecSelNeuron = [sAggNeuron.Violations1ms]<sParams.dblMaxViolations1ms & abs([sAggNeuron.NonStationarity])<sParams.dblMaxAbsNonStationarity;
    sSelNeuron = sAggNeuron(vecSelNeuron);
    fprintf('[%s] %d/%d clusters passed quality control\n',getTime,numel(sSelNeuron),numel(sAggNeuron));
else
    sSelNeuron = sAggNeuron;
end

%add to record
record.sAggNeuron = sAggNeuron; %all units in area
record.sSelNeuron = sSelNeuron; %units that pass QC

%% find stimulus block
indStim = [];
for intBlock = 1:numel(sAP.cellBlock)
    if strcmp(sAP.cellBlock{1,intBlock}.sAllDots.strStimSet,sParams.strStimulus)
        indStim = [indStim,intBlock]; %#ok<AGROW>
    end
end
if numel(indStim)>1 %ask user, if >1 block is available
    fprintf('[%s] Multiple "%s" blocks are available for %s: blocks %s\n',getTime,sParams.strStimulus,record.sessionid,num2str(indStim));
    indBlock = input('Select block to use: ');
else
    indBlock = indStim;
end
fprintf('[%s] Found "%s" stimulus block\n',getTime,sParams.strStimulus);

%add to record
record.sStimBlock = sAP.cellBlock{1,indBlock};
record.sStimuli = record.sStimBlock.sAllDots;
record.intScreenWidth_pix = record.sStimBlock.sStimParams.intScreenWidth_pix;

%% re-calculate stimulus speeds from bounding rects
%set up only for stimuli moving horizontally
fprintf('[%s] Re-computing stimulus speeds from bounding boxes\n',getTime);
for intStim = 1:record.sStimuli.intStimulusConditions
    stimID = record.sStimuli.stimID(intStim);
    vecBoundingRect = record.sStimuli.vecBoundingRect{intStim};
    intSpeed_ppf = abs(mean(diff(vecBoundingRect(1,:)'))); %pix/frame
    if mod(intSpeed_ppf,1) ~= 0 %check for decimals
        error('Check bounding rects, strange values...');
    end
    intSpeed_pps = intSpeed_ppf/record.sStimBlock.dblStimFrameDur; %pix/s
    record.sStimuli.vecSpeed_pix(intStim) = intSpeed_pps;
    record.sStimBlock.vecSpeed_pix(record.sStimBlock.vecStimID==stimID) = intSpeed_pps;
    %also correct speed in deg/s
    dblSpeed_deg = intSpeed_pps/record.sStimBlock.sStimParams.dblPixelsPerDeg; %deg/s
    record.sStimuli.vecSpeed_deg(intStim) = dblSpeed_deg;
    record.sStimBlock.vecSpeed_deg(record.sStimBlock.vecStimID==stimID) = dblSpeed_deg;
end

%% compute measures
vecStimID = record.sStimBlock.vecStimID; %for ease of use
vecStimOn = record.sStimBlock.vecStimOnTime;
vecStimOff = record.sStimBlock.vecStimOffTime;
vecSpeed_pix = record.sStimuli.vecSpeed_pix;

%0= rightwards (starting LEFT), 180= leftwards (starting RIGHT)
indLeft = record.sStimuli.stimID(record.sStimuli.vecDirection==0);
indRight = record.sStimuli.stimID(record.sStimuli.vecDirection==180);

%loop through clusters
measures = struct();
for intNeuron = 1:numel(sSelNeuron)
    fprintf('[%s] Computing cluster %d (%d/%d)\n',getTime,sSelNeuron(intNeuron).Cluster,intNeuron,numel(sSelNeuron));
    measures(intNeuron).intClu = sSelNeuron(intNeuron).Cluster;
    vecSpikeTimes = sSelNeuron(intNeuron).SpikeTimes;
    measures(intNeuron).dblRateSpontaneous = computeRateSpontaneous(vecSpikeTimes,vecStimOn,vecStimOff,sParams);


    vecDuration = zeros(1,length(record.sStimuli.stimID));
    vecNRepeats = zeros(1,length(record.sStimuli.stimID));
    for intStim = 1:record.sStimuli.intStimulusConditions
        stimID = record.sStimuli.stimID(intStim);
        indStims = find(vecStimID == stimID);
        
        vecEventStarts = vecStimOn(indStims);
        vecDuration(intStim) = mean(vecStimOff(indStims)-vecStimOn(indStims));
        vecNRepeats(intStim) = length(indStims);

        [dblZetaP,sZETA,sRate,~] = zetatest(vecSpikeTimes,...
            vecEventStarts,vecDuration(intStim)); % ~ is essential for correct output
        
        measures(intNeuron).dblZetaP(intStim) = dblZetaP;
        measures(intNeuron).cellT{intStim} = sZETA.vecSpikeT;  % spiketimes with extra 0 and end point
        measures(intNeuron).cellSpikeT{intStim} = measures(intNeuron).cellT{intStim}(2:end-1);  % spiketimes
        measures(intNeuron).vecNSpikes(intStim) = length(measures(intNeuron).cellSpikeT{intStim});

        if ~isempty(sRate)
            measures(intNeuron).cellRate{intStim} = sRate.vecRate;
            measures(intNeuron).vecPeakTime(intStim) = sRate.dblPeakTime;
            measures(intNeuron).vecPeakRate(intStim) = sRate.dblPeakRate;
        else
            measures(intNeuron).cellRate{intStim} = [];
            measures(intNeuron).vecPeakTime(intStim) = NaN;
            measures(intNeuron).vecPeakRate(intStim) = NaN;
        end
        measures(intNeuron).vecMeanRate(intStim) = (length(measures(intNeuron).cellSpikeT{intStim}))/vecDuration(intStim)/vecNRepeats(intStim);
    end %intStim
    
   
end %intNeuron

record.measures = measures;

% if verbose
%     RH_ResultsDotSpeeds(record);
% end
% 
% function dblRateSpontaneous = computeRateSpontaneous(vecSpikeTimes,vecStimOnSecs,vecStimOffSecs,sParams)
% %compute unit's spontaneous/baseline rate
% intCount = 0;
% dblPeriod = 0;
% for intTrial=1:length(vecStimOffSecs)-1
%     intCount = intCount + ...
%         length(find(vecSpikeTimes>(vecStimOffSecs(intTrial)+sParams.dblSecsFromPrevStimOff) & ...
%        vecSpikeTimes<vecStimOnSecs(intTrial+1)));
%     dblPeriod = dblPeriod + vecStimOnSecs(intTrial+1) - (vecStimOffSecs(intTrial)+sParams.dblSecsFromPrevStimOff);
% end
% dblRateSpontaneous = intCount / dblPeriod;
% if dblPeriod<5 
%     fprintf('Less than 5s to compute spontaneous rate.\n')
% end
% if intCount<10 
%     fprintf('Less than 10 spikes to compute spontaneous rate.\n')
% end