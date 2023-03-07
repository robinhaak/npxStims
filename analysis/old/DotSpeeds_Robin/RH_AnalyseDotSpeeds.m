function record = RH_AnalyseDotSpeeds(record,verbose)
%RH_ANALYSEDOTSPEEDS
%analyse responses to dot stimuli moving at different speeds, from left>right & right>left
%
%    for modeling, assumes:
%         tLeft = x / v + w /(2v) + DeltaT
%         tRight = - x / v + w /(2v) + DeltaT
%         where tLeft is response time when stimulus moves from left
%         where tRight is response time when stimulus moves from right
%         x is stimulus position (x=0 is center of the screen, w is width of the screen, positive is to the right)
%
%         from this follows:
%         DeltaT = (tLeft + tRight)/2 - w/(2v)
%         x = 1/2 v (tLeft - tRight)
%
%2022, Robin Haak & Alexander Heimel

if nargin<2 || isempty(verbose)
    verbose = true;
end

%% set parameters
sParams = struct;
sParams.strStimulus = 'dot_speeds';
sParams.strArea = 'superior colliculus';
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
%load acquipix output file
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
%set(gca,'YScale','log'); %zero-values disappear
title([record.subject ' - ' record.date]);
xlabel('Non-stationarity');
ylabel('Violations 1ms');
fixfig;

%select clusters that pass quality control
if sParams.boolUseOnlyGoodUnits
    vecSelNeuron = [sAggNeuron.Violations1ms]<sParams.dblMaxViolations1ms & abs([sAggNeuron.NonStationarity])<sParams.dblMaxAbsNonStationarity;
    sSelNeuron = sAggNeuron(vecSelNeuron);
    fprintf('[%s] %d/%d clusters passed quality control\n',getTime,numel(sSelNeuron),numel(sAggNeuron));
else
    sSelNeuron = sAggNeuron;
end

%add to record
record.sAggNeuron = sAggNeuron;
record.sSelNeuron = sSelNeuron;

%% find stimulus block
vecIsStim = [];
for intBlock = 1:numel(sAP.cellBlock)
    if strcmp(sAP.cellBlock{1,intBlock}.sAllDots.strStimSet,sParams.strStimulus)
        vecIsStim = [vecIsStim,intBlock]; %#ok<AGROW>
    end
end
if numel(vecIsStim)>1 %ask user, if >1 block is available
    fprintf('[%s] Multiple "%s" blocks are available for %s: blocks [%d]\n',getTime,sParams.strStimulus,record.sessionid,vecIsStim);
    indBlock = input('Select block to use: ');
else
    indBlock = vecIsStim;
end
record.sStimuli = sAP.cellBlock{1,indBlock};
fprintf('[%s] Found "%s" stimulus block\n',getTime,sParams.strStimulus)

%% re-calculate stimulus speeds from bounding rects
fprintf('[%s] Re-computing stimulus speeds\n',getTime);
for intStim = 1:record.sStimuli.sAllDots.intStimulusConditions
    vecBoundingRect = record.sStimuli.sAllDots.vecBoundingRect{intStim};
    intSpeed_ppf = abs(mean(diff(vecBoundingRect(2,:)'))); %pix/frame
    if intSpeed_ppf == 0
        intSpeed_ppf = abs(mean(diff(vecBoundingRect(1,:)'))); %pix/frame, alternative direction
    end
    if mod(intSpeed_ppf,1) ~= 0
        error('Check bounding rects, strange values...');
    end
    intSpeed_pps = intSpeed_ppf/record.sStimuli.dblStimFrameDur; %pix/s
    record.sStimuli.sAllDots.vecSpeed_pix(intStim) = intSpeed_pps;
end

for intTrial = 1:record.sStimuli.intTrialNum
    intStimID = record.sStimuli.vecStimID(intTrial);
    record.sStimuli.vecSpeed_pix(intTrial) = record.sStimuli.sAllDots.vecSpeed_pix(record.sStimuli.sAllDots.stimID==intStimID);
end

%% compute measures
record.intScreenWidth_pix = record.sStimuli.sStimParams.intScreenWidth_pix;

%for ease of use
vecStimID = record.sStimuli.vecStimID;
vecStimOn = record.sStimuli.vecStimOnTime;
vecStimOff = record.sStimuli.vecStimOffTime;
vecSpeed_pix = record.sStimuli.vecSpeed_pix;

%DeltaT calculation depends on stimulus specifics, check how they are organized
%0= rightwards, 180= leftwards
indRight = record.sStimuli.sAllDots.stimID(record.sStimuli.sAllDots.vecDirection==0);
indLeft = record.sStimuli.sAllDots.stimID(record.sStimuli.sAllDots.vecDirection==180);
if sum(diff(abs(record.sStimuli.sAllDots.vecSpeed_pix(indLeft)))<0)>0 || sum(diff(abs(record.sStimuli.sAllDots.vecSpeed_pix(indRight)))<0)>0
    error('Order of vecSpeed_pix in sAllDots is weird! Please check');
end

%loop through clusters
measures = struct();
for intNeuron = 1:numel(sSelNeuron)
    fprintf('[%s] Computing cluster %d (%d/%d)\n',getTime,sSelNeuron(intNeuron).Cluster,intNeuron,numel(sSelNeuron));
    measures(intNeuron).intClu = sSelNeuron(intNeuron).Cluster;
    vecSpikeTimes = sSelNeuron(intNeuron).SpikeTimes;
    measures(intNeuron).dblRateSpontaneous = computeRateSpontaneous(vecSpikeTimes,vecStimOn,vecStimOff,sParams);

    vecDuration = zeros(1,length(record.sStimuli.sAllDots.stimID));
    vecNRepeats = zeros(1,length(record.sStimuli.sAllDots.stimID));
    for intStim = 1:record.sStimuli.sAllDots.intStimulusConditions
        stimID = record.sStimuli.sAllDots.stimID(intStim);
        indStims = find(vecStimID == stimID);

        vecEventStarts = vecStimOn(indStims);
        vecDuration(intStim) = mean( vecStimOff(indStims) - vecStimOn(indStims));
        vecNRepeats(intStim) = length(indStims);

        [dblZetaP,sZETA,sRate,~] = zetatest(vecSpikeTimes,...
            vecEventStarts,vecDuration(intStim)); % ~ is essential for correct output

        measures(intNeuron).dblZetaP(intStim) = dblZetaP;
        measures(intNeuron).cellT{intStim} = sZETA.vecSpikeT;  % spiketimes with extra 0 and end point
        measures(intNeuron).cellSpikeT{intStim} = measures(intNeuron).cellT{intStim}(2:end-1);  % spiketimes
        measures(intNeuron).vecNSpikes(intStim) = length( measures(intNeuron).cellSpikeT{intStim} );

        if ~isempty(sRate)
            measures(intNeuron).cellRate{intStim} = sRate.vecRate;
            measures(intNeuron).vecPeakTime(intStim) = sRate.dblPeakTime;
            measures(intNeuron).vecPeakRate(intStim) = sRate.dblPeakRate;
        else
            measures(intNeuron).cellRate{intStim} = [];
            measures(intNeuron).vecPeakTime(intStim) = NaN;
            measures(intNeuron).vecPeakRate(intStim) = NaN;
        end
        measures(intNeuron).vecMeanRate(intStim) = (length(measures(intNeuron).cellSpikeT))/vecDuration(intStim)/vecNRepeats(intStim);
    end % stim i

    % DeltaT>0 response to stimulus in the past.
    % DeltaT<0 response to stimulus in the future.
    %
    % Calculate DeltaT based on PeakTime
    % Define
    %     TLeft = peakTime left stimulus
    %     TRight = peakTime right stimulus
    %
    %     xRF = x0 + v tPass
    %
    %         tLeft = (x - x0Left ) / v + DeltaT
    %         tRight = - (x - x0Right ) / v + DeltaT
    %         where tLeft is response time when stimulus moves from left
    %         where tRight is response time when stimulus moves from right
    %         x is stimulus position (x=0 is center of the screen, w is width of the screen, positive is to the right)
    %         for leftward stimulus, x0 = w/2 - stimulus_radius
    %
    %         from this follows:
    %         DeltaT = (tLeft + tRight)/2 - w/(2v)
    %         x = 1/2 v (tLeft - tRight)

    vecTLeft = measures(intNeuron).vecPeakTime(indLeft);
    vecTRight = measures(intNeuron).vecPeakTime(indRight);

    dblStimulusRadius_pix = abs(record.sStimuli.sAllDots.vecBoundingRect{1,1}(1,1)) / 2;
    dblStartLeft_pix = -record.intScreenWidth_pix / 2 - dblStimulusRadius_pix  ;
    dblStartRight_pix = record.intScreenWidth_pix / 2 + dblStimulusRadius_pix ;

    measures(intNeuron).vecPeakXRF_pix = vecSpeed_pix(indLeft) .* (vecTLeft - vecTRight) / 2;
    measures(intNeuron).vecPeakDeltaT = (vecTLeft + vecTRight)/2  - (dblStartRight_pix - dblStartLeft_pix) / 2; %



    % Calculate DeltaT based on all spikes
    % Assume rate(t) = rateSpontaneous + s(dir) g(t - DeltaT - tPass)
    %     with \int_-\inf^+\inf g(t) dt = 1,
    %     support only between (-DeltaT-tPass) and (T-DeltaT-tPass)
    %     and  \int_-\inf^+\inf t g(t) dt = 0 (satisfied if g is symmetric around 0)
    %     and tPass = (xRF - x0)/v, where xRF is RF position, x0 is start
    %     of stimulus and v is speed.
    % then s(dir) = nSpikes - T rateSpontaneous,
    %     with T time interval over which nSpikes is computed
    % and an integration with t'=t+DeltaT-tPass shows
    %     Expectation( \sum_i t(i) ) = 1/2 T^2 rateSpontaneous + s(dir)(DeltaT + tpass)
    %
    % Now define:
    %   SLeft = nSpikesLeft - T rateSpontaneous
    %   SRight = nSpikesRight - T rateSpontaneous
    %   MLeft = E(\sum_i t(i))_left -1/2 T^2 rateSpontaneous
    %   MRight = E(\sum_i t(i))_right -1/2 T^2 rateSpontaneous
    % Then:
    %   DeltaT = -h/v + (MLeft SRight + MRight SLeft) / (2 SLeft SRight)
    %   xRF = v (MLeft SRight - MRight SLeft) / (2 SLeft SRight)
    %       with h = screenWidth/2 and xRF relative to center
    %
    %   when only spikes are occuring at tLeft and tRight, this simplifies
    %   to earlier expression.

    dblH = record.intScreenWidth_pix/2;

    vecSLeft = measures(intNeuron).vecNSpikes(indLeft) - measures(intNeuron).dblRateSpontaneous  * vecDuration(indLeft) .* vecNRepeats(indLeft);
    vecSRight = measures(intNeuron).vecNSpikes(indRight) - measures(intNeuron).dblRateSpontaneous * vecDuration(indRight) .* vecNRepeats(indRight);
    vecMLeft =   cellfun(@sum,measures(intNeuron).cellSpikeT(indLeft))  - 1/2 *  measures(intNeuron).dblRateSpontaneous * vecDuration(indLeft).^2 .* vecNRepeats(indLeft);
    vecMRight = cellfun(@sum,measures(intNeuron).cellSpikeT(indRight)) - 1/2 *  measures(intNeuron).dblRateSpontaneous * vecDuration(indRight).^2 .* vecNRepeats(indRight);

    measures(intNeuron).vecMeanDeltaT = - dblH ./ vecSpeed_pix(indLeft)  + ...
        (vecMLeft .* vecSRight + vecMRight .* vecSLeft) ./ (2 * vecSLeft .* vecSRight )  ;

    measures(intNeuron).vecMeanXRF_pix = vecSpeed_pix(indLeft) .* ...
        (vecMLeft .* vecSRight - vecMRight .* vecSLeft) ./ (2 * vecSLeft .* vecSRight ) ;

    indNaN = abs(measures(intNeuron).vecMeanXRF_pix) > dblH;
    measures(intNeuron).vecMeanXRF_pix(indNaN) =  NaN;
    measures(intNeuron).vecMeanDeltaT(indNaN) = NaN;
end 

record.measures = measures;

if verbose
    RH_ResultsDotSpeeds(record);
end

function dblRateSpontaneous = computeRateSpontaneous(vecSpikeTimes,vecStimOnSecs,vecStimOffSecs,sParams)
%compute unit's spontaneous/baseline rate
intCount = 0;
dblPeriod = 0;
for i=1:length(vecStimOffSecs)-1
    intCount = intCount + ...
        length(find(vecSpikeTimes>(vecStimOffSecs(i)+sParams.dblSecsFromPrevStimOff) & ...
        vecSpikeTimes<vecStimOnSecs(i+1)));
    dblPeriod = dblPeriod + vecStimOnSecs(i+1) - (vecStimOffSecs(i)+sParams.dblSecsFromPrevStimOff);
end
dblRateSpontaneous = intCount / dblPeriod;
if dblPeriod<5
    fprintf('Less than 5s to compute spontaneous rate!\n')
end
if intCount<10
    fprintf('Less than 10 spikes to compute spontaneous rate!\n')
end
