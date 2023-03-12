function record = analyse_movingdots(record,db,verbose)
%ANALYSE_MOVINGDOTS analysis moving dots stimuli
%
%  RECORD = ANALYSE_MOVINGDOTS( RECORD, DB, VERBOSE )
%
%  Wrapper analysis responses to dots stimuli,
%  real analysis in done in stimulus specific scripts.
%
% 2022-2023, Robin Haak, Alexander Heimel

if nargin<2 || isempty(db)
    db = [];
end
if nargin<3 || isempty(verbose)
    verbose = true;
end
sParams = RH_defaultParameters();

% if verbose
%    record = RH_analyse_synchronization(record,verbose);
% end

%% Load data
% check which channels or clusters to analyse
vecClustersToAnalyze = get_channels2analyze( record );

strSessionPath = fullfile(sParams.strOutputPath,record.project,'Data_collection',record.dataset,record.subject,record.sessionid);
if ~exist(strSessionPath,'dir')
    logmsg([strSessionPath ' does not exist']);
    return
end

strLog = fullfile(strSessionPath,[record.sessionid '.mat']);
sVars = whos('-file',strLog);
if contains('sSynthData',{sVars.name}) %contains('sCluster',{sVars.name})
    % single unit data
    load(strLog,'sSynthData'); %rm sSynthData later
    
    cellExpType = cellfun(@(x) x.structEP.strExpType, sSynthData.cellStim, 'UniformOutput', false);
    ind = find(contains(cellExpType,record.stimulus));
    if isempty(ind)
        logmsg(['Could not find stimulus ' record.stimulus ' in '  strLog]);
        return
    end
    if length(ind)>1
        logmsg(['More than one stimulus ' record.stimulus ' in '  strLog]);
        logmsg('Using first stimulus');
        ind = 1;
    end
    
    structEP = sSynthData.cellStim{ind}.structEP;
    sCluster = sSynthData.sCluster;
    
    if isempty(vecClustersToAnalyze)
        vecClustersToAnalyze = unique([sCluster.IdxClust]);
    end
    intNumClusters = length(vecClustersToAnalyze);
    vecSpikeTimesOfCluster = cell(intNumClusters,1);
    for c = 1:intNumClusters
        indClus = find([sCluster.IdxClust]==vecClustersToAnalyze(c),1);
        if ~isempty(indClus)
            vecSpikeTimesOfCluster{c} = sCluster(indClus).SpikeTimes;
        end
    end
    
    vecStimOnTime = structEP.vecStimOnTime;
    vecStimOffTime = structEP.vecStimOffTime;
    %     sMeta = sSynthData.sMetaNI;
    %     dblStartT = str2double(sMeta.firstSample) / str2double(sMeta.niSampRate);
    %     vecStimOnTime = structEP.ActOnNI - dblStartT;
    %     vecStimOffTime = structEP.ActOffNI - dblStartT;
else
    % unsorted channels
    load(strLog,'structEP');
    load(fullfile(strSessionPath,'spikes.mat'),'vecSpikeCh','vecSpikeSecs');
    
    if isempty(vecClustersToAnalyze)
        vecClustersToAnalyze = unique(vecSpikeCh);
    end
    intNumClusters = length(vecClustersToAnalyze);
    vecSpikeTimesOfCluster = cell(intNumClusters,1);
    for c = 1:intNumClusters
        vecSpikeTimesOfCluster{c} = vecSpikeSecs(vecSpikeCh==vecClustersToAnalyze(c));
    end
    
    %to be added: re-alignment of times based on diode data
    vecStimOnTime = structEP.ActOnNI; %NI stream times
    vecStimOffTime = structEP.ActOffNI;
end

%% Compute measures
record.sStimuli = structEP.sAllDots;

vecStimStartX_pix = cellfun( @(x) mean([x(1,1),x(3,1)]), record.sStimuli.vecBoundingRect) - structEP.sStimParams.intScreenWidth_pix/2 ; % center of screen is 0
record.sStimuli.vecStimStartX_pix = vecStimStartX_pix;

logmsg('Assuming Left and right stimuli are in the same order for delta T calculation');
warning('off','zetatest:InsufficientSamples');

measures = struct([]);
indMeasures = 0;
dblIntertrialInterval =  structEP.sStimParams.dblSecsPreBlank + structEP.sStimParams.vecSecsPostBlank(1);
for c = 1:length(vecClustersToAnalyze) % over clusters or channels
    logmsg(['Computing cluster/channel index ' num2str(vecClustersToAnalyze(c)) ...
        ' (' num2str(c) ' of '  num2str(intNumClusters) '), ' num2str(length(vecSpikeTimesOfCluster{c})) ' spikes'])
    if length(vecSpikeTimesOfCluster{c}) < 1
        continue
    end
    
    indMeasures = indMeasures + 1;
    measure.intIndex = vecClustersToAnalyze(c);
    
    if exist('sCluster','var')
        ind = find([sCluster.IdxClust]==measure.intIndex);
        measure.dblDepth_um = sCluster(ind).Depth; %#ok<FNDSB>
    else
        measure.dblDepth_um = 1000 - c; %before: measure.depth
    end
    
    measure.dblRateSpontaneous = computeRateSpontaneous( vecSpikeTimesOfCluster{c}, vecStimOnTime,vecStimOffTime, sParams);
    measure.vecDuration = zeros(1,length(structEP.sAllDots.stimID));
    measure.vecNRepeats = zeros(1,length(structEP.sAllDots.stimID));
    for i = 1:length(structEP.sAllDots.stimID)
        stimID = structEP.sAllDots.stimID(i);
        indStims = find(structEP.vecStimID == stimID);
        
        vecEventStarts = vecStimOnTime(indStims);
        measure.vecDuration(i) = mean(vecStimOffTime(indStims) - vecStimOnTime(indStims));
        measure.vecNRepeats(i) = length(indStims);
        
        [ dblZetaP,sZETA] = zetatest(vecSpikeTimesOfCluster{c},...
            vecEventStarts,measure.vecDuration(i)+dblIntertrialInterval); % ~ is essential for correct output
        
        measure.vecZetaP(i) = dblZetaP;

        measure.cellSpikeTimes{i} =sZETA.vecSpikeT(2:end-1);  % spiketimes, ZETA has padded extra points at beginning and end
        measure.vecNSpikes(i) = length(measure.cellSpikeTimes{i} );
        
        [vecSpikeCount,vecEdges] = histcounts(measure.cellSpikeTimes{i},'BinWidth',sParams.dblBinWidth);
        vecCenters = (vecEdges(1:end-1)+vecEdges(2:end))/2;
        [dblPeak,indPeak] = max(vecSpikeCount);
        measure.cellSpikeCounts{i} = vecSpikeCount;
        measure.cellEdges{i} = vecEdges;
        measure.vecPeakTime(i) = vecCenters(indPeak);
        measure.vecPeakRate(i) = dblPeak / sParams.dblBinWidth / measure.vecNRepeats(i);
        measure.vecOnsetTime(i) = NaN;
        measure.vecOnsetTime(i) = computeOnsetFromSpikeCount( measure.cellSpikeTimes{i} );
        measure.vecMeanRate(i) = (length(measure.cellSpikeTimes{i}))/measure.vecDuration(i)/measure.vecNRepeats(i);
    end % stim i

    intNumStimuli = length(measure.vecZetaP);
    measure.vecResponsive = measure.vecZetaP<min(sParams.dblThresholdResponsiveZetaP,1/intNumStimuli);
    
    switch record.sStimuli.strStimSet
        case 'dot_speeds'
            measure = compute_dot_speeds_measures( measure, record );
        case 'flashing_dots'
            measure = compute_flashing_dots_measures( measure, record, structEP.sStimParams );
        case 'dot_diffhist'
            measure = compute_dot_diffhist_measures( measure, record, db );
        otherwise
            logmsg(['Analysis not implemented of ' record.sStimuli.strStimSet])
    end
    
    if sum(measure.vecResponsive)>2
        measure.boolResponsive = true;
    else
        measure.boolResponsive = false;
    end
    if ~measure.boolResponsive
        % reduce database size
        logmsg('Cluster not responsive. Removing PSTHs');
        measure.cellSpikeTimes = {};
        measure.cellSpikeCounts = {};
        measure.cellEdges = {};
    end
    
    if indMeasures == 1
        measures = measure;
    else
        measures(indMeasures) = measure;
    end
end % cluster c

if isfield(record,'measures') && ~isempty(record.measures) && length(measures)>=1
    % Merge with old measures
    sMeasuresOld = record.measures;
    vecIndicesOld = [sMeasuresOld.intIndex];
    vecIndicesNew = [measures.intIndex];
    vecOverlap = intersect(vecIndicesOld,vecIndicesNew);
    sMeasuresOld(ismember(vecIndicesOld,vecOverlap)) = [];
    try
        measures = [sMeasuresOld measures];
        vecIndices = [measures.intIndex];
        [~,ind] = sort(vecIndices);
        measures = measures(ind);
    catch me
        logmsg(['Cannot merge with old measures. Removing old measures. ' me.message]) 
    end
end

record.measures = measures;

% % Add dots results
% h_db = get_fighandle('Neuropixels database*');
% sUserData = get(h_db,'userdata');
% db = sUserData.db;
% record = analyse_add_patches_to_dots( record, db, true);

logmsg(['Analyzed ' recordfilter(record)]);
end

function measure = compute_dot_speeds_measures( measure, record )

sParams = RH_defaultParameters();

indLeft = find([record.sStimuli.vecDirection]==0);
indRight = find([record.sStimuli.vecDirection]==180);

vecSpeed_pix = record.sStimuli.vecSpeed_pix;


% DeltaT>0 response to stimulus in the past.
% DeltaT<0 response to stimulus in the future.
% Calculate DeltaT based on PeakTime
% Define
%     x0 as the starting position of the center of the stimulus
%     v the speed of the stimulus.
%     Then the center stimulus position at time t is x = x0 + v t
%     If xRF is the receptive field (RF) center, then the time tPass at
%     which the center passes the RF center is,
%         tPass = (xRF - x0)/v
%
% Compute
%     tLeft = peakTime left stimulus
%     tRight = peakTime right stimulus
% Assume tLeft = tPassLeft + DeltaT, and tRight = tPassRight + DeltaT
% Then
%         tLeft = (xRF - x0Left) / v + DeltaT
%         tRight = - (xRF - X0Right) / v + DeltaT
%
% from this follows:
%         DeltaT = (tLeft + tRight)/2 + (x0Left - x0Right)/(2v)
%         xRF = 1/2 v (tLeft - tRight)

vecTLeft = measure.vecPeakTime(indLeft);
vecTRight = measure.vecPeakTime(indRight);

vecStimStartX_pix = record.sStimuli.vecStimStartX_pix;

measure.vecPeakXRF_pix = vecSpeed_pix(indLeft) .* (vecTLeft - vecTRight) / 2;
measure.vecPeakDeltaT = (vecTLeft + vecTRight)/2  + (vecStimStartX_pix(indLeft)-vecStimStartX_pix(indRight)) ./ vecSpeed_pix(indLeft) / 2;

% Calculate DeltaT based on all spikes
% Assume rate(t) = rateSpontaneous + s(dir) g(t - DeltaT - tPass)
%     with \int_-\inf^+\inf g(t) dt = 1,
%     support only between (-DeltaT-tPass) and (T-DeltaT-tPass)
%     and  \int_-\inf^+\inf t g(t) dt = 0 (satisfied if g is symmetric around 0)
%     and tPass = (xRF - x0)/v, where xRF is RF position, x0 is start
%     of stimulus and v is speed.
% then from Expectation( nSpikes ) = \int_0^T rate(t) follows:
%     s(dir) = nSpikes - T rateSpontaneous,
%     with T time interval over which nSpikes is computed
% and an integration of t * rate(t) over 0 to T, by substituting t'=t+DeltaT-tPass shows
%     \int_0^T dt t rate(t) = 1/2 T^2 rateSpontaneous + s(dir)(DeltaT + tpass)
% but the Expectation( spike time of spikes between 0 to T ) * Expectation(nSpikes) = \int_0^T dt t rate(t)
% thus
%    E(spike time of spikes between 0 to T) * E(nSpikes) = 1/2 T^2 rateSpontaneous + s(dir)(DeltaT + tpass)
%
%    Note that E(spike time of spikes between 0 to T) * E(nSpikes) = E(\sum_i t(i))
% Now define:
%   SLeft = nSpikesLeft - T rateSpontaneous
%   SRight = nSpikesRight - T rateSpontaneous
%   MLeft = E(\sum_i t(i))_left -1/2 T^2 rateSpontaneous
%   MRight = E(\sum_i t(i))_right -1/2 T^2 rateSpontaneous
% Then:
%   MLeft / SLeft   = DeltaT + (xRF - x0Left)/v
%   MRight / SRight = DeltaT - (xRF - x0Right)/v
%
% =>
%   DeltaT = 1/2 (MLeft/SLeft + MRight/SRight + x0Left/v - x0Right/v)
%   xRF = 1/2 v (MLeft/SLeft + MRight/SRight) + 1/2 (x0Left + x0Right)
%
%   when only spikes are occuring at tLeft and tRight, this simplifies
%   to earlier expression.

%dblH = structEP.sStimParams.intScreenWidth_pix/2;
vecSLeft = measure.vecNSpikes(indLeft) -measure.dblRateSpontaneous  * measure.vecDuration(indLeft) .* measure.vecNRepeats(indLeft);
vecSRight = measure.vecNSpikes(indRight) -measure.dblRateSpontaneous * measure.vecDuration(indRight) .* measure.vecNRepeats(indRight);
vecMLeft =   cellfun(@sum,measure.cellSpikeTimes(indLeft))  - 1/2 * measure.dblRateSpontaneous * measure.vecDuration(indLeft).^2 .* measure.vecNRepeats(indLeft);
vecMRight = cellfun(@sum,measure.cellSpikeTimes(indRight)) - 1/2 * measure.dblRateSpontaneous * measure.vecDuration(indRight).^2 .* measure.vecNRepeats(indRight);

measure.vecMeanDeltaT = 1/2 * (vecMLeft ./ vecSLeft + vecMRight ./ vecSRight + ...
    vecStimStartX_pix(indLeft)./ vecSpeed_pix(indLeft) - vecStimStartX_pix(indRight)./ vecSpeed_pix(indRight) )   ;

measure.vecMeanXRF_pix = vecSpeed_pix(indLeft) .* ...
    (vecMLeft ./ vecSLeft - vecMRight ./ vecSRight) / 2  + ...
    1/2*(vecStimStartX_pix(indLeft) + vecStimStartX_pix(indRight));

indNaN = (measure.vecMeanXRF_pix > max(vecStimStartX_pix) | measure.vecMeanXRF_pix < min(vecStimStartX_pix));
measure.vecMeanXRF_pix(indNaN) =  NaN;
measure.vecMeanDeltaT(indNaN) = NaN;

% Comput fits on based of peak times for different velocities
% peakTime =  1/v (xRF - x0) + DeltaT
% we can fit xRF and DeltaT independently for each direction

vecLeft = false(1,length(vecSpeed_pix));
vecLeft(indLeft) = true;
vecRight = false(1,length(vecSpeed_pix));
vecRight(indRight) = true;
vecResponsive = measure.vecResponsive;

if sParams.boolOnlyUseMiddleRangeSpeeds
    vecMiddleRange = record.sStimuli.vecSpeed_deg>5 & record.sStimuli.vecSpeed_deg<100;
    vecResponsive = vecResponsive & vecMiddleRange;
end

if sum(vecLeft & vecResponsive)>1
    measure.lmLeft = compact(fitlm(1./vecSpeed_pix(vecLeft & vecResponsive),measure.vecPeakTime(vecLeft & vecResponsive)));
    measure.dblDeltaTLeft = measure.lmLeft.Coefficients.Estimate(1);
    measure.dblXRFLeft_pix = measure.lmLeft.Coefficients.Estimate(2) + mean(vecStimStartX_pix(indLeft));
    measure.lmLeftFromOnset = compact(fitlm(1./vecSpeed_pix(vecLeft & vecResponsive),measure.vecOnsetTime(vecLeft & vecResponsive)));
    measure.dblDeltaTLeftFromOnset = measure.lmLeftFromOnset.Coefficients.Estimate(1);
    measure.dblXRFLeftFromOnset_pix = measure.lmLeftFromOnset.Coefficients.Estimate(2) + mean(vecStimStartX_pix(indLeft));
else
    measure.lmLeft = [];
    measure.dblDeltaTLeft = NaN;
    measure.dblXRFLeft_pix = NaN;
    measure.lmLeftFromOnset = [];
    measure.dblDeltaTLeftFromOnset = NaN;
    measure.dblXRFLeftFromOnset_pix = NaN;
end
if sum(vecRight & vecResponsive)>1
    measure.lmRight = compact(fitlm(1./vecSpeed_pix(vecRight & vecResponsive),measure.vecPeakTime(vecRight & vecResponsive)));
    measure.dblDeltaTRight = measure.lmRight.Coefficients.Estimate(1);
    measure.dblXRFRight_pix = -measure.lmRight.Coefficients.Estimate(2) + mean(vecStimStartX_pix(indRight));
    measure.lmRightFromOnset = compact(fitlm(1./vecSpeed_pix(vecRight & vecResponsive),measure.vecOnsetTime(vecRight & vecResponsive)));
    measure.dblDeltaTRightFromOnset = measure.lmRightFromOnset.Coefficients.Estimate(1);
    measure.dblXRFRightFromOnset_pix = -measure.lmRightFromOnset.Coefficients.Estimate(2) + mean(vecStimStartX_pix(indRight));
else
    measure.lmRight = [];
    measure.dblDeltaTRight = NaN;
    measure.dblXRFRight_pix = NaN;
    measure.lmRightFromOnset = [];
    measure.dblDeltaTRightFromOnset = NaN;
    measure.dblXRFRightFromOnset_pix = NaN;
end
end

function measure = compute_dot_diffhist_measures( measure, record, db )

sMeasuresFlashingDots = get_related_measures( record, 'stimulus=flashing_dots', db, measure.intIndex );
measure.lmBefore = [];

if ~isempty(sMeasuresFlashingDots ) && isfield( sMeasuresFlashingDots,'dblXRFLeft_pix') && ~isempty( sMeasuresFlashingDots.dblXRFLeft_pix)
    ind = find(record.sStimuli.vecStimStartX_pix<sMeasuresFlashingDots.dblXRFLeft_pix  & measure.vecResponsive) ;
    vecStimPos_pix = timeToStimulusCenterPosition(measure.vecOnsetTime(ind),record.sStimuli,ind);
    if length(ind)>2
        measure.lmBefore = fitlm(record.sStimuli.vecStimStartX_pix(ind), vecStimPos_pix);
    end
end


end

function measure = compute_flashing_dots_measures( measure, record, sStimParams )

measure.vecResponsive = measure.vecResponsive & measure.vecPeakTime<0.250 & measure.vecPeakTime>0.030;

measure.dblXRFLeft_pix = NaN;
measure.dblXRFRight_pix = NaN;

indLeft = find(measure.vecResponsive(1:end-1) & measure.vecResponsive(2:end),1,'first');
if ~isempty(indLeft)
    measure.dblXRFLeft_pix = ...
        record.sStimuli.vecBoundingRect{indLeft}(1) - sStimParams.intScreenWidth_pix/2;
end
indRight = find(measure.vecResponsive(2:end) & measure.vecResponsive(1:end-1),1,'last')+1;
if ~isempty(indRight)
    measure.dblXRFRight_pix = ...
        record.sStimuli.vecBoundingRect{indRight}(3) - sStimParams.intScreenWidth_pix/2;
end

end
