function record = analyse_gratingpatches(record,verbose)
%ANALYSE_DOTSPEEDS analysis moving dots stimulus with different speeds
%
%  RECORD = ANALYSE_DOTSPEEDS( RECORD, VERBOSE )
%
%  Analysis responses to grating patches
%
% 2022-2023, Robin Haak, Alexander Heimel

if nargin<2 || isempty(verbose)
    verbose = true;
end

sParams = RH_defaultParameters();


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
if contains('sCluster',{sVars.name})
    % single unit data
    load(strLog,'sCluster','structEP','sSynthData'); %rm sSynthData later

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


% get grid data
vecUniqueRects = unique(structEP.vecDstRect','rows'); %unique dst rects
vecUniqueStims = 1:length(vecUniqueRects);
vecStimIdx = zeros(size(structEP.vecDstRect,2),1);
for intStim = 1:length(vecUniqueRects)
    vecStimIdx(ismember(structEP.vecDstRect',vecUniqueRects(intStim,:),'rows')) = vecUniqueStims(intStim);
end

vecX_pix = unique(vecUniqueRects(:,1))+(vecUniqueRects(1,3)-unique(vecUniqueRects(1,1)))/2;
vecY_pix = unique(vecUniqueRects(:,2))+(vecUniqueRects(1,4)-unique(vecUniqueRects(1,2)))/2;


record.sStimuli = [];
record.sStimuli.vecX_pix = vecX_pix;
record.sStimuli.vecY_pix = vecY_pix;


%% loop through data

measures = struct([]);
indMeasures = 0;

dblDuration = mode(structEP.vecTrialStimOffSecs - structEP.vecTrialStimOnSecs);

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


    vecSpikesCh = vecSpikeTimesOfCluster{c};
    vecRate = zeros(1,structEP.intTrialNum);
    for intTrial = 1:structEP.intTrialNum
        vecSpikeT = vecSpikesCh(vecSpikesCh>vecStimOnTime(intTrial) & vecSpikesCh<vecStimOffTime(intTrial));
        vecRate(intTrial) = numel(vecSpikeT)/(vecStimOffTime(intTrial)-vecStimOnTime(intTrial));

    end
    matAvgResp = NaN(numel(vecY_pix),numel(vecX_pix));

    for intLoc = vecUniqueStims
        matAvgResp(intLoc) = mean(vecRate(vecStimIdx==intLoc));
    end
    measure.matAvgResp = matAvgResp - measure.dblRateSpontaneous;



    % Analyse responses for peak location
    measure.dblPeakRate = NaN;
    measure.dblPeakTime = NaN;
    [measure.dblResponseMax,indLoc] = max(measure.matAvgResp(:));
    indTrials = find(vecStimIdx == indLoc);

    [measure.dblZetaP,sZETA,sRate,~] = zetatest(vecSpikesCh,...
        vecStimOnTime(indTrials),dblDuration); % ~ is essential for correct output
    if ~isempty(sRate)
        measure.dblPeakRate = sRate.dblPeakRate;
        measure.dblPeakTime = sRate.dblPeakTime;
    end


    if indMeasures == 1
        measures = measure;
    else
        measures(indMeasures) = measure;
    end

end % cluster c






record.measures = measures;

logmsg(['Analyzed ' recordfilter(record)]);
end
