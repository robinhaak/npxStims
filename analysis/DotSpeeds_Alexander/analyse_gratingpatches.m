function record = analyse_gratingpatches(record,db,verbose)
%ANALYSE_DOTSPEEDS analysis moving dots stimulus with different speeds
%
%  RECORD = ANALYSE_DOTSPEEDS( RECORD, DB=[], VERBOSE=true )
%
%  Analysis responses to grating patches
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
%     record = RH_analyse_synchronization(record,verbose);
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
        logmsg('Taking last stimulus')
        ind = ind(end);
    end
    
    structEP = sSynthData.cellStim{ind}.structEP; 
    sCluster = sSynthData.sCluster;
    
    % sCluster(1).SpikeTimes: 0.0210, +0.0515, +0.0665, ... 
    % 
    % structEP fields:
    %   vecStimOnTime : 3562.873, +1.0001, +1.00003, ... synced to IMEC spiketimes? 
    %   ActOnNI: 3682.12, +0.9976, +1.0022, +0.9818, % from NI?
    %   T0 = 119.285  (roughly difference between top 2), but the
    %     difference between vecStimOnTime and ActOnNI is not completely
    %     constance, sometimes a 30ms difference
    %   ActOnSecs: 691753.585, +1.0003, +1.0018, ... from PTB?
    %   ActStartSecs: 691753.314, +1.0001, +0.9999, ... from PTB?
    %   vecTrialStartSecs: 5, +1, +1, ... planned trial starts
    
    vecStimOnTime = structEP.vecStimOnTime;
    vecStimOffTime = structEP.vecStimOffTime;

%    vecStimOnTime = structEP.ActOnNI;
%    vecStimOffTime = structEP.ActOffNI;

    if isempty(vecClustersToAnalyze)
        vecClustersToAnalyze = unique([sCluster.IdxClust]);
    end
    intNumClusters = length(vecClustersToAnalyze);
    vecSpikeTimesOfCluster = cell(intNumClusters,1)  ;
    
    
    for c = 1:intNumClusters
        indClus = find([sCluster.IdxClust]==vecClustersToAnalyze(c),1);
        if ~isempty(indClus)
            vecSpikeTimesOfCluster{c} = sCluster(indClus).SpikeTimes ;
        end
    end

    
    %     sMeta = sSynthData.sMetaNI;
    %     dblStartT = str2double(sMeta.firstSample) / str2double(sMeta.niSampRate);
    %     vecStimOnTime = structEP.ActOnNI - dblStartT;
    %     vecStimOffTime = structEP.ActOffNI - dblStartT;


else
    % Unsorted channels
    load(strLog,'structEP');
    
    strSpikesFile = fullfile(strSessionPath,'spikes.mat');
    load(strSpikesFile,'vecSpikeCh','vecSpikeSecs');

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
    matPValue = NaN(numel(vecY_pix),numel(vecX_pix));

    for intLoc = vecUniqueStims
        vecRateForLoc = vecRate(vecStimIdx==intLoc);
        vecRateOtherLocs = vecRate(vecStimIdx~=intLoc);
        matAvgResp(intLoc) = mean(vecRateForLoc);
        % Uncorrected test
        matPValue(intLoc) = ranksum(vecRateForLoc,vecRateOtherLocs,'tail','right');
    end
    measure.matAvgResp = matAvgResp - measure.dblRateSpontaneous;
    measure.matPValue = matPValue;

    [measure.dblPValue,measure.boolResponsive] = myanova(vecRate,vecStimIdx);

    measure.matSignificant = zeros(size(matPValue));
    measure.dblXRFLeft_pix = NaN;
    measure.dblXRFRight_pix = NaN;
    if measure.boolResponsive 
        % Only get largest connect image
        matSignificant = matPValue < 0.025;
        sComponents = bwconncomp(matSignificant);
        vecSizes = cellfun(@numel, sComponents.PixelIdxList);
        [~, idx] = max(vecSizes);
        mask = zeros(size(matSignificant));
        mask(sComponents.PixelIdxList{idx}) = 1;
        measure.matSignificant = mask;


        vecXLeftBorder_pix = unique(vecUniqueRects(:,1));
        vecXRightBorder_pix = unique(vecUniqueRects(:,3));

        measure.dblXRFLeft_pix = vecXLeftBorder_pix(find(max(mask),1)) - structEP.sStimParams.intScreenWidth_pix/2;
        measure.dblXRFRight_pix = vecXRightBorder_pix(find(max(mask),1,'last')) - structEP.sStimParams.intScreenWidth_pix/2;

    end
    
    % Analyse responses for peak location
    measure.vecPeakLocationSpikeT = [];
    measure.dblPeakRate = NaN;
    measure.dblPeakTime = NaN;
    measure.dblOnsetTime = NaN;
    [measure.dblResponseMax,indLoc] = max(measure.matAvgResp(:));

    [measure.dblZetaP,sZETA,sRate,~] = zetatest(vecSpikesCh,...
        vecStimOnTime(vecStimIdx == indLoc),dblDuration); % ~ is essential for correct output
    
    measure.vecPeakLocationSpikeT = [];
    
    if ~isempty(sRate) && measure.boolResponsive
        measure.vecPeakLocationSpikeT = sZETA.vecSpikeT(2:end-1);
        measure.dblPeakRate = sRate.dblPeakRate;
        measure.dblPeakTime = sRate.dblPeakTime;

        if sParams.boolFitGaussian
            intNRepeats = sum(vecStimIdx == indLoc);

            [measure.dblPeakTime,measure.dblOnsetTime] = ...
                fitGaussianPSTH( measure.vecPeakLocationSpikeT, ...
                measure.dblRateSpontaneous*intNRepeats, ...
                sRate.dblPeakTime, false, true, verbose );
        else 
            logmsg('Computing onset time is not yet implemented.')
        end
    end

    if indMeasures == 1
        measures = measure;
    else
%         try
%             measures(indMeasures) = measure;
%         catch me
%             switch me.identifier
%                 case 'MATLAB:catenate:structFieldBad'
%                     logmsg('Measures structure has changed since previous analysis. Removing previous results.');
%                 otherwise
%                     errormsg(me.message);
%             end
%             measures = [];
            measures(indMeasures) = measure;
%         end
    end

end % cluster c


record.measures = measures;


% Add dots results
h_db = get_fighandle('Neuropixels database*');
sUserData = get(h_db,'userdata');
db = sUserData.db;
record = analyse_add_dots_to_patches( record, db, true);



logmsg(['Analyzed ' recordfilter(record)]);
end
