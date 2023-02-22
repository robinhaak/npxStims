function record = RH_analyse_synchronization(record,verbose)
%RH_ANALYSE_SYNCHRONIZATION checks time stamps and synchronization
%
% 2023,  Alexander Heimel

if nargin<2 || isempty(verbose)
    verbose = true;
end

measures = record.measures;

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

    cellExpType = cellfun(@(x) x.structEP.strExpType, sSynthData.cellStim, 'UniformOutput', false);
    ind = find(contains(cellExpType,record.stimulus));
    if isempty(ind)
        logmsg(['Could not find stimulus ' record.stimulus ' in '  strLog]);
        return
    end
    if length(ind)>1
        logmsg(['More than one stimulus ' record.stimulus ' in '  strLog]);
        logmsg('Using first one stimulus set.');
        ind = ind(1);
    end
    
    structEP = sSynthData.cellStim{ind}.structEP; 
    
    % sCluster(1).SpikeTimes: 0.0210, +0.0515, +0.0665, ... 
    % 
    % structEP fields:
    %   vecStimOnTime : 3562.873, +1.0001, +1.00003, ... synced to IMEC spiketimes? 
    %   ActOnNI: 3682.12, +0.9976, +1.0022, +0.9818, % from NI?
    %   T0 = 119.285  (roughly difference between top 2), but the
    %     difference between vecStimOnTime and ActOnNI is not completely
    %     constance, sometimes a 30ms difference
    %   ActStartSecs: 691753.314, +1.0001, +0.9999, ... First Flip time of PTB in trial (blank stim)
    %   ActOnSecs: 691753.585, +1.0003, +1.0018, ... First stimulus PTB,
    %        just before time from ActOnNI (see RH_GratingPatches)
    %
    % for GratingPatches
    %   vecTrialStartSecs: 5, +1, +1, ... planned trial starts
    
   
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

    
    %     sMeta = sSynthData.sMetaNI;
    %     dblStartT = str2double(sMeta.firstSample) / str2double(sMeta.niSampRate);
    %     vecStimOnTime = structEP.ActOnNI - dblStartT;
    %     vecStimOffTime = structEP.ActOffNI - dblStartT;


else
    % Unsorted channels
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

end


    %%
    figure('Name','Sync','NumberTitle','off' );
    subplot(2,2,1);
        
    lm = fitlm(structEP.ActOnSecs,structEP.ActOnNI);
    dblNIClockMultiplier = 1/lm.Coefficients.Estimate(2);
    logmsg(['dblNIClockMultiplier = ' num2str(dblNIClockMultiplier,'%0.8f')]);
    dblShiftActOnNIActOnSecs = structEP.ActOnNI(1) - structEP.ActOnSecs(1);
    plot(structEP.ActOnSecs - structEP.ActOnSecs(1) ,structEP.ActOnNI*dblNIClockMultiplier - structEP.ActOnSecs - dblShiftActOnNIActOnSecs,'.');
    title([num2str(dblNIClockMultiplier,'%.9f') ' x']);
    xlabel('\DeltaActOnSecs ');
    ylabel('\DeltaActOnNI*Mult - \DeltaActOnSecs');
    
    if isfield(structEP,'vecStimOnTime')
        subplot(2,2,2);
        dblShiftVecStimOnTimeActOnNI = structEP.vecStimOnTime(1) - structEP.ActOnNI(1);
        plot(structEP.ActOnNI-structEP.ActOnNI(1),structEP.vecStimOnTime - structEP.ActOnNI - dblShiftVecStimOnTimeActOnNI,'.' );
        title([num2str(dblShiftVecStimOnTimeActOnNI,'%.1f') ' s']);
        xlabel('\DeltaActOnNI');
        ylabel('\DeltaVecStimOnTime-\DeltaActOnNI');
    end
    
    subplot(2,2,3);
    dblShiftActStartSecsActOnSecs = structEP.ActStartSecs(1) - structEP.ActOnSecs(1);
    plot(structEP.ActStartSecs-structEP.ActOnSecs(1),structEP.ActStartSecs - structEP.ActOnSecs- dblShiftActStartSecsActOnSecs,'.' );
    title([num2str(dblShiftActStartSecsActOnSecs,'%.1f') ' s']);
    xlabel('\DeltaActStartSecs');
    ylabel('\DeltaActOnSecs-\DeltaActStartSecs');

    subplot(2,2,4)
    histogram(structEP.ActStartSecs - structEP.ActOnSecs- dblShiftActStartSecsActOnSecs,20)
    ylabel('Count');
    set(gca,'xtick',[-2/60 -1/60 0 1/60 2/60]);
    %%

record.measures = measures;

logmsg(['Analyzed ' recordfilter(record)]);
end
