
%GratingPatchesAnalysis
%
%online analysis for 'GratingPatches' using NI stream times
%you need a GPU to run this script!
%Robin Haak, Marh 2023
% 
% vecChs = flip(1:385); %channels to plot (1= bottom Ch)
% 
% %% query user for file names & locations
% %ap.bin en ap.meta
% [strAp,strPathIm] = uigetfile('*.ap.bin','Select imec .ap.bin file','MultiSelect','off');
% sMetaIm = DP_ReadMeta(fullpath(strPathIm,[strAp(1:end-4) '.meta']));
% 
% %structEP
% [strLog,strLogPath] = uigetfile('*.mat','Select trial-based log file','MultiSelect','off');
% load(fullpath(strLogPath,strLog)); % %#ok<ows');
% 
% %% detect spikes on each channel
% [vecSpikeCh,vecSpikeT,~] = DP_DetectSpikesInBinaryFile(fullpath(strPathIm,strAp),[],[],'int16'); %strClass='int16'
% vecSpikeSecs = double(vecSpikeT)/str2double(sMetaIm.imSampRate)+...
%     str2double(sMetaIm.firstSample)/str2double(sMetaIm.imSampRate); %convert to seconds+add offset
% intNumChs = 385; %str2num(sMetaIm.nSavedChans); %#ok<ST2NM>

%% get stimulus onset times
% %to be added: re-alignment of times based on diode data
% vecStimOnSecs = structEP.ActOnNI; %NI stream times
% vecStimOffSecs = structEP.ActOffNI;

sAP = sSynthData;
structEP = sAP.cellStim{1,2}.structEP;  
vecStimOnSecs = structEP.vecStimOnTime;
vecStimOffSecs = structEP.vecStimOffTime;

%% get grid data
vecUniqueRects = unique(structEP.vecDstRect','rows'); %unique dst rects
vecUniqueStims = 1:length(vecUniqueRects);
vecStimIdx = zeros(size(structEP.vecDstRect,2),1);
for intStim = 1:length(vecUniqueRects)
    vecStimIdx(ismember(structEP.vecDstRect',vecUniqueRects(intStim,:),'rows')) = vecUniqueStims(intStim);
end

vecX_pix = unique(vecUniqueRects(:,1))+(vecUniqueRects(1,3)-unique(vecUniqueRects(1,1)))/2;
vecY_pix = unique(vecUniqueRects(:,2))+(vecUniqueRects(1,4)-unique(vecUniqueRects(1,2)))/2;

%% loop through data
matAvgRespAll = NaN(numel(vecY_pix),numel(vecX_pix),intNumChs);
for intCh = 1:intNumChs
    vecSpikesCh = vecSpikeSecs(vecSpikeCh==intCh);
    vecRate = zeros(1,structEP.intTrialNum);
    for intTrial = 1:structEP.intTrialNum
        vecSpikeT = vecSpikesCh(vecSpikesCh>vecStimOnSecs(intTrial)&vecSpikesCh<vecStimOffSecs(intTrial));
        vecRate(intTrial) = numel(vecSpikeT)/(vecStimOffSecs(intTrial)-vecStimOnSecs(intTrial));

    end
    matAvgResp = NaN(numel(vecY_pix),numel(vecX_pix));

    for intLoc = vecUniqueStims
        matAvgResp(intLoc) = mean(vecRate(vecStimIdx==intLoc));
    end
    sParams.dblSecsFromPrevStimOff = 0.1; %s, for computing unit's baseline rate
    dblRateSpontaneous = computeRateSpontaneous(vecSpikeT,vecStimOnSecs,vecStimOffSecs,sParams);
    matAvgRespAll(:,:,intCh) = matAvgResp-dblRateSpontaneous;
end

%% plot data
%interpolate
vecX_pix_interp = linspace(vecX_pix(1),vecX_pix(end),16);
vecY_pix_interp = linspace(vecY_pix(1),vecY_pix(end),9);

%get colormap(s)
cellColorMaps = RH_ColorMaps;
%%
%loop through channels
set(0,'DefaultFigureWindowStyle','docked')
for intCh = vecChs
    matAvgRespAll_interp = matAvgRespAll(:,:,intCh);
    figure; hold on; 
    title(['Channel: ' num2str(intCh)]);
    imagesc(vecX_pix_interp,vecY_pix_interp,matAvgRespAll_interp);
    set(gca, 'YDir','reverse');
    colormap(cellColorMaps{2});
    cb=colorbar;
    cb.Label.String='spks/s';
    axis image
    fixfig;
end

%%
function dblRateSpontaneous = computeRateSpontaneous(vecSpikeTimes,vecStimOnSecs,vecStimOffSecs,sParams)
%compute unit's spontaneous/baseline rate
intCount = 0;
dblPeriod = 0;
for intTrial=1:length(vecStimOffSecs)-1
    intCount = intCount + ...
        length(find(vecSpikeTimes>(vecStimOffSecs(intTrial)+sParams.dblSecsFromPrevStimOff) & ...
        vecSpikeTimes<vecStimOnSecs(intTrial+1)));
    dblPeriod = dblPeriod + vecStimOnSecs(intTrial+1) - (vecStimOffSecs(intTrial)+sParams.dblSecsFromPrevStimOff);
end
dblRateSpontaneous = intCount / dblPeriod;
if dblPeriod<5
    fprintf('Less than 5s to compute spontaneous rate.\n')
end
if intCount<10
    fprintf('Less than 10 spikes to compute spontaneous rate.\n')
end
end
