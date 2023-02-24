

%RH_GratingPatchesAnalysis
%
%Robin Haak, 2023

%% query user for file names & locations
% %ap.bin en ap.meta
% [strAp,strPathIm] = uigetfile('*.ap.bin','Select imec .ap.bin file','MultiSelect','off');
% sMetaIm = DP_ReadMeta(fullpath(strPathIm,[strAp(1:end-4) '.meta']));
% 
% %structEP
% [strLog,strLogPath] = uigetfile('*.mat','Select trial-based log file','MultiSelect','off');
% load(fullpath(strLogPath,strLog)); % %#ok<ows');

%% detect spikes on each channel
% [vecSpikeCh,vecSpikeT,~] = DP_DetectSpikesInBinaryFile(fullpath(strPathIm,strAp),[],[],'int16'); %strClass='int16'
% vecSpikeSecs = double(vecSpikeT)/str2double(sMetaIm.imSampRate)+...
%     str2double(sMetaIm.firstSample)/str2double(sMetaIm.imSampRate); %convert to seconds+add offset
% intNumChs = str2num(sMetaIm.nSavedChans); %#ok<ST2NM>

%% get stimulus onset times
% %to be added: re-alignment of times based on diode data
% vecStimOnSecs = structEP.ActOnNI; %NI stream times
% vecStimOffSecs = structEP.ActOffNI;

structEP = sAP.cellBlock{1,1};
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
sCluster = sAP.sCluster;
intNumClu = numel(sCluster);

matAvgRespAll = NaN(numel(vecY_pix),numel(vecX_pix),intNumClu);
matAvgRespAllBlSub = NaN(numel(vecY_pix),numel(vecX_pix),intNumClu);

for intClu = 1:intNumClu
    vecRate = zeros(1,structEP.intTrialNum);
    vecBase = zeros(1,structEP.intTrialNum);
    vecSpikesCh = sAP.sCluster(intClu).SpikeTimes;  
    for intTrial = 1:structEP.intTrialNum
        vecSpikeT = vecSpikesCh(vecSpikesCh>vecStimOnSecs(intTrial)&vecSpikesCh<vecStimOffSecs(intTrial));
        vecRate(intTrial) = numel(vecSpikeT)/(vecStimOffSecs(intTrial)-vecStimOnSecs(intTrial));
        vecSpikeB = vecSpikesCh(vecSpikesCh>vecStimOnSecs(intTrial)-1&vecSpikesCh<vecStimOnSecs(intTrial));
        vecBase(intTrial) = numel(vecSpikeB)/1;
    end
    matAvgResp = NaN(numel(vecY_pix),numel(vecX_pix));
    matZeta = NaN(numel(vecY_pix),numel(vecX_pix));
%     matAvgRespBlSub = NaN(numel(vecY_pix),numel(vecX_pix));

    for intLoc = vecUniqueStims
        matAvgResp(intLoc) = mean(vecRate(vecStimIdx==intLoc));
%         [dblZetaP,sZeta] = zetatest(vecSpikesCh,vecStimOnSecs(vecStimIdx==intLoc),mean(vecStimOffSecs(vecStimIdx==intLoc)-vecStimOnSecs(vecStimIdx==intLoc)));
%         matZeta(intLoc) = dblZetaP;
%         matAvgRespBlSub(intLoc) = mean(vecRate(vecStimIdx==intLoc))-mean(vecBase(vecStimIdx==intLoc));
    end
    matAvgRespAll(:,:,intClu) = matAvgResp;
    sParams.dblSecsFromPrevStimOff = 0.1;
    dblRateSpontaneous = computeRateSpontaneous(vecSpikesCh,vecStimOnSecs,vecStimOffSecs,sParams);
    matAvgRespAllBlSub(:,:,intClu) = matAvgResp - dblRateSpontaneous; 
end

%% plot data
%interpolate
vecX_pix_interp = linspace(vecX_pix(1),vecX_pix(end),37);
vecY_pix_interp = linspace(vecY_pix(1),vecY_pix(end),21);

%get colormap(s)
cellColorMaps = RH_ColorMaps;
%%
%loop through channels
for intClu = 223:225 %numel(sCluster)
    matAvgRespAll_interp = interp2(matAvgRespAll(:,:,intClu),2);
    matAvgRespAllBlSub_interp = interp2(matAvgRespAllBlSub(:,:,intClu),2);
    figure;hold on;title(['Cluster: ' num2str(intClu)]);
%     subplot(2,1,1);
%     imagesc(vecX_pix_interp,vecY_pix_interp,matAvgRespAll_interp);
%     set(gca, 'YDir','reverse'); colormap(cellColorMaps{2});cb=colorbar;cb.Label.String='spks/s';
%     fixfig;
%     subplot(2,1,2);
   imagesc(vecX_pix_interp,vecY_pix_interp,matAvgRespAllBlSub_interp);
    %imagesc(vecX_pix,vecY_pix,matAvgRespAllBlSub(:,:,intClu));
    set(gca, 'YDir','reverse'); colormap(cellColorMaps{2});colorbar;cb=colorbar;cb.Label.String='spks/s';
    fixfig;
   % pause
end






