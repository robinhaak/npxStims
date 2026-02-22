% RH_GratingPatchesAnalysis_MultiShank_RF_NI_NoInterp
% Spike detection + RF mapping (NI times, no interpolation)
% Computes RF maps, COM-based sharpness, and RF SNR
%
% Robin Haak, 2026

clc; clear; close all;

%% 1. Select files
[strAp,strPathIm] = uigetfile('*.ap.bin','Select imec .ap.bin file','MultiSelect','off');
metaFile = fullfile(strPathIm, [strAp(1:end-4) '.meta']);

[strLog,strLogPath] = uigetfile('*.mat','Select trial-based log file','MultiSelect','off');
load(fullfile(strLogPath,strLog)); % expects variable containing structEP

%% 2. Detect spikes
fprintf('Detecting spikes...\n');
[vecSpikeCh, vecSpikeT, ~] = DP_DetectSpikesInBinaryFile(fullfile(strPathIm,strAp), [], [], 'int16');
vecSpikeCh = double(vecSpikeCh);

meta = ReadMeta(metaFile);
vecSpikeSecs = double(vecSpikeT)/str2double(meta.imSampRate) + ...
               str2double(meta.firstSample)/str2double(meta.imSampRate);

%% 3. Get shank/channel geometry
if isfield(meta,'snsGeomMap')
    [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = geomMapToGeom(meta);
elseif isfield(meta,'snsShankMap')
    [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = shankMapToGeom(meta);
else
    error('No geometry information found in metadata.');
end

shankInd = shankInd + 1;

% Optional: plot shank layout
plotSaved(xCoord, yCoord, shankInd, meta);

%% 4. Extract stimulus info
vecStimOnSecs  = structEP.ActOnNI;
vecStimOffSecs = structEP.ActOffNI;

%% 5. Get grid positions
vecUniqueRects = unique(structEP.vecDstRect','rows'); 
vecUniqueStims = 1:length(vecUniqueRects);

vecStimIdx = zeros(size(structEP.vecDstRect,2),1);
for intStim = 1:length(vecUniqueRects)
    vecStimIdx(ismember(structEP.vecDstRect', vecUniqueRects(intStim,:), 'rows')) = vecUniqueStims(intStim);
end

vecX_pix = unique(vecUniqueRects(:,1)) + (vecUniqueRects(1,3)-unique(vecUniqueRects(1,1)))/2;
vecY_pix = unique(vecUniqueRects(:,2)) + (vecUniqueRects(1,4)-unique(vecUniqueRects(1,2)))/2;

%% 6. Compute RF per channel
intNumChs = numel(shankInd);
nX = numel(vecX_pix);
nY = numel(vecY_pix);

matAvgRespAll = nan(nY, nX, intNumChs);

fprintf('Computing RF for %d channels...\n', intNumChs);

vecSpikeCh = double(vecSpikeCh);
validIdx = vecSpikeCh >= 1 & vecSpikeCh <= intNumChs;
vecSpikeCh = vecSpikeCh(validIdx);
vecSpikeSecs = vecSpikeSecs(validIdx);

stimDur = vecStimOffSecs - vecStimOnSecs;

for intCh = 1:intNumChs
    vecSpikesCh = vecSpikeSecs(vecSpikeCh == intCh);
    spikeCountsPerTrial = arrayfun(@(on,off) sum(vecSpikesCh >= on & vecSpikesCh < off), ...
                                   vecStimOnSecs, vecStimOffSecs);

    rfMat = nan(nY, nX);
    for stimIdx = 1:length(vecUniqueStims)
        trials = find(vecStimIdx == stimIdx);
        if isempty(trials), continue; end
        [yPos, xPos] = ind2sub([nY, nX], stimIdx);
        rfMat(yPos, xPos) = mean(spikeCountsPerTrial(trials) ./ stimDur(trials));
    end

    matAvgRespAll(:,:,intCh) = rfMat;
end

%% 7. Compute RF COM-based sharpness and RF SNR
vecRF_COMspread = nan(intNumChs,1); % COM-based sharpness
vecRF_SNR      = nan(intNumChs,1);  % RF SNR

for intCh = 1:intNumChs
    rfMat = matAvgRespAll(:,:,intCh);
    rfVals = rfMat(~isnan(rfMat));

    if numel(rfVals) < 2 || std(rfVals) == 0
        continue
    end

    % COM-based sharpness
    [xGrid, yGrid] = meshgrid(1:nX, 1:nY);
    rfMatFilled = rfMat; rfMatFilled(isnan(rfMatFilled)) = 0;
    rfNorm = rfMatFilled / sum(rfMatFilled(:));
    xCOM = sum(xGrid(:).*rfNorm(:));
    yCOM = sum(yGrid(:).*rfNorm(:));
    rfSpread = sum(rfNorm(:).* ((xGrid(:)-xCOM).^2 + (yGrid(:)-yCOM).^2));
    vecRF_COMspread(intCh) = 1 / (rfSpread + eps);

    % RF SNR
    vecRF_SNR(intCh) = (max(rfVals) - mean(rfVals)) / std(rfVals);
end

% Remove extreme outliers
vecRF_COMspread(vecRF_COMspread > prctile(vecRF_COMspread,95)) = prctile(vecRF_COMspread,95);
vecRF_SNR(vecRF_SNR > prctile(vecRF_SNR,95)) = prctile(vecRF_SNR,95);

% Ensure NaNs remain red
vecRF_COMspread(isnan(vecRF_COMspread)) = NaN;
vecRF_SNR(isnan(vecRF_SNR)) = NaN;

%% 8. Plot all shanks together
% --------------------
% COM-based RF sharpness
figure('Name','RF Sharpness (COM Spread)','Color','w'); hold on;
for sI = 1:nShank
    chIdx = find(shankInd == sI);
    xOffset = shankPitch*(sI-1);

    nanCh = isnan(vecRF_COMspread(chIdx));
    validCh = ~nanCh;

    scatter(xOffset + xCoord(chIdx(validCh)), ...
            yCoord(chIdx(validCh)), ...
            40, vecRF_COMspread(chIdx(validCh)), 'filled', 'square');
    scatter(xOffset + xCoord(chIdx(nanCh)), ...
            yCoord(chIdx(nanCh)), 40, 'r', 'filled', 'square');
end
colormap parula; colorbar; clim([0 max(vecRF_COMspread)]);
title('RF Sharpness (1 / COM Spread) — All Shanks');
xlabel('Probe X'); ylabel('Probe Y'); axis equal;

% --------------------
% RF SNR
figure('Name','RF SNR','Color','w'); hold on;
for sI = 1:nShank
    chIdx = find(shankInd == sI);
    xOffset = shankPitch*(sI-1);

    nanCh = isnan(vecRF_SNR(chIdx));
    validCh = ~nanCh;

    scatter(xOffset + xCoord(chIdx(validCh)), ...
            yCoord(chIdx(validCh)), ...
            40, vecRF_SNR(chIdx(validCh)), 'filled', 'square');
    scatter(xOffset + xCoord(chIdx(nanCh)), ...
            yCoord(chIdx(nanCh)), 40, 'r', 'filled', 'square');
end
colormap parula; colorbar; clim([0 max(vecRF_SNR)]);
title('RF SNR — All Shanks');
xlabel('Probe X'); ylabel('Probe Y'); axis equal;

%% 9. Plot example RF maps per shank (top N channels)
topN = 20;
for sh = 1:nShank
    channelsInShank = find(shankInd == sh);
    if isempty(channelsInShank), continue; end

    [~, sortIdx] = sort(max(matAvgRespAll(:,:,channelsInShank),[],[1 2]), 'descend');
    topChans = channelsInShank(sortIdx(1:min(topN,numel(sortIdx))));

    figure('Name',sprintf('RF Examples — Shank %d',sh),'Color','w'); hold on; %maxfig;
    nRows = 4; nCols = ceil(topN/nRows);
    for k = 1:numel(topChans)
        intCh = topChans(k);
        subplot(nRows,nCols,k); hold on;
        imagesc(vecX_pix, vecY_pix, squeeze(matAvgRespAll(:,:,intCh)));
        colormap parula; colorbar; title(colorbar,'Rate (sp/s)');
        xlabel('X (pix)'); ylabel('Y (pix)'); axis equal; axis tight;
        title(sprintf('Channel %d', intCh));
    end
    sgtitle(sprintf('RF Examples — Shank %d',sh));
end