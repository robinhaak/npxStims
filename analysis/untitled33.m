%% =========================================================
% RH_GratingPatchesAnalysis_MultiShank_RF_NI_NoInterp
% Spike detection + RF mapping (NI times, no interpolation)
% Uses structEP from loaded data
% =========================================================
%
% Robin Haak, 2026
clc; clear; close all;

%% 1. Select files
[apFile, apPath] = uigetfile('*.ap.bin','Select imec .ap.bin file','MultiSelect','off');
metaFile = fullfile(apPath, [apFile(1:end-4) '.meta']);

[logFile, logPath] = uigetfile('*.mat','Select trial-based log file','MultiSelect','off');
load(fullfile(logPath, logFile)); % expects variable containing structEP

%% 2. Detect spikes
fprintf('Detecting spikes...\n');
[vecSpikeCh, vecSpikeT, ~] = DP_DetectSpikesInBinaryFile(fullfile(apPath, apFile), [], [], 'int16');
vecSpikeCh = double(vecSpikeCh);

meta = ReadMeta(metaFile);
vecSpikeSecs = double(vecSpikeT) / str2double(meta.imSampRate) + ...
    str2double(meta.firstSample) / str2double(meta.imSampRate);

save('tempspikes.mat');
pause;

%% 3. Get shank/channel geometry
if isfield(meta, 'snsGeomMap')
    [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = geomMapToGeom(meta);
elseif isfield(meta, 'snsShankMap')
    [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = shankMapToGeom(meta);
else
    error('No geometry information found in metadata.');
end

% Optional: plot shank layout
plotSaved(xCoord, yCoord, shankInd, meta);

%% 4. Extract stimulus info
vecStimOnSecs  = structEP.ActOnNI;
vecStimOffSecs = structEP.ActOffNI;
stimDur = vecStimOffSecs - vecStimOnSecs;

%% 5. Get unique stimulus positions
vecUniqueRects = unique(structEP.vecDstRect','rows');
vecUniqueStims = 1:size(vecUniqueRects,1);

% map trials to stimulus indices
vecStimIdx = zeros(size(structEP.vecDstRect,2),1);
for stimIdx = 1:length(vecUniqueStims)
    vecStimIdx(ismember(structEP.vecDstRect', vecUniqueRects(stimIdx,:), 'rows')) = vecUniqueStims(stimIdx);
end

% compute center positions for each stimulus
vecX_pix = unique(vecUniqueRects(:,1)) + mean(vecUniqueRects(:,3)-vecUniqueRects(:,1))/2;
vecY_pix = unique(vecUniqueRects(:,2)) + mean(vecUniqueRects(:,4)-vecUniqueRects(:,2))/2;

%% 6. Compute RF per channel
intNumChs = numel(shankInd);
nX = numel(vecX_pix);
nY = numel(vecY_pix);
matAvgRespAll = nan(nY, nX, intNumChs);

% filter valid spike channels
validCh = vecSpikeCh >= 1 & vecSpikeCh <= intNumChs;
vecSpikeCh = vecSpikeCh(validCh);
vecSpikeSecs = vecSpikeSecs(validCh);

fprintf('Computing RF for %d channels...\n', intNumChs);

for ch = 1:intNumChs
    spikesCh = vecSpikeSecs(vecSpikeCh == ch);
    spikeCountsPerTrial = arrayfun(@(on,off) sum(spikesCh >= on & spikesCh < off), ...
        vecStimOnSecs, vecStimOffSecs);
    rfMat = nan(nY, nX);

    for stimIdx = 1:length(vecUniqueStims)
        trials = find(vecStimIdx == stimIdx);
        if isempty(trials), continue; end
        [yPos, xPos] = ind2sub([nY, nX], stimIdx);
        rfMat(yPos, xPos) = mean(spikeCountsPerTrial(trials) ./ stimDur(trials));
    end
    matAvgRespAll(:,:,ch) = rfMat;
end

%% 6b. Fit 2D Gaussian to each RF
rfCenterX_Gauss = nan(intNumChs,1);
rfCenterY_Gauss = nan(intNumChs,1);
rfSigmaX_Gauss  = nan(intNumChs,1);
rfSigmaY_Gauss  = nan(intNumChs,1);
rfAmp_Gauss     = nan(intNumChs,1);

[Xgrid, Ygrid] = meshgrid(vecX_pix, vecY_pix);
fprintf('Fitting 2D Gaussian to RF maps...\n');
for ch = 1:intNumChs
    rates = squeeze(matAvgRespAll(:,:,ch));
    if all(isnan(rates(:))), continue; end
    [cx, cy, sx, sy, amp] = fitGaussian2D_simple(rates, Xgrid, Ygrid);
    rfCenterX_Gauss(ch) = cx;
    rfCenterY_Gauss(ch) = cy;
    rfSigmaX_Gauss(ch) = sx;
    rfSigmaY_Gauss(ch) = sy;
    rfAmp_Gauss(ch)    = amp;
end

%% 7. Compute RF CoM, amplitude, and SNR
rfCenterX_CoM = nan(intNumChs,1);
rfCenterY_CoM = nan(intNumChs,1);
rfAmplitude   = nan(intNumChs,1);
rfSNR         = nan(intNumChs,1);

for ch = 1:intNumChs
    rates = squeeze(matAvgRespAll(:,:,ch));
    if all(isnan(rates(:))), continue; end
    rates(isnan(rates)) = 0;
    
    total = sum(rates(:));
    if total > 0
        [Ygrid, Xgrid] = meshgrid(vecY_pix, vecX_pix);
        Xgrid = Xgrid'; Ygrid = Ygrid'; % match dimensions
        rfCenterX_CoM(ch) = sum(rates(:) .* Xgrid(:)) / total;
        rfCenterY_CoM(ch) = sum(rates(:) .* Ygrid(:)) / total;
    end
    
    rfAmplitude(ch) = max(rates(:));
    
    mu = mean(rates(:),'omitnan');
    sigma = std(rates(:),'omitnan');
    if sigma>0
        rfSNR(ch) = mu/sigma;
    else
        rfSNR(ch) = NaN;
    end
end

%% 8. Visualization
shankIdx = shankInd(:) + 1;  % convert 0-based to 1-based
nShank = max(shankIdx);
shankColors = lines(nShank);
dx = mean(diff(vecX_pix)); dy = mean(diff(vecY_pix));
minSNR = 3;
validChans = find(~isnan(rfCenterX_CoM) & rfSNR >= minSNR);

%% 8a. Combined RF Hotspots + CoM + Gaussian
figure('Name','RF Hotspots + Centers','Color','w'); hold on;
for ch = validChans'
    rates = squeeze(matAvgRespAll(:,:,ch));
    if all(isnan(rates(:))), continue; end
    zMap = (rates - mean(rates(:),'omitnan')) / std(rates(:),'omitnan');
    mask = zMap >= 2;
    if ~any(mask(:)), continue; end
    
    [yIdx, xIdx] = find(mask);
    for k = 1:length(xIdx)
        x = vecX_pix(xIdx(k)) - dx/2;
        y = vecY_pix(yIdx(k)) - dy/2;
        intensity = min((zMap(yIdx(k), xIdx(k))-2)/3,1);
        rectangle('Position',[x,y,dx,dy], 'FaceColor',shankColors(shankIdx(ch),:)*(0.3+0.7*intensity), 'EdgeColor','none');
    end
end

% plot CoM
for ch = validChans'
    plot(rfCenterX_CoM(ch), rfCenterY_CoM(ch), 's', ...
        'MarkerEdgeColor','k', 'MarkerFaceColor', shankColors(shankIdx(ch),:), ...
        'MarkerSize', 4+min(rfSNR(ch),20)/20*8, 'LineWidth', 1.5);
end

% plot Gaussian centers
for ch = validChans'
    plot(rfCenterX_Gauss(ch), rfCenterY_Gauss(ch), 'o', ...
        'MarkerEdgeColor','k', 'MarkerFaceColor', shankColors(shankIdx(ch),:), ...
        'MarkerSize', 2+min(rfSNR(ch),20)/20*8, 'LineWidth', 1.5);
end

xlabel('RF X (px)'); ylabel('RF Y (px)');
title(sprintf('CoM (□) & Gaussian (○) Centers | SNR ≥ %.1f', minSNR));
axis equal; grid on;

%% 8b. Per-Shank Overview
for sh = 1:nShank
    chansShank = find(shankIdx==sh & ~isnan(rfCenterX_CoM) & rfSNR>=minSNR);
    if isempty(chansShank), continue; end
    figure('Name',sprintf('Shank %d RF Overview',sh),'Color','w'); hold on;
    
    for ch = chansShank'
        rates = squeeze(matAvgRespAll(:,:,ch));
        mu = mean(rates(:),'omitnan'); sigma = std(rates(:),'omitnan');
        if sigma==0, continue; end
        zMap = (rates - mu)/sigma;
        mask = zMap >= 2;
        [yIdx, xIdx] = find(mask);
        for k = 1:length(xIdx)
            x = vecX_pix(xIdx(k)) - dx/2;
            y = vecY_pix(yIdx(k)) - dy/2;
            intensity = min((zMap(yIdx(k), xIdx(k))-2)/3,1);
            rectangle('Position',[x,y,dx,dy],'FaceColor',shankColors(sh,:)*(0.3+0.7*intensity),'EdgeColor','none');
        end
    end
    
    % plot centers
    plot(rfCenterX_CoM(chansShank), rfCenterY_CoM(chansShank), 's', ...
        'MarkerEdgeColor','k','MarkerFaceColor',shankColors(sh,:), 'MarkerSize',8,'LineWidth',1.5);
    plot(rfCenterX_Gauss(chansShank), rfCenterY_Gauss(chansShank), 'o', ...
        'MarkerEdgeColor','k','MarkerFaceColor',shankColors(sh,:), 'MarkerSize',10,'LineWidth',1.5);
    
    xlabel('RF X (px)'); ylabel('RF Y (px)');
    title(sprintf('Shank %d RF Hotspots + Centers', sh));
    axis equal; grid on;
end

%% 9. Plot example RF maps per shank (top N channels)
topN = 12;
for sh = 1:nShank
    chansShank = find(shankIdx == sh);
    if isempty(chansShank), continue; end
    [~, sortIdx] = sort(max(matAvgRespAll(:,:,chansShank),[],[1 2]), 'descend');
    topChans = chansShank(sortIdx(1:min(topN,numel(sortIdx))));
    
    figure('Name',sprintf('RF Examples — Shank %d',sh),'Color','w');
    nRows = 3; nCols = ceil(topN/nRows);
    for k = 1:numel(topChans)
        ch = topChans(k);
        subplot(nRows,nCols,k); imagesc(vecX_pix, vecY_pix, squeeze(matAvgRespAll(:,:,ch)));
        colormap hot; colorbar;
        xlabel('X (pix)'); ylabel('Y (pix)'); axis equal; axis tight;
        title(sprintf('Channel %d', ch));
    end
end
