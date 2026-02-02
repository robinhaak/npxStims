%% 6. Compute RF per channel (Y × X matrices, no interpolation)

% Use number of physical channels from metadata, not max spike channel
intNumChs = numel(shankInd); 
nX = numel(vecX_pix);
nY = numel(vecY_pix);

matAvgRespAll = nan(nY, nX, intNumChs);

fprintf('Computing RF for %d channels...\n', intNumChs);

% Ensure spike channel indices are valid
vecSpikeCh = double(vecSpikeCh);
validIdx = vecSpikeCh >= 1 & vecSpikeCh <= intNumChs;
vecSpikeCh = vecSpikeCh(validIdx);
vecSpikeSecs = vecSpikeSecs(validIdx);

% Precompute stimulus duration
stimDur = vecStimOffSecs - vecStimOnSecs;

for intCh = 1:intNumChs
    % spikes for this channel
    vecSpikesCh = vecSpikeSecs(vecSpikeCh == intCh);
    
    % count spikes per trial
    spikeCountsPerTrial = arrayfun(@(on,off) sum(vecSpikesCh >= on & vecSpikesCh < off), ...
                                   vecStimOnSecs, vecStimOffSecs);

    % compute RF matrix
    rfMat = nan(nY, nX);
    for stimIdx = 1:length(vecUniqueStims)
        trials = find(vecStimIdx == stimIdx);
        if isempty(trials), continue; end
        
        % convert stimIdx to (y,x) index in grid
        [yPos, xPos] = ind2sub([nY, nX], stimIdx);
        rfMat(yPos, xPos) = mean(spikeCountsPerTrial(trials) ./ stimDur(trials));
    end
    
    matAvgRespAll(:,:,intCh) = rfMat;
end

%% 7. Improved RF center estimation

% Initialize output variables
rfCenterX_CoM = nan(1, intNumChs);      % Weighted center of mass
rfCenterY_CoM = nan(1, intNumChs);
rfPeakX = nan(1, intNumChs);            % Peak response location
rfPeakY = nan(1, intNumChs);
rfCenterX_Gauss = nan(1, intNumChs);    % Gaussian fit center
rfCenterY_Gauss = nan(1, intNumChs);
rfSigmaX = nan(1, intNumChs);           % Gaussian width
rfSigmaY = nan(1, intNumChs);
rfAmplitude = nan(1, intNumChs);        % Response amplitude
rfBaseline = nan(1, intNumChs);         % Baseline firing rate
rfArea = nan(1, intNumChs);             % RF area
rfSNR = nan(1, intNumChs);              % Signal-to-noise ratio

[Xgrid, Ygrid] = meshgrid(vecX_pix, vecY_pix);

fprintf('Estimating RF centers for %d channels...\n', intNumChs);

for intCh = 1:intNumChs
    rates = squeeze(matAvgRespAll(:,:,intCh));
    
    % Skip channels with no data
    if all(isnan(rates(:))), continue; end
    
    % Replace NaNs with zeros for computation
    rates(isnan(rates)) = 0;
    
    % Smooth to reduce noise (adjust sigma as needed)
    ratesSmooth = imgaussfilt(rates, 0.8, 'Padding', 'replicate');
    
    % Estimate baseline (median of lower 25% of responses)
    sortedRates = sort(ratesSmooth(:));
    baseline = median(sortedRates(1:ceil(0.25*length(sortedRates))));
    rfBaseline(intCh) = baseline;
    
    % Subtract baseline
    ratesNorm = ratesSmooth - baseline;
    ratesNorm(ratesNorm < 0) = 0;
    
    % Skip if no above-baseline responses
    if max(ratesNorm(:)) == 0, continue; end
    
    % ===== Method 1: Weighted Center of Mass (no thresholding) =====
    totalWeight = sum(ratesNorm(:));
    if totalWeight > 0
        rfCenterX_CoM(intCh) = sum(ratesNorm(:) .* Xgrid(:)) / totalWeight;
        rfCenterY_CoM(intCh) = sum(ratesNorm(:) .* Ygrid(:)) / totalWeight;
    end
    
    % ===== Method 2: Peak Location =====
    [maxRate, maxIdx] = max(ratesNorm(:));
    rfAmplitude(intCh) = maxRate;
    rfPeakX(intCh) = Xgrid(maxIdx);
    rfPeakY(intCh) = Ygrid(maxIdx);
    
    % ===== Method 3: 2D Gaussian Fit =====
    try
        [cX, cY, sX, sY, amp] = fitGaussian2D_simple(ratesNorm, Xgrid, Ygrid);
        if ~isnan(cX)
            rfCenterX_Gauss(intCh) = cX;
            rfCenterY_Gauss(intCh) = cY;
            rfSigmaX(intCh) = sX;
            rfSigmaY(intCh) = sY;
        end
    catch
        % Gaussian fit failed, keep as NaN
    end
    
    % ===== RF Quality Metrics =====
    % Area above half-max
    threshold = 0.5 * max(ratesNorm(:));
    rfArea(intCh) = sum(ratesNorm(:) >= threshold);
    
    % Signal-to-noise ratio
    noise = std(sortedRates(1:ceil(0.25*length(sortedRates))));
    if noise > 0
        rfSNR(intCh) = maxRate / noise;
    end
end

fprintf('RF centers: %d Gaussian fits, %d CoM estimates, %d total\n', ...
    sum(~isnan(rfCenterX_Gauss)), sum(~isnan(rfCenterX_CoM)), ...
    sum(~isnan(rfCenterX_CoM) | ~isnan(rfCenterX_Gauss)));

%% 8. Quality-weighted RF visualizations

shankIdx = shankInd + 1; % convert 0-based to 1-based
shankColors = lines(nShank);

% Quality threshold
minSNR = 3;
validChans = find(~isnan(rfCenterX_CoM) & rfSNR >= minSNR);

%% 8a. CoM Centers - Size by SNR, Color by Amplitude
figure('Name', 'RF Centers - CoM (Quality Weighted)', 'Color', 'w'); hold on;

for intCh = validChans'
    if isnan(rfCenterX_CoM(intCh)), continue; end
    
    % Size by SNR
    snrNorm = min(rfSNR(intCh), 20) / 20;
    mSize = 6 + snrNorm * 14;
    
    % Color intensity by amplitude (per shank)
    currentShank = shankIdx(intCh);
    shankmask = (shankIdx(:) == currentShank);
    localMaxAmp = max(rfAmplitude(shankmask & ~isnan(rfAmplitude)));
    if isempty(localMaxAmp) || localMaxAmp == 0
        ampNorm = 0.5;
    else
        ampNorm = rfAmplitude(intCh) / localMaxAmp;
    end
    
    % Blend color (darker = stronger response)
    baseColor = shankColors(currentShank,:);
    color = baseColor * (0.3 + 0.7 * ampNorm);
    
    plot(rfCenterX_CoM(intCh), rfCenterY_CoM(intCh), 'o', ...
        'MarkerEdgeColor', baseColor, ...
        'MarkerFaceColor', color, ...
        'MarkerSize', mSize, 'LineWidth', 1.5);
end

xlabel('RF X (px)'); ylabel('RF Y (px)');
title(sprintf('RF Centers - CoM (n=%d, SNR≥%.1f) | Size=SNR, Brightness=Amplitude', ...
    sum(~isnan(rfCenterX_CoM(validChans))), minSNR));
axis equal; grid on;

% Legend
legendEntries = arrayfun(@(s) sprintf('Shank %d', s), 1:nShank, 'UniformOutput', false);
legend(legendEntries, 'Location', 'best');

%% 8b. Gaussian Centers - Size by SNR, Color by Amplitude
validGauss = validChans(~isnan(rfCenterX_Gauss(validChans)));

figure('Name', 'RF Centers - Gaussian (Quality Weighted)', 'Color', 'w'); hold on;

for intCh = validGauss'
    % Size by SNR
    snrNorm = min(rfSNR(intCh), 20) / 20;
    mSize = 6 + snrNorm * 14;
    
    % Color intensity by amplitude (per shank)
    currentShank = shankIdx(intCh);
    shankmask = (shankIdx(:) == currentShank);
    localMaxAmp = max(rfAmplitude(shankmask & ~isnan(rfAmplitude)));
    if isempty(localMaxAmp) || localMaxAmp == 0
        ampNorm = 0.5;
    else
        ampNorm = rfAmplitude(intCh) / localMaxAmp;
    end
    
    % Blend color
    baseColor = shankColors(currentShank,:);
    color = baseColor * (0.3 + 0.7 * ampNorm);
    
    plot(rfCenterX_Gauss(intCh), rfCenterY_Gauss(intCh), 'o', ...
        'MarkerEdgeColor', baseColor, ...
        'MarkerFaceColor', color, ...
        'MarkerSize', mSize, 'LineWidth', 1.5);
    
    % Add Gaussian extent ellipse
    if ~isnan(rfSigmaX(intCh))
        theta = linspace(0, 2*pi, 50);
        ellipseX = rfCenterX_Gauss(intCh) + 2*rfSigmaX(intCh)*cos(theta);
        ellipseY = rfCenterY_Gauss(intCh) + 2*rfSigmaY(intCh)*sin(theta);
        plot(ellipseX, ellipseY, '-', 'Color', [baseColor, 0.3], 'LineWidth', 1);
    end
end

xlabel('RF X (px)'); ylabel('RF Y (px)');
title(sprintf('RF Centers - Gaussian (n=%d, SNR≥%.1f) | Size=SNR, Brightness=Amplitude', ...
    length(validGauss), minSNR));
axis equal; grid on;
legend(legendEntries, 'Location', 'best');

%% 8c. Comparison: CoM vs Gaussian
figure('Name', 'CoM vs Gaussian Comparison', 'Color', 'w'); hold on;

for intCh = validGauss'
    currentShank = shankIdx(intCh);
    
    % Plot CoM
    plot(rfCenterX_CoM(intCh), rfCenterY_CoM(intCh), 's', ...
        'MarkerEdgeColor', shankColors(currentShank,:), ...
        'MarkerFaceColor', 'none', 'MarkerSize', 8, 'LineWidth', 1.5);
    
    % Plot Gaussian
    plot(rfCenterX_Gauss(intCh), rfCenterY_Gauss(intCh), 'o', ...
        'MarkerEdgeColor', shankColors(currentShank,:), ...
        'MarkerFaceColor', shankColors(currentShank,:), ...
        'MarkerSize', 8, 'LineWidth', 1.5);
    
    % Connect with line
    plot([rfCenterX_CoM(intCh), rfCenterX_Gauss(intCh)], ...
         [rfCenterY_CoM(intCh), rfCenterY_Gauss(intCh)], '-', ...
         'Color', [shankColors(currentShank,:), 0.3], 'LineWidth', 1);
end

xlabel('RF X (px)'); ylabel('RF Y (px)');
title(sprintf('CoM (□) vs Gaussian (○) Centers (n=%d)', length(validGauss)));
axis equal; grid on;
legend(legendEntries, 'Location', 'best');

%% 8d. Combined View: CoM + Gaussian + Hotspots
figure('Name', 'RF Centers + Hotspots (Combined)', 'Color', 'w'); hold on;

% Get pixel dimensions
if numel(vecX_pix) > 1
    dx = mean(diff(vecX_pix));
else
    dx = 1;
end
if numel(vecY_pix) > 1
    dy = mean(diff(vecY_pix));
else
    dy = 1;
end

% Plot hotspots first (background)
for intCh = validChans'
    currentShank = shankIdx(intCh);
    rates = squeeze(matAvgRespAll(:,:,intCh));
    if all(isnan(rates(:))), continue; end
    
    mu = mean(rates(:), 'omitnan');
    sigma = std(rates(:), 'omitnan');
    if sigma == 0, continue; end
    zMap = (rates - mu) / sigma;
    
    mask = zMap >= 2;
    if ~any(mask(:)), continue; end
    
    [yIdx, xIdx] = find(mask);
    for k = 1:length(xIdx)
        x = vecX_pix(xIdx(k)) - dx/2;
        y = vecY_pix(yIdx(k)) - dy/2;
        
        rectangle('Position', [x, y, dx, dy], ...
            'FaceColor', shankColors(currentShank,:), ...
            'EdgeColor', 'none', 'FaceAlpha', 0.15);
    end
end

% Plot CoM centers
for intCh = validChans'
    if isnan(rfCenterX_CoM(intCh)), continue; end
    
    currentShank = shankIdx(intCh);
    snrNorm = min(rfSNR(intCh), 20) / 20;
    mSize = 6 + snrNorm * 10;
    
    plot(rfCenterX_CoM(intCh), rfCenterY_CoM(intCh), 's', ...
        'MarkerEdgeColor', shankColors(currentShank,:), ...
        'MarkerFaceColor', 'none', ...
        'MarkerSize', mSize, 'LineWidth', 2);
end

% Plot Gaussian centers on top
for intCh = validGauss'
    currentShank = shankIdx(intCh);
    snrNorm = min(rfSNR(intCh), 20) / 20;
    mSize = 8 + snrNorm * 12;
    
    plot(rfCenterX_Gauss(intCh), rfCenterY_Gauss(intCh), 'o', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', shankColors(currentShank,:), ...
        'MarkerSize', mSize, 'LineWidth', 2);
end

xlabel('RF X (px)'); ylabel('RF Y (px)');
title(sprintf('RF Analysis (SNR≥%.1f) | □=CoM, ○=Gauss, Shaded=Hotspots', minSNR));
axis equal; grid on;
legend(legendEntries, 'Location', 'best');

%% 8e. Per-Shank Detail Plots
for sh = 1:nShank
    channelsInShank = find(shankIdx == sh & ~isnan(rfCenterX_CoM) & rfSNR >= minSNR);
    if isempty(channelsInShank), continue; end
    
    figure('Name', sprintf('Shank %d - RF Details', sh), 'Color', 'w');
    
    % Subplot 1: Hotspots + CoM
    subplot(1,2,1); hold on;
    
    % Plot hotspots
    for intCh = channelsInShank'
        rates = squeeze(matAvgRespAll(:,:,intCh));
        if all(isnan(rates(:))), continue; end
        
        mu = mean(rates(:), 'omitnan');
        sigma = std(rates(:), 'omitnan');
        if sigma == 0, continue; end
        zMap = (rates - mu) / sigma;
        
        mask = zMap >= 2;
        if ~any(mask(:)), continue; end
        
        [yIdx, xIdx] = find(mask);
        for k = 1:length(xIdx)
            x = vecX_pix(xIdx(k)) - dx/2;
            y = vecY_pix(yIdx(k)) - dy/2;
            
            zVal = zMap(yIdx(k), xIdx(k));
            intensity = min((zVal - 2) / 3, 1);
            color = shankColors(sh,:) * (0.3 + 0.7 * intensity);
            
            rectangle('Position', [x, y, dx, dy], ...
                'FaceColor', color, 'EdgeColor', 'none');
        end
    end
    
    % Plot CoM centers
    for intCh = channelsInShank'
        if isnan(rfCenterX_CoM(intCh)), continue; end
        
        snrNorm = min(rfSNR(intCh), 20) / 20;
        mSize = 10 + snrNorm * 15;
        
        plot(rfCenterX_CoM(intCh), rfCenterY_CoM(intCh), 's', ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', shankColors(sh,:), ...
            'MarkerSize', mSize, 'LineWidth', 2);
        
        % Label high-SNR channels
        if rfSNR(intCh) >= 8
            text(rfCenterX_CoM(intCh), rfCenterY_CoM(intCh), sprintf('%d', intCh), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', 8, 'FontWeight', 'bold', 'Color', 'w');
        end
    end
    
    xlabel('RF X (px)'); ylabel('RF Y (px)');
    title(sprintf('CoM Centers (n=%d)', length(channelsInShank)));
    axis equal; grid on;
    
    % Subplot 2: Gaussian centers with ellipses
    subplot(1,2,2); hold on;
    
    gaussChans = channelsInShank(~isnan(rfCenterX_Gauss(channelsInShank)));
    
    % Plot Gaussian ellipses
    for intCh = gaussChans'
        if ~isnan(rfSigmaX(intCh))
            theta = linspace(0, 2*pi, 50);
            ellipseX = rfCenterX_Gauss(intCh) + 2*rfSigmaX(intCh)*cos(theta);
            ellipseY = rfCenterY_Gauss(intCh) + 2*rfSigmaY(intCh)*sin(theta);
            
            % Fill ellipse
            fill(ellipseX, ellipseY, shankColors(sh,:), ...
                'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
    end
    
    % Plot Gaussian centers
    for intCh = gaussChans'
        snrNorm = min(rfSNR(intCh), 20) / 20;
        mSize = 10 + snrNorm * 15;
        
        plot(rfCenterX_Gauss(intCh), rfCenterY_Gauss(intCh), 'o', ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', shankColors(sh,:), ...
            'MarkerSize', mSize, 'LineWidth', 2);
        
        % Label high-SNR channels
        if rfSNR(intCh) >= 8
            text(rfCenterX_Gauss(intCh), rfCenterY_Gauss(intCh), sprintf('%d', intCh), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', 8, 'FontWeight', 'bold', 'Color', 'w');
        end
    end
    
    xlabel('RF X (px)'); ylabel('RF Y (px)');
    title(sprintf('Gaussian Centers (n=%d)', length(gaussChans)));
    axis equal; grid on;
    
    sgtitle(sprintf('Shank %d: RF Analysis | Size=SNR', sh));
end

%% 9. Plot RF quality metrics

figure('Name', 'RF Quality Metrics', 'Color', 'w');

subplot(2,2,1); histogram(rfSNR(~isnan(rfSNR)), 20);
xlabel('SNR'); ylabel('Count'); title('Signal-to-Noise Ratio');

subplot(2,2,2); histogram(rfArea(~isnan(rfArea)), 20);
xlabel('Area (pixels)'); ylabel('Count'); title('RF Area (half-max)');

subplot(2,2,3); histogram(rfAmplitude(~isnan(rfAmplitude)), 20);
xlabel('Amplitude (sp/s)'); ylabel('Count'); title('Peak Response');

subplot(2,2,4); 
plot(rfAmplitude, rfSNR, 'ko', 'MarkerFaceColor', [0.5 0.5 0.5]);
xlabel('Peak Amplitude (sp/s)'); ylabel('SNR'); title('Amplitude vs SNR');
grid on;

%% 10. Plot example RF maps per shank (top N channels by SNR)

topN = 12;

for sh = 1:nShank
    channelsInShank = find(shankIdx == sh);
    if isempty(channelsInShank), continue; end

    % Pick top N channels by SNR
    [~, sortIdx] = sort(rfSNR(channelsInShank), 'descend', 'MissingPlacement', 'last');
    topChans = channelsInShank(sortIdx(1:min(topN, numel(sortIdx))));
    
    % Remove channels with NaN SNR
    topChans = topChans(~isnan(rfSNR(topChans)));
    if isempty(topChans), continue; end

    figure('Name', sprintf('RF Examples — Shank %d', sh), 'Color', 'w');
    nRows = 3; nCols = ceil(min(topN, numel(topChans)) / nRows);

    for k = 1:min(numel(topChans), topN)
        intCh = topChans(k);
        subplot(nRows, nCols, k); hold on;
        
        % Plot RF map
        rfMap = squeeze(matAvgRespAll(:,:,intCh));
        imagesc(vecX_pix, vecY_pix, rfMap);
        colormap hot; colorbar; 
        
        % Overlay RF center markers
        if ~isnan(rfCenterX_CoM(intCh))
            % CoM
            plot(rfCenterX_CoM(intCh), rfCenterY_CoM(intCh), 'cs', ...
                'MarkerSize', 10, 'LineWidth', 2);
        end
        
        if ~isnan(rfCenterX_Gauss(intCh))
            % Gaussian
            plot(rfCenterX_Gauss(intCh), rfCenterY_Gauss(intCh), 'c+', ...
                'MarkerSize', 12, 'LineWidth', 2);
        end
        
        % Peak location
        plot(rfPeakX(intCh), rfPeakY(intCh), 'gx', ...
            'MarkerSize', 10, 'LineWidth', 2);
        
        xlabel('X (pix)'); ylabel('Y (pix)'); 
        axis equal; axis tight;
        title(sprintf('Ch %d | SNR=%.1f', intCh, rfSNR(intCh)));
    end
    
    % Add legend to first subplot
    subplot(nRows, nCols, 1);
    legend({'', 'CoM', 'Gauss', 'Peak'}, 'TextColor', 'w', 'Location', 'best');
end

%% 11. Save results

save('RF_analysis_results.mat', 'matAvgRespAll', ...
    'rfCenterX_CoM', 'rfCenterY_CoM', ...
    'rfCenterX_Gauss', 'rfCenterY_Gauss', 'rfPeakX', 'rfPeakY', ...
    'rfSigmaX', 'rfSigmaY', 'rfAmplitude', 'rfBaseline', 'rfArea', 'rfSNR', ...
    'vecX_pix', 'vecY_pix', 'shankInd', 'xCoord', 'yCoord');

fprintf('Analysis complete. Results saved to RF_analysis_results.mat\n');
fprintf('Generated figures: 4 overview + %d per-shank + %d RF examples = %d total\n', ...
    nShank, nShank, 4 + 2*nShank);

%% ========== Helper Function ==========

function [centerX, centerY, sigmaX, sigmaY, amplitude] = fitGaussian2D_simple(rates, Xgrid, Ygrid)
    % Simple 2D Gaussian fitting (circular, no rotation)
    
    % Initial guess from data
    [maxVal, maxIdx] = max(rates(:));
    x0 = Xgrid(maxIdx);
    y0 = Ygrid(maxIdx);
    sigma0 = mean([range(Xgrid(:)), range(Ygrid(:))]) / 4;
    baseline = min(rates(:));
    
    % Define circular Gaussian model
    gauss2D = @(p, X, Y) p(4) + p(5) * exp(-((X - p(1)).^2 + (Y - p(2)).^2) / (2 * p(3)^2));
    
    % Initial parameters: [x0, y0, sigma, baseline, amplitude]
    params0 = [x0, y0, sigma0, baseline, maxVal - baseline];
    
    % Lower and upper bounds
    lb = [min(Xgrid(:)), min(Ygrid(:)), 0.5*sigma0, 0, 0];
    ub = [max(Xgrid(:)), max(Ygrid(:)), 2*sigma0, maxVal, 2*maxVal];
    
    % Fit using lsqcurvefit
    options = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxIterations', 100);
    
    try
        % Create coordinate arrays
        coords = [Xgrid(:), Ygrid(:)];
        ratesVec = rates(:);
        
        % Fit
        params = lsqcurvefit(@(p, coords) gauss2D(p, coords(:,1), coords(:,2)), ...
            params0, coords, ratesVec, lb, ub, options);
        
        % Extract parameters
        centerX = params(1);
        centerY = params(2);
        sigmaX = abs(params(3));
        sigmaY = abs(params(3)); % circular
        amplitude = params(5);
        
        % Quality check: center should be within grid bounds
        if centerX < min(Xgrid(:)) || centerX > max(Xgrid(:)) || ...
           centerY < min(Ygrid(:)) || centerY > max(Ygrid(:))
            centerX = NaN;
            centerY = NaN;
        end
        
    catch
        % If fitting fails, return NaN
        centerX = NaN;
        centerY = NaN;
        sigmaX = NaN;
        sigmaY = NaN;
        amplitude = NaN;
    end
end