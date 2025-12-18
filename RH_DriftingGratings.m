% RH_DriftingGratings
% Shows full-field sinusoidal drifting gratings
% 
% Robin Haak, 2025

%% suppress m-lint warnings
%#ok<*MCCD,*NASGU,*ASGLU,*CTCH>
clear; close all; Screen('CloseAll');

%% define variables
fprintf('Starting %s [%s]\n',mfilename,getTime);
boolUseSGL = false;
boolUseNI = false;
boolDebug = false;
dblPupilLightMultiplier = 1;
dblSyncLightMultiplier = 0.5;
strHostAddress = '192.87.11.8';
objDaqOut = [];
if exist('sExpMeta','var')
    if isfield(sExpMeta,'dblPupilLightMultiplier'),dblPupilLightMultiplier=sExpMeta.dblPupilLightMultiplier;end
    if isfield(sExpMeta,'dblSyncLightMultiplier'),dblSyncLightMultiplier=sExpMeta.dblSyncLightMultiplier;end
    if isfield(sExpMeta,'strHostAddress'),strHostAddress=sExpMeta.strHostAddress;end
    if isfield(sExpMeta,'objDaqOut'),objDaqOut=sExpMeta.objDaqOut;end
    if isfield(sExpMeta,'boolUseSGL'),boolUseSGL=sExpMeta.boolUseSGL;end
    if isfield(sExpMeta,'boolUseNI'),boolUseNI=sExpMeta.boolUseNI;end
else
    sExpMeta = [];
end

%% query user input for recording name
if exist('sStimParamsSettings','var') && isfield(sStimParamsSettings,'strRecording')
    strRecording = sStimParamsSettings.strRecording;
else
    strRecording = input('Recording name (e.g., MouseX): ', 's');
end

%% input params (use your defaults if nothing provided)
fprintf('Loading settings...\n');
if ~exist('sStimParamsSettings','var') || isempty(sStimParamsSettings) || ...
        ~(isfield(sStimParamsSettings,'strStimType') && (strcmpi(sStimParamsSettings.strStimType,'SquareGrating') || strcmpi(sStimParamsSettings.strStimType,'SineGrating')))
    sStimParamsSettings = struct;
    sStimParamsSettings.strStimType = 'GratingFullField';
    sStimParamsSettings.strOutputPath = 'C:\_Data\Exp';
    sStimParamsSettings.strTempObjectPath = 'C:\temp';
    sStimParamsSettings.dblScreenDistance_cm = 23;
    sStimParamsSettings.intCornerTrigger = 2;
    sStimParamsSettings.dblCornerSize = 1/30;
    sStimParamsSettings.dblScreenWidth_cm = 51;
    sStimParamsSettings.dblScreenHeight_cm = 29;
    sStimParamsSettings.dblScreenWidth_deg = 2 * atand(sStimParamsSettings.dblScreenWidth_cm / (2 * sStimParamsSettings.dblScreenDistance_cm));
    sStimParamsSettings.dblScreenHeight_deg = 2 * atand(sStimParamsSettings.dblScreenHeight_cm / (2 * sStimParamsSettings.dblScreenDistance_cm));
    sStimParamsSettings.intUseScreen = 2;
    sStimParamsSettings.dblBackground = 0.5;
    sStimParamsSettings.intBackground = round(mean(sStimParamsSettings.dblBackground)*255);
   
    %stimulus parameters
    sStimParamsSettings.dblSpatialFrequency_cd = 0.06; % cycles/deg default
    sStimParamsSettings.dblTemporalFrequency = 1; % Hz default, gives 25? deg/s
    sStimParamsSettings.dblSecsDuration = 2; % seconds per stimulus
    sStimParamsSettings.dblSecsInitialBlank = 5;
    sStimParamsSettings.dblSecsPreBlank = 0.35;
    sStimParamsSettings.dblSecsPostBlank = 0.15;
    sStimParamsSettings.dblSecsEndBlank = 5;
    sStimParamsSettings.intStimulusRepeats = 25;
    sStimParamsSettings.vecDirections = 0:45:315; % 8 directions
    sStimParamsSettings.str90Deg = '0 degrees is rightward motion; 90 degrees is upward motion';
    sStimParamsSettings.intUseDaqDevice = 1;
    sStimParamsSettings.intUseParPool = 0;
    sStimParamsSettings.intUseGPU = 1;
else
    % evaluate fields if they contain expressions
    cellFields = fieldnames(sStimParamsSettings);
    for intField=1:numel(cellFields)
        try
            sStimParamsSettings.(cellFields{intField}) = eval(sStimParamsSettings.(cellFields{intField}));
        catch
            sStimParamsSettings.(cellFields{intField}) = sStimParamsSettings.(cellFields{intField});
        end
    end
end
if boolDebug == 1
    intUseScreen = 0;
else
    intUseScreen = sStimParamsSettings.intUseScreen;
end
sStimParams = sStimParamsSettings;

%% set output locations for logs
strOutputPath = sStimParamsSettings.strOutputPath;
strTempObjectPath = sStimParamsSettings.strTempObjectPath;
strThisFilePath = mfilename('fullpath');
[strFilename,strLogDir,strTempDir,strTexDir] = RE_assertPaths(strOutputPath,strRecording,strTempObjectPath,strThisFilePath);
fprintf('Saving output in directory %s\n',strLogDir);

%% initialize connection with SpikeGLX
if boolUseSGL
    if exist('sExpMeta','var') && isfield(sExpMeta,'hSGL') && isfield(sExpMeta,'strRunName') && isfield(sExpMeta,'sParamsSGL')
        hSGL = sExpMeta.hSGL;
        strRunName = sExpMeta.strRunName;
        sParamsSGL = sExpMeta.sParamsSGL;
        intOutFlag = StartRecordingSGL(hSGL);
    else
        fprintf('Opening SpikeGLX connection & starting recording "%s" [%s]...\n',strRecording,getTime);
        [hSGL,strRunName,sParamsSGL] = InitSGL(strRecording,strHostAddress);
    end
    fprintf('SGL saving to "%s", matlab saving to "%s.mat" [%s]...\n',strRunName,strFilename,getTime);
    intStreamNI = 0;
    dblSampFreqNI = GetStreamSampleRate(hSGL, intStreamNI, strHostAddress);
    strDataDirSGL = GetDataDir(hSGL);
    jFileObj = java.io.File(strDataDirSGL);
    dblFreeGB = (jFileObj.getFreeSpace)/(1024^3);
    if dblFreeGB < 100,warning([mfilename ':LowDiskSpace'],'Low disk space available (%.0fGB) for Neuropixels data (dir: %s)',dblFreeGB,strDataDirSGL);end
else
    sParamsSGL = struct;
end

%% parallel / GPU initialization
if sStimParams.intUseParPool > 0 && isempty(gcp('nocreate'))
    parpool(sStimParams.intUseParPool * [1 1]);
end
if sStimParams.intUseGPU > 0
    objGPU = gpuDevice(sStimParams.intUseGPU);
end

%% NI I/O initialization
if boolUseNI
    fprintf('Connecting to National Instruments box...\n');
    boolDaqOutRunning = false;
    if exist('objDaqOut','var') && ~isempty(objDaqOut)
        try
            stop(objDaqOut);
            outputData1 = dblSyncLightMultiplier*cat(1,linspace(3, 3, 200)',linspace(0, 0, 50)');
            outputData2 = dblPupilLightMultiplier*linspace(3, 3, 250)';
            queueOutputData(objDaqOut,[outputData1 outputData2]);
            prepare(objDaqOut);
            pause(0.1);
            startBackground(objDaqOut)
            boolDaqOutRunning = true;
        catch
        end
    end
    if ~boolDaqOutRunning
        objDaqOut = openDaqOutput(sStimParamsSettings.intUseDaqDevice);
        stop(objDaqOut);
        outputData1 = dblSyncLightMultiplier*cat(1,linspace(3, 3, 200)',linspace(0, 0, 50)');
        outputData2 = dblPupilLightMultiplier*linspace(3, 3, 250)';
        queueOutputData(objDaqOut,[outputData1 outputData2]);
        prepare(objDaqOut);
        pause(0.1);
        startBackground(objDaqOut)
    end
end

try
    %% INITIALIZE SCREEN
    fprintf('Starting PsychToolBox extension...\n');
    AssertOpenGL;
    KbName('UnifyKeyNames');
    intOldVerbosity = Screen('Preference', 'Verbosity',1);
    if boolDebug == 1, vecInitRect = [0 0 640 640];else, vecInitRect = [];end
    try
        Screen('Preference', 'SkipSyncTests', 0);
        [ptrWindow,vecRect] = Screen('OpenWindow', intUseScreen, sStimParams.intBackground, vecInitRect);
    catch ME
        warning([mfilename ':ErrorPTB'],'Psychtoolbox error, attempting with sync test skip [msg: %s]',ME.message);
        Screen('Preference', 'SkipSyncTests', 1);
        [ptrWindow,vecRect] = Screen('OpenWindow', intUseScreen, sStimParams.intBackground, vecInitRect);
    end
    if exist('C:\Users\neuropixels\Desktop\GammaCorrection\gammaTable.mat','file')
        load("C:\Users\neuropixels\Desktop\GammaCorrection\gammaTable.mat");
        Screen('LoadNormalizedGammaTable', ptrWindow, gammaTable*[1 1 1]);
    end
    sStimParams.ptrWindow = ptrWindow;
    sStimParams.vecRect = vecRect;
    sStimParams.intScreenWidth_pix = vecRect(3)-vecRect(1);
    sStimParams.intScreenHeight_pix = vecRect(4)-vecRect(2);
    sStimParams.intWhite = WhiteIndex(sStimParams.ptrWindow);
    intOldPriority = 0;
    if boolDebug == 0
        intPriorityLevel=MaxPriority(ptrWindow);
        intOldPriority = Priority(intPriorityLevel);
    end
    dblStimFrameRate=Screen('FrameRate', ptrWindow);
    intStimFrameRate = round(dblStimFrameRate);
    dblStimFrameDur = mean(1/dblStimFrameRate);
    dblInterFlipInterval = Screen('GetFlipInterval', ptrWindow);
    if dblStimFrameDur/dblInterFlipInterval > 1.05 || dblStimFrameDur/dblInterFlipInterval < 0.95
        warning([mfilename ':InconsistentFlipDur'],'Something iffy with flip speed and monitor refresh rate detected; frame duration is %fs, while flip interval is %fs!',dblStimFrameDur,dblInterFlipInterval);
    end

    %% check esscape
    if CheckEsc(),error([mfilename ':EscapePressed'],'Esc pressed; exiting');end

    %% prepare stimuli
    % parameters
    numDirections = length(sStimParams.vecDirections);
    spatialFreq_cpd = sStimParams.dblSpatialFrequency_cd; % cycles/deg
    temporalFreq_hz = sStimParams.dblTemporalFrequency;   % Hz
    stimDuration = sStimParams.dblSecsDuration;           % seconds per stimulus
    initialBlank = sStimParams.dblSecsInitialBlank;
    preBlank = sStimParams.dblSecsPreBlank;
    postBlank = sStimParams.dblSecsPostBlank;
    endBlank = sStimParams.dblSecsEndBlank;
    repeats = sStimParams.intStimulusRepeats;

    %pixel conversion
    sStimParams.dblPixelsPerDeg = mean([(sStimParams.intScreenWidth_pix/sStimParams.dblScreenWidth_deg) ...
        (sStimParams.intScreenHeight_pix/sStimParams.dblScreenHeight_deg)]);
    pixelsPerCycle = sStimParams.dblPixelsPerDeg / spatialFreq_cpd; % pixels per grating cycle

    %build a texture slightly larger than screen (pad by a couple cycles to allow shifting)
    padCycles = 5; % number of cycles padding around screen to avoid edge artifacts
    texW = ceil(sStimParams.intScreenWidth_pix + ceil(padCycles * pixelsPerCycle));
    texH = ceil(sStimParams.intScreenHeight_pix + ceil(padCycles * pixelsPerCycle));
    %make even
    if mod(texW,2)~=0, texW = texW+1; end
    if mod(texH,2)~=0, texH = texH+1; end

    %create coordinate grids centered on zero (pixels)
    halfW = texW/2; halfH = texH/2;
    [Xpix,Ypix] = meshgrid( (1:texW) - halfW - 0.5, (1:texH) - halfH - 0.5 );

    % Pre-generate one base grating for phase=0 and rotate via sampling
    % (we build oriented gratings separately to avoid expensive rotation per frame)
    texCollection = zeros(numDirections,1);

    fprintf('Generating full-field gratings for %d directions (tex %dx%d, pixelsPerCycle=%.2f)...\n', numDirections, texW, texH, pixelsPerCycle);
    for d = 1:numDirections
        dirDeg = sStimParams.vecDirections(d); % 0=rightward
        theta_rad = deg2rad(mod(180-dirDeg,360));

        % compute coordinate along drift direction (pixels)
        coord_along = (Xpix * cos(theta_rad) + Ypix * sin(theta_rad)) / sStimParams.dblPixelsPerDeg; % degrees
        gr = 0.5 + 0.5 * sin(2*pi*spatialFreq_cpd .* coord_along); % base phase 0

        texIm = uint8(round(gr * 255));
        texCollection(d) = Screen('MakeTexture', ptrWindow, texIm);
    end
    fprintf('Texture generation complete.\n');

    
    %% build trial sequence: random order of directions with repeats
    vecDirSeq = repmat(1:numDirections, 1, repeats);
    vecDirSeq = vecDirSeq(randperm(length(vecDirSeq)));
    intTotalTrials = length(vecDirSeq);

    % prepare structEP (timing and logs)
    structEP = struct;
    structEP.strExpType = sStimParams.strStimType;
    structEP.intTrialNum = intTotalTrials;
    structEP.TrialNumber = nan(1,structEP.intTrialNum);
    structEP.dblStimFrameDur = dblStimFrameDur;
    structEP.dblInterFlipInterval = dblInterFlipInterval;
    structEP.ActOnSecs = nan(1,structEP.intTrialNum);
    structEP.ActOffSecs = nan(1,structEP.intTrialNum);
    structEP.ActStartSecs = nan(1,structEP.intTrialNum);
    structEP.ActStopSecs = nan(1,structEP.intTrialNum);
    structEP.ActOnNI = nan(1,structEP.intTrialNum);
    structEP.ActOffNI = nan(1,structEP.intTrialNum);
    structEP.vecDirection = nan(1,structEP.intTrialNum);
    structEP.vecTexturePhase = nan(1,structEP.intTrialNum);
   
    %compute trial timing in absolute seconds relative to experiment start
    trialDur = preBlank + stimDuration + postBlank;
    structEP.vecTrialStartSecs = initialBlank:trialDur:(initialBlank + trialDur*(intTotalTrials-1));
    structEP.vecTrialStimOnSecs = structEP.vecTrialStartSecs + preBlank;
    structEP.vecTrialStimOffSecs = structEP.vecTrialStimOnSecs + stimDuration;
    structEP.vecTrialEndSecs = structEP.vecTrialStimOffSecs + postBlank;

    %% present stimuli: wait for user
    opts = struct;
    opts.Default = 'Start';
    opts.Interpreter = 'tex';
    strAns = questdlg('Would you like to start the stimulation?', 'Start Stimulation', 'Start','Cancel',opts);
    if ~strcmp(strAns,opts.Default)
        error([mfilename ':RunCancelled'],'Cancelling');
    end

    % hide cursor
    %HideCursor(ptrWindow);

    % initial flips and blanking
    hTic = tic;
    dblLastFlip = Screen('Flip', ptrWindow);
    dblInitialFlip = dblLastFlip;
    structEP.strStartDate = getDate();
    structEP.strStartTime = getTime();

    % calculate photodiode rect if requested
    if sStimParams.intCornerTrigger > 0
        intCornerPix = floor(sStimParams.dblCornerSize * sStimParams.intScreenWidth_pix);
        if sStimParams.intCornerTrigger == 1
            vecDiodeRect = [0 0 intCornerPix intCornerPix];
        elseif sStimParams.intCornerTrigger == 2
            vecDiodeRect = [(vecRect(3)-intCornerPix) 0 vecRect(3) intCornerPix];
        elseif sStimParams.intCornerTrigger == 3
            vecDiodeRect = [0 (vecRect(4)-intCornerPix) intCornerPix vecRect(4)];
        elseif sStimParams.intCornerTrigger == 4
            vecDiodeRect = [(vecRect(3)-intCornerPix) (vecRect(4)-intCornerPix) vecRect(3) vecRect(4)];
        end
    else
        vecDiodeRect = [];
    end

    %% initial blank
    fprintf('Starting initial blank (dur=%.2fs) [%s]\n', initialBlank, getTime);
    dblInitialBlankDur = 0;
    while dblInitialBlankDur < (initialBlank - dblStimFrameDur)
        Screen('FillRect', ptrWindow, sStimParams.intBackground);
        dblLastFlip = Screen('Flip', ptrWindow, dblLastFlip + dblStimFrameDur/2);
        dblInitialBlankDur = dblLastFlip - dblInitialFlip;
    end

    %% MAIN TRIAL LOOP
    intThisTrial = 1;
    while intThisTrial <= intTotalTrials && ~CheckEsc()
        if CheckPause()
            fprintf('\n\n<strong>STIMULUS SCRIPT PAUSED!</strong> [%s]\n\n',getTime);
            WaitSecs(1);
            opts = struct; opts.Default = 'Restart'; opts.Interpreter = 'tex';
            strAns = questdlg('STIMULUS SCRIPT PAUSED! Would you like to restart the stimulation?', 'Restart Stimulation','Restart',opts);
            WaitSecs(1);
        end

        % trial start (background)
        Screen('FillRect', ptrWindow, sStimParams.intBackground);
        dblTrialStartFlip = Screen('Flip', ptrWindow);

        % DAQ fill if used (preserved)
        if boolUseNI
            stop(objDaqOut);
            outputData1 = dblSyncLightMultiplier*cat(1,linspace(3, 3, 200)',linspace(0, 0, 50)');
            outputData2 = dblPupilLightMultiplier*linspace(3, 3, 250)';
            queueOutputData(objDaqOut,[outputData1 outputData2]);
            prepare(objDaqOut);
        end

        % select direction index for this trial
        dirIdx = vecDirSeq(intThisTrial);
        currentDirection = sStimParams.vecDirections(dirIdx);

        % random phase for this trial (radians)
        thisPhase_rad = 2*pi*rand();

        % log direction and phase
        structEP.vecDirection(intThisTrial) = currentDirection;
        structEP.vecTexturePhase(intThisTrial) = thisPhase_rad;

        % timing references
        dblStartSecs = structEP.vecTrialStartSecs(intThisTrial);
        dblStimOnSecs = structEP.vecTrialStimOnSecs(intThisTrial);
        dblStimOffSecs = structEP.vecTrialStimOffSecs(intThisTrial);
        dblStimDurSecs = dblStimOffSecs - dblStimOnSecs;
        dblEndSecs = structEP.vecTrialEndSecs(intThisTrial);

        fprintf('%d/%d direction: %d deg, phase(rad)=%.2f [%s]\n', intThisTrial, intTotalTrials, currentDirection, thisPhase_rad, getTime);

        % pre-blank
        dblPreBlankDur = 0;
        Screen('FillRect', ptrWindow, sStimParams.intBackground);
        dblLastFlip = Screen('Flip', ptrWindow);
        while dblPreBlankDur < (dblStimOnSecs - dblStartSecs - dblStimFrameDur/2)
            Screen('FillRect', ptrWindow, sStimParams.intBackground);
            dblLastFlip = Screen('Flip', ptrWindow, dblLastFlip + dblStimFrameDur/2);
            dblPreBlankDur = dblLastFlip - dblTrialStartFlip;
        end

        % 250ms DAQ pulse at stim start
        if boolUseNI, startBackground(objDaqOut); end

        % start stimulus presentation
        refTime = tic;
        boolFirstFlip = false;

        % choose pre-generated texture for this direction
        texId = texCollection(dirIdx);

        % compute drift speed in pixels/s
        driftSpeed_deg_s = temporalFreq_hz / spatialFreq_cpd; % deg/sec
        driftSpeed_pix_s = driftSpeed_deg_s * sStimParams.dblPixelsPerDeg;

        % central sampling window in texture coordinates (centered)
        dstRect = [0 0 sStimParams.intScreenWidth_pix sStimParams.intScreenHeight_pix];
        texCenterX = texW/2;
        texCenterY = texH/2;
        halfScreenW = sStimParams.intScreenWidth_pix/2;
        halfScreenH = sStimParams.intScreenHeight_pix/2;

        while toc(refTime) < (dblStimDurSecs - 2*dblStimFrameDur)
            elapsed = toc(refTime);
            % shift magnitude along drift vector (wrap by pixelsPerCycle to keep phase correct)
            shift_pix = mod(driftSpeed_pix_s * elapsed, pixelsPerCycle);

            % convert scalar shift to x,y components (0deg = right)
            dx = -cosd(currentDirection) * shift_pix;
            dy = sind(currentDirection) * shift_pix;

            % src rectangle sampling (centered then shifted)
            srcX0 = texCenterX - halfScreenW + dx;
            srcY0 = texCenterY - halfScreenH + dy;
            srcRect = [srcX0, srcY0, srcX0 + sStimParams.intScreenWidth_pix, srcY0 + sStimParams.intScreenHeight_pix];

            % apply trial phase shift (phase -> pixel shift)
            phasePixShift = (thisPhase_rad / (2*pi)) * pixelsPerCycle;
            srcRect = srcRect + [phasePixShift, phasePixShift, phasePixShift, phasePixShift];

            % ensure srcRect stays within texture bounds (clamp)
            srcRect = max(srcRect, 1);
            srcRect(3) = min(srcRect(3), texW);
            srcRect(4) = min(srcRect(4), texH);
            % if clamping changed size, re-center to maintain width/height
            srcW = srcRect(3)-srcRect(1); srcH = srcRect(4)-srcRect(2);
            if abs(srcW - sStimParams.intScreenWidth_pix) > 0.5 || abs(srcH - sStimParams.intScreenHeight_pix) > 0.5
                srcRect = [texCenterX - halfScreenW, texCenterY - halfScreenH, texCenterX + halfScreenW, texCenterY + halfScreenH];
            end

            Screen('DrawTexture', ptrWindow, texId, srcRect, dstRect, 0);

            if sStimParams.intCornerTrigger > 0
                Screen('FillRect', ptrWindow, sStimParams.intWhite, vecDiodeRect);
            end
            Screen('DrawingFinished', ptrWindow);
            dblLastFlip = Screen('Flip', ptrWindow, dblLastFlip + dblInterFlipInterval/2);

            % log first flip timestamps
            if ~boolFirstFlip
                boolFirstFlip = true;
                if boolUseSGL
                    dblStimOnNI = GetStreamSampleCount(hSGL, intStreamNI, strHostAddress)/dblSampFreqNI;
                else
                    dblStimOnNI = nan;
                end
                dblStimOnFlip = dblLastFlip;
            end

        end % end stim draw loop

        % back to background
        Screen('FillRect', ptrWindow, sStimParams.intBackground);
        dblStimOffFlip = Screen('Flip', ptrWindow, dblLastFlip + dblStimFrameDur/2);

        if boolUseSGL
            dblStimOffNI = GetStreamSampleCount(hSGL, intStreamNI, strHostAddress)/dblSampFreqNI;
        else
            dblStimOffNI = nan;
        end

        % post-blank
        dblTrialDur = 0;
        Screen('FillRect', ptrWindow, sStimParams.intBackground);
        dblLastFlip = Screen('Flip', ptrWindow);
        while dblTrialDur < (dblEndSecs - dblStartSecs - 2*dblStimFrameDur)
            Screen('FillRect', ptrWindow, sStimParams.intBackground);
            dblLastFlip = Screen('Flip', ptrWindow, dblLastFlip + dblStimFrameDur/2);
            dblTrialDur = dblLastFlip - dblTrialStartFlip;
        end

        % log into structEP
        intStimNumber = intThisTrial;
        structEP.TrialNumber(intStimNumber) = intThisTrial;
        structEP.ActStartSecs(intStimNumber) = dblTrialStartFlip;
        structEP.ActOnSecs(intStimNumber) = dblStimOnFlip;
        structEP.ActOffSecs(intStimNumber) = dblStimOffFlip;
        structEP.ActStopSecs(intStimNumber) = dblLastFlip;
        structEP.ActOnNI(intStimNumber) = dblStimOnNI;
        structEP.ActOffNI(intStimNumber) = dblStimOffNI;

        % increment
        intThisTrial = intThisTrial + 1;
    end % end all trials

    % save results
    structEP.sStimParams = sStimParams;
    save(fullfile(strLogDir,strFilename), 'structEP','sParamsSGL');

    fprintf('Finished experiment & data saving at [%s], waiting for end blank (dur=%.2fs)\n',getTime,endBlank);
    % end blank
    dblEndBlankDur = 0;
    while dblEndBlankDur < endBlank
        Screen('FillRect', ptrWindow, sStimParams.intBackground);
        dblEndFlip = Screen('Flip', ptrWindow);
        dblEndBlankDur = dblEndFlip - dblLastFlip;
    end

    % cleanup textures
    for d = 1:numDirections
        if texCollection(d) ~= 0, Screen('Close', texCollection(d)); end
    end
    Screen('Close', ptrWindow);
    Screen('Close');
    Screen('CloseAll');
    ShowCursor;
    Priority(0);
    Screen('Preference', 'Verbosity', intOldVerbosity);

    if boolUseNI && ~(exist('sExpMeta','var') && isfield(sExpMeta,'objDaqOut'))
        try closeDaqOutput(objDaqOut); catch, end
    end

catch ME
    % handle escape gracefully
    if exist('ME','var') && strcmp(ME.identifier,[mfilename ':EscapePressed'])
        fprintf('\nEscape pressed at [%s], closing down and cleaning up...\n',getTime);
        structEP.sStimParams = sStimParams;
        save(fullfile(strLogDir,strFilename), 'structEP','sParamsSGL');
        if exist('texCollection','var')
            for d = 1:numel(texCollection)
                try Screen('Close', texCollection(d)); catch, end
            end
        end
        Screen('Close', ptrWindow);
        Screen('Close');
        Screen('CloseAll');
        ShowCursor;
        Priority(0);
        Screen('Preference', 'Verbosity', intOldVerbosity);
        if boolUseNI && ~(exist('sExpMeta','var') && isfield(sExpMeta,'objDaqOut'))
            try closeDaqOutput(objDaqOut); catch, end
        end
    else
        % generic error: try to save and cleanup
        fprintf('\n\n\nError occurred! Trying to save data and clean up...\n\n\n');
        try
            structEP.sStimParams = sStimParams;
            if ~exist('sParamsSGL','var'), sParamsSGL = []; end
            save(fullfile(strLogDir,strFilename), 'structEP','sParamsSGL');
        catch
        end
        if exist('texCollection','var')
            for d = 1:numel(texCollection)
                try Screen('Close', texCollection(d)); catch, end
            end
        end
        try
            Screen('Close'); Screen('CloseAll'); ShowCursor; Priority(0); Screen('Preference', 'Verbosity', intOldVerbosity); 
        catch
        end
        if boolUseNI && ~(exist('sExpMeta','var') && isfield(sExpMeta,'objDaqOut'))
            try closeDaqOutput(objDaqOut); catch, end
        end
        rethrow(ME);
    end
end
