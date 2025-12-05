% RH_RANDOMDOTKINEMATOGRAM
% Show a full-screen Random Dot Kinematogram
%
% Robin Haak, 2025

%% Suppress m-lint warnings
%#ok<*MCCD,*NASGU,*ASGLU,*CTCH,*AGROW>
clear; close all; Screen('CloseAll');

%% Define variables
fprintf('Starting %s [%s]\n',mfilename,getTime);
boolUseSGL = false; % Set to true if using SpikeGLX
boolUseNI = false;  % Set to true if using National Instruments DAQ
boolDebug = false;  % Set to true for debugging (smaller window, no sync tests)
dblPupilLightMultiplier = 1; % Strength of infrared LEDs
dblSyncLightMultiplier = 0.5;
strHostAddress = '192.87.11.8'; % '192.87.11.133'; % Default host address
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

%% Ask user for recording name
if exist('sStimParamsSettings','var') && isfield(sStimParamsSettings,'strRecording')
    strRecording = sStimParamsSettings.strRecording;
else
    strRecording = input('Recording name (e.g., MouseX): ', 's');
end

%% Parameters
fprintf('Loading settings...\n');
if ~exist('sStimParamsSettings','var') || isempty(sStimParamsSettings) || ~strcmpi(sStimParamsSettings.strStimType,'RandomDotKinematogram')
    sStimParamsSettings = struct;
    sStimParamsSettings.strStimType = 'RandomDotKinematogram';
    sStimParamsSettings.strOutputPath = 'C:\_Data\Exp'; % Appends date
    sStimParamsSettings.strTempObjectPath = 'C:\temp';
    % Visual space parameters
    sStimParamsSettings.dblScreenDistance_cm = 23; % cm; measured [~23]
    % Screen variables
    sStimParamsSettings.intCornerTrigger = 2; % Integer switch; 0=none,1=upper left, 2=upper right, 3=lower left, 4=lower right
    sStimParamsSettings.dblCornerSize = 1/30; % Fraction of screen width
    sStimParamsSettings.dblScreenWidth_cm = 51; % cm; measured [51]
    sStimParamsSettings.dblScreenHeight_cm = 29; % cm; measured [29]
    sStimParamsSettings.dblScreenWidth_deg = 2 * atand(sStimParamsSettings.dblScreenWidth_cm / (2 * sStimParamsSettings.dblScreenDistance_cm));
    sStimParamsSettings.dblScreenHeight_deg = 2 * atand(sStimParamsSettings.dblScreenHeight_cm / (2 * sStimParamsSettings.dblScreenDistance_cm));
    sStimParamsSettings.intUseScreen = 2; %1; % Screen to use
    sStimParamsSettings.dblBackground = 0.5; % Background intensity (dbl, [0 1])
    sStimParamsSettings.intBackground = round(sStimParamsSettings.dblBackground * 255);

    % Dot parameters (common to all sets)
    sStimParamsSettings.dblDotSize_deg = 2; % Dot size in visual degrees
    sStimParamsSettings.dblDotSpeed_degPerSec = 10; % Dot speed, deg/s
    sStimParamsSettings.dblCoveragePercent = 0.12; % Screen coverage by dots
    sStimParamsSettings.str90Deg = '0 degrees is rightward motion; 90 degrees is upward motion';

    % Stimulus set selection (more robust version)
    tuningSetLabel = 'Tuning Set (8 directions, high coherence)';
    coherenceSetLabel = 'Coherence Set (4 directions, 5 coherence levels, 1 lifetime level )';
    lifetimeSetLabel = 'Lifetime Set (4 directions, 1 coherence level, 5 lifetime levels)';

    stimLabels = {tuningSetLabel, coherenceSetLabel, lifetimeSetLabel};

    % Show list dialog with these options
    [intSelection, tfOK] = listdlg( ...
        'PromptString', 'Choose a Stimulus Set:', ...
        'SelectionMode', 'single', ...
        'ListString', stimLabels, ...
        'ListSize', [350 150]);

    if ~tfOK
        fprintf('No stimulus set selected. Exiting script.\n');
        return; % User canceled, exit gracefully
    end

    choice = stimLabels{intSelection};

    switch choice
        case tuningSetLabel
            sStimParamsSettings.vecDirections_deg = 0:45:315; % All 8 directions
            sStimParamsSettings.vecCoherenceLevels = 0.96; % Single high coherence
            sStimParamsSettings.vecDotLifetime_s = 1; % Dot lifetime in seconds
            fprintf('Selected: Tuning Set (8 directions, 1 coherence)\n');
        case coherenceSetLabel
            sStimParamsSettings.vecDirections_deg = 45:45:180; % Four directions
            sStimParamsSettings.vecCoherenceLevels = [0.06 0.12 0.24 0.48 0.96];
            sStimParamsSettings.vecDotLifetime_s = 1; % Dot lifetime in seconds
            fprintf('Selected: Coherence Set (4 directions, 5 coherences, 1 lifetime)\n');
        case lifetimeSetLabel
            sStimParamsSettings.vecDirections_deg = 45:45:180; % Four directions
            sStimParamsSettings.vecCoherenceLevels = 0.48;
            sStimParamsSettings.vecDotLifetime_s = linspace(0.1,1,5); % Dot lifetime in seconds
            fprintf('Selected: LifeTime Set (4 directions, 1 coherence, 5 lifetimes)\n');
        otherwise
            fprintf('Unexpected selection. Exiting.\n');
            return;
    end


    % General timing parameters (common), seconds
    sStimParamsSettings.dblSecsDuration = 2; %
    sStimParamsSettings.dblSecsInitialBlank = 5;
    sStimParamsSettings.dblSecsPreBlank = 0.5;
    sStimParamsSettings.dblSecsPostBlank = 0.5;
    sStimParamsSettings.dblSecsEndBlank = 5;

    % Control variables (common)
    sStimParamsSettings.intNumRepetitionsPerUniqueStim = 25; % Repetitions for each unique stimulus
    sStimParamsSettings.intNumRNGSeeds = 5; % Number of RNG seeds to use
    sStimParamsSettings.intUseDaqDevice = 1; % ID of DAQ device
    sStimParamsSettings.intUseParPool = 0; % Number of workers in parallel pool; [2]
    sStimParamsSettings.intUseGPU = 0; % Set to 1 if you have a GPU and Psychtoolbox supports it (for MakeTexture etc)
else
    % Evaluate and assign pre-defined values to structure if sStimParamsSettings already existed
    cellFields = fieldnames(sStimParamsSettings);
    for intField=1:numel(cellFields)
        try
            sStimParamsSettings.(cellFields{intField}) = eval(sStimParamsSettings.(cellFields{intField}));
        catch
            % Already evaluated or not a string that needs evaluation
        end
    end
end
if boolDebug == 1
    intUseScreen = 0;
else
    intUseScreen = sStimParamsSettings.intUseScreen;
end
sStimParams = sStimParamsSettings;

%% Set output locations for logs
strOutputPath = sStimParamsSettings.strOutputPath;
strTempObjectPath = sStimParamsSettings.strTempObjectPath; % This variable is no longer explicitly used for saving dot files
strThisFilePath = mfilename('fullpath');
[strFilename,strLogDir,strTempDir,~] = RE_assertPaths(strOutputPath,strRecording,strTempObjectPath,strThisFilePath); % strTempDir is still returned but not used for dot files
fprintf('Saving output in directory %s\n',strLogDir);

%% Initialize connection with SpikeGLX
if boolUseSGL
    % Check if data are supplied
    if exist('sExpMeta','var') && isfield(sExpMeta,'hSGL') && isfield(sExpMeta,'strRunName') && isfield(sExpMeta,'sParamsSGL')
        % Get data
        hSGL = sExpMeta.hSGL;
        strRunName = sExpMeta.strRunName;
        sParamsSGL = sExpMeta.sParamsSGL;
        % Start recording
        intOutFlag = StartRecordingSGL(hSGL);
    else
        % Start connection
        fprintf('Opening SpikeGLX connection & starting recording "%s" [%s]...\n',strRecording,getTime);
        [hSGL,strRunName,sParamsSGL] = InitSGL(strRecording,strHostAddress);
    end
    fprintf('SGL saving to "%s", matlab saving to "%s.mat" [%s]...\n',strRunName,strFilename,getTime);

    % Retrieve some parameters
    intStreamNI = 0; %-1;
    dblSampFreqNI = GetStreamSampleRate(hSGL, intStreamNI, strHostAddress);

    %% Check disk space available
    strDataDirSGL = GetDataDir(hSGL);
    jFileObj = java.io.File(strDataDirSGL);
    dblFreeGB = (jFileObj.getFreeSpace)/(1024^3);
    if dblFreeGB < 100,warning([mfilename ':LowDiskSpace'],'Low disk space available (%.0fGB) for Neuropixels data (dir: %s)',dblFreeGB,strDataDirSGL);end
else
    sParamsSGL = struct;
end

%% Initialize parallel pool && GPU
if sStimParams.intUseParPool > 0 && isempty(gcp('nocreate'))
    parpool(sStimParams.intUseParPool * [1 1]);
end
if sStimParams.intUseGPU > 0
    objGPU = gpuDevice(sStimParams.intUseGPU);
end
%% Initialize NI I/O box
if boolUseNI
    %initialize
    fprintf('Connecting to National Instruments box...\n');
    boolDaqOutRunning = false;
    if exist('objDaqOut','var') && ~isempty(objDaqOut)
        try
            %turns leds on
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
        %turns leds on
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
    %% INITALIZE SCREEN
    fprintf('Starting PsychToolBox extension...\n');
    % Open window
    AssertOpenGL;
    KbName('UnifyKeyNames');
    intOldVerbosity = Screen('Preference', 'Verbosity',1); %stop PTB spamming
    if boolDebug == 1, vecInitRect = [0 0 640 640];else, vecInitRect = [];end
    try
        Screen('Preference', 'SkipSyncTests', 0);
        [ptrWindow,vecRect] = Screen('OpenWindow', intUseScreen,sStimParams.intBackground,vecInitRect);
    catch ME
        warning([mfilename ':ErrorPTB'],'Psychtoolbox error, attempting with sync test skip [msg: %s]',ME.message);
        Screen('Preference', 'SkipSyncTests', 1);
        [ptrWindow,vecRect] = Screen('OpenWindow', intUseScreen,sStimParams.intBackground,vecInitRect);
    end
    % Calibrate monitor (if gammaTable.mat exists)
    if exist('C:\Users\neuropixels\Desktop\GammaCorrection\gammaTable.mat','file')
        load("C:\Users\neuropixels\Desktop\GammaCorrection\gammaTable.mat", 'gammaTable');
        Screen('LoadNormalizedGammaTable', ptrWindow, gammaTable*[1 1 1]);
    end
    % Window variables
    sStimParams.ptrWindow = ptrWindow;
    sStimParams.vecRect = vecRect;
    sStimParams.intScreenWidth_pix = vecRect(3)-vecRect(1);
    sStimParams.intScreenHeight_pix = vecRect(4)-vecRect(2);
    sStimParams.intWhite = WhiteIndex(sStimParams.ptrWindow);
    sStimParams.intBlack = BlackIndex(sStimParams.ptrWindow); % Define black index for dots
    % Calculate pixels per degree
    % Using the average of horizontal and vertical pixels per degree for consistency
    sStimParams.dblPixelsPerDeg = mean([(sStimParams.intScreenWidth_pix/sStimParams.dblScreenWidth_deg) ...
        (sStimParams.intScreenHeight_pix/sStimParams.dblScreenHeight_deg)]);

    %% MAXIMIZE PRIORITY
    intOldPriority = 0;
    if boolDebug == 0
        intPriorityLevel=MaxPriority(ptrWindow);
        intOldPriority = Priority(intPriorityLevel);
    end

    %% Get refresh rate
    dblStimFrameRate=Screen('FrameRate', ptrWindow);
    if dblStimFrameRate == 0, dblStimFrameRate = 60; end  % Fallback if refresh rate is unknown
    intStimFrameRate = round(dblStimFrameRate);
    dblStimFrameDur = mean(1/dblStimFrameRate);
    dblInterFlipInterval = Screen('GetFlipInterval', ptrWindow);
    if dblStimFrameDur/dblInterFlipInterval > 1.05 || dblStimFrameDur/dblInterFlipInterval < 0.95
        warning([mfilename ':InconsistentFlipDur'],sprintf('Something iffy with flip speed and monitor refresh rate detected; frame duration is %fs, while flip interval is %fs!',dblStimFrameDur,dblInterFlipInterval)); %#ok<SPWRN>
    end

    %% Check escape
    if CheckEsc(),error([mfilename ':EscapePressed'],'Esc pressed; exiting');end

    %% Prepare RDK stimulus parameters and precompute trajectories
    fprintf('Precomputing RDK dot trajectories... [%s]\n', getTime);
    % Get RDK parameters from sStimParams (now using degrees)
    dotSizeDeg = sStimParams.dblDotSize_deg;
    dotSpeedDegPerSec = sStimParams.dblDotSpeed_degPerSec;
    coveragePercent = sStimParams.dblCoveragePercent;
    % dotLifetimeSeconds = sStimParams.dblDotLifetime_s;
    coherenceLevels = sStimParams.vecCoherenceLevels;
    cohDirectionsDeg = sStimParams.vecDirections_deg;
    numRepetitionsPerUniqueStim = sStimParams.intNumRepetitionsPerUniqueStim;
    numSeeds = sStimParams.intNumRNGSeeds;
    % Convert degree-based parameters to pixels and frames
    dotSizePix = round(dotSizeDeg * sStimParams.dblPixelsPerDeg);
    dotSpeedPixPerFrame = dotSpeedDegPerSec * sStimParams.dblPixelsPerDeg * dblStimFrameDur; % Speed per frame in pixels
    % Calculate derived parameters
    dotRadius = dotSizePix / 2;
    dotArea = pi * (dotRadius)^2;
    screenArea = sStimParams.intScreenWidth_pix * sStimParams.intScreenHeight_pix;
    numDots = round((coveragePercent * screenArea) / dotArea);
    % Ensure numDots is at least 1 to avoid issues with empty arrays
    if numDots < 1, numDots = 1; warning('Calculated number of dots less than 1, setting to 1.'); end
    % dotLifetimeFrames = round(dotLifetimeSeconds / dblStimFrameDur);
    numFramesPerTrial = round(sStimParams.dblSecsDuration / dblStimFrameDur);
    % Generate initial RNG states for reproducibility
    initialRNGstates = cell(1, numSeeds);
    for s = 1:numSeeds
        rng(s); % Initialize RNG with a specific seed
        initialRNGstates{s} = rng; % Store the state
    end
    % Create unique stimulus combinations and the full trial list
    % [cohGrid, dirGrid] = meshgrid(coherenceLevels, cohDirectionsDeg);
    % uniqueStimCombos = [cohGrid(:), dirGrid(:)]; % [coherence, direction]
    [cohGrid, dirGrid, lifeGrid] = ndgrid(coherenceLevels, cohDirectionsDeg, sStimParams.vecDotLifetime_s);
    uniqueStimCombos = [cohGrid(:), dirGrid(:), lifeGrid(:)]; % [coherence, direction, lifetime]
    % trialList = []; % Stores [coherence, direction, seedIdx] for each trial
    % for i = 1:size(uniqueStimCombos, 1)
    %     coh = uniqueStimCombos(i, 1);
    %     dir = uniqueStimCombos(i, 2);
    %     assignedSeeds = [];
    %     % Assign seeds ensuring repetitions are spread across available seeds
    %     numFullSets = floor(numRepetitionsPerUniqueStim / numSeeds);
    %     for s = 1:numFullSets
    %         assignedSeeds = [assignedSeeds, 1:numSeeds];
    %     end
    %     remaining = mod(numRepetitionsPerUniqueStim, numSeeds);
    %     if remaining > 0
    %         assignedSeeds = [assignedSeeds, 1:remaining];
    %     end
    %     assignedSeeds = assignedSeeds(randperm(length(assignedSeeds))); % Randomize seed assignment
    %     for r = 1:numRepetitionsPerUniqueStim
    %         trialList = [trialList; coh, dir, assignedSeeds(r)];
    %     end
    % end
    trialList = []; % Will store [coherence, direction, lifetime, seedIdx]
    for i = 1:size(uniqueStimCombos, 1)
        coh = uniqueStimCombos(i, 1);
        dir = uniqueStimCombos(i, 2);
        lifetime = uniqueStimCombos(i, 3);

        % Assign seeds
        assignedSeeds = [];
        numFullSets = floor(numRepetitionsPerUniqueStim / numSeeds);
        for s = 1:numFullSets
            assignedSeeds = [assignedSeeds, 1:numSeeds];
        end
        remaining = mod(numRepetitionsPerUniqueStim, numSeeds);
        if remaining > 0
            assignedSeeds = [assignedSeeds, 1:remaining];
        end
        assignedSeeds = assignedSeeds(randperm(length(assignedSeeds))); % shuffle

        for r = 1:numRepetitionsPerUniqueStim
            trialList = [trialList; coh, dir, lifetime, assignedSeeds(r)];
        end
    end

    rng('shuffle'); % Shuffle the main RNG for actual trial order
    trialList = trialList(randperm(size(trialList,1)), :);

    % Precompute dot trajectories and store directly into a temporary cache.
    % We will then transfer relevant objects to structEP later.
    precomputedDotsCache = containers.Map; % Cache for quicker access during presentation
    dotColor = sStimParams.intBlack; % Use black dots
    backgroundColor = sStimParams.intBackground;
    uniqueKeysForPrecomputation = unique(trialList(:, [1 2 3 4]), 'rows'); % Get unique combinations of (coh, dir, lifetime, seed)
    for i = 1:size(uniqueKeysForPrecomputation,1)
        coh = uniqueKeysForPrecomputation(i,1);
        dirDeg = uniqueKeysForPrecomputation(i,2);
        lifetimeSeconds = uniqueKeysForPrecomputation(i,3);
        seedIdx = uniqueKeysForPrecomputation(i,4);

        dotLifetimeFrames = round(lifetimeSeconds / dblStimFrameDur);

        key = sprintf('coh_%.2f_dir_%.0f_life_%.2f_seed_%d', coh, dirDeg, lifetimeSeconds,seedIdx);
        rng(initialRNGstates{seedIdx}); % Reset RNG to the stored state for this specific seed
        dirRad = deg2rad(dirDeg);
        posX = rand(1, numDots) * sStimParams.intScreenWidth_pix;
        posY = rand(1, numDots) * sStimParams.intScreenHeight_pix;
        lifetimes = randi(dotLifetimeFrames, 1, numDots); % Initial lifetimes
        % Determine which dots are coherent/incoherent
        coherent = rand(1, numDots) < coh;
        incoherent = ~coherent;
        % Assign initial directions
        incohDirs = rand(1, numDots) * 2 * pi; % Random directions for incoherent dots
        dirX = cos(incohDirs);
        dirY = -sin(incohDirs); % Y-axis is inverted in PTB (positive is down)
        % Set coherent dot directions
        dirX(coherent) = cos(dirRad);
        dirY(coherent) = -sin(dirRad); % Y-axis is inverted in PTB
        currentFrameData = cell(1, numFramesPerTrial); % Cell array to hold positions for each frame
        for f = 1:numFramesPerTrial
            % Update positions
            posX = posX + dotSpeedPixPerFrame * dirX;
            posY = posY + dotSpeedPixPerFrame * dirY;
            % Wrap around screen boundaries
            posX(posX < 0) = posX(posX < 0) + sStimParams.intScreenWidth_pix;
            posX(posX >= sStimParams.intScreenWidth_pix) = posX(posX >= sStimParams.intScreenWidth_pix) - sStimParams.intScreenWidth_pix;
            posY(posY < 0) = posY(posY < 0) + sStimParams.intScreenHeight_pix;
            posY(posY >= sStimParams.intScreenHeight_pix) = posY(posY >= sStimParams.intScreenHeight_pix) - sStimParams.intScreenHeight_pix;
            % Update lifetimes and re-seed dead dots
            lifetimes = lifetimes - 1;
            dead = lifetimes <= 0;
            if any(dead)
                numDead = sum(dead);
                % Regenerate positions for dead dots
                posX(dead) = rand(1, numDead) * sStimParams.intScreenWidth_pix;
                posY(dead) = rand(1, numDead) * sStimParams.intScreenHeight_pix;
                lifetimes(dead) = dotLifetimeFrames; % Reset lifetime
                % Re-assign coherence and direction for new dots
                coherent(dead) = rand(1, numDead) < coh;
                incoherent(dead) = ~coherent(dead);
                incohIdx = find(dead & incoherent);
                newDirs = rand(1, length(incohIdx)) * 2 * pi;
                dirX(incohIdx) = cos(newDirs);
                dirY(incohIdx) = -sin(newDirs);
                cohIdx = find(dead & coherent);
                dirX(cohIdx) = cos(dirRad);
                dirY(cohIdx) = -sin(dirRad);
            end
            currentFrameData{f} = [posX; posY]; % Store positions for this frame (2xNumDots matrix)
        end
        precomputedDotsCache(key) = currentFrameData; % Store in cache
    end
    intTotalTrials = size(trialList, 1); % Total number of trials to run
    fprintf('Estimated total stimulus time: %.1f minutes\n', intTotalTrials * (sStimParams.dblSecsDuration + sStimParams.dblSecsPostBlank + sStimParams.dblSecsPreBlank) / 60);
    %% Build structEP
    % Stimulus timing info
    initialBlank = sStimParams.dblSecsInitialBlank;
    % Total trial duration including pre and post blanks (effectively ITI)
    trialDur = sStimParams.dblSecsPreBlank + sStimParams.dblSecsDuration + sStimParams.dblSecsPostBlank;
    endBlank = sStimParams.dblSecsEndBlank;
    totalLength = initialBlank + trialDur * intTotalTrials + endBlank;
    % build structEP, stim-based logs
    structEP = struct;
    structEP.strExpType = sStimParams.strStimType;
    structEP.intTrialNum = intTotalTrials;
    structEP.TrialNumber = nan(1,structEP.intTrialNum);
    structEP.dblStimFrameDur = dblStimFrameDur;
    structEP.dblInterFlipInterval = dblInterFlipInterval;
    % RDK Specific Fields
    structEP.vecCoherence = nan(1,structEP.intTrialNum);
    structEP.vecDirection_deg = nan(1,structEP.intTrialNum);
    structEP.intRNGSeed = nan(1,structEP.intTrialNum);
    structEP.vecDotSize_deg = repmat(sStimParams.dblDotSize_deg, [1 structEP.intTrialNum]);
    structEP.vecDotSpeed_degPerSec = repmat(sStimParams.dblDotSpeed_degPerSec, [1 structEP.intTrialNum]);
    structEP.vecCoveragePercent = repmat(sStimParams.dblCoveragePercent, [1 structEP.intTrialNum]);
    % structEP.dblDotLifetime_s = repmat(sStimParams.dblDotLifetime_s, [1 structEP.intTrialNum]);
    structEP.vecDotLifetime_s =nan(1,structEP.intTrialNum);
    structEP.intNumDots = repmat(numDots, [1 structEP.intTrialNum]); % Store actual number of dots
    % Removed: structEP.strDotLocationsFilePath = cell(1,structEP.intTrialNum); % New field to store path to saved dot locations
    structEP.stimObjects = cell(1, structEP.intTrialNum); % **NEW FIELD** to store precomputed dot locations
    % General timing log fields
    structEP.ActOnSecs = nan(1,structEP.intTrialNum);
    structEP.ActOffSecs = nan(1,structEP.intTrialNum);
    structEP.ActStartSecs = nan(1,structEP.intTrialNum);
    structEP.ActStopSecs = nan(1,structEP.intTrialNum);
    structEP.ActOnNI = nan(1,structEP.intTrialNum);
    structEP.ActOffNI = nan(1,structEP.intTrialNum);
    % Define dummy/placeholder fields for compatibility if necessary
    structEP.vecDstRect= nan(4,structEP.intTrialNum); % RDK typically full screen
    structEP.vecDirection = nan(1,structEP.intTrialNum); % Will store RDK direction
    structEP.vecTexture = nan(1,structEP.intTrialNum); % Not applicable for RDK
    % Pre-calculated expected timings
    structEP.vecTrialStartSecs = initialBlank:trialDur:(totalLength-endBlank-trialDur-eps); % Use eps to avoid floating point issues
    structEP.vecTrialStimOnSecs = structEP.vecTrialStartSecs + sStimParams.dblSecsPreBlank;
    structEP.vecTrialStimOffSecs = structEP.vecTrialStimOnSecs + sStimParams.dblSecsDuration;
    structEP.vecTrialEndSecs = structEP.vecTrialStimOffSecs + sStimParams.dblSecsPostBlank;

    %% PRESENT STIMULI
    % wait for signal
    opts = struct;
    opts.Default = 'Start';
    opts.Interpreter = 'tex';
    strAns = questdlg('Would you like to start the stimulation?', ...
        'Start Stimulation', ...
        'Start','Cancel',opts);
    if ~strcmp(strAns,opts.Default)
        error([mfilename ':RunCancelled'],'Cancelling');
    end
    % Hide cursor
    % HideCursor(ptrWindow);
    % Set timers
    hTic = tic;
    dblLastFlip = Screen('Flip', ptrWindow);
    dblInitialFlip = dblLastFlip;
    % Timestamp start
    structEP.strStartDate = getDate();
    structEP.strStartTime = getTime();

    %% Calculate diode rectangle location
    vecDiodeRect = []; % Initialize as empty
    if sStimParams.intCornerTrigger > 0
        % PTB uses [top,left, right, bottom] with origin at upper left
        intCornerPix = floor(sStimParams.dblCornerSize*sStimParams.intScreenWidth_pix);
        if sStimParams.intCornerTrigger == 1 %upper left
            vecDiodeRect =  [0 0 intCornerPix intCornerPix];
        elseif sStimParams.intCornerTrigger == 2 %upper right
            vecDiodeRect =  [(vecRect(3)-intCornerPix) 0 vecRect(3) intCornerPix];
        elseif sStimParams.intCornerTrigger == 3 %lower left
            vecDiodeRect =  [0 (vecRect(4) - intCornerPix) intCornerPix  vecRect(4)];
        elseif sStimParams.intCornerTrigger == 4 %lower right
            vecDiodeRect =  [(vecRect(3)-intCornerPix) (vecRect(4) - intCornerPix) vecRect(3)  vecRect(4)];
        end
    end
    %% Wait initial blanking
    fprintf('Starting initial blank (dur=%.2fs) [%s]\n',sStimParams.dblSecsInitialBlank,getTime);
    dblInitialBlankDur = 0;
    while dblInitialBlankDur < (sStimParams.dblSecsInitialBlank-dblStimFrameDur/2)
        Screen('FillRect',ptrWindow,sStimParams.intBackground);
        dblLastFlipFlip = Screen('Flip', ptrWindow,dblLastFlip + dblStimFrameDur/2);
        dblInitialBlankDur = dblLastFlipFlip - dblInitialFlip;
    end

    %% Draw stimuli on screen
    intThisTrial = 1;
    while intThisTrial <= intTotalTrials && ~CheckEsc()
        % Check whether script execution should be paused
        if CheckPause()
            fprintf('\n\n<strong>STIMULUS SCRIPT PAUSED!</strong> [%s]\n\n',getTime);
            WaitSecs(1);
            % Get user input
            opts = struct;
            opts.Default = 'Restart';
            opts.Interpreter = 'tex';
            strAns = questdlg('STIMULUS SCRIPT PAUSED! Would you like to restart the stimulation?', ...
                'Restart Stimulation', ...
                'Restart',opts);
            WaitSecs(1);
        end
        % Trial start; background
        Screen('FillRect',ptrWindow, sStimParams.intBackground);
        dblTrialStartFlip = Screen('Flip', ptrWindow); % Actual start of trial timing
        % Fill DAQ with data (for 250ms pulse if used)
        if boolUseNI
            stop(objDaqOut);
            outputData1 = dblSyncLightMultiplier*cat(1,linspace(3, 3, 200)',linspace(0, 0, 50)'); % e.g., for synclight
            outputData2 = dblPupilLightMultiplier*linspace(3, 3, 250)'; % e.g., for pupil light
            queueOutputData(objDaqOut,[outputData1 outputData2]);
            prepare(objDaqOut);
        end

        % Get RDK parameters for this trial
        coh = trialList(intThisTrial,1);
        dir = trialList(intThisTrial,2);
        life = trialList(intThisTrial,3);
        seed = trialList(intThisTrial,4);
        key = sprintf('coh_%.2f_dir_%.0f_life_%.2f_seed_%d', coh, dir, life, seed);

        % Retrieve precomputed frames from cache & store
        frames = precomputedDotsCache(key);
        structEP.stimObjects{intThisTrial} = frames;

        % Populate structEP for this trial
        structEP.vecCoherence(intThisTrial) = coh;
        structEP.vecDotLifetime_s(intThisTrial) = life;
        structEP.vecDirection_deg(intThisTrial) = dir;
        structEP.intRNGSeed(intThisTrial) = seed;
        structEP.vecDstRect(:,intThisTrial) = [0; 0; sStimParams.intScreenWidth_pix; sStimParams.intScreenHeight_pix]; % RDK is full screen
        structEP.vecDirection(intThisTrial) = dir; % Store the main direction
        structEP.vecTexture(intThisTrial) = nan; % Not applicable for RDK
        % Removed: structEP.strDotLocationsFilePath{intThisTrial} = dotLocFilePath; % No longer storing file path
        % Get *expected* timing for this trial from pre-calculated values
        expectedStimOnSecs = structEP.vecTrialStimOnSecs(intThisTrial);
        expectedStimOffSecs = structEP.vecTrialStimOffSecs(intThisTrial);
        expectedEndSecs = structEP.vecTrialEndSecs(intThisTrial); % This is the end of the *trial block*
        % Display in command window
        fprintf('%d/%d: Coherence=%.2f, Direction=%.0f deg, Lifetime=%.2f s Seed=%d [%s]\n', ...
            intThisTrial, intTotalTrials, coh, dir, life, seed, getTime);

        %% Wait pre-blanking
        % Target time for pre-blank to end is (expectedStimOnSecs - dblInitialFlip)
        % The actual start of the trial was dblTrialStartFlip.
        % So, we want the pre-blank to last for sStimParams.dblSecsPreBlank seconds.
        Screen('FillRect',ptrWindow,sStimParams.intBackground);
        dblLastFlip = Screen('Flip', ptrWindow); % First flip of pre-blank
        % Keep flipping until the expected stimulus onset time is reached, relative to dblTrialStartFlip
        while (GetSecs - dblTrialStartFlip) < (sStimParams.dblSecsPreBlank - dblStimFrameDur/2)
            Screen('FillRect',ptrWindow, sStimParams.intBackground);
            % Schedule flip for next frame
            dblLastFlip = Screen('Flip', ptrWindow, dblLastFlip + dblStimFrameDur/2);
        end

        %% 250ms pulse at stim start (if using NI DAQ)
        if boolUseNI,startBackground(objDaqOut);end % This triggers the pre-defined DAQ output

        %% Present RDK stimulus
        boolFirstFlip = false;
        for f = 1:numFramesPerTrial % Iterate through precomputed frames
            Screen('FillRect', sStimParams.ptrWindow, backgroundColor); % Clear screen with background color
            Screen('DrawDots', sStimParams.ptrWindow, frames{f}, dotSizePix, dotColor, [], 3); % Draw dots
            % Draw diode trigger (if enabled) - white during stimulus
            if ~isempty(vecDiodeRect)
                Screen('FillRect',sStimParams.ptrWindow,sStimParams.intWhite,vecDiodeRect);
            end
            Screen('DrawingFinished',sStimParams.ptrWindow); % Optimize drawing commands
            dblLastFlip = Screen('Flip',sStimParams.ptrWindow,dblLastFlip+dblInterFlipInterval/2); % Flip at desired interval
            % Send trigger for stim start on the very first flip of the stimulus presentation
            if ~boolFirstFlip
                boolFirstFlip = true;
                % Log NI timestamp
                if boolUseSGL
                    dblStimOnNI = GetStreamSampleCount(hSGL, intStreamNI, strHostAddress)/dblSampFreqNI;
                else
                    dblStimOnNI = nan;
                end
                % Log Psychtoolbox flip timestamp
                dblStimOnFlip = dblLastFlip;
            end
        end
        % Transition back to background after stimulus
        Screen('FillRect',ptrWindow, sStimParams.intBackground);

        dblStimOffFlip = Screen('Flip', ptrWindow, dblLastFlip+dblStimFrameDur/2); % Actual stim off flip
        % dblStimDur = dblStimOffFlip-dblStimOnFlip; % Actual stimulus duration (already logged from expected)
        % Log NI timestamp for stim off
        if boolUseSGL
            dblStimOffNI =  GetStreamSampleCount(hSGL, intStreamNI, strHostAddress)/dblSampFreqNI;
        else
            dblStimOffNI = nan;
        end
        %% Wait post-blanking (effectively ITI)
        % We want the post-blank to last for sStimParams.dblSecsPostBlank seconds.
        % It started after dblStimOffFlip.
        Screen('FillRect',ptrWindow, sStimParams.intBackground);
        % The first flip of the post-blank period already happened at dblStimOffFlip.
        % So, we start the timing from that point.
        while (GetSecs - dblStimOffFlip) < (sStimParams.dblSecsPostBlank - dblStimFrameDur/2)
            Screen('FillRect',ptrWindow, sStimParams.intBackground);
            % Schedule flip for next frame
            dblLastFlip = Screen('Flip', ptrWindow, dblLastFlip + dblStimFrameDur/2);
        end
        % Capture the actual end time of this trial block (end of post-blank)
        dblTrialEndFlip = dblLastFlip;
        % Log actual timing data into structEP
        structEP.TrialNumber(intThisTrial) = intThisTrial;
        structEP.ActStartSecs(intThisTrial) = dblTrialStartFlip; % Time when pre-blank started
        structEP.ActOnSecs(intThisTrial) = dblStimOnFlip; % Time when stimulus started
        structEP.ActOffSecs(intThisTrial) = dblStimOffFlip; % Time when stimulus ended
        structEP.ActStopSecs(intThisTrial) = dblTrialEndFlip; % Time when post-blank ended
        structEP.ActOnNI(intThisTrial) = dblStimOnNI;
        structEP.ActOffNI(intThisTrial) = dblStimOffNI;
        %%increment trial number
        intThisTrial = intThisTrial+1;
    end
    % Save data before final blank/cleanup
    structEP.sStimParams = sStimParams; % Store all stimulus parameters
    save(fullfile(strLogDir,strFilename), 'structEP','sParamsSGL'); % sParamsSGL might be empty if not using SGL
    % show summary
    fprintf('Finished experiment & data saving at [%s], waiting for end blank (dur=%.2fs)\n',getTime,sStimParams.dblSecsEndBlank)

    %% Wait end-blanking
    dblEndBlankStartFlip = dblLastFlip; % Start time for the final blank from the last flip of the last trial
    dblEndBlankDur = 0;
    while dblEndBlankDur < sStimParams.dblSecsEndBlank
        Screen('FillRect',ptrWindow, sStimParams.intBackground);
        dblCurrentFlip = Screen('Flip', ptrWindow); % Get current flip time
        dblEndBlankDur = dblCurrentFlip - dblEndBlankStartFlip; % Measure duration from the start of the end blank
    end
    % clean up
    fprintf('\nExperiment is finished at [%s], closing down and cleaning up...\n',getTime);
    Screen('Close',ptrWindow);
    Screen('Close');
    Screen('CloseAll');
    ShowCursor;
    Priority(0);
    Screen('Preference', 'Verbosity',intOldVerbosity);

    %% Close Daq IO
    if boolUseNI && ~(exist('sExpMeta','var') && isfield(sExpMeta,'objDaqOut'))
        try
            closeDaqOutput(objDaqOut);
        catch
            % Error closing DAQ, can ignore for clean exit if already closing
        end
    end
catch ME

    %% General error handling block
    ShowCursor; % Ensure cursor is visible
    Priority(0); % Reset priority
    Screen('CloseAll'); % Close all Psychtoolbox windows
    Screen('Preference', 'Verbosity',intOldVerbosity); % Restore verbosity

    % Attempt to save data if an error occurred
    fprintf('\n\n\nError occurred! Trying to save data and clean up...\n\n\n');
    if exist('structEP','var') % Check if structEP was initialized
        structEP.sStimParams = sStimParams; % Store all stimulus parameters
        if ~exist('sParamsSGL','var'), sParamsSGL = []; end % Ensure sParamsSGL exists for saving
        save(fullfile(strLogDir,strFilename), 'structEP','sParamsSGL');
        fprintf('Partial data saved to %s\n', fullfile(strLogDir,strFilename));
    else
        fprintf('structEP not initialized, no data to save.\n');
    end

    % Close DAQ if it was opened and not passed in from sExpMeta
    if boolUseNI && ~(exist('sExpMeta','var') && isfield(sExpMeta,'objDaqOut'))
        try
            closeDaqOutput(objDaqOut);
        catch
            % Error closing DAQ during error, often ignorable
        end
    end

    % Re-throw the error unless it was an intentional user exit
    if strcmp(ME.identifier,[mfilename ':EscapePressed'])
        fprintf('Experiment aborted by user (Escape key pressed).\n');
    else
        rethrow(ME); % Re-throw other errors for debugging
    end
end
