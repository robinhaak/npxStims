% RH_DotStimLR
% Dot motion stimulus, opposing directions L/R
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

% Check if sExpMeta is provided by a calling script/function
if exist('sExpMeta','var')
    if isfield(sExpMeta,'dblPupilLightMultiplier'),dblPupilLightMultiplier=sExpMeta.dblPupilLightMultiplier;end
    if isfield(sExpMeta,'dblSyncLightMultiplier'),dblSyncLightMultiplier=sExpMeta.dblSyncLightMultiplier;end
    if isfield(sExpMeta,'strHostAddress'),strHostAddress=sExpMeta.strHostAddress;end
    if isfield(sExpMeta,'objDaqOut'),objDaqOut=sExpMeta.objDaqOut;end
    if isfield(sExpMeta,'boolUseSGL'),boolUseSGL=sExpMeta.boolUseSGL;end
    if isfield(sExpMeta,'boolUseNI'),boolUseNI=sExpMeta.boolUseNI;end
else
    sExpMeta = []; % Initialize if not already existing
end

%% Ask user for recording name
if exist('sStimParamsSettings','var') && isfield(sStimParamsSettings,'strRecording')
    strRecording = sStimParamsSettings.strRecording;
else
    strRecording = input('Recording name (e.g., MouseX): ', 's');
end

%% Parameters
fprintf('Loading settings...\n');
% Check if sStimParamsSettings exists and matches the stimulus type.
% If not, or if it's the wrong type, initialize new settings.
if ~exist('sStimParamsSettings','var') || isempty(sStimParamsSettings) || ~strcmpi(sStimParamsSettings.strStimType,'LeftRightRDK')
    sStimParamsSettings = struct;
    sStimParamsSettings.strStimType = 'LeftRightRDK'; % Custom stimulus type
    sStimParamsSettings.strOutputPath = 'C:\_Data\Exp'; % Appends date
    sStimParamsSettings.strTempObjectPath = 'C:\temp'; % This variable is no longer explicitly used for saving dot files

    % Visual space parameters
    sStimParamsSettings.dblScreenDistance_cm = 23; % cm; measured [~23]

    % Screen variables
    sStimParamsSettings.intCornerTrigger = 2; % Integer switch; 0=none,1=upper left, 2=upper right, 3=lower left, 4=lower right
    sStimParamsSettings.dblCornerSize = 1/30; % Fraction of screen width
    sStimParamsSettings.dblScreenWidth_cm = 51; % cm; measured [51]
    sStimParamsSettings.dblScreenHeight_cm = 29; % cm; measured [29]
    sStimParamsSettings.dblScreenWidth_deg = 2 * atand(sStimParamsSettings.dblScreenWidth_cm / (2 * sStimParamsSettings.dblScreenDistance_cm));
    sStimParamsSettings.dblScreenHeight_deg = 2 * atand(sStimParamsSettings.dblScreenHeight_cm / (2 * sStimParamsSettings.dblScreenDistance_cm));
    sStimParamsSettings.intUseScreen = 2; % Screen to use
    sStimParamsSettings.dblBackground = 0.5; % Background intensity (dbl, [0 1])
    sStimParamsSettings.intBackground = round(sStimParamsSettings.dblBackground * 255);

    % Dot parameters (common to all sets)
    sStimParamsSettings.dblDotSize_deg = 2; % Dot size in visual degrees
    sStimParamsSettings.dblDotSpeed_degPerSec = 10; % Dot speed, deg/s
    sStimParamsSettings.dblCoveragePercent = 0.12; % Screen coverage by dots
    sStimParamsSettings.dblDotLifetime_s = 1; % Dot lifetime in seconds

    % CUSTOM STIMULUS SETTINGS FOR LEFT/RIGHT RDK
    sStimParamsSettings.vecDirections_deg = 0; % Use a dummy direction (e.g., 0 degrees for right). Coherence defines the actual left/right split.
    sStimParamsSettings.vecCoherenceLevels = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]; % Full range of coherence
    sStimParamsSettings.strMotionDefinition = 'Coherence 0: All dots Right; Coherence 1: All dots Left; 0.5: 50/50 split.';
    fprintf('Selected: Custom Left/Right RDK (Coherence 0=Right, 1=Left)\n');

    % General timing parameters (common), seconds
    sStimParamsSettings.dblSecsDuration = 2; % Duration of active stimulus presentation
    sStimParamsSettings.dblSecsInitialBlank = 5; % Blank screen at start of experiment
    sStimParamsSettings.dblSecsPreBlank = 0.5; % Blank screen before each stimulus presentation
    sStimParamsSettings.dblSecsPostBlank = 0.5; % Blank screen after each stimulus presentation (effectively ITI)
    sStimParamsSettings.dblSecsEndBlank = 5; % Blank screen at end of experiment

    % Control variables (common)
    sStimParamsSettings.intNumRepetitionsPerUniqueStim = 20; % Repetitions for each unique stimulus
    sStimParamsSettings.intNumRNGSeeds = 4; % Number of RNG seeds to use per unique stimulus
    sStimParamsSettings.intUseDaqDevice = 1; % ID of DAQ device
    sStimParamsSettings.intUseParPool = 0; % Number of workers in parallel pool; [2]
    sStimParamsSettings.intUseGPU = 0; % Set to 1 if you have a GPU and Psychtoolbox supports it (for MakeTexture etc)
else
    % If sStimParamsSettings already exists and is from this script,
    % evaluate and assign pre-defined values to structure.
    % This handles cases where sStimParamsSettings might be loaded from a file
    % where values are stored as strings that need evaluation (e.g., '2*atand(...)').
    cellFields = fieldnames(sStimParamsSettings);
    for intField=1:numel(cellFields)
        try
            % Check if the field value is a string and needs evaluation
            if ischar(sStimParamsSettings.(cellFields{intField}))
                sStimParamsSettings.(cellFields{intField}) = eval(sStimParamsSettings.(cellFields{intField}));
            end
        catch
            % Already evaluated or not a string that needs evaluation,
            % or contains an error in evaluation string.
            % For debugging, you might want to add: warning('Failed to evaluate field %s', cellFields{intField});
        end
    end
end

% Set screen for Psychtoolbox based on debug mode
if boolDebug == 1
    intUseScreen = 0; % Use primary screen for debugging (often smaller window)
else
    intUseScreen = sStimParamsSettings.intUseScreen; % Use specified screen for full experiment
end
sStimParams = sStimParamsSettings; % Copy settings to sStimParams for use in script

%% Set output locations for logs
strOutputPath = sStimParamsSettings.strOutputPath;
strTempObjectPath = sStimParamsSettings.strTempObjectPath;
strThisFilePath = mfilename('fullpath');
[strFilename,strLogDir,strTempDir,~] = RE_assertPaths(strOutputPath,strRecording,strTempObjectPath,strThisFilePath);
fprintf('Saving output in directory %s\n',strLogDir);

%% Initialize connection with SpikeGLX
if boolUseSGL
    % Check if data are supplied from a parent script/function
    if exist('sExpMeta','var') && isfield(sExpMeta,'hSGL') && isfield(sExpMeta,'strRunName') && isfield(sExpMeta,'sParamsSGL')
        % Get data
        hSGL = sExpMeta.hSGL;
        strRunName = sExpMeta.strRunName;
        sParamsSGL = sExpMeta.sParamsSGL;
        % Start recording
        intOutFlag = StartRecordingSGL(hSGL); 
    else
        % Start new connection if not provided
        fprintf('Opening SpikeGLX connection & starting recording "%s" [%s]...\n',strRecording,getTime);
        [hSGL,strRunName,sParamsSGL] = InitSGL(strRecording,strHostAddress);
    end
    fprintf('SGL saving to "%s", matlab saving to "%s.mat" [%s]...\n',strRunName,strFilename,getTime);
    % Retrieve some parameters
    intStreamNI = 0; % Or whatever stream index corresponds to NI data in SGL
    dblSampFreqNI = GetStreamSampleRate(hSGL, intStreamNI, strHostAddress);

    %% Check disk space available
    strDataDirSGL = GetDataDir(hSGL);
    jFileObj = java.io.File(strDataDirSGL);
    dblFreeGB = (jFileObj.getFreeSpace)/(1024^3);
    if dblFreeGB < 100,warning([mfilename ':LowDiskSpace'],'Low disk space available (%.0fGB) for Neuropixels data (dir: %s)',dblFreeGB,strDataDirSGL);end
else
    sParamsSGL = struct; % Initialize as empty struct if SGL not used
end

%% Initialize parallel pool && GPU
if sStimParams.intUseParPool > 0 && isempty(gcp('nocreate'))
    parpool(sStimParams.intUseParPool * [1 1]); % Start parallel pool
end
if sStimParams.intUseGPU > 0
    objGPU = gpuDevice(sStimParams.intUseGPU); % % Select GPU device
end

%% Initialize NI I/O box
if boolUseNI
    fprintf('Connecting to National Instruments box...\n');
    boolDaqOutRunning = false;
    if exist('objDaqOut','var') && ~isempty(objDaqOut)
        try
            % If DAQ object already exists and is valid, try to use it
            stop(objDaqOut);
            outputData1 = dblSyncLightMultiplier*cat(1,linspace(3, 3, 200)',linspace(0, 0, 50)'); % Sync pulse example
            outputData2 = dblPupilLightMultiplier*linspace(3, 3, 250)'; % Pupil light example
            queueOutputData(objDaqOut,[outputData1 outputData2]);
            prepare(objDaqOut);
            pause(0.1); % Short pause to ensure preparation
            startBackground(objDaqOut) % Start background output
            boolDaqOutRunning = true;
        catch
            % If existing object failed, it will be re-initialized
        end
    end
    if ~boolDaqOutRunning
        objDaqOut = openDaqOutput(sStimParamsSettings.intUseDaqDevice); % Custom function to open DAQ
        % Turn LEDs on as initial state
        stop(objDaqOut);
        outputData1 = dblSyncLightMultiplier*cat(1,linspace(3, 3, 200)',linspace(0, 0, 50)');
        outputData2 = dblPupilLightMultiplier*linspace(3, 3, 250)';
        queueOutputData(objDaqOut,[outputData1 outputData2]);
        prepare(objDaqOut);
        pause(0.1);
        startBackground(objDaqOut)
    end
end

try % Main experiment try-catch block for robust error handling
    %% INITALIZE SCREEN
    fprintf('Starting PsychToolBox extension...\n');
    % Open window
    AssertOpenGL; % Ensure OpenGL is supported
    KbName('UnifyKeyNames'); % Standardize keyboard names
    intOldVerbosity = Screen('Preference', 'Verbosity',1); % Suppress PTB spam
    if boolDebug == 1, vecInitRect = [0 0 640 640];else, vecInitRect = [];end % Debug window size
    try
        Screen('Preference', 'SkipSyncTests', 0); % Attempt strict sync tests first
        [ptrWindow,vecRect] = Screen('OpenWindow', intUseScreen,sStimParams.intBackground,vecInitRect);
    catch ME
        warning([mfilename ':ErrorPTB'],'Psychtoolbox error, attempting with sync test skip [msg: %s]',ME.message);
        Screen('Preference', 'SkipSyncTests', 1); % Fallback to skip sync tests
        [ptrWindow,vecRect] = Screen('OpenWindow', intUseScreen,sStimParams.intBackground,vecInitRect);
    end

    % Calibrate monitor (if gammaTable.mat exists)
    % IMPORTANT: Update this path to your gamma table location if you use one
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
    sStimParams.dblPixelsPerDeg = mean([(sStimParams.intScreenWidth_pix/sStimParams.dblScreenWidth_deg) ...
        (sStimParams.intScreenHeight_pix/sStimParams.dblScreenHeight_deg)]);

    %% MAXIMIZE PRIORITY
    intOldPriority = 0;
    if boolDebug == 0
        intPriorityLevel=MaxPriority(ptrWindow);
        intOldPriority = Priority(intPriorityLevel); % Elevate MATLAB process priority
    end

    %% Get refresh rate
    dblStimFrameRate=Screen('FrameRate', ptrWindow);
    if dblStimFrameRate == 0, dblStimFrameRate = 60; end  % Fallback if refresh rate is unknown
    intStimFrameRate = round(dblStimFrameRate); 
    dblStimFrameDur = mean(1/dblStimFrameRate); % Duration of one frame
    dblInterFlipInterval = Screen('GetFlipInterval', ptrWindow); % Measured inter-flip interval
    if dblStimFrameDur/dblInterFlipInterval > 1.05 || dblStimFrameDur/dblInterFlipInterval < 0.95
        warning([mfilename ':InconsistentFlipDur'],sprintf('Something iffy with flip speed and monitor refresh rate detected; frame duration is %fs, while flip interval is %fs!',dblStimFrameDur,dblInterFlipInterval)); %#ok<SPWRN>
    end

    %% Check escape (custom function)
    if CheckEsc(),error([mfilename ':EscapePressed'],'Esc pressed; exiting');end

    %% Prepare RDK stimulus parameters and precompute trajectories
    fprintf('Precomputing RDK dot trajectories (Left/Right only)... [%s]\n', getTime);

    % Get RDK parameters from sStimParams (now using degrees)
    dotSizeDeg = sStimParams.dblDotSize_deg;
    dotSpeedDegPerSec = sStimParams.dblDotSpeed_degPerSec;
    coveragePercent = sStimParams.dblCoveragePercent;
    dotLifetimeSeconds = sStimParams.dblDotLifetime_s;
    coherenceLevels = sStimParams.vecCoherenceLevels;
    % cohDirectionsDeg = sStimParams.vecDirections_deg; % No longer directly used for multiple directions

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
    dotLifetimeFrames = round(dotLifetimeSeconds / dblStimFrameDur);
    numFramesPerTrial = round(sStimParams.dblSecsDuration / dblStimFrameDur);

    % Generate initial RNG states for reproducibility
    initialRNGstates = cell(1, numSeeds);
    for s = 1:numSeeds
        rng(s); % Initialize RNG with a specific seed
        initialRNGstates{s} = rng; % Store the state
    end

    % Create unique stimulus combinations and the full trial list
    % For Left/Right RDK, unique combinations are just coherence levels and seeds
    % The 'direction' from sStimParams.vecDirections_deg (which is just 0)
    % will be included in the grid for consistency with structEP, but its value
    % will be ignored for actual motion calculation in this RDK type.
    [cohGrid, ~] = meshgrid(coherenceLevels, sStimParams.vecDirections_deg); % cohGrid takes all coherences, dirGrid will just repeat the 0.
    uniqueStimCombos = [cohGrid(:), repmat(sStimParams.vecDirections_deg(1),numel(cohGrid),1)]; % [coherence, dummy_direction]

    trialList = []; % Stores [coherence, dummy_direction, seedIdx] for each trial
    for i = 1:size(uniqueStimCombos, 1)
        coh_val = uniqueStimCombos(i, 1);
        dir_dummy = uniqueStimCombos(i, 2); % % This value will be 0 and largely ignored
        assignedSeeds = [];
        % Assign seeds ensuring repetitions are spread across available seeds
        numFullSets = floor(numRepetitionsPerUniqueStim / numSeeds);
        for s = 1:numFullSets
            assignedSeeds = [assignedSeeds, 1:numSeeds]; 
        end
        remaining = mod(numRepetitionsPerUniqueStim, numSeeds);
        if remaining > 0
            assignedSeeds = [assignedSeeds, 1:remaining]; 
        end
        assignedSeeds = assignedSeeds(randperm(length(assignedSeeds))); % Randomize seed assignment
        for r = 1:numRepetitionsPerUniqueStim
            trialList = [trialList; coh_val, dir_dummy, assignedSeeds(r)]; 
        end
    end
    rng('shuffle'); % Shuffle the main RNG for actual trial order
    trialList = trialList(randperm(size(trialList,1)), :); % Randomize trial order

    % Precompute dot trajectories and store directly into a temporary cache.
    precomputedDotsCache = containers.Map; % Cache for quicker access during presentation
    dotColor = sStimParams.intBlack; % Use black dots
    backgroundColor = sStimParams.intBackground;

    % Get unique combinations of (coh, dir, seed) that need precomputation
    uniqueKeysForPrecomputation = unique(trialList(:,1:3), 'rows');

    for i = 1:size(uniqueKeysForPrecomputation,1)
        coh_trial = uniqueKeysForPrecomputation(i,1);
        dirDeg_dummy = uniqueKeysForPrecomputation(i,2); % % This is the dummy direction (0)
        seedIdx = uniqueKeysForPrecomputation(i,3);
        key = sprintf('coh_%.2f_dir_%.0f_seed_%d', coh_trial, dirDeg_dummy, seedIdx);

        rng(initialRNGstates{seedIdx}); % Reset RNG to the stored state for this specific seed

        % Initialize dot positions
        posX = rand(1, numDots) * sStimParams.intScreenWidth_pix;
        posY = rand(1, numDots) * sStimParams.intScreenHeight_pix;
        lifetimes = randi(dotLifetimeFrames, 1, numDots); % Initial lifetimes

        % --- CUSTOM LEFT/RIGHT COHERENCE LOGIC ---
        % Coherence (coh_trial): 0 = all right, 1 = all left, 0.5 = 50/50
        % Proportion moving left: coh_trial
        % Proportion moving right: 1 - coh_trial

        % Determine initial directions for all dots based on coherence
        dotRandoms = rand(1, numDots); % Random number for each dot to determine its direction
        dirX = zeros(1, numDots);
        dirY = zeros(1, numDots); % Y-direction is always 0 for purely left/right motion

        % Assign left-moving dots (180 degrees)
        leftMovingIdx = dotRandoms < coh_trial;
        dirX(leftMovingIdx) = cos(deg2rad(180)); % -1 for left
        dirY(leftMovingIdx) = -sin(deg2rad(180)); % 0

        % Assign right-moving dots (0 degrees)
        rightMovingIdx = ~leftMovingIdx;
        dirX(rightMovingIdx) = cos(deg2rad(0)); % 1 for right
        dirY(rightMovingIdx) = -sin(deg2rad(0)); % 0
        % --- END CUSTOM COHERENCE LOGIC ---

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
            dead = lifetimes <= 0; % This is a 1 x numDots logical array

            if any(dead)
                numDead = sum(dead);
                % Regenerate positions for dead dots
                posX(dead) = rand(1, numDead) * sStimParams.intScreenWidth_pix;
                posY(dead) = rand(1, numDead) * sStimParams.intScreenHeight_pix;
                lifetimes(dead) = dotLifetimeFrames; % Reset lifetime

                % --- CUSTOM RE-SEEDING LOGIC (Corrected) ---
                % Re-assign directions for new dots based on the trial's coherence
                % Create an array of random numbers, one for each *dead* dot.
                dotRandoms_for_dead = rand(1, numDead);

                % Determine which of these *dead* dots will move left based on coherence
                temp_newLeftMoving = dotRandoms_for_dead < coh_trial; % This is 1 x numDead

                % Now, we need to apply these new directions back to the original dirX/dirY arrays.
                % Find the indices of the dead dots.
                dead_indices = find(dead);

                % Assign directions for newly re-seeded dots
                % For dead dots that should move left:
                left_moving_dead_indices = dead_indices(temp_newLeftMoving);
                dirX(left_moving_dead_indices) = cos(deg2rad(180)); % -1 for left
                dirY(left_moving_dead_indices) = -sin(deg2rad(180)); % 0

                % For dead dots that should move right:
                right_moving_dead_indices = dead_indices(~temp_newLeftMoving);
                dirX(right_moving_dead_indices) = cos(deg2rad(0)); % 1 for right
                dirY(right_moving_dead_indices) = -sin(deg2rad(0)); % 0
                % --- END CUSTOM RE-SEEDING LOGIC ---
            end
            currentFrameData{f} = [posX; posY]; % Store positions for this frame (2xNumDots matrix)
        end
        precomputedDotsCache(key) = currentFrameData; % Store in cache
    end
    intTotalTrials = size(trialList, 1); % Total number of trials to run
    fprintf('Estimated total stimulus time: %.1f minutes\n', intTotalTrials * (sStimParams.dblSecsDuration + sStimParams.dblSecsPostBlank + sStimParams.dblSecsPreBlank) / 60);

    %% Build structEP (Experiment Protocol structure for logging)
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

    % RDK Specific Fields (for Left/Right only)
    structEP.dblCoherence = nan(1,structEP.intTrialNum);
    structEP.dblDirection_deg = nan(1,structEP.intTrialNum); % This will always be 0 in this config
    structEP.intRNGSeed = nan(1,structEP.intTrialNum);
    structEP.dblDotSize_deg = repmat(sStimParams.dblDotSize_deg, [1 structEP.intTrialNum]);
    structEP.dblDotSpeed_degPerSec = repmat(sStimParams.dblDotSpeed_degPerSec, [1 structEP.intTrialNum]);
    structEP.dblCoveragePercent = repmat(sStimParams.dblCoveragePercent, [1 structEP.intTrialNum]);
    structEP.dblDotLifetime_s = repmat(sStimParams.dblDotLifetime_s, [1 structEP.intTrialNum]);
    structEP.intNumDots = repmat(numDots, [1 structEP.intTrialNum]); % Store actual number of dots
    structEP.stimObjects = cell(1, structEP.intTrialNum); % Stores precomputed dot locations

    % General timing log fields
    structEP.ActOnSecs = nan(1,structEP.intTrialNum);
    structEP.ActOffSecs = nan(1,structEP.intTrialNum);
    structEP.ActStartSecs = nan(1,structEP.intTrialNum);
    structEP.ActStopSecs = nan(1,structEP.intTrialNum);
    structEP.ActOnNI = nan(1,structEP.intTrialNum);
    structEP.ActOffNI = nan(1,structEP.intTrialNum);

    % Define dummy/placeholder fields for compatibility if necessary
    structEP.vecDstRect= nan(4,structEP.intTrialNum); % RDK typically full screen
    structEP.vecDirection = nan(1,structEP.intTrialNum); % Will store RDK direction (always 0 here)
    structEP.vecTexture = nan(1,structEP.intTrialNum); % Not applicable for RDK

    % Pre-calculated expected timings
    structEP.vecTrialStartSecs = initialBlank:trialDur:(totalLength-endBlank-eps); % Use eps to avoid floating point issues
    structEP.vecTrialStimOnSecs = structEP.vecTrialStartSecs + sStimParams.dblSecsPreBlank;
    structEP.vecTrialStimOffSecs = structEP.vecTrialStimOnSecs + sStimParams.dblSecsDuration;
    structEP.vecTrialEndSecs = structEP.vecTrialStimOffSecs + sStimParams.dblSecsPostBlank;

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

    %% PRESENT STIMULI
    % wait for signal from user
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
    % HideCursor(ptrWindow); % Uncomment if you want to hide the cursor during experiment

    % Set timers
    hTic = tic; % % Start a stopwatch timer
    dblLastFlip = Screen('Flip', ptrWindow); % Initial flip to get a base timestamp
    dblInitialFlip = dblLastFlip; % Store this as the experiment start time

    % Timestamp start of experiment in metadata
    structEP.strStartDate = getDate(); % Custom function to get date string
    structEP.strStartTime = getTime(); % Custom function to get time string

    %% Wait initial blanking
    fprintf('Starting initial blank (dur=%.2fs) [%s]\n',sStimParams.dblSecsInitialBlank,getTime);
    dblInitialBlankDur = 0;
    % Loop to display background until initial blank duration is met
    while dblInitialBlankDur < (sStimParams.dblSecsInitialBlank - dblStimFrameDur/2)
        Screen('FillRect',ptrWindow,sStimParams.intBackground); % Draw background
        dblLastFlipFlip = Screen('Flip', ptrWindow, dblLastFlip + dblStimFrameDur/2); % Flip, timed to previous flip
        dblInitialBlankDur = dblLastFlipFlip - dblInitialFlip; % Calculate elapsed time since experiment start
    end

    %% Draw stimuli on screen (Main Trial Loop)
    intThisTrial = 1;
    while intThisTrial <= intTotalTrials && ~CheckEsc() % Loop through all trials until complete or Esc pressed
        % Check whether script execution should be paused (custom function)
        if CheckPause()
            fprintf('\n\n<strong>STIMULUS SCRIPT PAUSED!</strong> [%s]\n\n',getTime);
            WaitSecs(1); % Short pause before showing dialog
            % Get user input to restart or cancel
            opts = struct;
            opts.Default = 'Restart';
            opts.Interpreter = 'tex';
            strAns = questdlg('STIMULUS SCRIPT PAUSED! Would you like to restart the stimulation?', ...
                'Restart Stimulation', ...
                'Restart',opts);
            WaitSecs(1); % Short pause after dialog
            if ~strcmp(strAns,opts.Default) % If user did not choose 'Restart'
                 error([mfilename ':RunCancelled'],'Stimulation paused and not restarted; exiting.');
            end
        end

        % Trial start; display background and record flip time
        Screen('FillRect',ptrWindow, sStimParams.intBackground);
        dblTrialStartFlip = Screen('Flip', ptrWindow); % Actual start of this trial block timing

        % Fill DAQ with data (for 250ms pulse if used, prepared for later trigger)
        if boolUseNI
            stop(objDaqOut); % Stop any ongoing DAQ output
            outputData1 = dblSyncLightMultiplier*cat(1,linspace(3, 3, 200)',linspace(0, 0, 50)'); % Example sync pulse (200ms high, 50ms low)
            outputData2 = dblPupilLightMultiplier*linspace(3, 3, 250)'; % Example pupil light signal (250ms high)
            queueOutputData(objDaqOut,[outputData1 outputData2]); % Queue data
            prepare(objDaqOut); % Prepare DAQ for output
        end

        % Get RDK parameters for this specific trial from trialList
        coh = trialList(intThisTrial,1);
        dir = trialList(intThisTrial,2); % This will be the dummy direction (0)
        seed = trialList(intThisTrial,3);
        key = sprintf('coh_%.2f_dir_%.0f_seed_%d', coh, dir, seed); % Construct key for cache lookup

        % Retrieve precomputed frames from cache & store in structEP for logging
        frames = precomputedDotsCache(key);
        structEP.stimObjects{intThisTrial} = frames; % Store the precomputed frames for this trial

        % Populate structEP for this trial's specific parameters
        structEP.dblCoherence(intThisTrial) = coh;
        structEP.dblDirection_deg(intThisTrial) = dir; % Log the dummy direction (0)
        structEP.intRNGSeed(intThisTrial) = seed;
        structEP.vecDstRect(:,intThisTrial) = [0; 0; sStimParams.intScreenWidth_pix; sStimParams.intScreenHeight_pix]; % RDK is full screen
        structEP.vecDirection(intThisTrial) = dir; % Store the main direction (0)
        structEP.vecTexture(intThisTrial) = nan; % Not applicable for RDK

        % Get *expected* timing for this trial from pre-calculated values
        expectedStimOnSecs = structEP.vecTrialStimOnSecs(intThisTrial); 
        expectedStimOffSecs = structEP.vecTrialStimOffSecs(intThisTrial); 
        expectedEndSecs = structEP.vecTrialEndSecs(intThisTrial); 

        % Display current trial info in command window
        fprintf('%d/%d: Coherence=%.2f, Seed=%d (Motion: Coherence %.0f is Left, %.0f is Right) [%s]\n', ...
            intThisTrial, intTotalTrials, coh, seed, coh*100, (1-coh)*100, getTime);

        %% Wait pre-blanking
        Screen('FillRect',ptrWindow,sStimParams.intBackground); % Ensure background is drawn
        dblLastFlip = Screen('Flip', ptrWindow); % First flip of pre-blank period
        % Keep flipping background until the pre-blank duration is met relative to dblTrialStartFlip
        while (GetSecs - dblTrialStartFlip) < (sStimParams.dblSecsPreBlank - dblStimFrameDur/2)
            Screen('FillRect',ptrWindow, sStimParams.intBackground);
            dblLastFlip = Screen('Flip', ptrWindow, dblLastFlip + dblStimFrameDur/2); % Schedule next flip
        end

        %% 250ms pulse at stim start (if using NI DAQ)
        if boolUseNI,startBackground(objDaqOut);end % This triggers the pre-defined DAQ output pattern

        %% Present RDK stimulus
        boolFirstFlip = false; % Flag to capture first flip of stimulus
        for f = 1:numFramesPerTrial % Iterate through precomputed frames for this trial
            Screen('FillRect', sStimParams.ptrWindow, backgroundColor); % Clear screen with background color
            Screen('DrawDots', sStimParams.ptrWindow, frames{f}, dotSizePix, dotColor, [], 3); % Draw dots at precomputed positions

            % Draw diode trigger (if enabled) - white during stimulus presentation
            if ~isempty(vecDiodeRect)
                Screen('FillRect',sStimParams.ptrWindow,sStimParams.intWhite,vecDiodeRect);
            end

            Screen('DrawingFinished',sStimParams.ptrWindow); % Optimize drawing commands
            dblLastFlip = Screen('Flip',sStimParams.ptrWindow,dblLastFlip+dblInterFlipInterval/2); % Flip at desired interval

            % Send trigger for stim start on the very first flip of the stimulus presentation
            if ~boolFirstFlip
                boolFirstFlip = true;
                % Log NI timestamp if using SpikeGLX
                if boolUseSGL
                    dblStimOnNI = GetStreamSampleCount(hSGL, intStreamNI, strHostAddress)/dblSampFreqNI;
                else
                    dblStimOnNI = nan;
                end
                % Log Psychtoolbox flip timestamp for stimulus onset
                dblStimOnFlip = dblLastFlip;
            end
        end

        % Transition back to background after stimulus
        Screen('FillRect',ptrWindow, sStimParams.intBackground);
        dblStimOffFlip = Screen('Flip', ptrWindow, dblLastFlip+dblStimFrameDur/2); % Actual stim off flip
        % Log NI timestamp for stim off if using SpikeGLX
        if boolUseSGL
            dblStimOffNI =  GetStreamSampleCount(hSGL, intStreamNI, strHostAddress)/dblSampFreqNI;
        else
            dblStimOffNI = nan;
        end

        %% Wait post-blanking (effectively ITI - Inter-Trial Interval)
        % Ensure background is drawn for the post-blank period
        Screen('FillRect',ptrWindow, sStimParams.intBackground);
        % Continue flipping background until the post-blank duration is met
        while (GetSecs - dblStimOffFlip) < (sStimParams.dblSecsPostBlank - dblStimFrameDur/2)
            Screen('FillRect',ptrWindow, sStimParams.intBackground);
            dblLastFlip = Screen('Flip', ptrWindow, dblLastFlip + dblStimFrameDur/2);
        end
        % Capture the actual end time of this trial block (end of post-blank)
        dblTrialEndFlip = dblLastFlip;

        % Log actual timing data into structEP for this trial
        structEP.TrialNumber(intThisTrial) = intThisTrial;
        structEP.ActStartSecs(intThisTrial) = dblTrialStartFlip; % Time when pre-blank started
        structEP.ActOnSecs(intThisTrial) = dblStimOnFlip; % Time when stimulus started
        structEP.ActOffSecs(intThisTrial) = dblStimOffFlip; % Time when stimulus ended
        structEP.ActStopSecs(intThisTrial) = dblTrialEndFlip; % Time when post-blank ended
        structEP.ActOnNI(intThisTrial) = dblStimOnNI; % NI timestamp for stim on
        structEP.ActOffNI(intThisTrial) = dblStimOffNI; % NI timestamp for stim off

        %% Increment trial number
        intThisTrial = intThisTrial+1;
    end

    % Save data before final blank/cleanup
    structEP.sStimParams = sStimParams; % Store all stimulus parameters in the log struct
    save(fullfile(strLogDir,strFilename), 'structEP','sParamsSGL'); % sParamsSGL might be empty if not using SGL
    fprintf('Finished experiment & data saving at [%s], waiting for end blank (dur=%.2fs)\n',getTime,sStimParams.dblSecsEndBlank)

    %% Wait end-blanking (final blank screen)
    dblEndBlankStartFlip = dblLastFlip; % Start time for the final blank from the last flip of the last trial
    dblEndBlankDur = 0;
    while dblEndBlankDur < sStimParams.dblSecsEndBlank
        Screen('FillRect',ptrWindow, sStimParams.intBackground);
        dblCurrentFlip = Screen('Flip', ptrWindow); % Get current flip time
        dblEndBlankDur = dblCurrentFlip - dblEndBlankStartFlip; % Measure duration from the start of the end blank
    end

    % Clean up Psychtoolbox and system settings
    fprintf('\nExperiment is finished at [%s], closing down and cleaning up...\n',getTime);
    Screen('Close',ptrWindow); % Close the specific window
    Screen('Close'); % Close any other open screens/textures
    Screen('CloseAll'); % Ensure all PTB resources are released
    ShowCursor; % Restore mouse cursor visibility
    Priority(0); % Reset MATLAB process priority to normal
    Screen('Preference', 'Verbosity',intOldVerbosity); % Restore PTB verbosity

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
