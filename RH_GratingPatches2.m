% Sequential drifting grating patches, one location at a time, no overlap
% Square patch (no circular mask)
% Random initial phase per patch (handled by random texture selection and drifting phase)
% Press ESC to exit

%% suppress m-lint warnings
%#ok<*MCCD,*NASGU,*ASGLU,*CTCH>
% clearvars -except sStimPresets sStimParamsSettings sExpMeta;
clear; close all; Screen('CloseAll');

%% define variables
fprintf('Starting %s [%s]\n',mfilename,getTime);
boolUseSGL = false;
boolUseNI = false;
boolDebug = false;
%defaults
dblPupilLightMultiplier = 1; %strength of infrared LEDs
dblSyncLightMultiplier = 0.5;
strHostAddress = '192.87.11.8'; %'192.87.11.133'; %default host address
objDaqOut = [];
if exist('sExpMeta','var')
    %expand structure
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

%% input params
fprintf('Loading settings...\n');
if ~exist('sStimParamsSettings','var') || isempty(sStimParamsSettings) || ~(strcmpi(sStimParamsSettings.strStimType,'SquareGrating') || ~strcmpi(sStimParamsSettings.strStimType,'SineGrating'))
    %general
    sStimParamsSettings = struct;
    sStimParamsSettings.strStimType = 'GratingPatches';
    sStimParamsSettings.strOutputPath = 'C:\_Data\Exp'; %appends date
    sStimParamsSettings.strTempObjectPath = 'C:\temp';%X:\JorritMontijn\ or F:\Data\Temp\
    %visual space parameters
    % 	sStimParamsSettings.dblSubjectPosX_cm = 0; % cm; relative to center of screen
    % 	sStimParamsSettings.dblSubjectPosY_cm = 0; % cm; relative to center of screen
    sStimParamsSettings.dblScreenDistance_cm = 23; % cm; measured [~23]
    %screen variables
    sStimParamsSettings.intCornerTrigger = 2; % integer switch; 0=none,1=upper left, 2=upper right, 3=lower left, 4=lower right
    sStimParamsSettings.dblCornerSize = 1/30; % fraction of screen width
    sStimParamsSettings.dblScreenWidth_cm = 51; % cm; measured [51]
    sStimParamsSettings.dblScreenHeight_cm = 29; % cm; measured [29]
    sStimParamsSettings.dblScreenWidth_deg = 2 * atand(sStimParamsSettings.dblScreenWidth_cm / (2 * sStimParamsSettings.dblScreenDistance_cm));
    sStimParamsSettings.dblScreenHeight_deg = 2 * atand(sStimParamsSettings.dblScreenHeight_cm / (2 * sStimParamsSettings.dblScreenDistance_cm));
    sStimParamsSettings.intUseScreen = 2; % 2; %which screen to use
    sStimParamsSettings.dblBackground = 0.5; %background intensity (dbl, [0 1])
    sStimParamsSettings.intBackground = round(mean(sStimParamsSettings.dblBackground)*255);
    %stimulus parameters
    sStimParamsSettings.dblStimulusSize_deg = 9; %6 %deg; (approximate) size of grating patches
    sStimParamsSettings.dblSpatialFrequency_cd = 0.06;%* %0.04; %spatial frequency in cycles per degree (changed from 1/StimSizeDeg to match your first script)
    sStimParamsSettings.dblTemporalFrequency = 3; %temporal frequency in cycles per second
    sStimParamsSettings.dblSecsDuration = 1; %s (changed from 0.5s to match your first script)
    sStimParamsSettings.dblSecsInitialBlank = 5; %s
    sStimParamsSettings.dblSecsPreBlank = 0.35; %0.25 %s
    sStimParamsSettings.dblSecsPostBlank = 0.15; %0.25; %0.5 %s
    sStimParamsSettings.dblSecsEndBlank = 5; %s
    sStimParamsSettings.intStimulusRepeats = 5; %per location (irrespective of direction or phase)
    sStimParamsSettings.vecDirections = 45:45:180; %drifting directions
    sStimParamsSettings.str90Deg = '0 degrees is rightward motion; 90 degrees is upward motion';
    %*for size= 9deg & tf= 3Hz, spf= ~0.11c/deg, speed= 27deg/s
    %control variables
    sStimParamsSettings.intUseDaqDevice = 1; %ID of DAQ device
    sStimParamsSettings.intUseParPool = 0; %number of workers in parallel pool; [2]
    sStimParamsSettings.intUseGPU = 1;
else
    % evaluate and assign pre-defined values to structure
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
    %check if data are supplied
    if exist('sExpMeta','var') && isfield(sExpMeta,'hSGL') && isfield(sExpMeta,'strRunName') && isfield(sExpMeta,'sParamsSGL')
        %get data
        hSGL = sExpMeta.hSGL;
        strRunName = sExpMeta.strRunName;
        sParamsSGL = sExpMeta.sParamsSGL;
        %start recording
        intOutFlag = StartRecordingSGL(hSGL);
    else
        %start connection
        fprintf('Opening SpikeGLX connection & starting recording "%s" [%s]...\n',strRecording,getTime);
        [hSGL,strRunName,sParamsSGL] = InitSGL(strRecording,strHostAddress);
    end
    fprintf('SGL saving to "%s", matlab saving to "%s.mat" [%s]...\n',strRunName,strFilename,getTime);
    %retrieve some parameters
    intStreamNI = 0; %-1;
%     dblSampFreqNI = GetSampleRate(hSGL, intStreamNI);
    dblSampFreqNI = GetStreamSampleRate(hSGL, intStreamNI, strHostAddress);
    
    %% check disk space available
    strDataDirSGL = GetDataDir(hSGL);
    jFileObj = java.io.File(strDataDirSGL);
    dblFreeGB = (jFileObj.getFreeSpace)/(1024^3);
    if dblFreeGB < 100,warning([mfilename ':LowDiskSpace'],'Low disk space available (%.0fGB) for Neuropixels data (dir: %s)',dblFreeGB,strDataDirSGL);end
else
    sParamsSGL = struct;
end

%% initialize parallel pool && gpu
if sStimParams.intUseParPool > 0 && isempty(gcp('nocreate'))
    parpool(sStimParams.intUseParPool * [1 1]);
end
if sStimParams.intUseGPU > 0
    objGPU = gpuDevice(sStimParams.intUseGPU);
end

%% initialize NI I/O box
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
    %open window
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
    %calibrate monitor
    if exist('C:\Users\neuropixels\Desktop\GammaCorrection\gammaTable.mat','file')
        load("C:\Users\neuropixels\Desktop\GammaCorrection\gammaTable.mat");
        Screen('LoadNormalizedGammaTable', ptrWindow, gammaTable*[1 1 1]);
    end
    %window variables
    sStimParams.ptrWindow = ptrWindow;
    sStimParams.vecRect = vecRect;
    sStimParams.intScreenWidth_pix = vecRect(3)-vecRect(1);
    sStimParams.intScreenHeight_pix = vecRect(4)-vecRect(2);
    sStimParams.intWhite = WhiteIndex(sStimParams.ptrWindow);
    %% MAXIMIZE PRIORITY
    intOldPriority = 0;
    if boolDebug == 0
        intPriorityLevel=MaxPriority(ptrWindow);
        intOldPriority = Priority(intPriorityLevel);
    end
    %% get refresh rate
    dblStimFrameRate=Screen('FrameRate', ptrWindow);
    intStimFrameRate = round(dblStimFrameRate);
    dblStimFrameDur = mean(1/dblStimFrameRate);
    dblInterFlipInterval = Screen('GetFlipInterval', ptrWindow);
    if dblStimFrameDur/dblInterFlipInterval > 1.05 || dblStimFrameDur/dblInterFlipInterval < 0.95
        warning([mfilename ':InconsistentFlipDur'],sprintf('Something iffy with flip speed and monitor refresh rate detected; frame duration is %fs, while flip interval is %fs!',dblStimFrameDur,dblInterFlipInterval)); %#ok<SPWRN>
    end
    %% check escape
    if CheckEsc(),error([mfilename ':EscapePressed'],'Esc pressed; exiting');end
    %% prepare stimuli
    %set additional parameters
    sStimParams.dblPixelsPerDeg = mean([(sStimParams.intScreenWidth_pix/sStimParams.dblScreenWidth_deg) ...
        (sStimParams.intScreenHeight_pix/sStimParams.dblScreenHeight_deg)]);
    patchSizePix = ceil(sStimParams.dblStimulusSize_deg*sStimParams.dblPixelsPerDeg); %pix, rounded up

    % Ensure patch size is even for texture generation convenience
    if mod(patchSizePix,2)~=0
        patchSizePix = patchSizePix+1; % Make it even
    end
    sStimParams.intStimulusSize_pix = patchSizePix; % Update sStimParams with the final pixel size

    % Grating parameters
    spatialFreq = sStimParams.dblSpatialFrequency_cd; % cycles per degree
    temporalFreq = sStimParams.dblTemporalFrequency; % Hz

    % Calculate the texture dimension needed to avoid artifacts during drift.
    % We need to ensure that when we shift the source rectangle, we don't
    % run out of grating pattern, especially for diagonal drifts.
    % A common strategy is to make the texture large enough to always contain
    % the patch, even after shifting, and for the pattern to wrap around.
    % To allow for arbitrary shifts, the texture needs to be at least
    % patchSizePix + max_shift_in_one_cycle.
    
    % The spatial period in pixels for a grating with this spatial frequency:
    pixelsPerCycle = 1 / (spatialFreq / sStimParams.dblPixelsPerDeg);
    
    % We want the texture to be larger than the patch
    textureDim = max(patchSizePix + ceil(pixelsPerCycle), ceil(patchSizePix * 5)); 
    
    % Make textureDim even
    if mod(textureDim,2)~=0
        textureDim = textureDim+1;
    end

    % Pre-generate textures for each direction and two phases
    numDirections = length(sStimParams.vecDirections);
    texCollection = cell(numDirections, 2);
    
    % Shift magnitude per frame in pixels (along the drift direction)
    % This represents the distance the grating pattern moves per frame.
    dblShiftPerFrame = temporalFreq * pixelsPerCycle * (1/dblStimFrameRate);

    fprintf('Pre-generating grating textures for %d directions... ', numDirections);
    for intDir = 1:numDirections
        direction = sStimParams.vecDirections(intDir);
        theta_rad = deg2rad(mod(180-direction,360)); % Convert drift direction to radians
        
        % Grating lines are perpendicular to drift direction.
        % So, if drift direction is 'theta_rad', grating lines are at 'theta_rad + pi/2'.
        % The spatial components (fx, fy) define the orientation of the grating lines.
        % For a grating, the phase changes along the direction perpendicular to the lines.
        % This is the direction of drift.
        
        % This means:
        % fx = spatialFreq * cos(theta_rad)
        % fy = spatialFreq * sin(theta_rad)
        % The phase of the grating will be 2 * pi * (fx * X_texture + fy * Y_texture)
        
        % Create meshgrids for the larger texture
        % The texture coordinates should be centered around zero for proper rotation
        % and to calculate phase.
        texture_dim_deg = textureDim / sStimParams.dblPixelsPerDeg;
        x_tex_deg = linspace(-texture_dim_deg/2, texture_dim_deg/2, textureDim);
        y_tex_deg = linspace(-texture_dim_deg/2, texture_dim_deg/2, textureDim);
        [TextureX, TextureY] = meshgrid(x_tex_deg, y_tex_deg);

        % Grating calculation for Phase 0
        % The 'direction' variable here is the drift direction, so fx and fy define
        % the spatial components that lead to drift in that direction.
        phase_0 = 2 * pi * (spatialFreq * (cos(theta_rad) * TextureX + sin(theta_rad) * TextureY));
        grating_0 = 0.5 + 0.5 * sin(phase_0);
        texCollection{intDir, 1} = Screen('MakeTexture', ptrWindow, uint8(grating_0 * 255));

        % Grating calculation for Phase pi (180 degrees shifted)
        phase_pi = 2 * pi * (spatialFreq * (cos(theta_rad) * TextureX + sin(theta_rad) * TextureY)) + pi;
        grating_pi = 0.5 + 0.5 * sin(phase_pi);
        texCollection{intDir, 2} = Screen('MakeTexture', ptrWindow, uint8(grating_pi * 255));
    end
    fprintf('Done.\n');

    %get grid data
    intPosX = floor(sStimParams.intScreenWidth_pix/sStimParams.intStimulusSize_pix);
    intPosY = floor(sStimParams.intScreenHeight_pix/sStimParams.intStimulusSize_pix);
    vecPosX = linspace(0,(intPosX-1)*sStimParams.intStimulusSize_pix,intPosX);
    vecPosY = linspace(0,(intPosY-1)*sStimParams.intStimulusSize_pix,intPosY);
    %adjust so that grid is centered
    vecPosX = vecPosX + ((sStimParams.intScreenWidth_pix-(max(vecPosX)+sStimParams.intStimulusSize_pix))/2);
    vecPosY = vecPosY + ((sStimParams.intScreenHeight_pix-(max(vecPosY)+sStimParams.intStimulusSize_pix))/2);
    %create grid of left top corner locations
    [matPosX, matPosY] = meshgrid(vecPosX,vecPosY);
    %randomize locations for stimulus presentation
    vecStimulusConditionIDs = 1:(length(vecPosX)*length(vecPosY));
    vecStimulusPresentationIDs = []; for intRepeat = 1:sStimParams.intStimulusRepeats, ...
            vecStimulusPresentationIDs = [vecStimulusPresentationIDs vecStimulusConditionIDs(randperm(length(vecStimulusConditionIDs)))]; end %#ok<AGROW>
    intTotalTrials = length(vecStimulusPresentationIDs);
    %create randomized direction vector (indices into vecDirections)
    vecDirectionIdx = randi(numDirections, 1, intTotalTrials); % Randomly select an index for direction

    %randomize stimulus phase (0 or pi)
    vecPhaseIdx = randi([1 2], 1, intTotalTrials); % 1 for phase 0, 2 for phase pi
    
    %% build structEP
    %stimulus timing info
    initialBlank = sStimParams.dblSecsInitialBlank;
    trialDur = sStimParams.dblSecsPreBlank+sStimParams.dblSecsDuration+sStimParams.dblSecsPostBlank;
    endBlank = sStimParams.dblSecsEndBlank;
    totalLength = initialBlank+trialDur*intTotalTrials+endBlank;
    %build structEP, stim-based logs
    structEP = struct;
    structEP.strExpType = sStimParams.strStimType;
    structEP.intTrialNum = intTotalTrials;
    structEP.TrialNumber = nan(1,structEP.intTrialNum);
    structEP.dblStimFrameDur = dblStimFrameDur;
    structEP.dblInterFlipInterval = dblInterFlipInterval;
    structEP.dblStimulusSize_deg = repmat(sStimParams.dblStimulusSize_deg,[1 structEP.intTrialNum]);
    structEP.intStimulusSize_pix = repmat(sStimParams.intStimulusSize_pix,[1 structEP.intTrialNum]);
    structEP.ActOnSecs = nan(1,structEP.intTrialNum);
    structEP.ActOffSecs = nan(1,structEP.intTrialNum);
    structEP.ActStartSecs = nan(1,structEP.intTrialNum);
    structEP.ActStopSecs = nan(1,structEP.intTrialNum);
    structEP.ActOnNI = nan(1,structEP.intTrialNum);
    structEP.ActOffNI = nan(1,structEP.intTrialNum);
    structEP.vecDstRect= nan(4,structEP.intTrialNum);
    structEP.vecDirection = nan(1,structEP.intTrialNum);
    structEP.vecTexturePhase = nan(1,structEP.intTrialNum); % Log which phase was used (0 or 1 for pi)
    structEP.vecTrialStartSecs = initialBlank:trialDur:(totalLength-endBlank-1);
    structEP.vecTrialStimOnSecs = structEP.vecTrialStartSecs+sStimParams.dblSecsPreBlank;
    structEP.vecTrialStimOffSecs = structEP.vecTrialStimOnSecs+sStimParams.dblSecsDuration;
    structEP.vecTrialEndSecs = structEP.vecTrialStimOffSecs+sStimParams.dblSecsPostBlank;
    %% PRESENT STIMULI
    %wait for signal
    opts = struct;
    opts.Default = 'Start';
    opts.Interpreter = 'tex';
    strAns = questdlg('Would you like to start the stimulation?', ...
        'Start Stimulation', ...
        'Start','Cancel',opts);
    if ~strcmp(strAns,opts.Default)
        error([mfilename ':RunCancelled'],'Cancelling');
    end
    %hide cursor
    %HideCursor(ptrWindow);
    %set timers
    hTic = tic;
    dblLastFlip = Screen('Flip', ptrWindow);
    dblInitialFlip = dblLastFlip;
    %timestamp start
    structEP.strStartDate = getDate();
    structEP.strStartTime = getTime();
    %% calculate diode rectangle location
    if sStimParams.intCornerTrigger > 0
        %PTB uses [top,left, right, bottom] with origin at upper left
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
    %% wait initial blanking
    fprintf('Starting initial blank (dur=%.2fs) [%s]\n',sStimParams.dblSecsInitialBlank,getTime);
    dblInitialBlankDur = 0;
    while dblInitialBlankDur < (sStimParams.dblSecsInitialBlank-dblStimFrameDur)
        %do nothing
        Screen('FillRect',ptrWindow,sStimParams.intBackground);
        dblLastFlipFlip = Screen('Flip', ptrWindow,dblLastFlip + dblStimFrameDur/2);
        dblInitialBlankDur = dblLastFlipFlip - dblInitialFlip;
    end
    %% draw stimuli on screen
    intThisTrial = 1;
    while intThisTrial <= intTotalTrials && ~CheckEsc()
        %check whether script execution should be paused
        if CheckPause()
            fprintf('\n\n<strong>STIMULUS SCRIPT PAUSED!</strong> [%s]\n\n',getTime);
            WaitSecs(1);
            %get user input
            opts = struct;
            opts.Default = 'Restart';
            opts.Interpreter = 'tex';
            strAns = questdlg('STIMULUS SCRIPT PAUSED! Would you like to restart the stimulation?', ...
                'Restart Stimulation', ...
                'Restart',opts);
            WaitSecs(1);
        end
        % trial start; background
        Screen('FillRect',ptrWindow, sStimParams.intBackground);
        dblTrialStartFlip = Screen('Flip', ptrWindow);
        %fill DAQ with data
        if boolUseNI
            stop(objDaqOut);
            outputData1 = dblSyncLightMultiplier*cat(1,linspace(3, 3, 200)',linspace(0, 0, 50)');
            outputData2 = dblPupilLightMultiplier*linspace(3, 3, 250)';
            queueOutputData(objDaqOut,[outputData1 outputData2]);
            prepare(objDaqOut);
        end
        %get stimulus position on screen
        vecDstRect = [matPosX(vecStimulusPresentationIDs(intThisTrial)),matPosY(vecStimulusPresentationIDs(intThisTrial)),...
            matPosX(vecStimulusPresentationIDs(intThisTrial))+sStimParams.intStimulusSize_pix,...
            matPosY(vecStimulusPresentationIDs(intThisTrial))+sStimParams.intStimulusSize_pix];
        
        % Get the actual drift direction for this trial
        currentDriftDirection = sStimParams.vecDirections(vecDirectionIdx(intThisTrial));
        
        % Get the pre-generated texture for this direction and phase
        intGratingTex = texCollection{vecDirectionIdx(intThisTrial), vecPhaseIdx(intThisTrial)};
        
        %add to structEP
        structEP.vecDstRect(:,intThisTrial) = vecDstRect;
        structEP.vecDirection(intThisTrial) = currentDriftDirection;
        structEP.vecTexturePhase(intThisTrial) = vecPhaseIdx(intThisTrial) - 1; % 0 for phase 0, 1 for phase pi
        
        %get timing
        dblStartSecs = structEP.vecTrialStartSecs(intThisTrial);
        dblStimOnSecs = structEP.vecTrialStimOnSecs(intThisTrial);
        dblStimOffSecs = structEP.vecTrialStimOffSecs(intThisTrial);
        dblStimDurSecs = dblStimOffSecs - dblStimOnSecs;
        dblEndSecs = structEP.vecTrialEndSecs(intThisTrial);
        %display in command window
        fprintf('%d/%d location: [%d %d %d %d], drift direction: %d, phase: %d [%s]\n',intThisTrial,intTotalTrials,vecDstRect(1),vecDstRect(2),vecDstRect(3),vecDstRect(4),currentDriftDirection,structEP.vecTexturePhase(intThisTrial),getTime);
        %% wait pre-blanking
        dblPreBlankDur = 0;
        Screen('FillRect',ptrWindow,sStimParams.intBackground);
        dblLastFlip = Screen('Flip', ptrWindow);
        dblDAQ_Dur = 0; %measured time to set NI DAQ switch
        while dblPreBlankDur < (dblStimOnSecs - dblStartSecs - dblDAQ_Dur - dblStimFrameDur/2)
            %do nothing
            Screen('FillRect',ptrWindow, sStimParams.intBackground);
            dblLastFlip = Screen('Flip', ptrWindow, dblLastFlip + dblStimFrameDur/2);
            dblPreBlankDur = dblLastFlip - dblTrialStartFlip;
        end
        %% 250ms pulse at stim start
        if boolUseNI,startBackground(objDaqOut);end
        %% present grating patch
        refTime = tic;
        boolFirstFlip = false;
        intFrame = 1;
        currentDriftDirection = mod(180-currentDriftDirection,360);
        while toc(refTime)<(sStimParams.dblSecsDuration-2*dblStimFrameDur)
            % The drift velocity in pixels per second
            % This is the magnitude of the velocity vector
            driftSpeed_pix_per_sec = (temporalFreq / spatialFreq) * sStimParams.dblPixelsPerDeg;
            
            % Shift per frame in pixels along the drift direction
            dblShiftPerFrame = driftSpeed_pix_per_sec * dblStimFrameDur;

            % Get the shift vector in pixels per frame along the orientation
            dx = cosd(currentDriftDirection) * dblShiftPerFrame;
            dy = sind(currentDriftDirection) * dblShiftPerFrame;

            % Compute accumulated shifts, wrapping around the grating pattern
            % We wrap around 'pixelsPerCycle' for the x-component and
            % 'pixelsPerCycle' for the y-component (as the grating repeats along both axes
            % when rotated and considering diagonal drift).
            % However, the texture itself is rectangular, with dimensions `textureDim`.
            % The shifts should wrap within the bounds of the texture, effectively
            % simulating an infinite grating.
            
            % It's more intuitive to think of the shifts as wrapping around
            % the 'pixelsPerCycle' *along the direction of the grating's internal X/Y axes*.
            % Since the grating is oriented, the relevant "cycle" for wrapping might be
            % projected onto the texture's x or y axis.
            % A simpler approach is to let the shifts grow, and the `mod` applies to the texture dimensions.
            
            % For pre-oriented textures, the actual pattern repeats along the *drift direction*.
            % So, the effective period for wrapping is `pixelsPerCycle` along that drift vector.
            % Let's use that as the wrapping modulus for the shift *magnitude*.
            
            % The actual shift in the texture's (x, y) coordinates
            accumulatedShift_x_pixels = (intFrame - 1) * dx;
            accumulatedShift_y_pixels = (intFrame - 1) * dy;

            % Apply modulo to wrap the shifts within a cycle.
            % Since `Screen('DrawTexture')` can handle subpixel shifts, we don't `round` yet.
            % The `mod` should be by the 'pixelsPerCycle' for the effective pattern repeat.
            % However, for `srcRect` shifts, if the texture contains multiple cycles,
            % the shift can be greater than one cycle. It should wrap within the
            % texture's bounds (which are multiples of the cycle for tiling).

            % Let's use the `textureDim` for wrapping, and ensure textureDim is a multiple of pixelsPerCycle
            % or sufficiently large for smooth looping.
            
            % The simplest is to ensure the texture contains *more than* patchSizePix
            % and then wrap the *effective source rectangle origin* around a conceptual larger space.
            % Psychtoolbox's DrawTexture is powerful. We specify the `srcRect` and it will draw it.
            % The *wrap-around* needs to happen at the cycle level.
            
            % Re-evaluating the most robust way:
            % If we pre-calculate the grating, its phase `phi = 2*pi*(fx*X + fy*Y)`.
            % To drift, this becomes `phi_drift = 2*pi*(fx*X + fy*Y - TemporalFreq * time)`.
            % The `- TemporalFreq * time` is our effective phase offset.
            % We need to convert this phase offset into a pixel offset in the texture's original (X,Y) space.
            % This conversion is tricky because `fx` and `fy` are components.

            % Let's stick to the principle: the pre-oriented texture is static.
            % To achieve drift, we shift the `srcRect` to 'sample' a different part of this static texture.
            % The shift is always *along the direction of motion*.
            % So, a shift in pixels along the *drift direction* `currentDriftDirection`.
            
            % This means if the grating drifts right (0 deg), we shift srcRect horizontally.
            % If it drifts up (270 deg), we shift srcRect vertically downwards in texture coords.
            % If it drifts diagonally, we shift srcRect diagonally.

            % Let's refine `dblShiftPerFrame` and `srcRect` based on this:
            % dblShiftPerFrame is the *magnitude* of movement per frame in pixels along the drift vector.

            % Current position within the drift cycle (in pixels along the drift vector)
            currentDriftPosition = mod((intFrame - 1) * dblShiftPerFrame, pixelsPerCycle);

            % Convert this drift position (a scalar value along the drift vector)
            % into x and y offsets for the texture's source rectangle.
            % These offsets are relative to the *center* of the texture,
            % because the grating was generated around the center of its meshgrid.

            % Calculate the (x,y) components of this shift
            shiftX_from_center = cosd(currentDriftDirection) * currentDriftPosition;
            shiftY_from_center = sind(currentDriftDirection) * currentDriftPosition;

            srcRect_origin_x = (textureDim / 2) - (sStimParams.intStimulusSize_pix / 2) + shiftX_from_center;
            srcRect_origin_y = (textureDim / 2) - (sStimParams.intStimulusSize_pix / 2) + shiftY_from_center;

            vecSrcRect = [srcRect_origin_x, srcRect_origin_y, ...
                          srcRect_origin_x + sStimParams.intStimulusSize_pix, ...
                          srcRect_origin_y + sStimParams.intStimulusSize_pix];
            
            % No `angle` parameter here, as the texture is already oriented.
            Screen('DrawTexture',ptrWindow,intGratingTex,vecSrcRect,vecDstRect,0); 
            Screen('FillRect',ptrWindow,sStimParams.intWhite,vecDiodeRect);
            Screen('DrawingFinished',ptrWindow);
            dblLastFlip = Screen('Flip',ptrWindow,dblLastFlip+dblInterFlipInterval/2);
            intFrame = intFrame+1;
            % send trigger for stim start
            if ~boolFirstFlip
                %set switch
                boolFirstFlip = true;
                %log NI timestamp
                if boolUseSGL
%                     dblStimOnNI = GetScanCount(hSGL, intStreamNI)/dblSampFreqNI;
					dblStimOnNI = GetStreamSampleCount(hSGL, intStreamNI, strHostAddress)/dblSampFreqNI;
                else
                    dblStimOnNI = nan;
                end
                %log flip
                dblStimOnFlip = dblLastFlip;
            end
        end
        %back to background
        Screen('FillRect',ptrWindow, sStimParams.intBackground);
        dblStimOffFlip = Screen('Flip', ptrWindow, dblLastFlip+dblStimFrameDur/2);
        dblStimDur = dblStimOffFlip-dblStimOnFlip;
        %log NI timestamp
        if boolUseSGL
% 			dblStimOffNI = GetScanCount(hSGL, intStreamNI)/dblSampFreqNI;
            dblStimOffNI =  GetStreamSampleCount(hSGL, intStreamNI, strHostAddress)/dblSampFreqNI;
        else
            dblStimOffNI = nan;
        end
        %% wait post-blanking
        dblTrialDur = 0;
        Screen('FillRect',ptrWindow, sStimParams.intBackground);
        dblLastFlip = Screen('Flip', ptrWindow);
        while dblTrialDur < (dblEndSecs - dblStartSecs - 2*dblStimFrameDur)
            %do nothing
            Screen('FillRect',ptrWindow, sStimParams.intBackground);
            dblLastFlip = Screen('Flip', ptrWindow, dblLastFlip + dblStimFrameDur/2);
            dblTrialDur = dblLastFlip - dblTrialStartFlip;
        end
        dblPostBlankDur = dblLastFlip-dblStimOffFlip;
        %new stim-based output
        intStimNumber = intThisTrial;
        structEP.TrialNumber(intStimNumber) = intThisTrial;
        structEP.ActStartSecs(intStimNumber) = dblTrialStartFlip;
        structEP.ActOnSecs(intStimNumber) = dblStimOnFlip;
        structEP.ActOffSecs(intStimNumber) = dblStimOffFlip;
        structEP.ActStopSecs(intStimNumber) = dblLastFlip;
        structEP.ActOnNI(intStimNumber) = dblStimOnNI;
        structEP.ActOffNI(intStimNumber) = dblStimOffNI;
        %%increment trial number
        intThisTrial = intThisTrial+1;
    end
    %save data
    structEP.sStimParams = sStimParams;
    save(fullfile(strLogDir,strFilename), 'structEP','sParamsSGL');
    %show summary
    fprintf('Finished experiment & data saving at [%s], waiting for end blank (dur=%.2fs)\n',getTime,sStimParams.dblSecsEndBlank)
    %% wait end-blanking
    dblEndBlankDur = 0;
    while dblEndBlankDur < sStimParams.dblSecsEndBlank
        %do nothing
        Screen('FillRect',ptrWindow, sStimParams.intBackground);
        dblEndFlip = Screen('Flip', ptrWindow);
        dblEndBlankDur = dblEndFlip - dblLastFlip;
    end
    %clean up
    fprintf('\nExperiment is finished at [%s], closing down and cleaning up...\n',getTime);
    % Close all textures to free up memory
    for intDir = 1:numDirections
        Screen('Close', texCollection{intDir, 1});
        Screen('Close', texCollection{intDir, 2});
    end
    Screen('Close',ptrWindow);
    Screen('Close');
    Screen('CloseAll');
    ShowCursor;
    Priority(0);
    Screen('Preference', 'Verbosity',intOldVerbosity);
    %% close Daq IO
    if boolUseNI && ~(exist('sExpMeta','var') && isfield(sExpMeta,'objDaqOut'))
        try
            closeDaqOutput(objDaqOut);
        catch
        end
    end
catch ME
    %% check if escape
    if strcmp(ME.identifier,[mfilename ':EscapePressed'])
        fprintf('\nEscape pressed at [%s], closing down and cleaning up...\n',getTime);
        %save data
        structEP.sStimParams = sStimParams;
        save(fullfile(strLogDir,strFilename), 'structEP','sParamsSGL');
        %clean up
        fprintf('\nExperiment is finished at [%s], closing down and cleaning up...\n',getTime);
        % Close all textures to free up memory on error too
        if exist('texCollection', 'var')
            for intDir = 1:numDirections
                if ~isempty(texCollection{intDir,1}) && IsTexture(texCollection{intDir,1}), Screen('Close', texCollection{intDir, 1}); end
                if ~isempty(texCollection{intDir,2}) && IsTexture(texCollection{intDir,2}), Screen('Close', texCollection{intDir, 2}); end
            end
        end
        Screen('Close',ptrWindow);
        Screen('Close');
        Screen('CloseAll');
        ShowCursor;
        Priority(0);
        Screen('Preference', 'Verbosity',intOldVerbosity);
        %% close Daq IO
        if boolUseNI && ~(exist('sExpMeta','var') && isfield(sExpMeta,'objDaqOut'))
            try
                closeDaqOutput(objDaqOut);
            catch
            end
        end
    else
        %% catch me and throw me
        fprintf('\n\n\nError occurred! Trying to save data and clean up...\n\n\n');
        %save data
        structEP.sStimParams = sStimParams;
        if ~exist('sParamsSGL','var'),sParamsSGL=[];end
        save(fullfile(strLogDir,strFilename), 'structEP','sParamsSGL');
        %% catch me and throw me
        % Close all textures to free up memory on error too
        if exist('texCollection', 'var')
            for intDir = 1:numDirections
                if ~isempty(texCollection{intDir,1}) && IsTexture(texCollection{intDir,1}), Screen('Close', texCollection{intDir, 1}); end
                if ~isempty(texCollection{intDir,2}) && IsTexture(texCollection{intDir,2}), Screen('Close', texCollection{intDir, 2}); end
            end
        end
        Screen('Close');
        Screen('CloseAll');
        ShowCursor;
        Priority(0);
        Screen('Preference', 'Verbosity',intOldVerbosity);
        %% close Daq IO
        if boolUseNI && ~(exist('sExpMeta','var') && isfield(sExpMeta,'objDaqOut'))
            try
                closeDaqOutput(objDaqOut);
            catch
            end
        end
        %% show error
        rethrow(ME);
    end
end