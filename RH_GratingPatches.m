
%GRATINGSPATCHES
%show patches of (drifting) grating at different locations on the screen the screen
%
%Robin Haak, updated March 6th 2023

%% suppress m-lint warnings
%#ok<*MCCD,*NASGU,*ASGLU,*CTCH>
% clearvars -except sStimPresets sStimParamsSettings sExpMeta;
clear; close all; Screen('CloseAll');

%% define variables
fprintf('Starting %s [%s]\n',mfilename,getTime);
boolUseSGL = true;
boolUseNI = true;
boolDebug = false;

%defaults
dblPupilLightMultiplier = 1; %strength of infrared LEDs
dblSyncLightMultiplier = 0.5;
strHostAddress = '192.87.11.133'; %default host address
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
    sStimParamsSettings.strTempObjectPath = 'X:\JorritMontijn\';%X:\JorritMontijn\ or F:\Data\Temp\

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
    sStimParamsSettings.dblSpatialFrequency_cd = 1/sStimParamsSettings.dblStimulusSize_deg;%* %0.04; %spatial frequency in cycles per degree
    sStimParamsSettings.dblTemporalFrequency = 3; %temporal frequency in cycles per second
    sStimParamsSettings.dblSecsDuration = 0.5; %s
    sStimParamsSettings.dblSecsInitialBlank = 5; %s
    sStimParamsSettings.dblSecsPreBlank = 0.25; %s
    sStimParamsSettings.dblSecsPostBlank = 0.25; %0.5 %s
    sStimParamsSettings.dblSecsEndBlank = 5; %s
    sStimParamsSettings.intStimulusRepeats = 5; %per location (irrespective of direction or phase)
    sStimParamsSettings.vecDirections = [0 180]; %[0 90 180 270]; %drifting directions (only [0 90 180 270], otherwise it gets funky)
    sStimParamsSettings.str90Deg = '0 degrees is rightward motion; 90 degrees is downward motion';
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
    intStreamNI = -1;
    dblSampFreqNI = GetSampleRate(hSGL, intStreamNI);

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
    %sStimParams.dblPixelsPerDeg = sStimParams.intScreenWidth_pix/sStimParams.dblScreenWidth_deg;
    sStimParams.dblPixelsPerDeg = mean([(sStimParams.intScreenWidth_pix/sStimParams.dblScreenWidth_deg) ...
        (sStimParams.intScreenHeight_pix/sStimParams.dblScreenHeight_deg)]); %added this line(and commented the one before Feb 10th '23)
    sStimParams.intStimulusSize_pix = ceil(sStimParams.dblStimulusSize_deg*sStimParams.dblPixelsPerDeg); %pix, rounded up
    if mod(sStimParams.intStimulusSize_pix,2)==0
        intFullTexSize_pix = sStimParams.intStimulusSize_pix*2;
    else
        sStimParams.intStimulusSize_pix = sStimParams.intStimulusSize_pix-1;
        intFullTexSize_pix = (sStimParams.intStimulusSize_pix-1)*2;
    end
    sStimParams.dblSpatialFrequency_cp = sStimParams.dblSpatialFrequency_cd/sStimParams.dblPixelsPerDeg; %cycles/pix
    intPixelsPerCycle = ceil(1/sStimParams.dblSpatialFrequency_cp); %pix/cycle, rounded up
    dblSpatialFrequency_rad = sStimParams.dblSpatialFrequency_cp*2*pi; %radians

    %create grating textures
    vecX = meshgrid(-intFullTexSize_pix:intFullTexSize_pix+intPixelsPerCycle,1);
    matGrating = sStimParams.intBackground+(sStimParams.intWhite-sStimParams.intBackground)*cos(dblSpatialFrequency_rad*vecX);
    intGratingTex_1 = Screen('MakeTexture',ptrWindow,matGrating);
    intGratingTex_2 = Screen('MakeTexture',ptrWindow,sStimParams.intWhite-matGrating);

    %create mask
    matMask = ones(2*intFullTexSize_pix+1,2*intFullTexSize_pix+1,2)*sStimParams.intBackground;
    [matX,matY] = meshgrid(-intFullTexSize_pix:intFullTexSize_pix,-intFullTexSize_pix:intFullTexSize_pix);
    matMask(:,:,2) = ~buildCircularCosineRampPix(size(matX,1),sStimParams.intStimulusSize_pix,2);
    intMaskTex=Screen('MakeTexture',ptrWindow,matMask);

    %calculate shift per frame
    intPixelsPerCycle = 1/sStimParams.dblSpatialFrequency_cp; %pix/cycle, re-calculate w/o ceil() to prevent errors later on
    dblShiftPerFrame = sStimParams.dblTemporalFrequency*intPixelsPerCycle*(1/dblStimFrameRate);

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


    %create randomized direction vector
    vecDirection = repmat(sStimParams.vecDirections,[1 floor(intTotalTrials/length(sStimParams.vecDirections))]);
    vecDirection = [vecDirection vecDirection(1:(intTotalTrials-length(vecDirection)))];
    vecDirection = vecDirection(randperm(length(vecDirection)));

    %randomize stimulus texture
    vecTexture = repmat([intGratingTex_1 intGratingTex_2],[1 intTotalTrials/2]);
    vecTexture = vecTexture(randperm(length(vecTexture)));

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
    structEP.vecTexture = nan(1,structEP.intTrialNum);
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

        %so that 0 = rightward motion
        dblAngle = vecDirection(intThisTrial)+180;
        if dblAngle >= 360, dblAngle = dblAngle-360; end

        %get texture
        intGratingTex = vecTexture(intThisTrial);

        %add to structEP
        structEP.vecDstRect(:,intThisTrial) = vecDstRect;
        structEP.vecDirection(intThisTrial) = vecDirection(intThisTrial);
        if vecTexture(intThisTrial) == intGratingTex_1
            structEP.vecTexture(intThisTrial) = 0;
        elseif vecTexture(intThisTrial) == intGratingTex_2
            structEP.vecTexture(intThisTrial) = 0;
        end

        %get timing
        dblStartSecs = structEP.vecTrialStartSecs(intThisTrial);
        dblStimOnSecs = structEP.vecTrialStimOnSecs(intThisTrial);
        dblStimOffSecs = structEP.vecTrialStimOffSecs(intThisTrial);
        dblStimDurSecs = dblStimOffSecs - dblStimOnSecs;
        dblEndSecs = structEP.vecTrialEndSecs(intThisTrial);

        %display in command window
        fprintf('%d/%d location: [%d %d %d %d], direction: %d [%s]\n',intThisTrial,intTotalTrials,vecDstRect(1),vecDstRect(2),vecDstRect(3),vecDstRect(4),dblAngle,getTime);

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
        refTimeLocal = tic;
        while toc(refTime)<(sStimParams.dblSecsDuration-2*dblStimFrameDur)
            dblOffsetX = mod((intFrame-1)*dblShiftPerFrame,intPixelsPerCycle);
            vecSrcRect=[dblOffsetX 0 dblOffsetX + sStimParams.intStimulusSize_pix sStimParams.intStimulusSize_pix];
            Screen('DrawTexture',ptrWindow,intGratingTex,vecSrcRect,vecDstRect,dblAngle);
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
                    dblStimOnNI = GetScanCount(hSGL, intStreamNI)/dblSampFreqNI;
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
            dblStimOffNI = GetScanCount(hSGL, intStreamNI)/dblSampFreqNI;
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
        Screen('Close');
        Screen('CloseAll');
        ShowCursor;
        Priority(0);
        Screen('Preference', 'Verbosity',intOldVerbosity);

        %% close Daq IO
        if boolUseNI && ~(exist('sExpMeta','var') && isfield(sExpMeta,'objDaqOut'))
            try
                closeDaqOutput(objDaqOut);ds
            catch
            end
        end

        %% show error
        rethrow(ME);
    end
end
%end
