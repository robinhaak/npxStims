%RH_MovingDots

%% suppress m-lint warnings
%#ok<*MCCD,*NASGU,*ASGLU,*CTCH>
clearvars -except sStimPresets sStimParamsSettings sExpMeta;

%% define variables
fprintf('Starting %s [%s]\n',mfilename,getTime);
intStimSet = 1; % 1=dot_grid (for RF mapping), 2=dot_variations, 3=dot_speeds, 4=dot_reversal
boolUseSGL = false;
boolUseNI = false;
boolDebug = false;
if exist('sExpMeta','var')
    %defaults
    dblPupilLightMultiplier = 1; %strength of infrared LEDs
    dblSyncLightMultiplier = 0.5;
    strHostAddress = '192.87.11.133'; %default host address
    objDaqOut = [];

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
if ~exist('sStimParamsSettings','var') || isempty(sStimParamsSettings) || ~(strcmpi(sStimParamsSettings.strStimType,'MovingDots'))
    %general
    sStimParamsSettings = struct;
    sStimParamsSettings.strStimType = 'MovingDots';
    sStimParamsSettings.strOutputPath = 'C:\_Data\Exp'; %appends date
    sStimParamsSettings.strTempObjectPath = 'X:\JorritMontijn\';%X:\JorritMontijn\ or F:\Data\Temp\

    %visual space parameters
    sStimParamsSettings.dblSubjectPosX_cm = 0; % cm; relative to center of screen
    sStimParamsSettings.dblSubjectPosY_cm = 0; % cm; relative to center of screen, -3.5
    sStimParamsSettings.dblScreenDistance_cm = 14; % cm; measured, 17
    sStimParamsSettings.vecUseMask = 0; %[1] if mask to emulate retinal-space, [0] use screen-space

    %receptive field size & location parameters (values are in pixels and not dva)
    sStimParamsSettings.intRespPosX_pix = 0; % pix; x screen pos. of reponse zone
    sStimParamsSettings.intRespPosY_pix = 0; % pix; y screen pos.
    sStimParamsSettings.intRespSize_pix = 0; % pix; diameter of the response zone

    %screen variables
    sStimParamsSettings.intCornerTrigger = 2; % integer switch; 0=none,1=upper left, 2=upper right, 3=lower left, 4=lower right
    sStimParamsSettings.dblCornerSize = 1/30; % fraction of screen width
    sStimParamsSettings.dblScreenWidth_cm = 51; % cm; measured [51]
    sStimParamsSettings.dblScreenHeight_cm = 29; % cm; measured [29]
    sStimParamsSettings.dblScreenWidth_deg = 2*atand(sStimParamsSettings.dblScreenWidth_cm/(2*sStimParamsSettings.dblScreenDistance_cm));
    sStimParamsSettings.dblScreenHeight_deg = 2*atand(sStimParamsSettings.dblScreenHeight_cm/(2*sStimParamsSettings.dblScreenDistance_cm));
    sStimParamsSettings.intUseScreen = 2; %which screen to use

    %stimulus control variables
    sStimParamsSettings.intUseDaqDevice = 1; %ID of DAQ device
    sStimParamsSettings.intUseParPool = 0; %number of workers in parallel pool; [2]
    sStimParamsSettings.intUseGPU = 0;
    sStimParamsSettings.intAntiAlias = 0;
    sStimParamsSettings.str90Deg = '0 degrees is rightward motion; 90 degrees is downward motion';
    sStimParamsSettings.dblBackground = 0.5; %background intensity (dbl, [0 1])
    sStimParamsSettings.intBackground = round(mean(sStimParamsSettings.dblBackground)*255);
    sStimParamsSettings.dblSecsInitialBlank = 5; % s
    sStimParamsSettings.dblSecsPreBlank = 0.25; % s
    sStimParamsSettings.vecSecsPostBlank = [0.75 1.25]; % s, random within range
    sStimParamsSettings.dblSecsEndBlank = 5; %s
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

%keep stuff the same as other Acquipix-based scripts
sStimParams = sStimParamsSettings;
if boolDebug == 1
    intUseScreen = 0;
else
    intUseScreen = sStimParamsSettings.intUseScreen;
end

%% set output locations for logs
strOutputPath = sStimParamsSettings.strOutputPath;
strTempObjectPath = sStimParamsSettings.strTempObjectPath;
strThisFilePath = mfilename('fullpath');
[strFilename,strLogDir,strTempDir,strTexDir] = RE_assertPaths(strOutputPath,strRecording,strTempObjectPath,strThisFilePath);
fprintf('Saving output in directory %s; loading textures from %s\n',strLogDir,strTexDir); %no textures are loaded for this script

%% query user for stimulus set & location of the response zone/receptive field

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

    %% check disk space available
    strDataDirSGL = GetDataDir(hSGL);
    jFileObj = java.io.File(strDataDirSGL);
    dblFreeGB = (jFileObj.getFreeSpace)/(1024^3);
    if dblFreeGB < 100,warning([mfilename ':LowDiskSpace'],'Low disk space available (%.0fGB) for Neuropixels data (dir: %s)',dblFreeGB,strDataDirSGL);end

    %% prepare real-time stream
    %retrieve some parameters
    dblSampFreqNI = GetSampleRate(hSGL, intStreamNI);
    dblNI2V = (sParamsSGL.niAiRangeMax) / (double(intmax('int16'))/2);

    %assign data
    sStream = SC_populateStreamCoreStructure([]);
    sStream.NI2V = dblNI2V;
    sStream.intStreamNI = intStreamNI;
    sStream.intStimOnsetChanNI=intStimOnsetChanNI;
    sStream.intRunningChanNI=intRunningChanNI;
    sStream.dblSampFreqNI=dblSampFreqNI;
    sStream.hSGL = hSGL;
else
    strHostAddress = '';
    hSGL = [];
    sParamsSGL = struct;
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
        objDaqOut = openDaqOutput(sStimParams.intUseDaqDevice);
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
    if boolDebug == 1,vecInitRect = [0 0 640 640];else,vecInitRect = [];end
    try
        Screen('Preference', 'SkipSyncTests', 0);
        [ptrWindow,vecRect] = Screen('OpenWindow', intUseScreen,sStimParams.intBackground,vecInitRect);
    catch ME
        warning([mfilename ':ErrorPTB'],'Psychtoolbox error, attempting with sync test skip [msg: %s]',ME.message);
        Screen('Preference', 'SkipSyncTests', 1);
        [ptrWindow,vecRect] = Screen('OpenWindow', intUseScreen,sStimParams.intBackground,vecInitRect);
    end
    %window variables
    sStimParams.ptrWindow = ptrWindow;
    sStimParams.vecRect = vecRect;
    sStimParams.intScreenWidth_pix = vecRect(3)-vecRect(1);
    sStimParams.intScreenHeight_pix = vecRect(4)-vecRect(2);
    sStimParams.dblPixelsPerDeg = mean([(sStimParams.intScreenWidth_pix/sStimParams.dblScreenWidth_deg) ...
        (sStimParams.intScreenHeight_pix/sStimParams.dblScreenHeight_deg)]);
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
    sStimParams.intStimFrameRate = intStimFrameRate;
    dblStimFrameDur = mean(1/dblStimFrameRate);
    dblInterFlipInterval = Screen('GetFlipInterval', ptrWindow);
    if dblStimFrameDur/dblInterFlipInterval > 1.05 || dblStimFrameDur/dblInterFlipInterval < 0.95
        warning([mfilename ':InconsistentFlipDur'],...
            sprintf('Something iffy with flip speed and monitor refresh rate detected; frame duration is %fs, while flip interval is %fs!',...
            dblStimFrameDur,dblInterFlipInterval)); %#ok<SPWRN>
    end

    %% check escape
    if CheckEsc(),error([mfilename ':EscapePressed'],'Esc pressed; exiting');end

    %% create presentation vectors
    %get info for each individual stimulus condition
    [sAllDots,sDotParams] = RH_CreateDotTrajectories(intStimSet,sStimParams);

    %build structEP (trial-based output)
    structEP = struct;
    structEP.intTrialNum = length(sAllDots.stimID)*sDotParams.intRepsPerCondition;
    structEP.TrialNum = nan(1,structEP.intTrialNum);
    structEP.dblStimFrameRate = dblStimFrameRate;
    structEP.intStimFrameRate = intStimFrameRate;
    structEP.vecStimID = [];
    for intRep = 1:sDotParams.intRepsPerCondition
        structEP.vecStimID = [structEP.vecStimID sAllDots.stimID(randperm(length(sAllDots.stimID)))];
    end
    for intTrial = 1:structEP.intTrialNum
        structEP.vecBoundingRect{intTrial} = sAllDots.vecBoundingRect{sAllDots.stimID==structEP.vecStimID(intTrial)};
        structEP.vecColor{intTrial} = sAllDots.vecColor{sAllDots.stimID==structEP.vecStimID(intTrial)};
        structEP.vecDirection(intTrial) = sAllDots.vecDirection(sAllDots.stimID==structEP.vecStimID(intTrial));
        structEP.vecSpeed_deg(intTrial) = sAllDots.vecSpeed_deg(sAllDots.stimID==structEP.vecStimID(intTrial));
        structEP.vecSpeed_pix(intTrial) = sAllDots.vecSpeed_pix(sAllDots.stimID==structEP.vecStimID(intTrial));
        structEP.vecStimName{intTrial} = sAllDots.vecStimName{sAllDots.stimID==structEP.vecStimID(intTrial)};
    end
    structEP.ActOnSecs = nan(1,structEP.intTrialNum);
    structEP.ActOffSecs = nan(1,structEP.intTrialNum);
    structEP.ActStartSecs = nan(1,structEP.intTrialNum);
    structEP.ActStopSecs = nan(1,structEP.intTrialNum);
    structEP.ActOnNI = nan(1,structEP.intTrialNum);
    structEP.ActOffNI = nan(1,structEP.intTrialNum);

    %% PRESENT STIMULI
    %wait for user input
    opts = struct;
    opts.Default = 'Start';
    opts.Interpreter = 'tex';
    strAns = questdlg('Would you like to start the stimulation?', ...
        'Start Stimulation', ...
        'Start','Cancel',opts);
    if ~strcmp(strAns,opts.Default)
        error([mfilename ':RunCancelled'],'Cancelling');
    end

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
    while intThisTrial <= structEP.intTrialNum && ~CheckEsc()

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

        %get trial-specific parameters
        strStimName = structEP.vecStimName{intThisTrial};
        vecBoundingRect = structEP.vecBoundingRect{intThisTrial};
        intFrames = length(vecBoundingRect);
        dblDirection = structEP.vecDirection(intThisTrial);
        intSpeed_ppf = round(structEP.vecSpeed_pix(intThisTrial)/sStimParams.intStimFrameRate);
        vecColor = structEP.vecColor{intThisTrial};
        dblSecsPostBlank = sStimParams.vecSecsPostBlank(1)+(sStimParams.vecSecsPostBlank(2)-sStimParams.vecSecsPostBlank(1))*rand();

        %display in command window
        fprintf('%d/%d "%s" dot, dir=%d speed=%d (PBdur=%.2fs) [%s]\n',intThisTrial,structEP.intTrialNum,strStimName,dblDirection,...
            intSpeed_ppf,dblSecsPostBlank,getTime);

        %% wait pre-blanking
        dblPreBlankDur = 0;
        Screen('FillRect',ptrWindow,sStimParams.intBackground);
        dblLastFlip = Screen('Flip', ptrWindow);
        dblDAQ_Dur = 0; %measured time to set NI DAQ switch
        while dblPreBlankDur < (sStimParams.dblSecsPreBlank-dblDAQ_Dur-dblStimFrameDur)
            %do nothing
            Screen('FillRect',ptrWindow,sStimParams.intBackground);
            dblLastFlip = Screen('Flip',ptrWindow, dblLastFlip+dblStimFrameDur/2);
            dblPreBlankDur = dblLastFlip-dblTrialStartFlip;
        end

        %% 250ms pulse at stim start
        if boolUseNI,startBackground(objDaqOut);end

        %% present moving dot
        refTime = tic;
        boolFirstFlip = false;
        intFrame = 1;
        refTimeLocal = tic;

        while intFrame < intFrames
            Screen('FillOval',ptrWindow,vecColor(intFrame),vecBoundingRect(:,intFrame)); %dot
            Screen('FillRect',ptrWindow,sStimParams.intWhite,vecDiodeRect); %diode
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
        dblPostBlankDur = 0;
        Screen('FillRect',ptrWindow,sStimParams.intBackground);
        dblLastFlip = Screen('Flip',ptrWindow);
        dblStimOffFlip = dblLastFlip;
        while dblPostBlankDur < (sStimParams.dblSecsPreBlank-dblStimFrameDur)
            %do nothing
            Screen('FillRect',ptrWindow,sStimParams.intBackground);
            dblLastFlip = Screen('Flip',ptrWindow,dblLastFlip+dblStimFrameDur/2);
            dblPostBlankDur = dblLastFlip-dblStimOffFlip;
        end

        % new stim-based output
        intStimNumber = intThisTrial;
        structEP.TrialNumber(intStimNumber) = intThisTrial;
        structEP.ActStartSecs(intStimNumber) = dblTrialStartFlip;
        structEP.ActOnSecs(intStimNumber) = dblStimOnFlip;
        structEP.ActOffSecs(intStimNumber) = dblStimOffFlip;
        structEP.ActStopSecs(intStimNumber) = dblLastFlip;
        structEP.ActOnNI(intStimNumber) = dblStimOnNI;
        structEP.ActOffNI(intStimNumber) = dblStimOffNI;

        %increment trial
        intThisTrial = intThisTrial+1;
    end

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
    if strcmp(ME.identifier,'RH_MovingDots:EscapePressed')
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
                closeDaqOutput(objDaqOut);
            catch
            end
        end

        %% show error
        rethrow(ME);
    end
end

