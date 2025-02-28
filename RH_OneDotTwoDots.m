% STRUCTEP = RH_ONEDOTTWODOTS()
%
% Robin Haak
%
% history:
% 31 Aug 2023
%   - created by Robin Haak
% 25 Sep 2023
%   - minor changes to code in prep for exp
% 30 May 2024
%	- added three 'hidden' stim sets (5-7)
% 28 Feb 2025
%   - changed to work with new spikeGLX version


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
strHostAddress = '192.87.11.84';%'192.87.11.133'; %default host address
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

%% get user input for recording name
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
    sStimParamsSettings.strStimType = 'OneDotTwoDots';
    sStimParamsSettings.strOutputPath = 'C:\_Data\Exp'; %appends date
    sStimParamsSettings.strTempObjectPath = 'C:\temp';%'X:\JorritMontijn\';%X:\JorritMontijn\ or F:\Data\Temp\

    %visual space parameters
    %     sStimParamsSettings.dblSubjectPosX_cm = 0; % cm; relative to center of screen
    %     sStimParamsSettings.dblSubjectPosY_cm = 0; % cm; relative to center of screen
    sStimParamsSettings.dblScreenDistance_cm = 15; % cm; measured [23]

    %receptive field size & location parameters (values are in pixels and not dva)
    sStimParamsSettings.intRfPosX_pix = nan; % pix; x screen pos. of RF
    sStimParamsSettings.intRfPosY_pix = nan; % pix; y screen pos.
    sStimParamsSettings.intRfSize_pix = nan; % pix; diameter

    %screen variables
    sStimParamsSettings.intCornerTrigger = 2; % integer switch; 0=none,1=upper left, 2=upper right, 3=lower left, 4=lower right
    sStimParamsSettings.dblCornerSize = 1/30; % fraction of screen width
    sStimParamsSettings.dblScreenWidth_cm = 39.6;%70.15; %51; % cm; measured [51]
    sStimParamsSettings.dblScreenHeight_cm = 70.15; %39.6; %29; % cm; measured [29]
    sStimParamsSettings.dblScreenWidth_deg = 2*atand(sStimParamsSettings.dblScreenWidth_cm/(2*sStimParamsSettings.dblScreenDistance_cm));
    sStimParamsSettings.dblScreenHeight_deg = 2*atand(sStimParamsSettings.dblScreenHeight_cm/(2*sStimParamsSettings.dblScreenDistance_cm));
    sStimParamsSettings.intUseScreen = 2; %which screen to use

    %stimulus control variables
    sStimParamsSettings.intReps = 5;
    sStimParamsSettings.intUseDaqDevice = 1; %ID of DAQ device
    sStimParamsSettings.intUseParPool = 0; %number of workers in parallel pool; [2]
    sStimParamsSettings.intUseGPU = 0;
    sStimParamsSettings.dblBackground = 0.5; %background intensity (dbl, [0 1])
    sStimParamsSettings.intBackground = round(mean(sStimParamsSettings.dblBackground)*255);
    sStimParamsSettings.dblSecsInitialBlank = 5; % s
    sStimParamsSettings.dblSecsPreBlank = 1;
    sStimParamsSettings.dblSecsPostBlank = 1.7; %1.5; % total ITI is ~2.5s
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
if boolDebug
    intUseScreen = 0;
else
    intUseScreen = sStimParamsSettings.intUseScreen;
end

%% set output locations for logs
strOutputPath = sStimParamsSettings.strOutputPath;
strTempObjectPath = sStimParamsSettings.strTempObjectPath;
strThisFilePath = mfilename('fullpath');
[strFilename,strLogDir,strTempDir,strTexDir] = RE_assertPaths(strOutputPath,strRecording,strTempObjectPath,strThisFilePath);
fprintf('Saving output in directory %s;\n',strLogDir); %no textures are loaded for this script

%% query user for stimulus set & approx. location of the receptive field
%stimulus set
fprintf(['\n--<strong>Select stimulus set</strong>--\nAvailable sets:\n' ...
    '(1)"inside-outside"\n(2)"inside <strong>control</strong>"\n(3)"outside-outside"\n(4)"outside <strong>control</strong>"...\n\n']);
intStimSet  = input('set: ');
%5 and 6 are 'hidden'sets
if sum(intStimSet==1:7)<1,error('Non-existent stim set');end
sStimParams.intStimSet = intStimSet;

%receptive field center
fprintf('\n--<strong>Receptive field center</strong>--\n');
sStimParams.intRfPosX_pix = round(input('intRfPosX_pix= ')); % pix; we only need x-coordinate

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
    fprintf('\nStarting PsychToolBox extension...\n');
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

    %calibrate monitor
	if exist('C:\Users\neuropixels\Desktop\GammaCorrection\gammaTable_20230926.mat','file')
		load("C:\Users\neuropixels\Desktop\GammaCorrection\gammaTable_20230926.mat");
		Screen('LoadNormalizedGammaTable', ptrWindow, gammaTable*[1 1 1]);
		fprintf('Calibrated monitor!\n')
		sStimParams.boolCalibMonitor = true;
	else
		Screen('LoadNormalizedGammaTable', ptrWindow, linspace(0,1,256)'*[1 1 1]);
		sStimParams.boolCalibMonitor = false;
	end
		
    %window variables
	sStimParams.ptrWindow = ptrWindow;
	sStimParams.vecRect = vecRect;
	sStimParams.intScreenWidth_pix = vecRect(3)-vecRect(1);
	sStimParams.intScreenHeight_pix = vecRect(4)-vecRect(2);
	
	%     sStimParams.dblPixelsPerDeg = mean([(sStimParams.intScreenWidth_pix/sStimParams.dblScreenWidth_deg) ...
	%         (sStimParams.intScreenHeight_pix/sStimParams.dblScreenHeight_deg)]);
	
% 	sStimParams.dblPixelsPerDeg = sStimParams.intScreenHeight_pix/sStimParams.dblScreenHeight_deg;
	sStimParams.dblPixelsPerDeg = sStimParams.intScreenWidth_pix/sStimParams.dblScreenWidth_deg;

	sStimParams.intWhite = WhiteIndex(sStimParams.ptrWindow);

    %% MAXIMIZE PRIORITY
    intOldPriority = 0;
    if ~boolDebug
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
    [sDots,sStimParams] = RH_CreateOneDotTwoDots(intStimSet,sStimParams);

    %build structEP
    structEP = struct;
    structEP.strExpType = sStimParams.strStimType;
    structEP.strStimSet = sDots.strStimSet;
    structEP.dblStimFrameDur = dblStimFrameDur;
    structEP.dblInterFlipInterval = dblInterFlipInterval;

    vecStimID = [repmat(sDots.vecBlockStructure,[1 sStimParams.intReps]) sDots.vecBlockStructure(1)];
    structEP.vecStimID = vecStimID;
    structEP.intTrialNum = length(structEP.vecStimID);
    structEP.TrialNum = nan(1,structEP.intTrialNum);
    for intTrial = 1:structEP.intTrialNum
        structEP.cellBoundingRect{intTrial} = sDots.cellBoundingRect{vecStimID(intTrial)};
        structEP.cellColor{intTrial} = sDots.cellColor{vecStimID(intTrial)};
        structEP.vecSpeed_deg(intTrial) = sDots.vecSpeed_deg(vecStimID(intTrial));
        structEP.vecSpeed_pix(intTrial) = sDots.vecSpeed_pix(vecStimID(intTrial));
    end
    %pre-allocate
    structEP.ActOnSecs = nan(1,structEP.intTrialNum);
    structEP.ActOffSecs = nan(1,structEP.intTrialNum);
    structEP.ActStartSecs = nan(1,structEP.intTrialNum);
    structEP.ActStopSecs = nan(1,structEP.intTrialNum);
    structEP.ActOnNI = nan(1,structEP.intTrialNum);
    structEP.ActOffNI = nan(1,structEP.intTrialNum);

    %% PRESENT STIMULI
    %wait for user input
    sOptions = struct;
    sOptions.Default = 'Start';
    sOptions.Interpreter = 'tex';
    strAns = questdlg('Would you like to start the stimulation?', ...
        'Start Stimulation', ...
        'Start','Cancel',sOptions);
    if ~strcmp(strAns,sOptions.Default)
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
    while dblInitialBlankDur < (sStimParams.dblSecsInitialBlank-dblStimFrameDur/2)
        %do nothing
        Screen('FillRect',ptrWindow,sStimParams.intBackground);
        dblLastFlipFlip = Screen('Flip', ptrWindow,dblLastFlip + dblStimFrameDur/2);
        dblInitialBlankDur = dblLastFlipFlip - dblInitialFlip;
    end

    %% draw stimuli on screen
    intThisTrial = 1;
    while intThisTrial <= structEP.intTrialNum && ~CheckEsc()

        %check whether script execution should be paused
        if CheckPause()
            fprintf('\n\n<strong>STIMULUS SCRIPT PAUSED!</strong> [%s]\n\n',getTime);
            WaitSecs(1);
            %get user input
            sOptions = struct;
            sOptions.Default = 'Restart';
            sOptions.Interpreter = 'tex';
            strAns = questdlg('STIMULUS SCRIPT PAUSED! Would you like to restart the stimulation?', ...
                'Restart Stimulation', ...
                'Restart',sOptions);
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

        %get trial-specific parameters
        vecThisBoundingRect = structEP.cellBoundingRect{intThisTrial};
        intNumFrames = length(vecThisBoundingRect);
        vecThisColor = structEP.cellColor{intThisTrial};

        %display in command window
        fprintf('[%s] Showing stim %d/%d\n',getTime,intThisTrial,structEP.intTrialNum);


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

        %% present stimulus
        refTime = tic;
        boolFirstFlip = false;
        intFrame = 1;
        while intFrame < intNumFrames
            Screen('FillOval',ptrWindow,vecThisColor(intFrame),vecThisBoundingRect(:,intFrame)); %dot
            Screen('FillRect',ptrWindow,sStimParams.intWhite,vecDiodeRect); %diode
            dblLastFlip = Screen('Flip',ptrWindow,dblLastFlip+dblInterFlipInterval/2);
            %send trigger for stim start
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
            %increment frame
            intFrame = intFrame+1;
        end

        %back to background
        Screen('FillRect',ptrWindow, sStimParams.intBackground);
        dblLastFlip = Screen('Flip', ptrWindow, dblLastFlip+dblStimFrameDur/2);

        %log NI timestamp
        if boolUseSGL
            % 			dblStimOffNI = GetScanCount(hSGL, intStreamNI)/dblSampFreqNI;
            dblStimOffNI =  GetStreamSampleCount(hSGL, intStreamNI, strHostAddress)/dblSampFreqNI;
        else
            dblStimOffNI = nan;
        end

        %get stim duration
        dblStimOffFlip = dblLastFlip;
        dblStimDur = dblStimOffFlip-dblStimOnFlip;

        %% wait post-blanking
        dblPostBlankDur = 0;
        Screen('FillRect',ptrWindow,sStimParams.intBackground);
        dblLastFlip = Screen('Flip',ptrWindow);
        while dblPostBlankDur < (sStimParams.dblSecsPostBlank-dblStimFrameDur*2)
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

    %% end
    %save data
    structEP.sStimParams = sStimParams;
    structEP.sDots = sDots;
    save(fullfile(strLogDir,strFilename), 'structEP','sParamsSGL');

    %show trial summary
    fprintf('Finished experiment & data saving at [%s], waiting for end blank (dur=%.3fs)\n',getTime,sStimParams.dblSecsEndBlank)

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
        structEP.sDots = sDots;
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
        structEP.sDots = sDots;
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
