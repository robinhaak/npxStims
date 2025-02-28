
%function structEP = RH_OptoStim

%% suppress m-lint warnings
%#ok<*MCCD,*NASGU,*ASGLU,*CTCH>
clear; %clearvars -except sStimPresets sStimParamsSettings sExpMeta;

%% define variables
fprintf('Starting %s [%s]\n',mfilename,getTime);
intStimSet = 1;
boolUseSGL = true;
boolUseNI = true;
boolDebug = false;
boolGammaCorrectScreen = true;
%defaults
dblPupilLightMultiplier = 1; %strength of infrared LEDs
dblLaserTriggerMultiplier = 0.5;
strHostAddress = '192.87.10.228';%'192.87.11.133'; %default host address
objDaqOut = [];
if exist('sExpMeta','var')

    %expand structure
    if isfield(sExpMeta,'dblPupilLightMultiplier'),dblPupilLightMultiplier=sExpMeta.dblPupilLightMultiplier;end
    if isfield(sExpMeta,'dblLaserTriggerMultiplier'),dblLaserTriggerMultiplier=sExpMeta.dblLaserTriggerMultiplier;end
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
if ~exist('sStimParamsSettings','var') || isempty(sStimParamsSettings) || ~(strcmpi(sStimParamsSettings.strStimType,'OptoStim'))
    %general
    sStimParamsSettings = struct;
    sStimParamsSettings.strStimType = 'OptoStim';
    sStimParamsSettings.strOutputPath = 'C:\_Data\Exp'; %appends date
    sStimParamsSettings.strTempObjectPath = 'C:\temp'; %'X:\JorritMontijn\';%X:\JorritMontijn\ or F:\Data\Temp\

    %NB, all we show is a grey screen so this not super relevant
    %visual space parameters
    sStimParamsSettings.dblSubjectPosX_cm = 0; % cm; relative to center of screen
    sStimParamsSettings.dblSubjectPosY_cm = 0; %-2.5; % cm; relative to center of screen, -3.5
    sStimParamsSettings.dblScreenDistance_cm = 24; %17; % cm; measured, 14

    %screen variables
    sStimParamsSettings.dblScreenWidth_cm = 51; % cm; measured [51]
    sStimParamsSettings.dblScreenHeight_cm = 29; % cm; measured [29]
    sStimParamsSettings.dblScreenWidth_deg = 2 * atand(sStimParamsSettings.dblScreenWidth_cm / (2 * sStimParamsSettings.dblScreenDistance_cm));
    sStimParamsSettings.dblScreenHeight_deg = 2 * atand(sStimParamsSettings.dblScreenHeight_cm / (2 * sStimParamsSettings.dblScreenDistance_cm));
    sStimParamsSettings.intUseScreen = 2; %which screen to use

    %stimulus control variables
    sStimParamsSettings.intUseDaqDevice = 1; %ID of DAQ device
    sStimParamsSettings.intUseParPool = 0; %number of workers in parallel pool; [2]
    sStimParamsSettings.intUseGPU = 0;
    sStimParamsSettings.vecBackgrounds = 0.5; %background intensity (dbl, [0 1])
    sStimParamsSettings.intBackground = round(mean(sStimParamsSettings.vecBackgrounds)*255);

    sStimParamsSettings.vecPulseDuration = 2e-2; %[2e-3 5e-3 1e-2 2e-2 5e-2]; %s
    sStimParamsSettings.dblITI = 1; %s, inter-pulse interval
    sStimParamsSettings.intRepeats = 100; %reps per pulse duration

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
fprintf('Saving output in directory %s; loading textures from %s\n',strLogDir,strTexDir);

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

%% create presentation vectors/build structEP
structEP = struct;
structEP.strExpType = 'OptoStim';
structEP.intNumRepeats = sStimParams.intRepeats;
structEP.vecPulseDuration = sort(sStimParams.vecPulseDuration);
structEP.intStimTypes = length(structEP.vecPulseDuration);
structEP.dblITI = sStimParams.dblITI;
structEP.dblSecsBlankAtStart = 3;
structEP.dblSecsBlankAtEnd = 3;

%randomize presentation order
structEP.vecTrialStimTypes = [];
for intRep = 1:structEP.intNumRepeats
    structEP.vecTrialStimTypes = [structEP.vecTrialStimTypes randperm(structEP.intStimTypes)];
end
structEP.intTrialNum = numel(structEP.vecTrialStimTypes);

%create pulse data
dblSamplingRate = 1000;
cellPulseVolt = cell(1,structEP.intStimTypes);
for intType = 1:structEP.intStimTypes
	intLength = structEP.vecPulseDuration(intType)*dblSamplingRate;
	if intLength < 100 %min 100 scans to output buffer
		intExtraSamp = 100-intLength;
		cellPulseVolt{intType} = dblLaserTriggerMultiplier*cat(1,linspace(3, 3, intLength)',linspace(0, 0, intExtraSamp)');
	else
		cellPulseVolt{intType} = dblLaserTriggerMultiplier*cat(1,linspace(3, 3, intLength)',linspace(0, 0, 1)');
	end
end
structEP.PulseVolt = cellPulseVolt;

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
            outputData1 = dblLaserTriggerMultiplier*cat(1,linspace(3, 3, 200)',linspace(0, 0, 50)');
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
        outputData1 = dblLaserTriggerMultiplier*cat(1,linspace(3, 3, 200)',linspace(0, 0, 50)');
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
    if boolDebug == 1, vecInitRect = [0 0 640 640];else vecInitRect = [];end
    try
        Screen('Preference', 'SkipSyncTests', 0);
        [ptrWindow,vecRect] = Screen('OpenWindow', intUseScreen,sStimParams.intBackground,vecInitRect);
    catch ME
        warning([mfilename ':ErrorPTB'],'Psychtoolbox error, attempting with sync test skip [msg: %s]',ME.message);
        Screen('Preference', 'SkipSyncTests', 1);
        [ptrWindow,vecRect] = Screen('OpenWindow', intUseScreen,sStimParams.intBackground,vecInitRect);
    end

    %calibrate monitor
    if boolGammaCorrectScreen && exist('C:\Users\neuropixels\Desktop\GammaCorrection\gammaTable20230916.mat','file')
        load("C:\Users\neuropixels\Desktop\GammaCorrection\gammaTable_20230926.mat");
        Screen('LoadNormalizedGammaTable', ptrWindow, gammaTable*[1 1 1]);
    end

    %window variables
    sStimParams.ptrWindow = ptrWindow;
    sStimParams.vecRect = vecRect;
    sStimParams.intScreenWidth_pix = vecRect(3)-vecRect(1);
    sStimParams.intScreenHeight_pix = vecRect(4)-vecRect(2);

    %% MAXIMIZE PRIORITY
    intOldPriority = 0;
    if boolDebug == 0
        intPriorityLevel=MaxPriority(ptrWindow);
        intOldPriority = Priority(intPriorityLevel);
    end

    %% get refresh rate, most likely obsolete but will leave it in just in case
    dblStimFrameRate=Screen('FrameRate', ptrWindow);
    intStimFrameRate = round(dblStimFrameRate);
    dblStimFrameDur = mean(1/dblStimFrameRate);
    dblInterFlipInterval = Screen('GetFlipInterval', ptrWindow);
    if dblStimFrameDur/dblInterFlipInterval > 1.05 || dblStimFrameDur/dblInterFlipInterval < 0.95
        warning([mfilename ':InconsistentFlipDur'],sprintf('Something iffy with flip speed and monitor refresh rate detected; frame duration is %fs, while flip interval is %fs!',dblStimFrameDur,dblInterFlipInterval));
    end

    %% check escape
    if CheckEsc(),error([mfilename ':EscapePressed'],'Esc pressed; exiting');end

    %% PRESENT STIMULI
    %stim-based logs
    structEP.intStimNumber = structEP.intTrialNum;
    structEP.TrialNumber = nan(1,structEP.intStimNumber);
    structEP.ActStimType = nan(1,structEP.intStimNumber);
    structEP.ActOnSecs = nan(1,structEP.intStimNumber);
    structEP.ActOffSecs = nan(1,structEP.intStimNumber);
    structEP.ActOnNI = nan(1,structEP.intStimNumber);
    structEP.ActOffNI = nan(1,structEP.intStimNumber);
    structEP.dblStimFrameDur = dblStimFrameDur;

    %show trial summary
    fprintf('Finished preparation at [%s], will present %d repetitions of %d stimuli\n\n   Waiting for start signal\n',...
        getTime,structEP.intNumRepeats,structEP.intStimTypes);

    %wait for signal
    opts = struct;
    opts.Default = 'Start';
    opts.Interpreter = 'tex';
    strAns = questdlg('Would you like to start the stimulation?', ...r
        'Start Stimulation', ...
        'Start','Cancel',opts);
    if ~strcmp(strAns,opts.Default)
        error([mfilename ':RunCancelled'],'Cancelling');
    end

    %set timers
    refTime = tic;
    dblLastFlip = Screen('Flip', ptrWindow);
    dblInitialFlip = dblLastFlip;

    %timestamp start
    structEP.strStartDate = getDate();
    structEP.strStartTime = getTime();

    %% wait initial-blanking
    fprintf('Starting initial blank (dur=%.3fs) [%s]\n',structEP.dblSecsBlankAtStart,getTime);
    dblInitialBlankDur = 0;
    while dblInitialBlankDur < (structEP.dblSecsBlankAtStart - dblStimFrameDur)
        %do nothing
        Screen('FillRect',ptrWindow, sStimParams.intBackground);
        dblLastFlipFlip = Screen('Flip', ptrWindow,dblLastFlip + dblStimFrameDur/2);
        dblInitialBlankDur = dblLastFlipFlip - dblInitialFlip;
    end

    for intThisTrial = 1:structEP.intTrialNum
        hTicTrial = tic;
        %% check escape
        if CheckEsc(),error([mfilename ':EscapePressed'],'Esc pressed; exiting');end

        %% prep trial
        %trial start
        Screen('FillRect',ptrWindow, sStimParams.intBackground);
        dblTrialStartFlip = Screen('Flip', ptrWindow);

        %retrieve stimulus info
        intStimType = structEP.vecTrialStimTypes(intThisTrial);
        vecPulseOut1 = structEP.PulseVolt{intStimType}(:);
        vecPulseOut2 = dblPupilLightMultiplier*linspace(3, 3, length(vecPulseOut1));

        %fill DAQ with data
        if boolUseNI
            stop(objDaqOut);
            outputData1 = vecPulseOut1(:); %dblLaserTriggerMultiplier*cat(1,linspace(3, 3, 200)',linspace(0, 0, 50)');
            outputData2 = vecPulseOut2(:); %dblPupilLightMultiplier*linspace(3, 3, 250)';
            queueOutputData(objDaqOut,[outputData1 outputData2]);
            prepare(objDaqOut);
        end

        %wait
        while toc(hTicTrial) < (structEP.dblITI*0.9)
            pause((structEP.dblITI - toc(hTicTrial))*0.3);
        end
        while toc(hTicTrial) < structEP.dblITI
            %do nothing
        end

        %start stimulus
        fprintf('Stim [%d/%d] started at %.3fs; %.1fs ms pulse duration\n',intThisTrial,structEP.intTrialNum, toc(refTime),structEP.vecPulseDuration(intStimType)*1000);
        if boolUseNI,startBackground(objDaqOut);end
        if boolUseSGL %log timestamp
            dblStimOnNI = GetScanCount(hSGL, intStreamNI)/dblSampFreqNI;
            dblStimOnSecs = toc(refTime);
        else
            dblStimOnNI = nan;
            dblStimOnSecs = toc(refTime);
        end

        %wait
        if boolUseNI && sStimParamsSettings.intUseDaqDevice > 0
            dblTimeout = (length(vecPulseOut1)/dblSamplingRate)*1.5;
            wait(objDaqOut,dblTimeout);
        end

        %log NI timestamp
        if boolUseSGL
            dblStimOffNI = GetScanCount(hSGL, intStreamNI)/dblSampFreqNI;
            dblStimOffSecs = toc(refTime);
        else
            dblStimOffNI = nan;
            dblStimOffSecs = toc(refTime);
        end

        %new stim-based output
        intStimNumber = intThisTrial;
        structEP.TrialNumber(intStimNumber) = intThisTrial;
        structEP.ActStimType(intStimNumber) = intStimType;
        structEP.PulseDurTrial(intStimNumber) = structEP.vecPulseDuration(intStimType);
        structEP.ActOnSecs(intStimNumber) = dblStimOnSecs;
        structEP.ActOffSecs(intStimNumber) = dblStimOffSecs;
        structEP.ActOnNI(intStimNumber) = dblStimOnNI;
        structEP.ActOffNI(intStimNumber) = dblStimOffNI;


    end

    %% save data
    %save data
    structEP.sStimParams = sStimParams;
    save(fullfile(strLogDir,strFilename), 'structEP','sParamsSGL');

    %show trial summary
    fprintf('Finished experiment & data saving at [%s], waiting for end blank (dur=%.3fs)\n',getTime,structEP.dblSecsBlankAtEnd);

    %% wait end-blanking
    dblEndBlankDur = 0;
    while dblEndBlankDur < structEP.dblSecsBlankAtEnd
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
    if strcmp(ME.identifier,'RunDriftingGratings:EscapePressed')
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
%end
