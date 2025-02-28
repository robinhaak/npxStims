%RH_RunLocomotionFeedback
%December 2023, Robin Haak

%% suppress m-lint warnings
%#ok<*MCCD,*NASGU,*ASGLU,*CTCH>
%clearvars -except sStimPresets sStimParamsSettings sExpMeta;
clear; close all; Screen('CloseAll');

%% set default variables and debug switches
% intStimSet = 1;
boolUseSGL = true;
boolUseNI = true;
boolDebug = false;

%% set channel ids
%niCh0=onset
%niCh1=running
%niCh2=sync
intStimOnsetChanNI=0;
intRunningChanNI=1;
intStreamNI=-1;

%% define variables
%defaults
dblPupilLightMultiplier = 1; %strength of infrared LEDs
dblSyncLightMultiplier = 0.5;
strHostAddress = '192.87.11.84'; %'192.87.10.228'; %'192.87.11.133'; %default host address
objDaqOut = [];

fprintf('Starting %s [%s]\n',mfilename,getTime);
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

%% input parameters
fprintf('Loading settings...\n');
if ~exist('sStimParamsSettings','var') || isempty(sStimParamsSettings)
	%general
	sStimParamsSettings = struct;
	sStimParamsSettings.strStimType = 'LocomotionFeedback';
	sStimParamsSettings.strOutputPath = 'C:\_Data\Exp'; %appends date
	sStimParamsSettings.strTempObjectPath = 'C:\temp'; %'X:\JorritMontijn\';%X:\JorritMontijn\ or F:\Data\Temp\
	
	%visual space parameters
	% 	sStimParamsSettings.dblSubjectPosX_cm = 0; % cm; relative to center of screen
	% 	sStimParamsSettings.dblSubjectPosY_cm = 0; % cm; relative to center of screen
	sStimParamsSettings.dblScreenDistance_cm = 31; %35; %43; % cm; measured, 14
	
	%screen variables
	sStimParamsSettings.intCornerTrigger = 1; % integer switch; 0=none,1=upper left, 2=upper right, 3=lower left, 4=lower right
	sStimParamsSettings.dblCornerSize = 1/50; % fraction of screen width
	sStimParamsSettings.dblScreenWidth_cm = 70.15; % cm; measured [51]
	sStimParamsSettings.dblScreenHeight_cm = 39.6; % cm; measured [29]
	sStimParamsSettings.dblScreenWidth_deg = 2 * atand(sStimParamsSettings.dblScreenWidth_cm / (2 * sStimParamsSettings.dblScreenDistance_cm));
	sStimParamsSettings.dblScreenHeight_deg = 2 * atand(sStimParamsSettings.dblScreenHeight_cm / (2 * sStimParamsSettings.dblScreenDistance_cm));
	sStimParamsSettings.intUseScreen = 2; %which screen to use
	sStimParamsSettings.dblBackground = 0.5; %background intensity (dbl, [0 1])
	sStimParamsSettings.intBackground = round(mean(sStimParamsSettings.dblBackground)*255);
	
	%stimulus parameters
	sStimParamsSettings.dblSize_deg = 4; %deg
	sStimParamsSettings.dblSpeed_deg = 20; %deg/s
	sStimParamsSettings.intTrials = 40;
	sStimParamsSettings.intTrialsActive = sStimParamsSettings.intTrials/2; %first n trials coupled to locomotion
	sStimParamsSettings.dblMaxTimeActive_min = 35; %min
	sStimParamsSettings.dblTrialInterval = 60; %s
	sStimParamsSettings.dblTrialIntervalVar = 15; %s
	sStimParamsSettings.dblRunThreshold = 0.1; %m/s
	sStimParamsSettings.dblRunPeriod_s = 2; %s
	
	%stimulus control variables
	sStimParamsSettings.intUseDaqDevice = 1; %ID of DAQ device
	sStimParamsSettings.intUseParPool = 2; %number of workers in parallel pool; [2]
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

%RunExperiment will not load values from below here
%keeps stuff the same as other scripts
sStimParams = sStimParamsSettings;
if boolDebug == 1
	sStimParams.intUseScreen = 0;
end

%% set output location for the logs
strOutputPath = sStimParamsSettings.strOutputPath;
strTempObjectPath = sStimParamsSettings.strTempObjectPath;
strThisFilePath = mfilename('fullpath');
[strFilename,strLogDir,strTempDir,strTexDir] = RE_assertPaths(strOutputPath,strRecording,strTempObjectPath,strThisFilePath);
fprintf('Saving output in directory %s;\n',strLogDir); %no textures are loaded for this script

%% query user for the stimulus set & approx. location of the receptive field
fprintf(['\n--<strong>Select stimulus set--</strong>\nAvailable sets:\n'...
	'(1)"one dot"\n(2)"two dots"\n']);
intStimSet = input('set: ');
if sum(intStimSet==[1 2])<1,error('Non-existent stim set');end
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
	intStreamNI = -1;
	dblSampFreqNI = GetSampleRate(hSGL, intStreamNI);
	
	%% check disk space available
	strDataDirSGL = GetDataDir(hSGL);
	jFileObj = java.io.File(strDataDirSGL);
	dblFreeGB = (jFileObj.getFreeSpace)/(1024^3);
	if dblFreeGB < 100,warning([mfilename ':LowDiskSpace'],'Low disk space available (%.0fGB) for Neuropixels data (dir: %s)',dblFreeGB,strDataDirSGL);end
	
	%% prepare real-time stream
	%retrieve some parameters
	dblSampFreqNI = GetSampleRate(hSGL, intStreamNI);
	dblNI2V = (sParamsSGL.niAiRangeMax)/(double(intmax('int16'))/2);
	
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

%% initialize parallel pool && gpu
if sStimParams.intUseGPU > 0
	objGPU = gpuDevice(sStimParams.intUseGPU);
end

%% initialize NI I/O box
% if boolUseNI
% 	%initialize
% 	fprintf('Connecting to National Instruments box...\n');
% 	boolDaqOutRunning = false;
% 	if exist('objDaqOut','var') && ~isempty(objDaqOut)
% 		try
% 			%turns leds on
% 			stop(objDaqOut);
% 			outputData1 = dblSyncLightMultiplier*cat(1,linspace(3, 3, 200)',linspace(0, 0, 50)');
% 			outputData2 = dblPupilLightMultiplier*linspace(3, 3, 250)';
% 			queueOutputData(objDaqOut,[outputData1 outputData2]);
% 			prepare(objDaqOut);
% 			pause(0.1);
% 			startBackground(objDaqOut)
% 			boolDaqOutRunning = true;
% 		catch
% 		end
% 	end
% 	if ~boolDaqOutRunning
% 		objDaqOut = openDaqOutput(sStimParams.intUseDaqDevice);
% 		%turns leds on
% 		stop(objDaqOut);
% 		outputData1 = dblSyncLightMultiplier*cat(1,linspace(3, 3, 200)',linspace(0, 0, 50)');
% 		outputData2 = dblPupilLightMultiplier*linspace(3, 3, 250)';
% 		queueOutputData(objDaqOut,[outputData1 outputData2]);
% 		prepare(objDaqOut);
% 		pause(0.1);
% 		startBackground(objDaqOut)
% 	end
% end
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

%% start experiment
try
	%% pre-allocate & set parameters
	%parameters
	intNumTrials = sStimParams.intTrials;
	intTrialsActive = sStimParams.intTrialsActive;
	dblMaxTimeActive_s = sStimParams.dblMaxTimeActive_min*60;
	dblMinSampleTime = 2e-3;
	dblRunThreshold = sStimParams.dblRunThreshold;
	dblSampFreq = 300; %Hz, measured(?)
	vecRunThresholdCrossed = zeros(round(dblSampFreq)*round(sStimParams.dblRunPeriod_s),1);
	intRunThresholdPos = 0;
	vecTrialInterval = randi([-sStimParams.dblTrialIntervalVar sStimParams.dblTrialIntervalVar],[intNumTrials+1 1]) + ...
		sStimParams.dblTrialInterval;
	dblApproxStimDur = (sStimParams.dblScreenHeight_deg+sStimParams.dblSize_deg)/sStimParams.dblSpeed_deg;
	
	vecIsActive = zeros(sStimParams.intTrials,1);
	
	%fixed stuff
	sStimParams.strHostAddress=strHostAddress;
	mmapSignal = InitMemMap('dataswitch',[0 0]);
	mmapParams = InitMemMap('sStimParams',sStimParams);
	clear mmapParams;
	intStimNumber = 0;
	dblDaqRefillDur = 0.5; %must be less than the stimulus duration, and cannot be less than ~300ms
	boolMustRefillDaq = false;
	if boolUseNI
		%refill daq
		stop(objDaqOut);
		outputData1 = dblSyncLightMultiplier*cat(1,linspace(3, 3, 200)',linspace(0, 0, 50)');
		outputData2 = dblPupilLightMultiplier*linspace(3, 3, 250)';
		queueOutputData(objDaqOut,[outputData1 outputData2]);
		prepare(objDaqOut);
	end
	
    pause(2);
	
	%% wait for other matlab to join memory map
	fprintf('Preparation complete. Waiting for PTB matlab to join the memory map...\n');
	while mmapSignal.Data(2) == 0
		pause(0.1);
	end
	
	%% run until escape button is pressed
	dblLastStim = 0; %-inf;
	% 	intSample = 0; %for pilot (to be deleted later)
	% 	vecSampleTimes = nan(1e5,1);
	hTic = tic;
	dblTrialInterval = vecTrialInterval(1);
	
	while intStimNumber <= intNumTrials && ~CheckEsc()
		% 		intSample = intSample + 1;
		% 		vecSampleTimes(intSample) = toc(hTic);
		
		%check whether script execution should be paused
		if CheckPause()
			fprintf('\n\n<strong>STIMULUS SCRIPT PAUSED!</strong> [%s]\n\n',getTime);
			pauseTic = tic;
			WaitSecs(1);
			%get user input
			opts = struct;
			opts.Default = 'Restart';
			opts.Interpreter = 'tex';
			strAns = questdlg('STIMULUS SCRIPT PAUSED! Would you like to restart the stimulation?', ...
				'Restart Stimulation', ...
				'Restart',opts);
			WaitSecs(1);
			dblLastStim = dblLastStim + toc(pauseTic);
		end
		
		%break if all stimuli are shown
		if intStimNumber == intNumTrials && toc(hTic) > (dblLastStim + dblTrialInterval) %assuming that ITI > stim duration
			break
		end
		
		if toc(hTic) < dblMaxTimeActive_s && intStimNumber < intTrialsActive
			%read running speed
			if boolUseSGL
				[vecSamples,vecRunningData,vecSyncData,sStream]=readRunningSpeed(sStream,intRunningChanNI,intStimOnsetChanNI);
				sRunSpeed = PP_GetRunSpeed(vecRunningData,dblSampFreqNI);
				vecTimestamps_s = sRunSpeed.vecOutT + sStream.dblEphysTimeNI;
				vecTraversed_m = sRunSpeed.vecTraversed_m;
				vecRunningSpeed_mps = sRunSpeed.vecSpeed_mps; %this is filtered, so it might not be accurate for short data vectors
			else
				vecTimestamps_s = toc(hTic);
				vecRunningSpeed_mps = 2;
			end
			
			% 			%keep running speed for online analysis (rough), pilot (to be deleted later)
			% 			vecApproxRunTimes_s(intSample) = vecTimestamps_s(end); %#ok<*SAGROW>
			% 			vecApproxRunSpeed_mps(intSample) = vecRunningSpeed_mps(end);
			
			%do some sort of calculation here to determine whether the mouse crossed the threshold
			dblRunSpeed = mean(vecRunningSpeed_mps);
			if intRunThresholdPos < length(vecRunThresholdCrossed)
				intRunThresholdPos = intRunThresholdPos + 1;
			else
				intRunThresholdPos = 1;
			end
			if dblRunSpeed > dblRunThreshold
				vecRunThresholdCrossed(intRunThresholdPos) = 1;
			else
				vecRunThresholdCrossed(intRunThresholdPos) = 0;
			end
			boolRunningThresholdCrossed = sum(vecRunThresholdCrossed) > 0.95 * length(vecRunThresholdCrossed);
			boolIsActive = 1;
		else
			boolRunningThresholdCrossed = 1;
			boolIsActive = 0;
		end
		
		%check if we should refill the DAQ
		if boolUseNI && boolMustRefillDaq && (toc(hTic) > (dblDaqRefillDur + dblLastStim + dblTrialInterval))
			stop(objDaqOut);
			outputData1 = dblSyncLightMultiplier*cat(1,linspace(3, 3, 200)',linspace(0, 0, 50)');
			outputData2 = dblPupilLightMultiplier*linspace(3, 3, 250)';
			queueOutputData(objDaqOut,[outputData1 outputData2]);
			prepare(objDaqOut);
			boolMustRefillDaq = false;
		end
		
		%check if running speed is high enough
		if intStimNumber < intNumTrials && boolRunningThresholdCrossed && (toc(hTic) > (dblLastStim + dblApproxStimDur + dblTrialInterval))
			%% increment trial & log NI timestamp
			intStimNumber = intStimNumber + 1;
			dblLastStim = toc(hTic);
			vecIsActive(intStimNumber) = boolIsActive;
			fprintf('Stim %d started at %s (run speed was %.3f)\n',intStimNumber,getTime,dblRunSpeed);
			
			%% flash eye-tracker synchronization LED
			if boolUseNI
				startBackground(objDaqOut);
				boolMustRefillDaq = true;
			end
			
			%% start stimulus
			mmapSignal.Data(1) = intStimNumber;
			% 			mmapSignal.Data(2) = intStimType;
			
			%% get approximate timestamp for start
			if ~isempty(hSGL)
				dblApproxStimOnNI = GetScanCount(hSGL, intStreamNI)/dblSampFreqNI;
			else
				dblApproxStimOnNI = toc(hTic);
			end
			
			dblTrialInterval = vecTrialInterval(intStimNumber+1);
		end
		
	end
	%signal end
	mmapSignal.Data(1) = -1;
	% 	mmapSignal.Data(2) = -1;
	
	%% wait for other matlab to join memory map
	fprintf('Experiment complete. Waiting for PTB matlab to send data...\n');
	while ~all(mmapSignal.Data == -2)
		pause(0.1);
	end
	fprintf('Data received! Sending exit signal to PTB matlab and & saving data.\n');
	
	%% retrieve trial data
	mmapData = JoinMemMap('structEP','struct');
	structEP = mmapData.Data;
	
	%signal retrieval
	mmapSignal.Data(1) = -3;
	mmapSignal.Data(2) = -3;
	
	%% save data
	%save data
	structEP.vecIsActive = vecIsActive;
	structEP.sStimParams = sStimParams;
	% 	structEP.vecApproxRunTimes_s = vecApproxRunTimes_s;
	% 	structEP.vecApproxRunSpeed_mps = vecApproxRunSpeed_mps;
	save(fullfile(strLogDir,strFilename), 'structEP','sParamsSGL');
	
	%clean up
	fprintf('\nExperiment is finished at [%s], closing down and cleaning up...\n',getTime);
	
	%% close Daq IO
	if boolUseNI && ~(exist('sExpMeta','var') && isfield(sExpMeta,'objDaqOut'))
		try
			closeDaqOutput(objDaqOut);
		catch
		end
	end
	
catch ME
	%% catch me and throw me
	fprintf('\n\n\nError occurred! Trying to save data and clean up...\n\n\n');
	
	%save data
	structEP.vecIsActive = vecIsActive;
	structEP.sStimParams = sStimParams;
	% 	structEP.vecApproxRunTimes_s = vecApproxRunTimes_s;
	% 	structEP.vecApproxRunSpeed_mps = vecApproxRunSpeed_mps;
	if ~exist('sParamsSGL','var'),sParamsSGL=[];end
	save(fullfile(strLogDir,strFilename), 'structEP','sParamsSGL');
	
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
