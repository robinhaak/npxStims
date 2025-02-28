%RH_runAsyncStim()
%December 2023, Robin Haak

%% switches
boolDebug = false;
boolUseSGL = true; %always true!

%% start memory maps
try
	%% data switch
	mmapSignal = JoinMemMap('dataswitch');
	
	%% load stim params
	mmapParams = JoinMemMap('sStimParams','struct');
	sStimParams = mmapParams.Data;
	if ~isstruct(sStimParams) && isscalar(sStimParams) && sStimParams == 0
		error([mfilename ':DataMapNotInitialized'],'Data transfer failed. Did you start the other matlab first?');
	end
	strHostAddress = sStimParams.strHostAddress;
catch
	error([mfilename ':DataMapNotInitialized'],'Memory mapping and/or data transfer failed. Did you start the other matlab first?');
end

%% connect to spikeglx
if boolDebug == 1
	hSGL = [];
else
	hSGL = SpikeGL(strHostAddress);
end

intStreamNI = -1;
dblSampFreqNI = GetSampleRate(hSGL, intStreamNI);

%% start PTB
try
	fprintf('Starting PsychToolBox extension...\n');
	%% open window
	AssertOpenGL;
	KbName('UnifyKeyNames');
	intOldVerbosity = Screen('Preference', 'Verbosity',1); %stop PTB spamming
	if boolDebug == 1, vecInitRect = [0 0 640 640];else; vecInitRect = [];end
	try
		Screen('Preference', 'SkipSyncTests', 0);
		[ptrWindow,vecRect] = Screen('OpenWindow', sStimParams.intUseScreen,sStimParams.intBackground,vecInitRect);
	catch ME
		warning([mfilename ':ErrorPTB'],'Psychtoolbox error, attempting with sync test skip [msg: %s]',ME.message);
		Screen('Preference', 'SkipSyncTests', 1);
		[ptrWindow,vecRect] = Screen('OpenWindow', sStimParams.intUseScreen,sStimParams.intBackground,vecInitRect);
	end
	
	%correct gamma
	if exist('C:\Users\neuropixels\Desktop\GammaCorrection\gammaTable_20230926.mat','file')
		load('C:\Users\neuropixels\Desktop\GammaCorrection\gammaTable_20230926.mat');
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
	
	sStimParams.dblPixelsPerDeg = sStimParams.intScreenHeight_pix/sStimParams.dblScreenHeight_deg;
	sStimParams.intWhite = WhiteIndex(sStimParams.ptrWindow);
	
	%% MAXIMIZE PRIORITY
	intOldPriority = 0; %#ok<*NASGU>
	if boolDebug == 0
		intPriorityLevel=MaxPriority(ptrWindow);
		intOldPriority = Priority(intPriorityLevel);
	end
	
	%% get monitor refresh rate
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
	
	%% add parameters to stim struct
	intStimSet = sStimParams.intStimSet;
	[sDots,sStimParams] = RH_CreateDotsFromAbove(intStimSet,sStimParams);
	
	%build structEP
	structEP = struct;
	structEP.sStimParams = sStimParams;
	structEP.strExpType = sStimParams.strStimType;
	structEP.strStimSet = intStimSet;
	structEP.dblStimFrameDur = dblStimFrameDur;
	structEP.dblInterFlipInterval = dblInterFlipInterval;
	if intStimSet==1 %one location
		vecStimID = repmat(sDots.stimID(1),[1 sStimParams.intTrials]);
	elseif intStimSet==2
		intTrialsNonAct = sStimParams.intTrials-sStimParams.intTrialsActive;
		vecThisStimID = repmat(sDots.stimID,[1 sStimParams.intTrialsActive/2]);
		vecThisStimID = vecThisStimID(randperm(numel(vecThisStimID)));
		vecThisStimID_NonAct = repmat(sDots.stimID,[1 intTrialsNonAct/2]);
		vecThisStimID_NonAct = vecThisStimID_NonAct(randperm(numel(vecThisStimID_NonAct)));
		vecStimID = [vecThisStimID vecThisStimID_NonAct];
	end
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
% 	structEP.ActOnSecs = nan(1,structEP.intTrialNum);
% 	structEP.ActOffSecs = nan(1,structEP.intTrialNum);
% 	structEP.ActStartSecs = nan(1,structEP.intTrialNum);
% 	structEP.ActStopSecs = nan(1,structEP.intTrialNum);
	structEP.ActOnNI = nan(1,structEP.intTrialNum);
	structEP.ActOffNI = nan(1,structEP.intTrialNum);
	
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
	
	%% Start stimulation
	% trial start; background
	Screen('FillRect',ptrWindow, sStimParams.intBackground);
	dblLastFlip = Screen('Flip', ptrWindow);
	
	%% run until we get the signal to stop
	intStimNumber = mmapSignal.Data(1);
	%send signal we're ready to start
	fprintf('Preparation complete. Sending go-ahead signal!\n');
	mmapSignal.Data(2) = -1;
	intTrialCounter = 0;
	while intStimNumber ~= -1
		%check if we need to show a new stimulus
		intStimNumber = mmapSignal.Data(1);
		intStimType = mmapSignal.Data(2);
		if intStimNumber > 0
			%% set counter
			intTrialCounter = intTrialCounter + 1;
			intThisTrial = intTrialCounter;
			
			%get trial-specific parameters
			vecThisBoundingRect = structEP.cellBoundingRect{intThisTrial};
			intNumFrames = length(vecThisBoundingRect);
			vecThisColor = structEP.cellColor{intThisTrial};
			
			%display in command window
			fprintf('[%s] Showing stim %d/%d\n',getTime,intThisTrial,structEP.intTrialNum);
			
			%% prsent stimulus
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
						dblStimOnNI = GetScanCount(hSGL, intStreamNI)/dblSampFreqNI;
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
				dblStimOffNI = GetScanCount(hSGL, intStreamNI)/dblSampFreqNI;
			else
				dblStimOffNI = nan; %#ok<UNRCH>
			end
			
			%get stim duration
			dblStimOffFlip = dblLastFlip;
			dblStimDur = dblStimOffFlip-dblStimOnFlip;
			
			%% save data
			%save object
% 			sObject = sStimParams;
% 			sObject.dblStimOnNI = dblStimOnNI;
% 			sObject.dblStimOffNI = dblStimOffNI;
% 			sObject.intStimType = 1;%  intStimType;
% 			save(fullfile(sStimParams.strTempObjectPath,['Object',num2str(intStimNr),'.mat']),'sObject');
% 			
			structEP.TrialNumber(intThisTrial) = intThisTrial;
			structEP.ActStimType(intThisTrial) = 1; %intStimType
			structEP.ActOnNI(intThisTrial) = dblStimOnNI;
			structEP.ActOffNI(intThisTrial) = dblStimOffNI;
			
			%% reset signal
			mmapSignal.Data=[0 0];
		end
		%pause to avoid cpu overload
		pause(0.01);
	end
	
	%% export all data
	mmapData = InitMemMap('structEP',structEP);
	clear mmapData;
	
	%signal we're done
	mmapSignal.Data(1) = -2;
	mmapSignal.Data(2) = -2;
	
	%% close PTB
	Screen('Close',ptrWindow);
	Screen('CloseAll');
	ShowCursor;
	Priority(0);
	Screen('Preference', 'Verbosity',intOldVerbosity);
	
	%% wait until data has been received
	while mmapSignal.Data(1) ~= -3
		pause(0.1);
	end
	
catch ME
	%% catch me and throw me
	Screen('Close');
	Screen('CloseAll');
	ShowCursor;
	Priority(0);
	Screen('Preference', 'Verbosity',intOldVerbosity);
	%% show error
	rethrow(ME);
end
% end