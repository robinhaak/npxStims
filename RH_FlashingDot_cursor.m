
%RH_FLASHINGDOT_CURSOR
%Robin Haak, last update 7 November 2022

clear; close all; Screen('CloseAll');
fprintf('Starting %s [%s]\n',mfilename,getTime);

%% set debug switch
boolDebug  = false;

%% input parameters
sStimParams = struct;

%screen variables
sStimParams.intUseScreen = 2; %which screen to use
sStimParams.dblBackground = 0.5; %background color (dbl, [0 1]); 0 = black, 1 = white
sStimParams.intBackground = round(mean(sStimParams.dblBackground)*255);

%stimulus parameters
sStimParams.dblFlashRate = 3; %Hz; initial flash rate
sStimParams.dblSize = 100; %pix; initial stimulus size
sStimParams.dblStimulus = 0; %stimulus color (dbl, [0 1])
sStimParams.intStimulus = round(mean(sStimParams.dblStimulus)*255);
sStimParams.dblStimulusAlternate = 1; %stimulus color (dbl, [0 1])
sStimParams.intStimulusAlternate = round(mean(sStimParams.dblStimulusAlternate)*255);

if boolDebug == true
    sStimParams.intUseScreen = 0;
end

%% start PTB
try
    fprintf('Starting PsychToolBox extension...\n');
    %% open window
    AssertOpenGL;
    KbName('UnifyKeyNames');
    intOldVerbosity = Screen('Preference', 'Verbosity',1); %stop PTB spamming
    if boolDebug == 1, vecInitRect = [0 0 640 640]; else; vecInitRect = []; end
    try
        Screen('Preference', 'SkipSyncTests', 0);
        [ptrWindow,vecRect] = Screen('OpenWindow', sStimParams.intUseScreen,sStimParams.intBackground,vecInitRect);
    catch ME
        warning([mfilename ':ErrorPTB'],'Psychtoolbox error, attempting with sync test skip [msg: %s]',ME.message);
        Screen('Preference', 'SkipSyncTests', 1);
        [ptrWindow,vecRect] = Screen('OpenWindow', sStimParams.intUseScreen,sStimParams.intBackground,vecInitRect);
    end

    %calibrate monitor
    if exist('C:\Users\neuropixels\Desktop\GammaCorrection\gammaTable.mat','file')
        load("C:\Users\neuropixels\Desktop\GammaCorrection\gammaTable.mat");
        Screen('LoadNormalizedGammaTable', ptrWindow, gammaTable*[1 1 1]);
    end

    %window variables
    sStimParams.intScreenWidth_pix = vecRect(3)-vecRect(1);
    sStimParams.intScreenHeight_pix = vecRect(4)-vecRect(2);

    %% MAXIMIZE PRIORITY
    intOldPriority = 0; %#ok<*NASGU>
    if boolDebug == 0
        intPriorityLevel=MaxPriority(ptrWindow);
        intOldPriority = Priority(intPriorityLevel);
    end

    %% get monitor refresh rate
    dblStimFrameRate = Screen('FrameRate', ptrWindow);
    intStimFrameRate = round(dblStimFrameRate);
    dblStimFrameDur = mean(1/dblStimFrameRate);
    dblInterFlipInterval = Screen('GetFlipInterval', ptrWindow);
    if dblStimFrameDur/dblInterFlipInterval > 1.05 || dblStimFrameDur/dblInterFlipInterval < 0.95
        warning([mfilename ':InconsistentFlipDur'],sprintf('Something iffy with flip speed and monitor refresh rate detected; frame duration is %fs, while flip interval is %fs!',dblStimFrameDur,dblInterFlipInterval)); %#ok<SPWRN>
    end

    %% start stimulation
    %set initial position and size
    SetMouse(sStimParams.intScreenWidth_pix/2,sStimParams.intScreenHeight_pix/2,ptrWindow);%set the cursor to the top left of the screen to start with
    dblSize = sStimParams.dblSize; %pix
    intFlashRate = round(dblStimFrameRate/(sStimParams.dblFlashRate)); %frames

    %start stimulation
    HideCursor(ptrWindow);
    dblLastFlip = Screen('Flip',ptrWindow);
    intFrame = 0;
    hTic = tic;
    while ~CheckEsc()

        intFrame = intFrame + 1;

        %get current cursor position and button presses
        [dblPosX, dblPosY, vecButtons] = GetMouse(ptrWindow);
        dblPosX = min(dblPosX, sStimParams.intScreenWidth_pix);
        dblPosY = min(dblPosY, sStimParams.intScreenHeight_pix);

        %change stimulus size based on ui
        if vecButtons(1) == 1
            dblSize = dblSize - 1;
            if dblSize < 0, dblSize = 0; end %min
        elseif vecButtons(3) == 1
            dblSize = dblSize + 1;
            if dblSize > sStimParams.intScreenWidth_pix, dblSize = sStimParams.intScreenWidth_pix; end %max
        end

        %set color
        if intFrame <= intFlashRate
            intStimColor = sStimParams.intStimulus;
        elseif intFrame > intFlashRate
            intStimColor = sStimParams.intStimulusAlternate;
            if intFrame == 2*intFlashRate
                intFrame = 0;
            end
        end

        %draw and flip
        vecBoundingRect = [dblPosX-dblSize/2,dblPosY-dblSize/2,dblPosX+dblSize/2,dblPosY+dblSize/2];
        Screen('FillOval', ptrWindow,intStimColor,vecBoundingRect);
        dblLastFlip = Screen('Flip', ptrWindow, dblLastFlip+0.5/intStimFrameRate);
        fprintf('Position: %d (x), %d (y); Size: %d\n',dblPosX,dblPosY,dblSize);
    end

    %close
    Screen('Close',ptrWindow);
    Screen('Close');
    Screen('CloseAll');
    ShowCursor;
    Priority(0);
    Screen('Preference','Verbosity',intOldVerbosity);
catch

    %close
    Screen('Close',ptrWindow);
    Screen('Close');
    Screen('CloseAll');
    ShowCursor;
    Priority(0);
    Screen('Preference','Verbosity',intOldVerbosity);
end