% RH_FLASHINGDOT_CURSOR
% Robin Haak, 2022-2025
%
% Displays a flashing dot stimulus at the mouse cursor position.The dot
% alternates color at a set rate, mouse buttons control size.

clear; close all; Screen('CloseAll');
boolDebug  = false; %set debug switch

%% INPUT MENU for parameters
prompt = {
    'Screen number to use (e.g. 0, 1):', ...
    'Background color (0=black to 1=white):', ...
    'Corner trigger (1=UL, 2=UR, 3=LL, 4=LR):', ...
    'Dot flash rate (Hz):', ...
    'Initial dot size (pixels):', ...
    ['Color mode:' newline ...
    '1 = Alternating black/white' newline ...
    '2 = Flashing black' newline ...
    '3 = Flashing white'], ...
    'Dot expansion speed:'
};
dlgtitle = 'Set Flashing Dot Parameters';
dims = [1 60];
definput = {'2', '0.5', '2', '2', '100', '1', '2'};

answer = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(answer)
    error('User cancelled parameter selection.');
end

%parse inputs
sStimParams = struct;
sStimParams.intUseScreen = str2double(answer{1});
sStimParams.dblBackground = str2double(answer{2});
sStimParams.intCornerTrigger = str2double(answer{3});
sStimParams.dblFlashRate = str2double(answer{4});
sStimParams.dblSize = str2double(answer{5});
colorMode = str2double(answer{6});
sStimParams.dblExpansionSpeed = str2double(answer{7});

%validate inputs
if isnan(sStimParams.intUseScreen)
    error('Invalid screen number.');
end
if isnan(sStimParams.dblBackground) || sStimParams.dblBackground < 0 || sStimParams.dblBackground > 1
    error('Background color must be between 0 and 1.');
end
if isnan(sStimParams.intCornerTrigger) || ~ismember(sStimParams.intCornerTrigger, 1:4)
    error('Corner trigger must be 1, 2, 3, or 4.');
end
if isnan(sStimParams.dblFlashRate) || sStimParams.dblFlashRate < 0
    error('Flash rate must be >= 0.');
end
if isnan(sStimParams.dblSize) || sStimParams.dblSize <= 0
    error('Initial dot size must be positive.');
end
if isnan(colorMode) || ~ismember(colorMode, [1 2 3])
    error('Color mode must be 1, 2, or 3.');
end

sStimParams.intBackground = round(sStimParams.dblBackground * 255);

switch colorMode
    case 1  % alternating black and white
        sStimParams.dblStimulus = 0;          % black
        sStimParams.dblStimulusAlternate = 1; % white
    case 2  % flashing black
        sStimParams.dblStimulus = 0;                         % black
        sStimParams.dblStimulusAlternate = sStimParams.dblBackground; % background
    case 3  % flashing white
        sStimParams.dblStimulus = 1;                         % white
        sStimParams.dblStimulusAlternate = sStimParams.dblBackground; % background
end

if isnan(sStimParams.dblExpansionSpeed) || sStimParams.dblExpansionSpeed <= 0
    error('Expansion speed must be a positive number.');
end

sStimParams.intStimulus = round(sStimParams.dblStimulus * 255);
sStimParams.intStimulusAlternate = round(sStimParams.dblStimulusAlternate * 255);

sStimParams.dblCornerSize = 1/30; % fraction of screen width


% %screen variables
% sStimParams.intUseScreen = 1; %which screen to use
% sStimParams.dblBackground = 0.5; %background color (dbl, [0 1]); 0 = black, 1 = white
% sStimParams.intBackground = round(mean(sStimParams.dblBackground)*255);
% sStimParams.intCornerTrigger = 1; % 1; %2;
% sStimParams.dblCornerSize = 1/30; % fraction of screen width
%
% %stimulus parameters
% sStimParams.dblFlashRate = 2; %Hz; initial flash rate
% sStimParams.dblSize = 100; %pix; initial stimulus size
% sStimParams.dblStimulus = 0; %stimulus color (dbl, [0 1])
% sStimParams.intStimulus = round(mean(sStimParams.dblStimulus)*255);
% sStimParams.dblStimulusAlternate = 1; %0.5; %stimulus color (dbl, [0 1])
% sStimParams.intStimulusAlternate = round(mean(sStimParams.dblStimulusAlternate)*255);

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
    sStimParams.intWhite = WhiteIndex(ptrWindow);

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

    %% start stimulation
    %set initial position and size
    SetMouse(sStimParams.intScreenWidth_pix/2,sStimParams.intScreenHeight_pix/2,ptrWindow);%set the cursor to the top left of the screen to start with
    dblSize = sStimParams.dblSize; %pix
    intFlashRate = round(dblStimFrameRate/(sStimParams.dblFlashRate)); %frames

    %     HideCursor(ptrWindow);
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
            dblSize = dblSize - sStimParams.dblExpansionSpeed;
            if dblSize < 0, dblSize = 0; end
        elseif vecButtons(3) == 1
            dblSize = dblSize + sStimParams.dblExpansionSpeed;
            if dblSize > sStimParams.intScreenWidth_pix, dblSize = sStimParams.intScreenWidth_pix; end
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
        if intStimColor == sStimParams.intStimulus
            Screen('FillRect',ptrWindow,sStimParams.intWhite,vecDiodeRect); %diode
        end
        dblLastFlip = Screen('Flip', ptrWindow, dblLastFlip+0.5/intStimFrameRate);
        % fprintf('Position: %d (x), %d (y); Size: %d\n',dblPosX,dblPosY,dblSize);
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