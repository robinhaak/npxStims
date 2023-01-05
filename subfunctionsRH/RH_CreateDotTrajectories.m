function  [sAllDots, sStimParams] = RH_CreateDotTrajectories(intStimSet,sStimParams)
%I'm aware that this has become a monstrosity
%Robin Haak, last update: 22 November 2022

%% standard dot parameters
sStimParams.dblSize_deg = 3; % deg; approx. stimulus size in dva
sStimParams.intSize_pix = round(sStimParams.dblSize_deg*sStimParams.dblPixelsPerDeg);
sStimParams.dblColor = 0; %color (dbl, [0 1])
sStimParams.intColor = round(mean(sStimParams.dblColor)*255);
sStimParams.vecDotSpeeds_deg = 15; %30; %approx. speed in deg/s
sStimParams.vecDotSpeeds_pix = sStimParams.vecDotSpeeds_deg*sStimParams.dblPixelsPerDeg; %pixels/s
sStimParams.vecDotSpeeds_ppf = round(sStimParams.vecDotSpeeds_pix/sStimParams.intStimFrameRate);
sStimParams.vecDirections = [0 180]; %0 is rightward & 90 is downward motion, for now only 0 and 180 are available (with the exception of 'dot_grid')
sStimParams.intReps = 15;

if intStimSet == 1
    %% (1) 'dot_grid' - dots moving along multiple horizontal & vertical trajectories
    if length(sStimParams.vecDotSpeeds_pix)>1 %just to be sure
        error('vecDotSpeeds_pix should not be >1!');
    end
	%set-specific parameters
	sStimParams.intReps = 10;
	sStimParams.vecDotSpeeds_deg = 30; %approx. speed in deg/s
	sStimParams.vecDotSpeeds_pix = sStimParams.vecDotSpeeds_deg*sStimParams.dblPixelsPerDeg; %pixels/s
	sStimParams.vecDotSpeeds_ppf = round(sStimParams.vecDotSpeeds_pix/sStimParams.intStimFrameRate);
	sStimParams.vecSecsPostBlank = [0.25 0.25];

    %get number of trajectories
    sStimParams.intTrajSpacing_pix = round(6*sStimParams.dblPixelsPerDeg); %spacing is ~6 dva (as in Beltramo && Scanziani, Science 2019)
    intTrajX = round((sStimParams.intScreenWidth_pix-sStimParams.intTrajSpacing_pix)/sStimParams.intTrajSpacing_pix)*2; %2 directions
    intTrajY = round((sStimParams.intScreenHeight_pix-sStimParams.intTrajSpacing_pix)/sStimParams.intTrajSpacing_pix)*2;

    %get x (for vertical movement) and y (for horizontal movement) coordinates, center grid on screen
    vecCoordsX = (0:intTrajX/2-1)*sStimParams.intTrajSpacing_pix; vecCoordsX = vecCoordsX+(sStimParams.intScreenWidth_pix-max(vecCoordsX))/2;
    vecCoordsY = (0:intTrajY/2-1)*sStimParams.intTrajSpacing_pix; vecCoordsY = vecCoordsY+(sStimParams.intScreenHeight_pix-max(vecCoordsY))/2;

    %create 'sAllDots' struct containing parameters for each trajectory
    sAllDots = struct;
    sAllDots.strStimSet = 'dot_grid';
    intStimIdx = 1;
    for intStim = 1:intTrajY/2 %left-right (0)
        sAllDots.stimID(intStimIdx) = intStimIdx;
        vecBoundingRect = [];
        vecBoundingRect(1,:) = (0-sStimParams.intSize_pix):sStimParams.vecDotSpeeds_ppf:(sStimParams.intScreenWidth_pix+sStimParams.vecDotSpeeds_ppf);
        vecBoundingRect(2,:) = repmat(vecCoordsY(intStim)-sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        vecBoundingRect(3,:) = vecBoundingRect(1,:)+sStimParams.intSize_pix;
        vecBoundingRect(4,:) = repmat(vecCoordsY(intStim)+sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        sAllDots.vecBoundingRect{intStimIdx} = vecBoundingRect;
        sAllDots.vecDirection(intStimIdx) = 0;
        sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(vecBoundingRect(1,:)));
        intStimIdx = intStimIdx+1;
    end
    for intStim = 1:intTrajX/2 %up-down (90)
        sAllDots.stimID(intStimIdx) = intStimIdx;
        vecBoundingRect = [];
        vecBoundingRect(2,:) = (0-sStimParams.intSize_pix):sStimParams.vecDotSpeeds_ppf:(sStimParams.intScreenHeight_pix+sStimParams.vecDotSpeeds_ppf);
        vecBoundingRect(1,:) = repmat(vecCoordsX(intStim)-sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        vecBoundingRect(3,:) = repmat(vecCoordsX(intStim)+sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        vecBoundingRect(4,:) = vecBoundingRect(2,:)+sStimParams.intSize_pix;
        sAllDots.vecBoundingRect{intStimIdx} = vecBoundingRect;
        sAllDots.vecDirection(intStimIdx) = 90;
        sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(vecBoundingRect(1,:)));
        intStimIdx = intStimIdx+1;
    end
    for intStim = 1:intTrajY/2 %right-left (180)
        sAllDots.stimID(intStimIdx) = intStimIdx;
        vecBoundingRect = [];
        vecBoundingRect(1,:) = -(-sStimParams.intScreenWidth_pix:sStimParams.vecDotSpeeds_ppf:(0+sStimParams.intSize_pix+sStimParams.vecDotSpeeds_ppf));
        vecBoundingRect(2,:) = repmat(vecCoordsY(intStim)-sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        vecBoundingRect(3,:) = vecBoundingRect(1,:)+sStimParams.intSize_pix;
        vecBoundingRect(4,:) = repmat(vecCoordsY(intStim)+sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        sAllDots.vecBoundingRect{intStimIdx} = vecBoundingRect;
        sAllDots.vecDirection(intStimIdx) = 180;
        sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(vecBoundingRect(1,:)));
        intStimIdx = intStimIdx+1;
    end
    for intStim = 1:intTrajX/2 %down-up (270)
        sAllDots.stimID(intStimIdx) = intStimIdx;
        vecBoundingRect = [];
        vecBoundingRect(2,:) = -(-sStimParams.intScreenHeight_pix:sStimParams.vecDotSpeeds_ppf:(0+sStimParams.intSize_pix+sStimParams.vecDotSpeeds_ppf));
        vecBoundingRect(1,:) = repmat(vecCoordsX(intStim)-sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        vecBoundingRect(3,:) = repmat(vecCoordsX(intStim)+sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        vecBoundingRect(4,:) = vecBoundingRect(2,:)+sStimParams.intSize_pix;
        sAllDots.vecBoundingRect{intStimIdx} = vecBoundingRect;
        sAllDots.vecDirection(intStimIdx) = 270;
        sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(vecBoundingRect(1,:)));
        intStimIdx = intStimIdx+1;
    end
    sAllDots.intStimulusConditions = intStimIdx-1;
    sAllDots.vecStimName = repmat({'classic'},size(sAllDots.vecBoundingRect));
    sAllDots.vecSpeed_deg = repmat(sStimParams.vecDotSpeeds_deg,size(sAllDots.vecBoundingRect));
    sAllDots.vecSpeed_pix = repmat(sStimParams.vecDotSpeeds_pix,size(sAllDots.vecBoundingRect));
    sAllDots.vecReversalFrame = NaN(size(sAllDots.vecBoundingRect));

elseif intStimSet == 2
    %% (2) 'dot_variations' - different varieties of moving dots
    if length(sStimParams.vecDotSpeeds_pix)>1 %just to be sure
        error('vecDotSpeeds_pix should not be >1!');
    end

    %input desired stimulus conditions
    vecStimConditions = {'classic','appear','disappear','offset'}; %available: 'classic', 'appear', 'disappear', 'shuffle', 'offset'

    %create base trajectories
    if sum(sStimParams.vecDirections==0)>0 %left-right (0)
        vecBoundingRect = [];
        vecBoundingRect(1,:) = (0-sStimParams.intSize_pix):sStimParams.vecDotSpeeds_ppf:(sStimParams.intScreenWidth_pix+sStimParams.vecDotSpeeds_ppf);
        vecBoundingRect(2,:) = repmat(sStimParams.intRespPosY_pix-sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        vecBoundingRect(3,:) = vecBoundingRect(1,:)+sStimParams.intSize_pix;
        vecBoundingRect(4,:) = repmat(sStimParams.intRespPosY_pix+sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        vecLeftRight = vecBoundingRect;
    end
    if sum(sStimParams.vecDirections==180)>0 %right-left (180)
        vecBoundingRect = [];
        vecBoundingRect(1,:) = -(-sStimParams.intScreenWidth_pix:sStimParams.vecDotSpeeds_ppf:(0+sStimParams.intSize_pix+sStimParams.vecDotSpeeds_ppf));
        vecBoundingRect(2,:) = repmat(sStimParams.intRespPosY_pix-sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        vecBoundingRect(3,:) = vecBoundingRect(1,:)+sStimParams.intSize_pix;
        vecBoundingRect(4,:) = repmat(sStimParams.intRespPosY_pix+sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        vecRightLeft = vecBoundingRect;
    end

    %create 'sAllDots' struct containing parameters for each trajectory
    sAllDots = struct;
    sAllDots.strStimSet = 'dot_variations';
    intStimIdx = 1;
    vecRespBordersX = [sStimParams.intRespPosX_pix-sStimParams.intRespSize_pix/2 sStimParams.intRespPosX_pix+sStimParams.intRespSize_pix/2];
    %'classic', continuous movement
    if sum(strcmp(vecStimConditions(:),'classic'))>0
        if sum(sStimParams.vecDirections==0)>0 %left-right (0)
            sAllDots.stimID(intStimIdx) = intStimIdx;
            sAllDots.vecBoundingRect{intStimIdx} = vecLeftRight;
            sAllDots.vecDirection(intStimIdx) = 0;
            sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(sAllDots.vecBoundingRect{intStimIdx}(1,:)));
            sAllDots.vecStimName{intStimIdx} = 'classic';
            intStimIdx = intStimIdx+1;
        end
        if sum(sStimParams.vecDirections==180)>0  %right-left (180)
            sAllDots.stimID(intStimIdx) = intStimIdx;
            sAllDots.vecBoundingRect{intStimIdx} = vecRightLeft;
            sAllDots.vecDirection(intStimIdx) = 180;
            sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(sAllDots.vecBoundingRect{intStimIdx}(1,:)));
            sAllDots.vecStimName{intStimIdx} = 'classic';
            intStimIdx = intStimIdx+1;
        end
    end
    %'appear', stimulus appears just before estimated response zone(timing info same as for 'classic')
    if sum(strcmp(vecStimConditions(:),'appear'))>0
        if sum(sStimParams.vecDirections==0)>0 %left-right (0)
            sAllDots.stimID(intStimIdx) = intStimIdx;
            sAllDots.vecBoundingRect{intStimIdx} = vecLeftRight;
            sAllDots.vecDirection(intStimIdx) = 0;
            vecColor= repmat(sStimParams.dblColor,size(sAllDots.vecBoundingRect{intStimIdx}(1,:)));
            vecColor(vecLeftRight(3,:)<vecRespBordersX(1)) = sStimParams.intBackground; %(3,:)=leading edge
            sAllDots.vecColor{intStimIdx} = vecColor;
            sAllDots.vecStimName{intStimIdx} = 'appear';
            intStimIdx = intStimIdx+1;
        end
        if sum(sStimParams.vecDirections==180)>0 %right-left (180)
            sAllDots.stimID(intStimIdx) = intStimIdx;
            sAllDots.vecBoundingRect{intStimIdx} = vecRightLeft;
            sAllDots.vecDirection(intStimIdx) = 180;
            vecColor= repmat(sStimParams.dblColor,size(sAllDots.vecBoundingRect{intStimIdx}(1,:)));
            vecColor(vecRightLeft(1,:)>vecRespBordersX(2)) = sStimParams.intBackground; %(1,:)=leading edge
            sAllDots.vecColor{intStimIdx} = vecColor;
            sAllDots.vecStimName{intStimIdx} = 'appear';
            intStimIdx = intStimIdx+1;
        end
    end
    %'disappear', stimulus disappears at location of estimated response zone
    if sum(strcmp(vecStimConditions(:),'disappear'))>0
        sAllDots.stimID(intStimIdx) = intStimIdx;
        vecRespBordersX = [sStimParams.intRespPosX_pix-sStimParams.intRespSize_pix/2 sStimParams.intRespPosX_pix+sStimParams.intRespSize_pix/2];
        if sum(sStimParams.vecDirections==0)>0 %left-right (0)
            sAllDots.vecBoundingRect{intStimIdx} = vecLeftRight;
            sAllDots.vecDirection(intStimIdx) = 0;
            vecColor= repmat(sStimParams.dblColor,size(sAllDots.vecBoundingRect{intStimIdx}(1,:)));
            vecColor(vecLeftRight(3,:)>vecRespBordersX(1)) = sStimParams.intBackground; %(3,:)=leading edge
            sAllDots.vecColor{intStimIdx} = vecColor;
            sAllDots.vecStimName{intStimIdx} = 'disappear';
            intStimIdx = intStimIdx+1;
        end
        if sum(sStimParams.vecDirections==180)>0 %right-left (180)
            sAllDots.stimID(intStimIdx) = intStimIdx;
            sAllDots.vecBoundingRect{intStimIdx} = vecRightLeft;
            sAllDots.vecDirection(intStimIdx) = 180;
            vecColor= repmat(sStimParams.dblColor,size(sAllDots.vecBoundingRect{intStimIdx}(1,:)));
            vecColor(vecRightLeft(1,:)<vecRespBordersX(2)) = sStimParams.intBackground; %(1,:)=leading edge, (3,:)=trailing edge
            sAllDots.vecColor{intStimIdx} = vecColor;
            sAllDots.vecStimName{intStimIdx} = 'disappear';
            intStimIdx = intStimIdx+1;
        end
    end
    %'shuffle', order of dot position along trajectory before response field is temporally shuffled
    if sum(strcmp(vecStimConditions(:),'shuffle'))> 0
        vecRespBordersX = [sStimParams.intRespPosX_pix-sStimParams.intRespSize_pix/2 sStimParams.intRespPosX_pix+sStimParams.intRespSize_pix/2];
        if sum(sStimParams.vecDirections==0)>0 %left-right (0)
            sAllDots.stimID(intStimIdx) = intStimIdx;
            vecBoundingRect = vecLeftRight;
            vecNewOrder = Shuffle(find(vecLeftRight(3,:)<vecRespBordersX(1))); %(3,:)=leading edge
            vecBoundingRect(:,1:length(vecNewOrder)) = vecBoundingRect(:,vecNewOrder);
            sAllDots.vecBoundingRect{intStimIdx} = vecBoundingRect;
            sAllDots.vecDirection(intStimIdx) = 0;
            sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(sAllDots.vecBoundingRect{intStimIdx}(1,:)));
            sAllDots.vecStimName{intStimIdx} = 'shuffle';
            intStimIdx = intStimIdx+1;
        end
        if sum(sStimParams.vecDirections==180)>0 %right-left (180)
            sAllDots.stimID(intStimIdx) = intStimIdx;
            vecBoundingRect = vecRightLeft;
            vecNewOrder = Shuffle(find(vecRightLeft(1,:)>vecRespBordersX(2))); %(1,:)=leading edge
            vecBoundingRect(:,1:length(vecNewOrder)) = vecBoundingRect(:,vecNewOrder);
            sAllDots.vecBoundingRect{intStimIdx} = vecBoundingRect;
            sAllDots.vecDirection(intStimIdx) = 180;
            sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(sAllDots.vecBoundingRect{intStimIdx}(1,:)));
            sAllDots.vecStimName{intStimIdx} = 'shuffle';
            intStimIdx = intStimIdx+1;
        end
    end
    %'offset', dot trajectory before response zone is spatially offset
    intSpatialOffset_pix = round(12*sStimParams.dblPixelsPerDeg); % ~12 dva
    if sum(strcmp(vecStimConditions(:),'offset'))> 0
        vecRespBordersX = [sStimParams.intRespPosX_pix-sStimParams.intRespSize_pix/2 sStimParams.intRespPosX_pix+sStimParams.intRespSize_pix/2];
        if (sStimParams.intScreenHeight_pix-sStimParams.intRespPosY_pix) <= sStimParams.intScreenHeight_pix/2 %check if offset should be below or above
            intOffset = -intSpatialOffset_pix;
        else
            intOffset = intSpatialOffset_pix;
        end
        if sum(sStimParams.vecDirections==0)>0 %left-right (0)
            sAllDots.stimID(intStimIdx) = intStimIdx;
            vecBoundingRect = vecLeftRight;
            vecBoundingRect([2 4],vecLeftRight(3,:)<vecRespBordersX(1)) = vecBoundingRect([2 4],vecLeftRight(3,:)<vecRespBordersX(1))+intOffset; %(3,:)=leading edge
            sAllDots.vecBoundingRect{intStimIdx} = vecBoundingRect;
            sAllDots.vecDirection(intStimIdx) = 0;
            sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(sAllDots.vecBoundingRect{intStimIdx}(1,:)));
            sAllDots.vecStimName{intStimIdx} = 'offset';
            intStimIdx = intStimIdx+1;
        end
        if sum(sStimParams.vecDirections==180)>0 %right-left (180)
            sAllDots.stimID(intStimIdx) = intStimIdx;
            vecBoundingRect = vecRightLeft;
            vecBoundingRect([2 4],vecRightLeft(1,:)>vecRespBordersX(2)) = vecBoundingRect([2 4],vecRightLeft(1,:)>vecRespBordersX(2))+intOffset; %(1,:)=leading edge
            sAllDots.vecBoundingRect{intStimIdx} = vecBoundingRect;
            sAllDots.vecDirection(intStimIdx) = 180;
            sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(sAllDots.vecBoundingRect{intStimIdx}(1,:)));
            sAllDots.vecStimName{intStimIdx} = 'offset';
            intStimIdx = intStimIdx+1;
        end
    end
    sAllDots.intStimulusConditions = intStimIdx-1;
    sAllDots.vecSpeed_deg = repmat(sStimParams.vecDotSpeeds_deg,size(sAllDots.vecBoundingRect));
    sAllDots.vecSpeed_pix = repmat(sStimParams.vecDotSpeeds_pix,size(sAllDots.vecBoundingRect));
    sAllDots.vecReversalFrame = NaN(size(sAllDots.vecBoundingRect));

elseif intStimSet == 3
    %% (3) 'dot_speeds' - dots moving continuously at different speeds
    %set speeds
    sStimParams.vecDotSpeeds_deg = [3 6 15 30 60 120]; % deg/s
    sStimParams.vecDotSpeeds_pix = sStimParams.vecDotSpeeds_deg*sStimParams.dblPixelsPerDeg; % pixels/s
    sStimParams.vecDotSpeeds_ppf = round(sStimParams.vecDotSpeeds_pix/sStimParams.intStimFrameRate); % pixels/frame
    if sum(sStimParams.vecDotSpeeds_ppf == 0) > 0
        error('vecDotSpeeds_ppf cannot be 0')
    end

    %get number of trajectories
    intTrajY = sum(sStimParams.vecDirections==[0 180])*length(sStimParams.vecDotSpeeds_ppf);

    %create 'sAllDots' struct containing parameters for each trajectory
    sAllDots = struct;
    sAllDots.strStimSet = 'dot_speeds';
    intStimIdx = 1;
    for intStim = 1:intTrajY/sum(sStimParams.vecDirections==[0 180]) %left-right (0)
        sAllDots.stimID(intStimIdx) = intStimIdx;
        vecBoundingRect = [];
        vecBoundingRect(1,:) = (0-sStimParams.intSize_pix):sStimParams.vecDotSpeeds_ppf(intStim):(sStimParams.intScreenWidth_pix+sStimParams.vecDotSpeeds_ppf(intStim));
        vecBoundingRect(2,:) = repmat(sStimParams.intRespPosY_pix-sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        vecBoundingRect(3,:) = vecBoundingRect(1,:)+sStimParams.intSize_pix;
        vecBoundingRect(4,:) = repmat(sStimParams.intRespPosY_pix+sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        sAllDots.vecBoundingRect{intStimIdx} = vecBoundingRect;
        sAllDots.vecDirection(intStimIdx) = 0;
        sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(vecBoundingRect(1,:)));
        sAllDots.vecSpeed_deg(intStimIdx) = sStimParams.vecDotSpeeds_deg(intStim);
        sAllDots.vecSpeed_pix(intStimIdx) = sStimParams.vecDotSpeeds_pix(intStim);
        intStimIdx = intStimIdx+1;
    end
    for intStim = 1:intTrajY/sum(sStimParams.vecDirections==[0 180]) %right-left (180)
        sAllDots.stimID(intStimIdx) = intStimIdx;
        vecBoundingRect = [];
        vecBoundingRect(1,:) = -(-sStimParams.intScreenWidth_pix:sStimParams.vecDotSpeeds_ppf(intStim):(0+sStimParams.intSize_pix+sStimParams.vecDotSpeeds_ppf(intStim)));
        vecBoundingRect(2,:) = repmat(sStimParams.intRespPosY_pix-sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        vecBoundingRect(3,:) = vecBoundingRect(1,:)+sStimParams.intSize_pix;
        vecBoundingRect(4,:) = repmat(sStimParams.intRespPosY_pix+sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        sAllDots.vecBoundingRect{intStimIdx} = vecBoundingRect;
        sAllDots.vecDirection(intStimIdx) = 180;
        sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(vecBoundingRect(1,:)));
        sAllDots.vecSpeed_deg(intStimIdx) = sStimParams.vecDotSpeeds_deg(intStim);
        sAllDots.vecSpeed_pix(intStimIdx) = sStimParams.vecDotSpeeds_pix(intStim);
        intStimIdx = intStimIdx+1;
    end
    sAllDots.intStimulusConditions = intStimIdx-1;
    sAllDots.vecStimName = repmat({'classic'},size(sAllDots.vecBoundingRect));
    sAllDots.vecReversalFrame = NaN(size(sAllDots.vecBoundingRect));

elseif intStimSet == 4
    %% (4) 'dot_reversal' - dots reversing direction of motion
    if length(sStimParams.vecDotSpeeds_pix)>1 %just to be sure
        error('vecDotSpeeds_pix should not be >1!');
    end

    %input desired stimulus conditions
    vecStimConditions = {'classic','reversal_1','reversal_2','reversal_3'}; %available: 'classic', 'reversal_1', 'reversal_2', 'reversal_3'

    %create base trajectories
    if sum(sStimParams.vecDirections==0)>0 %left-right (0)
        vecBoundingRect = [];
        vecBoundingRect(1,:) = (0-sStimParams.intSize_pix):sStimParams.vecDotSpeeds_ppf:(sStimParams.intScreenWidth_pix+sStimParams.vecDotSpeeds_ppf);
        vecBoundingRect(2,:) = repmat(sStimParams.intRespPosY_pix-sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        vecBoundingRect(3,:) = vecBoundingRect(1,:)+sStimParams.intSize_pix;
        vecBoundingRect(4,:) = repmat(sStimParams.intRespPosY_pix+sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        vecLeftRight = vecBoundingRect;
    end
    if sum(sStimParams.vecDirections==180)>0 %right-left (180)
        vecBoundingRect = [];
        vecBoundingRect(1,:) = -(-sStimParams.intScreenWidth_pix:sStimParams.vecDotSpeeds_ppf:(0+sStimParams.intSize_pix+sStimParams.vecDotSpeeds_ppf));
        vecBoundingRect(2,:) = repmat(sStimParams.intRespPosY_pix-sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        vecBoundingRect(3,:) = vecBoundingRect(1,:)+sStimParams.intSize_pix;
        vecBoundingRect(4,:) = repmat(sStimParams.intRespPosY_pix+sStimParams.intSize_pix/2,size(vecBoundingRect(1,:)));
        vecRightLeft = vecBoundingRect;
    end

    %create 'sAllDots' struct containing parameters for each trajectory
    sAllDots = struct;
    sAllDots.strStimSet = 'dot_reversal';
    intStimIdx = 1;
    vecRespBordersX = [sStimParams.intRespPosX_pix-sStimParams.intRespSize_pix/2 sStimParams.intRespPosX_pix+sStimParams.intRespSize_pix/2];
    %'classic', continuous movement
    if sum(strcmp(vecStimConditions(:),'classic'))>0
        if sum(sStimParams.vecDirections==0)>0 %left-right (0)
            sAllDots.stimID(intStimIdx) = intStimIdx;
            sAllDots.vecBoundingRect{intStimIdx} = vecLeftRight;
            sAllDots.vecReversalFrame(intStimIdx) = NaN;
            sAllDots.vecDirection(intStimIdx) = 0;
            sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(sAllDots.vecBoundingRect{intStimIdx}(1,:)));
            sAllDots.vecStimName{intStimIdx} = 'classic';
            intStimIdx = intStimIdx+1;
        end
        if sum(sStimParams.vecDirections==180)>0  %right-left (180)
            sAllDots.stimID(intStimIdx) = intStimIdx;
            sAllDots.vecBoundingRect{intStimIdx} = vecRightLeft;
            sAllDots.vecReversalFrame(intStimIdx) = NaN;
            sAllDots.vecDirection(intStimIdx) = 180;
            sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(sAllDots.vecBoundingRect{intStimIdx}(1,:)));
            sAllDots.vecStimName{intStimIdx} = 'classic';
            intStimIdx = intStimIdx+1;
        end
    end
    %'reversal_1', dot motion direction reverses just before entering response field
    if sum(strcmp(vecStimConditions(:),'reversal_1'))>0
        if sum(sStimParams.vecDirections==0)>0 %left-right (0), initial direction
            sAllDots.stimID(intStimIdx) = intStimIdx;
            vecBoundingRect = vecLeftRight(:,vecLeftRight(3,:)<vecRespBordersX(1)+sStimParams.intSize_pix/2); %(3,:)=leading edge
            sAllDots.vecBoundingRect{intStimIdx} = [vecBoundingRect flip(vecBoundingRect(:,1:end-1),2)];
            sAllDots.vecReversalFrame(intStimIdx) = size(vecBoundingRect,2)+1;
            sAllDots.vecDirection(intStimIdx) = 0;
            sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(sAllDots.vecBoundingRect{intStimIdx}(1,:)));
            sAllDots.vecStimName{intStimIdx} = 'reversal_1';
            intStimIdx = intStimIdx+1;
        end
        if sum(sStimParams.vecDirections==180)>0 %right-left (180)
            sAllDots.stimID(intStimIdx) = intStimIdx;
            vecBoundingRect = vecRightLeft(:,vecRightLeft(1,:)>vecRespBordersX(2)-sStimParams.intSize_pix/2); %(1,:)=leading edge
            sAllDots.vecBoundingRect{intStimIdx} = [vecBoundingRect flip(vecBoundingRect(:,1:end-1),2)];
            sAllDots.vecReversalFrame(intStimIdx) = size(vecBoundingRect,2)+1;
            sAllDots.vecDirection(intStimIdx) = 180;
            sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(sAllDots.vecBoundingRect{intStimIdx}(1,:)));
            sAllDots.vecStimName{intStimIdx} = 'reversal_1';
            intStimIdx = intStimIdx+1;
        end
    end
    %'reversal_2', motion direction reverses in the middle of estimated response field
    if sum(strcmp(vecStimConditions(:),'reversal_2'))>0
        if sum(sStimParams.vecDirections==0)>0 %left-right (0), initial direction
            sAllDots.stimID(intStimIdx) = intStimIdx;
            vecBoundingRect = vecLeftRight(:,vecLeftRight(3,:)<sStimParams.intRespPosX_pix+sStimParams.intSize_pix/2); %(3,:)=leading edge
            sAllDots.vecBoundingRect{intStimIdx} = [vecBoundingRect flip(vecBoundingRect(:,1:end-1),2)];
            sAllDots.vecReversalFrame(intStimIdx) = size(vecBoundingRect,2)+1;
            sAllDots.vecDirection(intStimIdx) = 0;
            sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(sAllDots.vecBoundingRect{intStimIdx}(1,:)));
            sAllDots.vecStimName{intStimIdx} = 'reversal_2';
            intStimIdx = intStimIdx+1;
        end
        if sum(sStimParams.vecDirections==180)>0 %right-left (180)
            sAllDots.stimID(intStimIdx) = intStimIdx;
            vecBoundingRect = vecRightLeft(:,vecRightLeft(1,:)>sStimParams.intRespPosX_pix-sStimParams.intSize_pix/2); %(1,:)=leading edge
            sAllDots.vecBoundingRect{intStimIdx} = [vecBoundingRect flip(vecBoundingRect(:,1:end-1),2)];
            sAllDots.vecReversalFrame(intStimIdx) = size(vecBoundingRect,2)+1;
            sAllDots.vecDirection(intStimIdx) = 180;
            sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(sAllDots.vecBoundingRect{intStimIdx}(1,:)));
            sAllDots.vecStimName{intStimIdx} = 'reversal_2';
            intStimIdx = intStimIdx+1;
        end
    end
    %'reversal_3', motion direction reverses when leaving response field
    if sum(strcmp(vecStimConditions(:),'reversal_3'))>0
        if sum(sStimParams.vecDirections==0)>0 %left-right (0), initial direction
            sAllDots.stimID(intStimIdx) = intStimIdx;
            vecBoundingRect = vecLeftRight(:,vecLeftRight(1,:)<vecRespBordersX(2)-sStimParams.intSize_pix/2); %(3,:)=trailing edge
            sAllDots.vecBoundingRect{intStimIdx} = [vecBoundingRect flip(vecBoundingRect(:,1:end-1),2)];
            sAllDots.vecReversalFrame(intStimIdx) = size(vecBoundingRect,2)+1;
            sAllDots.vecDirection(intStimIdx) = 0;
            sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(sAllDots.vecBoundingRect{intStimIdx}(1,:)));
            sAllDots.vecStimName{intStimIdx} = 'reversal_3';
            intStimIdx = intStimIdx+1;
        end
        if sum(sStimParams.vecDirections==180)>0 %right-left (180)
            sAllDots.stimID(intStimIdx) = intStimIdx;
            vecBoundingRect = vecRightLeft(:,vecRightLeft(3,:)>vecRespBordersX(1)+sStimParams.intSize_pix/2); %(3,:)=leading edge
            sAllDots.vecBoundingRect{intStimIdx} = [vecBoundingRect flip(vecBoundingRect(:,1:end-1),2)];
            sAllDots.vecReversalFrame(intStimIdx) = size(vecBoundingRect,2)+1;
            sAllDots.vecDirection(intStimIdx) = 180;
            sAllDots.vecColor{intStimIdx} = repmat(sStimParams.dblColor,size(sAllDots.vecBoundingRect{intStimIdx}(1,:)));
            sAllDots.vecStimName{intStimIdx} = 'reversal_3';
            intStimIdx = intStimIdx+1;
        end
    end
    sAllDots.intStimulusConditions = intStimIdx-1;
    sAllDots.vecSpeed_deg = repmat(sStimParams.vecDotSpeeds_deg,size(sAllDots.vecBoundingRect));
    sAllDots.vecSpeed_pix = repmat(sStimParams.vecDotSpeeds_pix,size(sAllDots.vecBoundingRect));

    %% add more sets here
    %elseif intStimSet == 5
end