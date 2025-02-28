function  [sDots,sParams] = RH_CreateOneDotTwoDots(intSet,sParams)
% generate stimulus trajectories for RH_OneDotTwoDots
% there are 4 available stim sets:
%   (1) "inside-outside" : 1x in RF, 15x outside
%   (2) "inside control" : 1x in RF, 15x NO stim
%   (3) "outside-outside" : 1x outside RF, 15x outside RF on the other side
%   (4) "outside control": 1x outside RF, 15x NO stim
%
%	EXTRA:
%	(5) "inside only"
%	(6) "inside-outside closer"
%	(7) "inside control closer"
%
% history:
% 31 Aug 2023
%   - created by Robin Haak
% 30 May 2024
%	- added stim sets 5-7

%% default params
sParams.dblSize_deg = 15; %deg, approx. stimulus size
sParams.intSize_pix = round(sParams.dblSize_deg*sParams.dblPixelsPerDeg);
if mod(sParams.intSize_pix,2),sParams.intSize_pix = sParams.intSize_pix-1;end %avoid decimal pix vals
sParams.dblColor = 0; %[0 1]
sParams.intColor = round(mean(sParams.dblColor)*255);
sParams.dblSpeed_deg = 60; %deg/s, approx. stim speed
sParams.dblSpeed_pix = sParams.dblSpeed_deg*sParams.dblPixelsPerDeg; %pix/s
sParams.intSpeed_ppf = round(sParams.dblSpeed_pix/sParams.intStimFrameRate); %pix/frame
sParams.dblDotSep_deg = 65; %deg, separation between two dots
sParams.intDotSep_pix = round(sParams.dblDotSep_deg*sParams.dblPixelsPerDeg); %pix

% intPixFromEdge = round(10*sParams.dblPixelsPerDeg); %pix, from screen edge
intPixFromEdge = round(5*sParams.dblPixelsPerDeg); %pix, from screen edge

if intSet==1 || intSet==2
    %% inside-outside / inside control
    %optimize dot locations
    intLoc1 = sParams.intRfPosX_pix;
    intSize = sParams.intSize_pix;
    intScreenWidth = sParams.intScreenWidth_pix;
    intScreenHeight = sParams.intScreenHeight_pix;
    if intScreenWidth-intLoc1<intScreenWidth/2, intLoc2 = intLoc1-sParams.intDotSep_pix;
    else, intLoc2 = intLoc1+sParams.intDotSep_pix; end

    %check if stims are on screen
    boolOnScreen1 = intLoc1+intSize/2<intScreenWidth-intPixFromEdge & intLoc1-intSize/2>intPixFromEdge;
    if ~boolOnScreen1 %if 'inside' stim is not on the screen, there is nothing we can do
        if intLoc1+intSize/2>intScreenWidth || intLoc1-intSize/2<0
            error('The "inside" stimulus falls off the screen, move screen!');
        else
            fprintf('\n\n <strong>"INSIDE" STIMULUS IS VERY CLOSE TO THE BORDER, consider moving the screen!</strong>\n\n');
            pause(2);
        end
    end

    boolOnScreen2 = intLoc2+intSize/2<intScreenWidth-intPixFromEdge & intLoc2-intSize/2>intPixFromEdge;
    if ~boolOnScreen2 %if the 'outside' stim falls of the screen, propose new locations
        if intLoc2<intLoc1, intAltLoc2 = intPixFromEdge+intSize/2;
        else, intAltLoc2 = intScreenWidth-(intPixFromEdge+intSize/2); end

        %plot alternative location
        f = figure; hold on;
        rectangle('Position',[0 0 intScreenWidth intScreenHeight],'FaceColor',[0.5 0.5 0.5]);
        viscircles([intLoc1 intScreenHeight/2;intLoc2 intScreenHeight/2],[sParams.intSize_pix/2 sParams.intSize_pix/2],'Color','b');
        viscircles([intAltLoc2 intScreenHeight/2],sParams.intSize_pix/2,'Color','r');
        s = scatter(sParams.intRfPosX_pix,intScreenHeight/2,'kx');
        legend(s,'RF center'); xlabel('pix'); ylabel('pix');
        axis('auto xy');axis image;fixfig

        %get UI
        dblNewSep_deg = abs(diff([intLoc1 intAltLoc2]))/sParams.dblPixelsPerDeg;
        strBox = ['Would you like to adjust the stim location to B>R (new sep= ' num2str(dblNewSep_deg) ' deg)?'];
        sOptions.Default = 'Yes';
        sOptions.Interpreter = 'tex';
        strAns = questdlg(strBox,'','Yes','No',sOptions);
        switch strAns
            case 'Yes'
                intLoc2 = intAltLoc2;
                fprintf('\nSuccesfully adjusted stimulus location!\n');
            case 'No'
                return
        end
    close(f);
    end

    %create stimulus trajectories
    vecBoundingRect1 = []; %inside
    vecBoundingRect1(2,:) = -intSize:sParams.intSpeed_ppf:(intScreenHeight+sParams.intSpeed_ppf);
    vecBoundingRect1(1,:) = repmat(intLoc1-intSize/2,size(vecBoundingRect1(2,:)));
    vecBoundingRect1(3,:) = repmat(intLoc1+intSize/2,size(vecBoundingRect1(2,:)));
    vecBoundingRect1(4,:) = vecBoundingRect1(2,:)+intSize;

    vecBoundingRect2 = []; %outside
    vecBoundingRect2(2,:) = -intSize:sParams.intSpeed_ppf:(intScreenHeight+sParams.intSpeed_ppf);
    vecBoundingRect2(1,:) = repmat(intLoc2-intSize/2,size(vecBoundingRect1(2,:)));
    vecBoundingRect2(3,:) = repmat(intLoc2+intSize/2,size(vecBoundingRect1(2,:)));
    vecBoundingRect2(4,:) = vecBoundingRect1(2,:)+intSize;

    %create sDots
    sDots = struct;
    if intSet==1,sDots.strStimSet = 'inside-outside';
    elseif intSet==2,sDots.strStimSet = 'inside control';end

    sDots.stimID(1) = 1; %inside
    sDots.cellBoundingRect{1} = vecBoundingRect1;
    sDots.cellColor{1} = repmat(sParams.intColor,size(vecBoundingRect1(1,:)));
    sDots.vecSpeed_deg(1) = sParams.dblSpeed_deg;
    sDots.vecSpeed_pix(1) = sParams.dblSpeed_pix;

    sDots.stimID(2) = 2; %outside
    sDots.cellBoundingRect{2} = vecBoundingRect2;
    if intSet==1,sDots.cellColor{2} = repmat(sParams.intColor,size(vecBoundingRect2(1,:)));
    elseif intSet==2,sDots.cellColor{2} = repmat(sParams.intBackground,size(vecBoundingRect2(1,:)));end
    sDots.vecSpeed_deg(2) = sParams.dblSpeed_deg;
    sDots.vecSpeed_pix(2) = sParams.dblSpeed_pix;

    sDots.intStimulusConditions = 2;
    sDots.vecBlockStructure = [1 repmat(2,[1 15])];

elseif intSet==3 || intSet==4
    %% outside-outside / outside control
    %optimize dot locations
    if mod(sParams.intDotSep_pix,2),intDotSep=sParams.intDotSep_pix-1;
    else,intDotSep = sParams.intDotSep;end
    intSize = sParams.intSize_pix;
    intScreenWidth = sParams.intScreenWidth_pix;
    intScreenHeight = sParams.intScreenHeight_pix;
    intLoc1 = sParams.intRfPosX_pix-intDotSep/2; %left
    intLoc2 = sParams.intRfPosX_pix+intDotSep/2; %right

    %check if stim are on the screen
    boolOnScreen1 =  intLoc1-intSize/2>intPixFromEdge;
    boolOnScreen2 = intLoc2+intSize/2<intScreenWidth-intPixFromEdge;
    hTic = tic;
    while ~boolOnScreen1 || ~boolOnScreen2
        if ~boolOnScreen1
            intLoc1 = intLoc1+1;
            intLoc2 = intLoc2+1;
        end
        if ~boolOnScreen2
            intLoc2 = intLoc2-1;
            intLoc1 = intLoc1-1;
        end
        boolOnScreen1 =  intLoc1-intSize/2>intPixFromEdge;
        boolOnScreen2 = intLoc2+intSize/2<intScreenWidth-intPixFromEdge;
        if toc(hTic) > 5
            error('cannot optimize stimulus locations');
        end
    end

    %set intLoc1 to be the stim closest to the RF
    if abs(sParams.intRfPosX_pix-intLoc1) > abs(sParams.intRfPosX_pix-intLoc2)
        intNewLoc1 = intLoc2;
        intLoc2 = intLoc1;
        intLoc1 = intNewLoc1;
    end

    %create stimulus trajectories
    vecBoundingRect1 = []; %inside
    vecBoundingRect1(2,:) = -intSize:sParams.intSpeed_ppf:(intScreenHeight+sParams.intSpeed_ppf);
    vecBoundingRect1(1,:) = repmat(intLoc1-intSize/2,size(vecBoundingRect1(2,:)));
    vecBoundingRect1(3,:) = repmat(intLoc1+intSize/2,size(vecBoundingRect1(2,:)));
    vecBoundingRect1(4,:) = vecBoundingRect1(2,:)+intSize;

    vecBoundingRect2 = []; %outside
    vecBoundingRect2(2,:) = -intSize:sParams.intSpeed_ppf:(intScreenHeight+sParams.intSpeed_ppf);
    vecBoundingRect2(1,:) = repmat(intLoc2-intSize/2,size(vecBoundingRect1(2,:)));
    vecBoundingRect2(3,:) = repmat(intLoc2+intSize/2,size(vecBoundingRect1(2,:)));
    vecBoundingRect2(4,:) = vecBoundingRect1(2,:)+intSize;

    %create sDots
    sDots = struct;
    if intSet==3,sDots.strStimSet = 'outside-outside';
    elseif intSet==4,sDots.strStimSet = 'outside control';end

    sDots.stimID(1) = 1; %outside 1
    sDots.cellBoundingRect{1} = vecBoundingRect1;
    sDots.cellColor{1} = repmat(sParams.intColor,size(vecBoundingRect1(1,:)));
    sDots.vecSpeed_deg(1) = sParams.dblSpeed_deg;
    sDots.vecSpeed_pix(1) = sParams.dblSpeed_pix;

    sDots.stimID(2) = 2; %outside 2
    sDots.cellBoundingRect{2} = vecBoundingRect2;
    if intSet==3,sDots.cellColor{2} = repmat(sParams.intColor,size(vecBoundingRect2(1,:)));
    elseif intSet==4,sDots.cellColor{2} = repmat(sParams.intBackground,size(vecBoundingRect2(1,:)));end
    sDots.vecSpeed_deg(2) = sParams.dblSpeed_deg;
    sDots.vecSpeed_pix(2) = sParams.dblSpeed_pix;

    sDots.intStimulusConditions = 2;
    sDots.vecBlockStructure = [1 repmat(2,[1 15])];
	
elseif intSet==5
	%% inside only
	%optimize dot locations
	intLoc1 = sParams.intRfPosX_pix;
	intSize = sParams.intSize_pix;
	intScreenWidth = sParams.intScreenWidth_pix;
	intScreenHeight = sParams.intScreenHeight_pix;
	
	%check if stims are on screen
	boolOnScreen1 = intLoc1+intSize/2<intScreenWidth-intPixFromEdge & intLoc1-intSize/2>intPixFromEdge;
	if ~boolOnScreen1 %if 'inside' stim is not on the screen, there is nothing we can do
		if intLoc1+intSize/2>intScreenWidth || intLoc1-intSize/2<0
			error('The "inside" stimulus falls off the screen, move screen!');
		else
			fprintf('\n\n <strong>"INSIDE" STIMULUS IS VERY CLOSE TO THE BORDER, consider moving the screen!</strong>\n\n');
			pause(2);
		end
	end
	
	%create stimulus trajectories
	vecBoundingRect1 = []; %inside
	vecBoundingRect1(2,:) = -intSize:sParams.intSpeed_ppf:(intScreenHeight+sParams.intSpeed_ppf);
	vecBoundingRect1(1,:) = repmat(intLoc1-intSize/2,size(vecBoundingRect1(2,:)));
	vecBoundingRect1(3,:) = repmat(intLoc1+intSize/2,size(vecBoundingRect1(2,:)));
	vecBoundingRect1(4,:) = vecBoundingRect1(2,:)+intSize;
	
	%create sDots
	sDots = struct;
	sDots.strStimSet ='inside only';
	
	sDots.stimID(1) = 1; %inside
	sDots.cellBoundingRect{1} = vecBoundingRect1;
	sDots.cellColor{1} = repmat(sParams.intColor,size(vecBoundingRect1(1,:)));
	sDots.vecSpeed_deg(1) = sParams.dblSpeed_deg;
	sDots.vecSpeed_pix(1) = sParams.dblSpeed_pix;
	sDots.intStimulusConditions = 1;
	sDots.vecBlockStructure = [1]; %#ok<NBRAK>
	
	sParams.intReps = 30;
	
elseif intSet==6 || intSet==7
	%% inside-outside closer / inside control closer
	%Alexander's recomputed measurements
	sParams.dblSize_deg = 11.3; %15; %deg, approx. stimulus size
	sParams.intSize_pix = round(sParams.dblSize_deg*sParams.dblPixelsPerDeg);
	if mod(sParams.intSize_pix,2),sParams.intSize_pix = sParams.intSize_pix-1;end %avoid decimal pix vals
	
	sParams.dblSpeed_deg = 45; %60; %deg/s, approx. stim speed
	sParams.dblSpeed_pix = sParams.dblSpeed_deg*sParams.dblPixelsPerDeg; %pix/s
	sParams.intSpeed_ppf = round(sParams.dblSpeed_pix/sParams.intStimFrameRate); %pix/frame
	
	sParams.dblDotSep_deg = 58; %65; %deg, separation between two dots
	sParams.intDotSep_pix = round(sParams.dblDotSep_deg*sParams.dblPixelsPerDeg); %pix
	
	%optimize dot locations
	intLoc1 = sParams.intRfPosX_pix;
    intSize = sParams.intSize_pix;
    intScreenWidth = sParams.intScreenWidth_pix;
    intScreenHeight = sParams.intScreenHeight_pix;
    if intScreenWidth-intLoc1<intScreenWidth/2, intLoc2 = intLoc1-sParams.intDotSep_pix;
    else, intLoc2 = intLoc1+sParams.intDotSep_pix; end

    %check if stims are on screen
    boolOnScreen1 = intLoc1+intSize/2<intScreenWidth-intPixFromEdge & intLoc1-intSize/2>intPixFromEdge;
    if ~boolOnScreen1 %if 'inside' stim is not on the screen, there is nothing we can do
        if intLoc1+intSize/2>intScreenWidth || intLoc1-intSize/2<0
            error('The "inside" stimulus falls off the screen, move screen!');
        else
            fprintf('\n\n <strong>"INSIDE" STIMULUS IS VERY CLOSE TO THE BORDER, consider moving the screen!</strong>\n\n');
            pause(2);
        end
    end

    boolOnScreen2 = intLoc2+intSize/2<intScreenWidth-intPixFromEdge & intLoc2-intSize/2>intPixFromEdge;
    if ~boolOnScreen2 %if the 'outside' stim falls of the screen, propose new locations
        if intLoc2<intLoc1, intAltLoc2 = intPixFromEdge+intSize/2;
        else, intAltLoc2 = intScreenWidth-(intPixFromEdge+intSize/2); end

        %plot alternative location
        f = figure; hold on;
        rectangle('Position',[0 0 intScreenWidth intScreenHeight],'FaceColor',[0.5 0.5 0.5]);
        viscircles([intLoc1 intScreenHeight/2;intLoc2 intScreenHeight/2],[sParams.intSize_pix/2 sParams.intSize_pix/2],'Color','b');
        viscircles([intAltLoc2 intScreenHeight/2],sParams.intSize_pix/2,'Color','r');
        s = scatter(sParams.intRfPosX_pix,intScreenHeight/2,'kx');
        legend(s,'RF center'); xlabel('pix'); ylabel('pix');
        axis('auto xy');axis image;fixfig

        %get UI
        dblNewSep_deg = abs(diff([intLoc1 intAltLoc2]))/sParams.dblPixelsPerDeg;
        strBox = ['Would you like to adjust the stim location to B>R (new sep= ' num2str(dblNewSep_deg) ' deg)?'];
        sOptions.Default = 'Yes';
        sOptions.Interpreter = 'tex';
        strAns = questdlg(strBox,'','Yes','No',sOptions);
        switch strAns
            case 'Yes'
                intLoc2 = intAltLoc2;
                fprintf('\nSuccesfully adjusted stimulus location!\n');
            case 'No'
                return
        end
    close(f);
    end

    %create stimulus trajectories
    vecBoundingRect1 = []; %inside
    vecBoundingRect1(2,:) = -intSize:sParams.intSpeed_ppf:(intScreenHeight+sParams.intSpeed_ppf);
    vecBoundingRect1(1,:) = repmat(intLoc1-intSize/2,size(vecBoundingRect1(2,:)));
    vecBoundingRect1(3,:) = repmat(intLoc1+intSize/2,size(vecBoundingRect1(2,:)));
    vecBoundingRect1(4,:) = vecBoundingRect1(2,:)+intSize;

    vecBoundingRect2 = []; %outside
    vecBoundingRect2(2,:) = -intSize:sParams.intSpeed_ppf:(intScreenHeight+sParams.intSpeed_ppf);
    vecBoundingRect2(1,:) = repmat(intLoc2-intSize/2,size(vecBoundingRect1(2,:)));
    vecBoundingRect2(3,:) = repmat(intLoc2+intSize/2,size(vecBoundingRect1(2,:)));
    vecBoundingRect2(4,:) = vecBoundingRect1(2,:)+intSize;

    %create sDots
    sDots = struct;
    if intSet==6,sDots.strStimSet = 'inside-outside closer';
    elseif intSet==7,sDots.strStimSet = 'inside control closer';end

    sDots.stimID(1) = 1; %inside
    sDots.cellBoundingRect{1} = vecBoundingRect1;
    sDots.cellColor{1} = repmat(sParams.intColor,size(vecBoundingRect1(1,:)));
    sDots.vecSpeed_deg(1) = sParams.dblSpeed_deg;
    sDots.vecSpeed_pix(1) = sParams.dblSpeed_pix;

    sDots.stimID(2) = 2; %outside
    sDots.cellBoundingRect{2} = vecBoundingRect2;
    if intSet==6,sDots.cellColor{2} = repmat(sParams.intColor,size(vecBoundingRect2(1,:)));
    elseif intSet==7,sDots.cellColor{2} = repmat(sParams.intBackground,size(vecBoundingRect2(1,:)));end
    sDots.vecSpeed_deg(2) = sParams.dblSpeed_deg;
    sDots.vecSpeed_pix(2) = sParams.dblSpeed_pix;

    sDots.intStimulusConditions = 2;
    sDots.vecBlockStructure = [1 repmat(2,[1 15])];
	
	sParams.dblSecsPostBlank = 0.5;
	
end
end

