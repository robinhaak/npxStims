function  [sDots,sParams] = RH_CreateDotsFromAbove(intSet,sParams)
% generate stimulus trajectories for
% RH_runLocomotionFeedack/RH_runAsyncStim
% there are 2 available stim set:
%   (1) "one dot"
%   (2) "two_dots" (x degrees apart), one through RF
%
% history:
% December 2023
%   - created by Robin Haak

%% default params
sParams.intSize_pix = round(sParams.dblSize_deg*sParams.dblPixelsPerDeg);
if mod(sParams.intSize_pix,2),sParams.intSize_pix = sParams.intSize_pix-1;end %avoid decimal pix vals
sParams.dblColor = 0; %[0 1]
sParams.intColor = round(mean(sParams.dblColor)*255);
sParams.dblSpeed_pix = sParams.dblSpeed_deg*sParams.dblPixelsPerDeg; %pix/s
sParams.intSpeed_ppf = round(sParams.dblSpeed_pix/sParams.intStimFrameRate); %pix/frame
if sParams.intSpeed_ppf==0
	sParams.intSpeed_ppf =1;
end
sParams.dblDotSep_deg = 30; %deg, separation between two dots
sParams.intDotSep_pix = round(sParams.dblDotSep_deg*sParams.dblPixelsPerDeg); %pix
% intPixFromEdge = round(10*sParams.dblPixelsPerDeg); %pix, from screen edge

%% create stimulus trajectories
%optimize dot locations
intLoc1 = sParams.intRfPosX_pix;
intSize = sParams.intSize_pix;
intScreenWidth = sParams.intScreenWidth_pix;
intScreenHeight = sParams.intScreenHeight_pix;

%if possible, location 2 is in the same hemifield as loc. 1
%nb, these calculations are specific to screen orientation
intLoc2_pos = intLoc1+sParams.intDotSep_pix;
intLoc2_alt = intLoc1-sParams.intDotSep_pix;

%give possible locations
f = figure; hold on;
rectangle('Position',[0 0 intScreenWidth intScreenHeight],'FaceColor',[0.5 0.5 0.5]);
viscircles([intLoc1 intScreenHeight/2],intSize/2,'Color','b'); %loc1
viscircles([intLoc2_pos intScreenHeight/2; intLoc2_alt intScreenHeight/2],[intSize/2; intSize/2],'Color','r'); %pos. locs2
s = scatter(sParams.intRfPosX_pix,intScreenHeight/2,'kx');
legend(s,'RF center'); xlabel('pix'); ylabel('pix');
axis('auto xy');axis image;fixfig

%get UI
fprintf('Stimulus locations:\n - Loc1: %d\n - Loc2_pos1: %d\n - Loc2_pos2: %d\n',intLoc1,intLoc2_pos,intLoc2_alt);
fprintf('What is/are the final stim location(s)?\n')
if intSet==1
    %note that if you want to show a stimulus outside the estimated RF,]
    %you can fill out one of the Loc2 vals
    intLoc1 = input('Location 1: ');
elseif intSet==2
    intLoc1 = input('Location 1: ');
    intLoc2 = input('Location 2: ');
end
close(f);

%create stimulus trajectories
vecBoundingRect1 = []; %1
vecBoundingRect1(2,:) = -intSize:sParams.intSpeed_ppf:(intScreenHeight+sParams.intSpeed_ppf);
vecBoundingRect1(1,:) = repmat(intLoc1-intSize/2,size(vecBoundingRect1(2,:)));
vecBoundingRect1(3,:) = repmat(intLoc1+intSize/2,size(vecBoundingRect1(2,:)));
vecBoundingRect1(4,:) = vecBoundingRect1(2,:)+intSize;

%create sDots
sDots = struct;
if intSet==1,sDots.strStimSet = 'one_loc';
elseif intSet==2,sDots.strStimSet = 'two_loc';end

sDots.stimID(1) = 1;
sDots.cellBoundingRect{1} = vecBoundingRect1;
sDots.cellColor{1} = repmat(sParams.intColor,size(vecBoundingRect1(1,:)));
sDots.vecSpeed_deg(1) = sParams.dblSpeed_deg;
sDots.vecSpeed_pix(1) = sParams.dblSpeed_pix;

if intSet==2
    vecBoundingRect2 = []; %2
    vecBoundingRect2(2,:) = -intSize:sParams.intSpeed_ppf:(intScreenHeight+sParams.intSpeed_ppf);
    vecBoundingRect2(1,:) = repmat(intLoc2-intSize/2,size(vecBoundingRect1(2,:)));
    vecBoundingRect2(3,:) = repmat(intLoc2+intSize/2,size(vecBoundingRect1(2,:)));
    vecBoundingRect2(4,:) = vecBoundingRect1(2,:)+intSize;

    sDots.stimID(2) = 2;
    sDots.cellBoundingRect{2} = vecBoundingRect2;
    sDots.cellColor{2} = repmat(sParams.intColor,size(vecBoundingRect2(1,:)));
    sDots.vecSpeed_deg(2) = sParams.dblSpeed_deg;
    sDots.vecSpeed_pix(2) = sParams.dblSpeed_pix;
end
