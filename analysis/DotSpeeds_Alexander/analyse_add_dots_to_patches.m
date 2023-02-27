function record = analyse_add_dots_to_patches( record, db, verbose)
%ANALYSE_ADD_DOTS_TO_PATCHES add results from moving dots to grating patches record
%
%   record = analyse_add_dots_to_patches( record, db, verbose)
%
% 2023, Alexander Heimel

if nargin<3 || isempty(verbose)
    verbose = true;
end
if nargin<2 || isempty(db)
    h_db = get_fighandle('Neuropixels database*');
    sUserData = get(h_db,'userdata');
    db = sUserData.db;
end

if ~strcmp(record.stimulus,'GratingPatches')
    return
end


crit = ['dataset=' record.dataset ',subject=' record.subject ... 
    ',sessionid=' record.sessionid ',stimulus=MovingDots'];
ind = find_record(db,crit);
if isempty(ind) && verbose
    logmsg(['No matching dots record for ' recordfilter(record)]);
    return
end
if length(ind)>1
    logmsg(['More than one matching dots record for ' recordfilter(record)]);
    logmsg('Using first record');
    ind = ind(1);
end
dotsrecord = db(ind);

measures = record.measures;
dotsmeasures = dotsrecord.measures;
dotsindex = [dotsmeasures.intIndex];
for i = 1:length(measures)
    ind = find(dotsindex == measures(i).intIndex);
    measures(i).dblXRFLeftFromOnsetFromMovingDots_pix = dotsmeasures(ind).dblXRFLeftFromOnset_pix;
    measures(i).dblXRFRightFromOnsetFromMovingDots_pix = dotsmeasures(ind).dblXRFRightFromOnset_pix;
    measures(i).dblXRFLeftFromMovingDots_pix = dotsmeasures(ind).dblXRFLeft_pix;
    measures(i).dblXRFRightFromMovingDots_pix = dotsmeasures(ind).dblXRFRight_pix;
    measures(i).dblXRFRightFromMovingDots_pix = dotsmeasures(ind).dblXRFRight_pix;
    measures(i).dblDeltaTLeftFromMovingDots = dotsmeasures(ind).dblDeltaTLeft;
    measures(i).dblDeltaTRightFromMovingDots = dotsmeasures(ind).dblDeltaTRight;
    measures(i).dblDeltaTLeftFromOnsetFromMovingDots = dotsmeasures(ind).dblDeltaTLeftFromOnset;
    measures(i).dblDeltaTRightFromOnsetFromMovingDots = dotsmeasures(ind).dblDeltaTRightFromOnset;
end
record.measures = measures;
logmsg(['Add dots results to ' recordfilter(record)]);
