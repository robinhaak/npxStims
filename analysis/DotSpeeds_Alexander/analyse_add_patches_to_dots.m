function record = analyse_add_patches_to_dots( record, db, verbose)
%ANALYSE_ADD_PATCHES_TO_DOTS add results from grating patches to dot speeds record
%
%   record = analyse_add_patches_to_dots( record, db, verbose)
%
% 2023, Alexander Heimel

logmsg('DEPRECATED. SHOULD NOW BE DONE IN ANALYSE_MOVINGDOTS AND RESULTS_MOVINGDOTS')


if nargin<3 || isempty(verbose)
    verbose = true;
end
if nargin<2 || isempty(db)
    h_db = get_fighandle('Neuropixels database*');
    sUserData = get(h_db,'userdata');
    db = sUserData.db;
end

if ~strcmp(record.stimulus,'MovingDots')
    return
end


crit = ['dataset=' record.dataset ',subject=' record.subject ... 
    ',sessionid=' record.sessionid ',stimulus=GratingPatches'];
ind = find_record(db,crit);
if isempty(ind) && verbose
    logmsg(['No matching grating patches record for ' recordfilter(record)]);
    return
end
if length(ind)>1
    logmsg(['More than one matching grating patches record for ' recordfilter(record)]);
    logmsg('Using first record');
    ind = ind(1);
end
patchesrecord = db(ind);

measures = record.measures;
patchesmeasures = patchesrecord.measures;
patchesindex = [patchesmeasures.intIndex];
for i = 1:length(measures)
    ind = find(patchesindex == measures(i).intIndex);
    measures(i).dblXRFLeftFromGratingPatches_pix = patchesmeasures(ind).dblXRFLeft_pix;
    measures(i).dblXRFRightFromGratingPatches_pix = patchesmeasures(ind).dblXRFRight_pix;
    measures(i).dblPeakTimeFromGratingPatches = patchesmeasures(ind).dblPeakTime;
    measures(i).dblOnsetTimeFromGratingPatches = patchesmeasures(ind).dblOnsetTime;
end
record.measures = measures;
logmsg(['Add patches results to ' recordfilter(record)]);
