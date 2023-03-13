function sRelatedRecord = get_related_record( record, strCrit, db )
%GET_RELATED_RECORD retrieves related record from experiment database
% 
% 2023, Alexander Heimel

sRelatedRecord = [];
if nargin<3 || isempty(db)
    h_db = get_fighandle('Neuropixels database*');
    if ~isempty(h_db)
        sUserData = get(h_db,'userdata');
        db = sUserData.db;
    else
        logmsg('Cannot retrieve database');
        return
    end
end

strCrit = ['dataset=' record.dataset ',subject=' record.subject ...
    ',date=' record.date ',' strCrit];
sRelatedRecord = db(find_record(db,strCrit));
end
