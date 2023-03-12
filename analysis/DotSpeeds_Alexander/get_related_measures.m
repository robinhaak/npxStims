function sRelatedMeasures = get_related_measures( record, strCrit, db, vecIndex)
%GET_RELATED_MEASURES retrieves measures from same units from related record
%
% 2023, Alexander Heimel

if nargin<3 || isempty(db)
    db = [];
end
if nargin<4 || isempty(vecIndex)
    vecIndex = [];
end

sRelatedMeasures = [];
sRelatedRecord = get_related_record( record, strCrit, db );
if ~isempty(sRelatedRecord)
    sRelatedMeasures = sRelatedRecord.measures;
    if ~isempty(vecIndex)
        sRelatedMeasures = sRelatedMeasures(ismember([sRelatedMeasures.intIndex],vecIndex));
    end
end
end