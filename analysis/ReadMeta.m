function meta = ReadMeta(metaFilePath)
    % metaFilePath: full path to .meta file
    fid = fopen(metaFilePath, 'r');
    if fid < 0
        error('Cannot open meta file: %s', metaFilePath);
    end
    C = textscan(fid, '%[^=] = %[^\r\n]');
    fclose(fid);
    
    meta = struct();
    for i = 1:length(C{1})
        tag = strtrim(C{1}{i});
        if tag(1) == '~'
            tag = tag(2:end); % remove leading ~
        end
        meta.(tag) = strtrim(C{2}{i});
    end
end
