function [ filename, filepath, ext] = parse_filename( filename )
%parse_filename chops the filename in only two parts, fills in filepath if
%needed

[filepath, file, ext] = fileparts(filename);
filename = [file, ext];

if isempty(filepath) || strcmp(filepath,'.')
    filepath = pwd;
end

if strcmp(filepath,'..')
    hd = pwd;
    cd ..
    filepath = pwd;
    cd(hd);
end

end

