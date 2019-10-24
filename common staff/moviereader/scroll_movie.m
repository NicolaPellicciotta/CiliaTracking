function [] = scroll_movie(filename)

if nargin < 1
    [filename, pathfile] = uigetfile('*.movie');
    filename = fullfile(pathfile,filename);
end
movieobj = moviereader(filename);
videofig(movieobj.NumberOfFrames, @redraw, movieobj, round(movieobj.FrameRate));
redraw(movieobj,1);

end