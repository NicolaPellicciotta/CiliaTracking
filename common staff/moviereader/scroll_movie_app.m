function [] = scroll_movie_app()

    [filename, pathfile] = uigetfile('*.movie');
    movieobj = moviereader(fullfile(pathfile,filename));
    videofig(movieobj.NumberOfFrames, @redraw, movieobj, round(movieobj.FrameRate));
    redraw(movieobj,1);

end