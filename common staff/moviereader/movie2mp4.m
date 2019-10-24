function [] = movie2mp4()
%movie2mp4 Interactive program to save a .movie file to an .mp4 video

[filename, pathname] = uigetfile('*.movie');

mo = moviereader(fullfile(pathname,filename));

choice = input(['This movie is ',num2str(mo.NumberOfFrames),' frames long. Do you want to convert all the frames? If yes, press [1], otherwise press any other key...']);

switch choice
    case{1}
        minFrame = 1;
        maxFrame = mo.NumberOfFrames;
    otherwise
        minFrame = input('Type first frame to convert: ');
        maxFrame = input('Type last frame to convert: ');
end

outputName = fullfile(pathname,filename);
outputName(end-5:end) = [];
outputName = [outputName,'.mp4'];

wo = VideoWriter(outputName,'MPEG-4');

if mo.FrameRate <= 100
    wo.FrameRate = mo.FrameRate;
else
    wo.FrameRate = 100;
    warning('Source Framerate too high!! set to 100fps in output.');
end

if mo.height > 1088 || mo.width > 1920
    warning('Max resolution is 1088x1920, the movie will be cropped');
end

out_height = mo.height;
out_width = mo.width;

quality = input('Type desired quality of video (higher quality equals bigger filesize). Has to be a number between 1 and 100: ');

quality = abs(round(quality));

if quality > 100
    warning('Max quality is 100, setting it now!');
    quality = 100;
end

wo.Quality = quality;

open(wo);

map = colormap(gray(256));
close(gcf);

choice = input('To rescale contrast globally press [1] (extremely slow for big movies), otherwise press any other key (video might flicker a bit but better contrast and faster)');

switch choice
    case {1}
        
        try
            fs = mo.read([minFrame, maxFrame]);
        catch
            warning('frames out of bound, loading the entire video');
            fs = mo.read;
        end
        
        fs = fs(1:out_height,1:out_width,:);
        
        sz = size(fs);
        fs = reshape(fs,[sz(1)*sz(3), sz(2)]);
        fs = imadjust(fs);
        fs = reshape(fs,sz);
        
        
        for f = 1:sz(3);
            IM = gray2ind(fs(1:out_height,1:out_width,f),256);
            frame = im2frame(flipud(IM),map);
            wo.writeVideo(frame);
        end
        
    otherwise
        
        
        for f = minFrame:maxFrame
            
            try
                IM = mo.read(f);
            catch
                return
            end
            
            IM2 = imadjust(IM);
            IM3 =  gray2ind(IM2,256);
            frame = im2frame(flipud(IM3),map);
            wo.writeVideo(frame);
            
        end
end

close(wo);

end
