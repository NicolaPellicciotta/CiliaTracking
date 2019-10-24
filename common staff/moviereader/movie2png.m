function [] = movie2png()
%movie2png Interactive program to save a .movie file to a series of pngs

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

outputDirName = fullfile(pathname,filename);
outputDirName(end-5:end) = [];
outputFileName = filename(1:end-6);

choice = input('To rescale contrast globally press [1] (extremely slow for big movies), otherwise press any other key (video might flicker a bit but better contrast and faster)');

switch choice
    case {1}
        
        try
            fs = mo.read([minFrame, maxFrame]);
        catch
            warning('frames out of bound, loading the entire video');
            fs = mo.read;
        end
        
        sz = size(fs);
        fs = reshape(fs,[sz(1)*sz(3), sz(2)]);
        fs = imadjust(fs);
        fs = reshape(fs,sz);
        
        mkdir(outputDirName);
        
        for f = 1:sz(3);
            filename = [ outputFileName, sprintf('_%.6d.png',f+minFrame-1) ];
            imwrite(flipud(fs(:,:,f)), fullfile(outputDirName,filename), 'png');
        end
        
    otherwise
        
        mkdir(outputDirName);
        
        for f = minFrame:maxFrame
            
            try
                IM = mo.read(f);
            catch
                return
            end
            
            IM2 = imadjust(IM);
            filename = [ outputFileName, sprintf('_%.6d.png', f) ];
            imwrite(flipud(IM2), fullfile(outputDirName,filename), 'png');
            
        end
end

end
