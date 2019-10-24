function redraw(movieobj,frame)

hf = gcf;
hf.Name = movieobj.Filename;

IM = movieobj.read(frame);
im = IM;
% imshow(imadjust(im,stretchlim(im,0.0001)),[])
% imshow(im,[])
imshow(imadjust(im),[])
set(gca,'units','normalized');
text(0,0,['frame ',num2str(frame)],'Color','r','Units','pixels','VerticalAlignment','bottom');

end