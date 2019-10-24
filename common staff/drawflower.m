clear all; close all; 

MatlabVersion = ver('Matlab');
numMatlabVersion = str2double(MatlabVersion.Version);

theta=[0:.01:2*pi];
r_vals=cos(5*theta);
% r_vals=cos(7*theta);
x_vals=r_vals.*cos(theta);
y_vals=r_vals.*sin(theta);

stem_y=[-2:.01:0]; %...and then generate a 'stem' for it.
stem_x=-.05.*stem_y.^2;

if numMatlabVersion >= 8.4
    
%     disp('Matlab r2014b or more recent detected.');
    
    figure
    hold on
    xlim([-1 1]);
    ylim([-2 1]);
    
    % stem
    h_stem = animatedline(stem_x(1),stem_y(1),'Color','g','LineWidth',2);
    for i=1:numel(stem_x)
        addpoints(h_stem,stem_x(i),stem_y(i));
        pause(0.01)
    end
    
    % centre
    hs = scatter(0,0,50,'y','filled');
    
    % flower
    h_rose = animatedline(x_vals(1), y_vals(1),'Color','r','LineWidth',2);
    uistack(hs, 'top');
    for i=1:numel(x_vals)
        addpoints(h_rose,x_vals(i),y_vals(i));
        pause(0.01)
    end
    
else
    
%     disp('Matlab r2014a or less recent detected.');
    
    figure
    shg
    hold on
    xlim([-1 1]);
    ylim([-2 1]);
    
    % stem
    h_stem = plot(stem_x(1),stem_y(1),'Color','g','LineWidth',2);
    for i=2:numel(stem_x)
        
        xd = get(h_stem,'xdata');
        yd = get(h_stem,'ydata');
        set(h_stem,'xdata',[xd(:) ; stem_x(i)],'ydata',[yd(:) ; stem_y(i)]) ;
        pause(0.01)
    end
    
    % centre
    hs = scatter(0,0,50,'y','filled'); %This creates a center for the flower.
    
    % flower
    h_rose = plot(x_vals(1), y_vals(1),'Color','r','LineWidth',2);
    uistack(hs, 'top');
    for i=1:numel(x_vals)
        
        xd = get(h_rose,'xdata');
        yd = get(h_rose,'ydata');
        set(h_rose,'xdata',[xd(:) ; x_vals(i)],'ydata',[yd(:) ; y_vals(i)]) ;
        pause(0.01)
    end
    
end


