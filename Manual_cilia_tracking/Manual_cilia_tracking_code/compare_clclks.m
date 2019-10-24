function [duplicate_matrix, duplicate_array] = compare_clclks(clclk)
%compare_clclks compares clclks by plotting the raw clicking points
%returns a duplicate_matrix in which entry (i,j) is true if clclk(i) ==
%clclk(j), false otherwise

figure;
hold on

cmap = lines(numel(clclk));
markarr = '^odsvp><h*+x.';

for cc = 1:numel(clclk)
    
    % find clicked frames
    ind_clicked_frames = [clclk(cc).clicked_frames.frame_number];
    
    % concatenate all x and y data
    x_plot = cat(1, clclk(cc).points(ind_clicked_frames).cilium_x);
    y_plot = cat(1, clclk(cc).points(ind_clicked_frames).cilium_y);
    
    % plot
    hp(cc) = plot(x_plot, y_plot, markarr(cc),...
        'color',cmap(cc,:),...
        'markersize',6);
    hp(cc).DisplayName = [num2str(cc),' --- ',clclk(cc).savename] ;
    
    % save quantities for automatic recognition of duplicate
    S(cc).ind_clicked_frames = ind_clicked_frames;
    S(cc).x_plot = x_plot;
    S(cc).y_plot = y_plot;
    S(cc).savename = clclk(cc).savename;
    
end % for


% for loop that identifies if any clclk is a duplicate
duplicate_matrix = false(numel(S));
duplicate_array = zeros(numel(S),1);

for c = 1:numel(S)
    for r = 1:numel(S)
        
        if numel(S(c).ind_clicked_frames) == numel(S(r).ind_clicked_frames)
            if numel(S(c).x_plot) == numel(S(r).x_plot) && numel(S(c).y_plot) == numel(S(r).y_plot)
                
                if all(S(c).x_plot == S(r).x_plot) && all(S(c).y_plot == S(r).y_plot)
                    duplicate_matrix(r,c) = true;
                    if r > c
                        disp(['File ', num2str(c) ,' == file ', num2str(r)])
                        duplicate_array(c) = r;
                    end
                end
            end
        end
        
    end % for r
end % for c


end