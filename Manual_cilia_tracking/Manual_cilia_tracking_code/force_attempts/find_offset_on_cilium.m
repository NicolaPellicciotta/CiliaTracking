function [ offset, cilium_to_offset ] = find_offset_on_cilium( cil, flag_plot )
% find_offset_on_cilium finds how to register couples of cilia so that they better overlap.
%   Two different cilia can have been clicked starting from a different
%   base point. Here I try and chop off a tiny bit from the beginning of
%   each to see if it improves the correspondence between points spaced 1um
%   along the cilia as taken on the two different cilia

%% input check

if nargin < 2 || isempty(flag_plot)
    flag_plot = false;
end

if numel(cil) ~= 2
    error('cil hase to have only 2 elements');
end


%% constants and initialisations

best_offset = [1; 1];
shortest_distance = zeros(2,1); % this is the (shortest, as a function of offset) average distance that separates a couple of (integer um) corresponding points on the two cilia

max_of = floor(numel(cil(1).tt)/10); %assuming not more than a 10% error on clicking beginning

% only look at interdistances within the first 3 um, otherwise the tip of 
% the cilium weighs too much and the tangential sliding becomes severely 
% underestimated
arclength_cilium_for_offset_um = 3; 


%% find best offset

for i=1:2 % movie offset on cilium i, keep j still still
    
    % when i=1, j=2 and viceversa
    if i == 1, j = 2; elseif i == 2,  j = 1; end
    
    % initialise average distance array
    adtvo = nan(max_of ,1);
    
    % plot stuff
    if flag_plot
        hf(i) = figure;
        clf;
        
        ha2 = axes('Position',[0.55 0.1 0.4 0.8]);
        hold on;
        box on;
        ha2.XLim = cil(i).tt_um([1,end/5]);
        xlabel(ha2,'offset');
        ylabel(ha2,'av distance')
        
        ha1 = axes('Position',[0.05 0.1 0.4 0.8]);
        hold on
        box on
        plot(ha1,cil(i).xx, cil(i).yy);
        plot(ha1,cil(j).xx, cil(j).yy);
        plot(ha1,cil(j).xr, cil(j).yr, '.r');
        axis image
        hold on
    end
    
    % change offset
    for of = 1:max_of
        
        % find the maximum amount of (integer um) points that would fit on
        % both the portion of cilium i and cilium j
        min_arclength_um = min( diff(cil(i).tt_um([of,end])), cil(j).tt_um(end) );
        
        % selcts the points to "match" (to find the distances with). use
        % only the first 3 um (to prevent tip to weigh too much)
        um_points_to_find = (0:1:min(arclength_cilium_for_offset_um,min_arclength_um))';
        
        % apply offset on cilium 1
        temp_tt_um = cil(i).tt_um - cil(i).tt_um(of);
        
        % find points at integer um distances from offset, and their
        % coordinates
        temp_idx_umm = knnsearch(temp_tt_um, um_points_to_find);
        temp_tr_um = temp_tt_um(temp_idx_umm);
        temp_xr = cil(i).xx(temp_idx_umm);
        temp_yr = cil(i).yy(temp_idx_umm);
        
        
        % take the first few points on cil(j).*r so they match in number the
        % points on cilium 1
        
        adtvo(of) = mean(hypot(temp_xr - cil(j).xr(1:numel(temp_tr_um)),...
            temp_yr - cil(j).yr(1:numel(temp_tr_um))));
        
        
        % plot bits
        if flag_plot
            if of > 1
                delete(hpof);
                delete(hpumm);
                delete(hpjp);
                delete(had);
            end
            
            hpumm = plot(ha1,temp_xr, temp_yr,'b.');
            hpof = plot(ha1,cil(i).xx(of), cil(i).yy(of),'g.');
            
            for ii=1:numel(temp_tr_um)
                hpjp(ii) = plot(ha1,[temp_xr(ii),cil(j).xr(ii)],[temp_yr(ii),cil(j).yr(ii)],'g');
            end
            
            had = plot(ha2,cil(i).tt_um(1:max_of),adtvo,'.b');
            
            drawnow
            pause(0.01);
        end
        
    end
    
    % find the shortest distance and corresponding offset for this ciilum
    [shortest_distance(i), best_offset(i)] = min(adtvo);
    
end

[~, cilium_to_offset] = min(shortest_distance);
offset = best_offset(cilium_to_offset);

% delete figures when user is ready
if flag_plot
    pause;
    delete(hf);
end

end

