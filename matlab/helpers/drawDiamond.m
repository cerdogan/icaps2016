%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% @name Circle.m
%% @author Can Erdogan
%% @date Nov 15, 2012
%% @brief Draws the diamond (l1 norm) that would fit in a circle 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = drawDiamond (pos, radius, color) 

    if(nargin < 3), color = 'k'; end;
    
    p1 = pos + [radius, 0];
    p2 = pos + [0, radius];
    p3 = pos + [-radius, 0];
    p4 = pos + [0, -radius];
    ps = [p1;p2;p3;p4;p1];
    plot(ps(:,1), ps(:,2), ['',color]); hold on;
end
