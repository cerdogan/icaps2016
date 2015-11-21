function [res, counter, coll, dist, p3] = genericTest (p1, p2, rLim, distFunc, solvePointFunc, ...
        drawFunc, checkCollisions, drawScene, drawOpts, counter)

    global p0 p4 r0 r1 r2 r3 r4 op0;
    coll = false;
    p3 = [0,0];
    res = false;
    
    % Check if the two gears can be connected
    dist = distFunc(p4, p2);
    if(dist > rLim), return; end;
      
    % Compute collision distances
    distObs1 = distFunc(op0, p1);
    distObs2 = distFunc(op0, p2);
    dist02 = distFunc(p0, p2);
    dist24 = distFunc(p4, p2);
    
    % Compute the location of new gear and compute collision for it
    p3 = solvePointFunc(p2,p4,r2+r3,r3+r4);
    distObs3 = distFunc(op0, p3);
    
    % Check for collisions
		assert(false); % Is not this always false if I don't want to check for collisions
    coll = checkCollisions && (distObs1 > (r1+0.1)) && ...
      (distObs2 > (r2+0.1)) && (distObs3 > (r3+0.1)) && ...
      (dist02 > r2+r0) && (dist24 > r2+r4);
    if(coll == false), return; end;
      
    % Visualize the result
    if(drawScene == 0)
        figure(1);
        subplot(2,1,drawOpts{1}); 
        plot(p1(1),p2(1),['o',drawOpts{2}],'MarkerSize',drawOpts{3}); 
        hold on; 
    elseif(1 && drawScene == 1)
        clf;
        drawFunc(op0, 0.1, 'k');
        drawFunc(p1, r1, 'r');
        drawFunc(p2, r2, 'r');
        drawFunc(p3, r3, 'r');
        drawFunc(p0, r0, 'k');
        drawFunc(p4, r4, 'k');
         axis([-4, 12, -4, 10]);  
        %return;
        %pause;
    end
    
    % Increment the counter
    counter = counter + 1;
    res = true;
end
