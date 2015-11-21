%%
% @file simplicial.m
% @author Can Erdogan
% @date 2015-10-10
% @brief Implements simplicial approximation.

% ------------------------------------------------------------------------
function [bPs, radii] = simplicial () 

  % Generate the symbolic constraints
  syms x1 x2
  global xs ineqs lims
  xs = [x1,x2];
  [ineqs, lims] = example2;
  global colors;
  colors = ['r', 'g', 'c', 'y'];
  global hConvex hCircle hMids hExtend
  hMids = []; hConvex = []; hCircle = []; hExtend = [];

  % Plot the constraints
  clf
  x = lims(1):0.1:lims(2);
  for i = 1 : numel(ineqs)
    ineq = ineqs(i);
    vals = coeffs(ineq,x2);
    plot(x,subs(-vals(1)/vals(2),x)); hold on;
  end
  axis(lims); axis equal
  pause(0.01);
  
  % -----------------------------------------------------------------------
  % Generate minimal simplex
  
  % Setup the functions for faster evaluation
  global ineqFs;
  ineqFs = {};
  for i = 1 : numel(ineqs), ineqFs{i} = matlabFunction(ineqs(i)); end;
  
  % Find a feasible sample 
  sourcePoint = [(lims(2) - lims(1))*rand + lims(1), (lims(4) - lims(3))*rand + lims(3)]';
  fun = @(x) norm(sourcePoint - x);
%   solution = fmincon(fun, sourcePoint, [], [], [], [], [-10, -60], [10 60], @nonlcons);
  solution = [9.395, -20.62]';
  plot(solution(1), solution(2), 'ko', 'MarkerSize', 7, 'LineWidth', 4);
  pause(0.01);
  
  % Choose a random line, find intersections with other inequalities,
  % find the closest one to the solution which we assume is in the
  % feasible space, collect three such points.
  bPs = [];
  dimensions = 2;
  while(true) 

    % Check if collected sufficient points
    if(size(bPs,1) > dimensions), break; end;
    
    % Try to find a new point to add to the boundary
    [success, val] = tryNewBoundaryPoint (solution, bPs, ineqs);
    if(not(success)), continue; end;
    
    % Add the point
    bPs(end+1,:) = val(:);
    fprintf('Added new point to initial simplex\n');
    pause
    %pause(0.01);
  end

  k = convhull(bPs(:,1), bPs(:,2));
  plot(bPs(k,1),bPs(k,2),['r-'],bPs(:,1), bPs(:,2),'b*','LineWidth', 2);
  
  fprintf('Initial hull: \n');
  bPs 
  return;
  
  % -----------------------------------------------------------------------
  % Extend the convex hull
  
  boundaryDists = [];
  hulls = {};
  global HULL;
  for HULL = 1 : 1
    
    % Create the secondary hulls
    if(HULL > 1)
      
      bPs
      
      % Compute the edge of the previous hull that has the largest distance
      % to the boundary
      % Alternatively, we can choose the edge with the largest hypersphere
      maxDist = -1;
      bestEdgeIdx = -1;
      maxMinClosestPoint = [];
      for i = 1 : (size(k,1)-1)
        
        % Generate the intersecting line from the midpoint to outer surface
        p1 = bPs(k(i),:); p2 = bPs(k(i+1),:);
        v21 = (p2 - p1) / norm(p2 - p1);
        vPerp = -[-v21(2), v21(1)];
        pm = (p1 + p2) / 2.0;
        pm2 = pm' + 100 * vPerp';
        
        % Find the intersection with each inequality
        minDist = 1000000;
        minClosestPoint = [];
        for ineqIdx = 1 : numel(ineqs)
          
          [res, val] = lineSegmentIntersect(ineqs(ineqIdx), pm', pm2);
          if(not(res)), fprintf('No solution\n'); continue; end;
          dist = norm(val-pm);
          if(dist < minDist)
            minDist = dist;
            minClosestPoint = val; 
          end
        end
        assert(minDist ~= 1000000);
        
        % Save the minimum distance
        boundaryDists(end+1,:) = [minDist, HULL-1, i, minClosestPoint];
        
        % Update the edge with the largest distance
        if(minDist > maxDist)
          maxDist = minDist;
          maxMinClosestPoint = minClosestPoint;
          bestEdgeIdx = i;
        end
      end
      
      % Choose the best edge
      boundaryDists = sortrows(boundaryDists, -1);
      bestEdgeIdx = boundaryDists(1,3);
      maxMinClosestPoint = boundaryDists(1,4:5);
      temp = hulls{boundaryDists(1,2)};
      boundaryDists(1,:) = [];
      prevBPs = temp{1};
      prevK = temp{2};
      
      % Remove the inequality of previous hull if necessary
      if(HULL > 2), ineqs(end) = []; end;
      
      % Add the other half-space of the chosen inequality 
      p1 = prevBPs(prevK(bestEdgeIdx),:); p2 = prevBPs(prevK(bestEdgeIdx+1),:);
      m = (p2(2) - p1(2)) / (p2(1) - p1(1));
      b = p2(2) - m * p2(1);
      newIneq = m * x1 + b - x2;
      meanPrevBPs = mean(prevBPs);
      temp = subs(newIneq, [x1,x2], meanPrevBPs);
      if(not(isa(temp, 'double'))), temp = eval(temp); end;
      if(temp < 0)
        fprintf('Flipped new inequality\n');
        newIneq = -newIneq;
      end;
      ineqs(end+1) = newIneq;
      
      % Initialize boundary points with the vertices of this edge
      fprintf('Found the new boundary point\n');
      bPs = [p1; p2; maxMinClosestPoint]
      %pause
      plot(bPs(:,1), bPs(:,2), 'rs'); hold on;
      k = convhull(bPs(:,1), bPs(:,2));
      colors = ['r', 'g', 'c', 'y'];
      hullColor = colors(mod(HULL-1,4)+1);
      plot(bPs(k,1),bPs(k,2),[hullColor,'-'],bPs(:,1), bPs(:,2),'b*','LineWidth', 1 + HULL);
      %pause
    end
    
    % Add points iteratively
    radii = [];
    for ITERS = 1 : 5
      pause
      [bPs, k, radius] = addPoint(bPs, ineqs);
      radii(end+1) = radius;
      if(numel(radii) > 1 && abs(radii(end-1) - radii(end)) < 1e-3), break; end;
    end
    cleanUp;
    hulls{HULL} = {bPs,k};
    hullColor = colors(mod(HULL-1,4)+1);
    plot(bPs(k,1),bPs(k,2),[hullColor,'-'],bPs(:,1), bPs(:,2),'b*','LineWidth', 1 + HULL);
  end
end

% ------------------------------------------------------------------------
function [] = plotConstraint (ineqs, index, color)  
  global xs;
  x2 = xs(2);
  ineq = ineqs(index);
  x = -10:0.1:10;
  vals = coeffs(ineq,x2);
  plot(x,subs(-vals(1)/vals(2),x), ['',color]); hold on; 
end

% ------------------------------------------------------------------------
function [success, val] = tryNewBoundaryPoint (solution, bPs, ineqs)

  success = false;
  global xs lims;
  x1 = xs(1); x2 = xs(2);
  
  % Find random line
  temp = [(lims(2) - lims(1))*rand + lims(1), (lims(4) - lims(3))*rand + lims(3)]';
%   temp = [-8,-16]';
  randDir = (temp - solution);
  randDir = randDir / norm(randDir);
  sample = solution + 60 * randDir;
  h = plot([sample(1); solution(1)], [sample(2); solution(2)], '-o');

  % Intersect the random line with all the constraint and pick closest
  minDist = 1000000;
  minClosestPoint = [];
  for ineqIdx = 1 : numel(ineqs)
    [res, val] = lineSegmentIntersect(ineqs(ineqIdx), solution, sample);
    if(not(res)), continue; end;
    dist = norm(val'-solution);
    if(dist < minDist)
      minDist = dist;
      minClosestPoint = val'; 
    end
  end
  assert(minDist ~= 1000000);
  val = minClosestPoint;
  
  % Check if the point is within the feasible region
  for j = 1 : numel(ineqs)
    value = subs(ineqs(j), [xs(1),xs(2)], [val(1), val(2)]);
    if(value > 0)
      fprintf('Neighbor segment goes directly into infeasible space\n');
      delete(h);
      return
    end
  end

  % Check if the point is too close to the others
  if(size(bPs,1) > 0)
    closestDist = min(sqrt(sum(abs(bPs - repmat(val',size(bPs,1),1)).^2,2)));
    if(closestDist < 0.5)
      fprintf('Too close to others: %f\n', closestDist);
      delete(h);
      return; 
    end
  end;

  % Check if the intersection can be connected with a line segment with the 
  % points in the convex hull without getting out of the feasible space
  for pIdx = 1 : size(bPs,1)
    
    % Check if the ray from the solution to the hull point moves
    % towards the unfeasible space
    % TODO Can be sped up by keeping track of which equation to point 
    % belongs to
    p = transpose(bPs(pIdx,:));
    v = (p - val) / norm(p - val);
    startPoint = val + v * 1e-2;
    for j = 1 : numel(ineqs)
      value = subs(ineqs(j), [xs(1),xs(2)], [startPoint(1), startPoint(2)]);
      if(value > 0)
        fprintf('Neighbor segment goes directly into infeasible space');
        delete(h);
        return
      end
    end
        
    % Check if the line to the hull point intersects the inequalities
    endPoint = p - v * 1e-2;
    for j = 1 : numel(ineqs)
      [res, ~] = lineSegmentIntersect(ineqs(j), startPoint, endPoint);
      if(res)
        fprintf('A part of the convex hull is in the infeasible region\n');
        delete(h);
        return; 
      else
        fprintf('\tInequality %d is clear: %d\n', j, res);
      end;
    end
  end

  plot(val(1), val(2), 'cs', 'MarkerSize', 2, 'LineWidth', 3); hold on;

%   keyboard
  success = true;
end

% ------------------------------------------------------------------------
function [success, val] = lineSegmentIntersect (mainEq, solution, sample)

  % Setup
  success = false;
  val = [];
  global xs;
  x1 = xs(1); x2 = xs(2);
  
  % Compute the line's values
  m = (sample(2) - solution(2)) / (sample(1) - solution(1));
  b = sample(2) - m * sample(1);
  lineEq = m * x1 - x2 + b;
  
  % Get the intersection points
  val = coeffs(lineEq, x2);
  temp1 = subs(mainEq, x2, -val(1)/val(2));
  sols_x2 = eval(solve(temp1));
  sols_y2 = (subs(-val(1)/val(2), x1, sols_x2));
  if(not(isa(sols_y2,'double'))), sols_y2 = eval(sols_y2); end;
  ps = [sols_x2, sols_y2];
  %plot(ps(:,1), ps(:,2), '-go'); hold on;

  % See if one of the intersections is in between the two points
  v = -(solution - sample) / norm(solution - sample);
  projs = (ps - repmat(solution', size(ps,1), 1)) * v;
  dist = norm(sample' - solution');
  res = abs(imag(projs)) < 1e-3 & real(projs) > 1e-3 & real(projs) < dist+1e-3;
  %if(sum(res) > 2), error('may need to choose the minimal distance'); end;
  if(sum(res) > 0)
    minProj = min(projs(res));
    minVal = min(projs(res));
    val = ps(abs(projs  - minVal) < 1e-3,:);
    plot(val(1), val(2), 'ro'); hold on;
%     if(res(1) == 1), val = (ps(1,:)); 
%     else val = (ps(2,:)); end;
  else
   % keyboard
    fprintf('Point not in between the two points\n');
    return;
  end
  
  assert(isreal(val));
  success = true;
end

% ------------------------------------------------------------------------
function [] = cleanUp ()
  global hConvex hCircle hMids hExtend
  delete(hMids);
  delete(hCircle);
  delete(hExtend);
  delete(hConvex);
  hMids = []; hConvex = []; hCircle = []; hExtend = [];  
end

% ------------------------------------------------------------------------
function [bPs, k, radius] = addPoint (bPs, ineqs)

  global xs lims
  x1 = xs(1); x2 = xs(2);
  global hConvex hCircle hMids hExtend
  cleanUp();
  global HULL colors;
  hullColor = colors(mod(HULL-1,4)+1);

  % B. Generate the convex hull
  hConvex = plot(bPs(:,1), bPs(:,2), 'rs'); hold on;
  k = convhull(bPs(:,1), bPs(:,2));
  hConvex=cat(1,hConvex, plot(bPs(k,1),bPs(k,2),[hullColor,'-'],bPs(:,1), bPs(:,2),'b*','LineWidth', 1 + HULL));
  pause(0.01);
  
  % C. Inscribe the largest hypersphere
  % To do so, solve a linear problem where the distances to each line
  % segment is some radius r. 

  % Generate the inequalities for the convex hull
  A = [];
  b = [];
  for i = 1 : (size(k,1) - 1)
    p1 = bPs(k(i),:); p2 = bPs(k(i+1),:);
    v21 = (p2 - p1) / norm(p2 - p1);
    hi = -[-v21(2), v21(1)];
    ci = dot(hi,p1);
    A(end+1,:) = hi;
    b(end+1,:) = ci;
  end
  if(dot(mean(bPs), A(2,:)) > b(2)), A = -A; b = -b; end;
  
  % Solve for the center of the circle
  A(:,end+1) = -ones(size(A,1),1);
  options=optimset('Display','none');
  [res, ~, exitflag] = linprog([0,0,1], A, b,[],[],[-10,-60,-60],[10,60,60], [], options);
  assert(exitflag == 1);
  center = res(1:2)'; 
  radius = -res(3);
  hCircle = drawCircle(center, radius);
  %axis([-10 10 -60 60]); axis square
  axis(lims); axis equal
  %keyboard
  
  % Determine the tangent "edge"s
  tangents = abs(A*res-b) < 1e-4;
  
  % Find the largest "edge"
  A(:,end) = [];
  A = A';
  A(end+1,:) = zeros(1,size(A,2));
  n = 2;
  A(:,end+1) = zeros(n+1,1);
%   largestEdge = -1;
%   largestEdgeRadius = -1;
  radii = [];
  for i = 1 : (size(k,1) - 1)

    % Check if the edge is tangent
    if(tangents((i)) ~= 1), continue; end;
    
    % Generate the equality constraint
    for j = 1 : (size(k,1) - 1)
      if(i == j)
        A(end,(j)) = 0;
        continue;
      end
      h1 = A(1:n,(i));
      h2 = A(1:n,(j));
      costh = dot(h1,h2) / (norm(h1) * norm(h2));
      th = acos(costh);
      A(end,(j)) = sin(th);
    end
    A(:,end) = [-A(1:n,(i)); 0];
    c = zeros(n+1,1);
    c(end) = 1;

    % Solve for the largest radius that can be fit
    Ap = A;
    A = real(Ap);
    coeffs_ = [b;-b((i))];
    [vals, ~, res] = linprog(coeffs_, [], [], A, c, zeros(size(k,1),1), [], [], options);
    assert(res == 1); 
    rad = dot(vals, coeffs_);
    radii(end+1, :) = [rad, i]; %#ok<AGROW>
    assert(rad > 0);

    % Find the largest radius
%     if(rad > largestEdgeRadius)
%       largestEdgeRadius = rad;
%       largestEdge = i;
%     end

    % Plot for fun
    p1 = bPs(k(i),:); p2 = bPs(k(i+1),:);
    v21 = (p2 - p1) / norm(p2 - p1);
    pc = p1 + v21 * rad;
    assert(abs((norm(p2-p1) / rad) - 2) < 1e-4);
    hMids(end+1) = plot(pc(1), pc(2), 'ok', 'MarkerSize', 8, 'LineWidth', 2); hold on; %#ok<AGROW>
    %keyboard
  end

  % Choose the largest one that is not on the boundary
  largestEdge = -1;
  largestEdgeRadius = -1;
  sortedRadii = sortrows(radii, -1);
  for idx = 1 : size(sortedRadii, 1)
    
    % Generate the center
    i = sortedRadii(idx,2);
    p1 = bPs(k(i),:); p2 = bPs(k(i+1), :);
    v21 = (p2 - p1) / norm(p2 - p1);
    pc = p1 + v21 * sortedRadii(idx,1);
    
    % Check if the center of the largest vertex is on the boundary 
    % If so, moving out would violate so make that the new point.
    onBoundary = false;
    for eqIdx = 1 : numel(ineqs)
      value = subs(ineqs(eqIdx), [xs(1),xs(2)], [pc(1), pc(2)]);
      if(abs(value) < 1e-3)
        fprintf('On boundary of equation %d\n', eqIdx);
        onBoundary = true;
        break;
      end
    end  
    
    % Set the largest radii
    if(not(onBoundary))
      largestEdge = i;
      largestEdgeRadius = sortedRadii(idx, 1);
      break;
    end
  end
  assert(largestEdge ~= -1);
  
  % Plot the largest one
  p1 = bPs(k(largestEdge),:); p2 = bPs(k(largestEdge+1),:);
  v21 = (p2 - p1) / norm(p2 - p1);
  pc = p1 + v21 * largestEdgeRadius;
  hExtend(end+1) = plot(pc(1), pc(2), 'mo', 'MarkerSize', 6, 'LineWidth', 2); hold on;
  %pause
  
  % Move out of the center of the largest face 
  vPerp = [-v21(2), v21(1)];
  if(dot(mean(bPs)-pc, vPerp) > 0), vPerp = -vPerp; end;
  pOut = pc + vPerp * 60;
  hExtend(end+1) = plot([pOut(1); pc(1)], [pOut(2); pc(2)], 'c-', 'MarkerSize', 6, 'LineWidth', 2); hold on;

  % Find the new point
  m = (pOut(2) - pc(2)) / (pOut(1) - pc(1));
  b = pOut(2) - m * pOut(1);
  lineEq = m * x1 - x2 + b;
  minDist = 1000000;
  for i = 1 : numel(ineqs)
    mainEq = ineqs(i);
    noInters = false;
    if(false)
      counter = 0;
      while(true)
        tic
        temp = subs(lineEq, x1, val(1)/val(2))
        [val, ~, flag] = fsolve(@lineSegmentIntersection, [20 * (rand-0.5); 120 * (rand-0.5)]', options);
        toc
        if(flag == 1), fprintf('Got solution!\n'); break; end;
        counter = counter+1;
        if(counter > 3), noInters = true; break; end;
      end
    else 
      
      % Get the intersection points
      val = coeffs(lineEq, x2);
      temp1 = subs(mainEq, x2, -val(1)/val(2));
      sols_x2 = eval(solve(temp1));
      sols_y2 = (subs(-val(1)/val(2), x1, sols_x2));
      if(not(isa(sols_y2,'double'))), sols_y2 = eval(sols_y2); end;
      ps = [sols_x2, sols_y2];
      %plot(ps(:,1), ps(:,2), '-go'); hold on;
      
      % See if one of the intersections is in between the two points
      v = (pOut - pc) / norm(pOut - pc);
      projs = (ps - repmat(pc, size(ps,1), 1)) * v';
      dist = norm(pc - pOut);
      res = projs > -1e-3 & projs < dist+1e-3;
      if(sum(res) > 0)
        if(res(1) == 1), val = ps(1,:); 
        else val = ps(2,:); end;
      else
        noInters = true;
      end
    end
    val = val';
%     counter = 0;
%     while(true)
%       [val, ~, flag] = fsolve(@lineSegmentIntersection, [20 * (rand-0.5), 120 * (rand-0.5)]', options);
%       if(flag == 1), break; end;
%       counter = counter+1;
%       if(counter > 3), noInters = true; break; end;
%     end
    if(noInters), continue; end;
    hExtend(end+1) = plot(val(1), val(2), 'r*'); hold on; %#ok<AGROW>
    dist = norm(val - pc');
    if(dist > 1e-4 && dist < minDist)
      minDist = dist;
      newP = val;
    end
  end
  
  % Check if the new point is going to cause three points to be linear
  % in the convex hull. If so, remove the one in the middle.
  for i = 1 : (size(k,1) - 1)
    
    % Check if the three points are colinear
    p1 = bPs(k(i),:); p2 = bPs(k(i+1),:);
    v21 = (p2 - p1) / norm(p2 - p1);
    perpV = -[-v21(2), v21(1)];
    perpDist = dot((newP' - p1), perpV);
    if(abs(perpDist) > 1e-4), continue; end;
      
    warning('May need to update k');
    % Make sure P1x < P2x.
    flipped = false;
    if(p1(1) > p2(1))
      flipped = true;
      temp = p1; p1 = p2; p2 = temp;
    end
      
    % Three cases: PnX < P1X < P2X, P1X < PnX < P2X, or P1X < P2X < PnX
    if(newP(1) < p1(1))
      if(flipped), bPs(k(i+1), :) = []; 
      else bPs(k(i), :) = []; end;
      break;
    elseif(newP(1) > p2(1))
      if(flipped), bPs(k(i), :) = []; 
      else bPs(k(i+1), :) = []; end;
      break;
    else
      assert(false && 'This case should have been not be possible due to convexity'); 
    end
    
  end
  
  % Add the new point
  bPs(end+1,:) = newP;
  
  % B. Generate the convex hull
  axis(lims); 
  k = convhull(bPs(:,1), bPs(:,2));
  %'Done'
  
end

% ------------------------------------------------------------------------
function [c, ceq] = nonlcons (x)
  ceq = [];
  global ineqFs;
  c = zeros(numel(ineqFs),1);
  for i = 1 : numel(ineqFs)
    f = ineqFs{i};
    c(i) = f(x(1), x(2));
  end
end
% ------------------------------------------------------------------------