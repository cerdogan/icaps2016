%%
% @file simplicial.m
% @author Can Erdogan
% @date 2015-10-10
% @brief Implements simplicial approximation.
%
% ------------------------------------------------------------------------
%
% Algorithm summary:
% Input: n constraints C in d-dimensions. 
%
% 1- H = generateInitialHull(C): Find d+1 points where each point (1) lives 
% in the feasible space F of the constraints and the line segments 
% connecting each point lives entirely in F.
%
% 2- tH = fitSphereToHull(H): Fit the largest hypersphere to the convex 
% hull and return which edges in H are tangent to it.
%
% 3- [m_tH, sc] = findLargestEdge(tH): Fit hypersheres to each edge in tH and
% return the edge with the maximum radius, m_th, and the center of the
% hypersphere it has.
%
% 4- pNew = extendOut (m_tH, sc, H, C): Along the outer normal of m_th,
% starting from sc, find a point that is in F (hopefully on boundary of C),
% such that any edge between this point and those in the hull H lives
% entirely in F. This is the new point to the hull.
%
% 5- Repeat from 2-4 as necessary.
%
% ------------------------------------------------------------------------
% 
% Helper functions:
%
% a- lineSegmentIntersect(): Given a constraint and a line, returns the 
% intersection point assuming general coordinates.
%
% b- checkLineSegmentsFeasible (p, H, C): Checks if the line segments
% from point p to any vertex that define hull H live in the feasible
% space defined by the constraints C.
%
% c- newP = tryNewBoundaryPoint (p0, H, C): Attempt to find a point in the 
% feasible space that can be used to connect to the other points in the 
% hull with all the line segments in % the feasible space. The random point 
% is generated by finding a random direction to shoot out from from the 
% central point, p0, and finding its projection to the constraints' boundary.
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
function [] = simplicial () 

  % Constants 
  vis = false;
  randomSeed = randi(1024,1);
  fprintf('Random seed: %d\n', randomSeed);
  rng(randomSeed);
  
  % Generate the symbolic constraints
  [C, Cf, lims] = example2;
  if(vis), plotConstraints(C, lims); end

  % Generate minimal simplex
  [H, k] = generateInitialHull(C, Cf, lims, vis);
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
        for ineqIdx = 1 : numel(C)
          
          [res, val] = lineSegmentIntersect(C(ineqIdx), pm', pm2);
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
      if(HULL > 2), C(end) = []; end;
      
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
      C(end+1) = newIneq;
      
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
      [bPs, k, radius] = addPoint(bPs, C);
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
function [H, k] = generateInitialHull (C, Cf, lims, vis) 

  % Find a feasible sample 
  global Cf_ optimOffset;
  Cf_ = Cf;
  optimOffset = 4;
  sourcePoint = [(lims(2) - lims(1))*rand + lims(1), (lims(4) - lims(3))*rand + lims(3)]';
  fun = @(x) norm(sourcePoint - x);
  p0 = fmincon(fun, sourcePoint, [], [], [], [], [-10, -60], [10 60], @nonlcons);
  % solution = [8.395, -20.62]';
  optimOffset = 0;

  % Visualize
  if(vis)
    plot(p0(1), p0(2), 'ko', 'MarkerSize', 7, 'LineWidth', 4);
    pause(0.01);
  end
  
  % Choose a random line, find intersections with other inequalities,
  % find the closest one to the solution which we assume is in the
  % feasible space, collect three such points.
  H = [];
  dimensions = 2;
  while(true) 

    % Check if collected sufficient points
    if(size(H,1) > dimensions), break; end;
    
    % Try to find a new point to add to the boundary
    [success, val] = addPointForHull0 (p0, H, C, lims, vis);
    if(not(success)), continue; end;
    
    % Add the point
    H(end+1,:) = val(:);
    fprintf('Added new point to initial simplex\n');
  end

  % Compute the convex hull
  k = convhull(H(:,1), H(:,2));
  
  if(vis)
    plot(H(k,1),H(k,2),['r-'],H(:,1), H(:,2),'b*','LineWidth', 2);
  end
  
end

% ------------------------------------------------------------------------
function [success, newPoint] = addPointForHull0 (p0, H, C, lims, vis)

  success = false;
  
  % Find random line to shoot out from the central point
  temp = [(lims(2) - lims(1))*rand + lims(1), (lims(4) - lims(3))*rand + lims(3)]';
  randDir = (temp - p0);
  randDir = randDir / norm(randDir);
  sample = p0 + 60 * randDir;
  
  % Visualize the line
  if(vis)
    h = plot([sample(1); p0(1)], [sample(2); p0(2)], '-o', 'LineWidth', 2);
  end

  % Intersect the random line with all the constraint and pick closest
  minDist = 1000000;
  minClosestPoint = [];
  for ineqIdx = 1 : numel(C)
    [res, point] = lineSegmentIntersect(C(ineqIdx), p0, sample);
    if(not(res)), continue; end;
    dist = norm(point'-p0);
    if(dist < minDist)
      minDist = dist;
      minClosestPoint = point'; 
    end
  end
  assert(minDist ~= 1000000);
  newPoint = minClosestPoint;
  
  % Check if the point is within the feasible region
  syms x1 x2;
  for j = 1 : numel(C)
    value = subs(C(j), [x1 x2], [newPoint(1), newPoint(2)]);
    if(value > 0)
      fprintf('Neighbor segment goes directly into infeasible space\n');
      if(vis), delete(h); end;
      return
    end
  end

  % Check if the point is too close to the others
  if(size(H,1) > 0)
    closestDist = min(sqrt(sum(abs(H - repmat(newPoint',size(H,1),1)).^2,2)));
    if(closestDist < 2)
      fprintf('Too close to others: %f\n', closestDist);
      if(vis), delete(h); end;
      return; 
    end
  end;

  % Check if the line segments that connect the point to the hull points
  % live entirely in the feasible space
  if(not(checkLineSegmentsFeasible(newPoint, H, C)))
    fprintf('Connecting line segments not in feasible space.\n');
    if(vis), delete(h); end;
    return; 
  end;
  
  % Check if the new point could be generated as a combination of the
  % others (collinearity) by checking the condition number of the points
  % and bounding it 
  if(size(H,1) > 1)
    eigs = svd([H; newPoint']);
    fprintf('Condition number: %f\n', abs(eigs(1) / eigs(end)));
    if(abs(eigs(1) / eigs(end)) > 5)
      fprintf('Colinearity.\n');
      if(vis), delete(h); end;
      return; 
    end
  end
    
  if(vis), plot(newPoint(1), newPoint(2), 'cs', 'MarkerSize', 2, 'LineWidth', 3); hold on; end;
  success = true;
end

% ------------------------------------------------------------------------
function bool = checkLineSegmentsFeasible (p, H, C) 
  
  bool = false;
  
  % Check if the intersection can be connected with a line segment with the 
  % points in the convex hull without getting out of the feasible space
  for pIdx = 1 : size(H,1)
    
    % Get edge information
    q = transpose(H(pIdx,:));
    v = (q - p) / norm(q - p);
    startPoint = p + v * 1e-2;
    endPoint = q - v * 1e-2;
    
    % Check if the start and end points are in the feasible space
    if(sum(nonlcons(startPoint) < 0) < 4)
      fprintf('Edge start outside feasible space.\n');
      return; 
    end;
    if(sum(nonlcons(endPoint) < 0) < 4)
      fprintf('Edge start outside feasible space.\n');
      return; 
    end;
  
    % Check if the line to the hull point intersects the inequalities
    for j = 1 : numel(C)
      [res, ~] = lineSegmentIntersect(C(j), startPoint, endPoint);
      if(res)
        fprintf('A part of the convex hull is in the infeasible region\n');
        return; 
      else
        fprintf('\tInequality %d is clear: %d\n', j, res);
      end;
    end
    
  end

  bool = true;
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
    val = real(ps(abs(projs  - minVal) < 1e-3,:));
    %plot(val(1), val(2), 'ro'); hold on;
%     if(res(1) == 1), val = (ps(1,:)); 
%     else val = (ps(2,:)); end;
  else
   % keyboard
    fprintf('Point not in between the two points\n');
    return;
  end
  
  success = true;
end

% ------------------------------------------------------------------------
function bestEdgeIdx = findLargestEdge () 

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
function [] = plotConstraints (C, lims)

  % Plot the constraints
  clf
  syms x2;
  x = lims(1):0.1:lims(2);
  for i = 1 : numel(C)
    ineq = C(i);
    vals = coeffs(ineq,x2);
    plot(x,subs(-vals(1)/vals(2),x)); hold on;
  end
  axis(lims); axis equal
  
end

% ------------------------------------------------------------------------
function [c, ceq] = nonlcons (x)
  ceq = [];
  global Cf_ optimOffset;
  c = zeros(numel(Cf_),1);
  for i = 1 : numel(Cf_)
    f = Cf_{i};
    c(i) = f(x(1), x(2)) + optimOffset;
  end
end
% ------------------------------------------------------------------------
