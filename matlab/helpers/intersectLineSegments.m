% -------------------------------------------------------------------------
function res = intersectLineSegments (xa, ya, xb, yb, xc, yc, xd, yd, off)

  % Check if the two lines are parallel
  res = [0,0,0];
  if ((abs(xb - xa) < 1e-4) && (abs(xd - xc) < 1e-4)), return;
  elseif ((abs(yb - ya) < 1e-4) && (abs(yd - yc) < 1e-4)), return;
  elseif ((abs(xb - xa) > 1e-4) && (abs(xd - xc) > 1e-4))
    m1 = (yb - ya) / (xb - xa); 
    m2 = (yd - yc) / (xd - xc); 
    if (abs(m1 - m2) < 1e-4), return; end
  end

  % Find the intersection of the two lines
  if (abs(xb - xa) < 1e-4)        % First line is vertical
    m2 = (yd - yc) / (xd - xc); 
    b2 = yc - m2 * xc;
    x = xa;
    y = m2 * xa + b2;
  elseif (abs(xd - xc) < 1e-4)    % Second line is vertical
    m1 = (yb - ya) / (xb - xa);
    b1 = ya - m1 * xa;
    x = xc;
    y = m1 * xc + b1;
  else                            % Normal case
    m1 = (yb - ya) / (xb - xa); 
    b1 = ya - m1 * xa;
    m2 = (yd - yc) / (xd - xc); 
    b2 = yc - m2 * xc;
    x = (b2 - b1) / (m1 - m2); 
    y = m1 * x + b1;
  end
                                                                                                          
  % Check if intersection is within the first line segment
  d1 = sqrt((xb - xa) * (xb - xa) + (yb - ya) * (yb - ya));
  v1x = (xb - xa) / d1;
  v1y = (yb - ya) / d1;
  proj1 = (x - xa) * v1x + (y - ya) * v1y;
  if((proj1 < off-1e-3) || (proj1 > d1 - off + 1e-3)), return; end;

  % Check if intersection is within the second line segment
  d2 = sqrt((xd - xc) * (xd - xc) + (yd - yc) * (yd - yc));
  v2x = (xd - xc) / d2;
  v2y = (yd - yc) / d2;
  proj2 = (x - xc) * v2x + (y - yc) * v2y;
  if (proj2 < off-1e-3 || proj2 > d2 - off + 1e-3), return; end;

  % If no parallels and intersection is in segments, return point
  res = [x,y,1];

end