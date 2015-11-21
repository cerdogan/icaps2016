bPs = [
   -0.1012    0.0102
    1.9036    3.6235
    4.4309   20.3668
   -3.3798   11.4230
    3.4571   11.9514];
  
% Generate the convex hull
clf
plot(bPs(:,1), bPs(:,2), 'rs'); hold on;
k = convhull(bPs(:,1), bPs(:,2));
plot(bPs(k,1),bPs(k,2),'r-',bPs(:,1), bPs(:,2),'b*','LineWidth', 2);
axis([min(bPs(:,1)), max(bPs(:,1)), min(bPs(:,2)), max(bPs(:,2))]); axis equal

% Generate the inequalities for the convex hull
A = [];
b = [];
for i = 1 : (size(k,1) - 1)
  p1 = bPs(k(i),:); p2 = bPs(k(i+1),:);
  v21 = (p2 - p1) / norm(p2 - p1);
  hi = -[-v21(2), v21(1)];
  ci = dot(hi,p1);
  pm = (p2 + p1) / 2.0;
  pe = pm + hi;
  ps = [pm; pe];
  plot(ps(:,1), ps(:,2), 'k-', 'LineWidth', 2); hold on;
  ps = [p1; p2];
  plot(ps(:,1), ps(:,2), 'k-', 'LineWidth', 3); hold on;
  pause
  axis equal  
  A(end+1,:) = hi;
  b(end+1,:) = ci;
end
if(dot(mean(bPs), A(2,:)) > b(2)), A = -A; b = -b; 'flipped', end;

% Find the largest "edge"
n = 2;
for i = 4 : 4 %(size(k,1) - 1)
  
  % Accumulate the constraint for each edge
  thi = [];
  for j = 1 : (size(k,1) - 1)
      if(i == j), continue; end
      h1 = A((i), 1:n);
      h2 = A((j), 1:n);
      costh = dot(h1,h2) / (norm(h1) * norm(h2));
      th = acos(costh);
      [i,j,h1,h2,th/pi*180]
      thi(end+1) = sin(th);
  end
  thi
end