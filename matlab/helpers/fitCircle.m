clf;
clear all;
bPs = [
         -2.64858623313084          19.4056550674766
          3.12613230789056          30.2272967935628
         -1.45725246313283          2.12358474130671
          3.95803228940697          15.6660196039882];
        
          
bPs =[

        -0.282198929381144         0.079636235743862
         -2.70124580696508          19.1950167721397
         -1.23057451993691          25.0777019202524
         -3.50423627068714          12.2796718407993
         ];
       
bPs = [        -0.101238231389933        0.0102491794949704
          1.90355852037285          3.62353504048406
          4.43093935781676          20.3667764073504
         -3.37978865798638          11.4229713726534];
       
% B. Generate the convex hull
plot(bPs(:,1), bPs(:,2), 'rs'); hold on;
k = convhull(bPs(:,1), bPs(:,2));
plot(bPs(k,1),bPs(k,2),'r-',bPs(:,1), bPs(:,2),'b*','LineWidth', 2)
pause(0.01);

% Generate the inequalities for the convex hull
A = [];
b = [];
for i = 1 : size(k,1) - 1
  p1 = bPs(k(i),:); p2 = bPs(k(i+1),:);
  v21 = (p2 - p1) / norm(p2 - p1);
  hi = -[-v21(2), v21(1)];
  pm = (p2 + p1) / 2.0;
  pe = pm + hi;
  ps = [pm; pe];
  plot(ps(:,1), ps(:,2), 'k-'); hold on;
  ps = [p1; p2];
  plot(ps(:,1), ps(:,2), 'k-', 'LineWidth', 3); hold on;
  axis equal  
  ci = dot(hi,p1);
  A(end+1,:) = hi;
  b(end+1,:) = ci;
end
mp = mean(bPs);
plot(mp(1), mp(2), 'ro'); hold on;
if(dot(mean(bPs), A(2,:)) > b(2)), A = -A; b = -b; 'flipped', end;

% Solve for the center of the circle
A(:,end+1) = -ones(size(A,1),1);
[res, ~, exitflag] = linprog([0,0,1], A, b,[],[],[-10,-60,-60],[10,60,60]);
exitflag
center = res(1:2)'
radius = res(3)
drawCircle(center, radius);
axis([min(bPs(:,1)), max(bPs(:,1)), min(bPs(:,2)), max(bPs(:,2))]); axis equal