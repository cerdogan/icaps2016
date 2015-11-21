% @name solveL2Distance.m
% @author Can Erdogan
% @date 2015-10-01
% @brief Given two points p0, p2 and distances d01 and d12, find p1 such 
% that |p0p1|^2 = d01 and |p1p2|^2 = d12. Distances in L2 norm.

function p1 = solveL2Distance (p0, p2, rA, rB)

xA = p0(1); yA = p0(2); 
xB = p2(1); yB = p2(2); 
d = sqrt((xB-xA).^2 + (yB-yA).^2);
K = 0.25 * sqrt(((rA+rB).^2-d.^2)*(d.^2-(rA-rB).^2));
x1 = (1/2)*(xB+xA) + (1/2)*(xB-xA)*(rA.^2-rB.^2)/(d.^2) - 2*(yB-yA)*K/(d.^2);
y1 = (1/2)*(yB+yA) + (1/2)*(yB-yA)*(rA.^2-rB.^2)/(d.^2) - -2*(xB-xA)*K/(d.^2);
p1 = [x1,y1];
%x2 = (1/2)*(xB+xA) + (1/2)*(xB-xA)*(rA.^2-rB.^2)/(d.^2) + 2*(yB-yA)*K/(d.^2);
%y2 = (1/2)*(yB+yA) + (1/2)*(yB-yA)*(rA.^2-rB.^2)/(d.^2) + -2*(xB-xA)*K/(d.^2);
%p2 = [x2,y2];
%fprintf('K: %f', 2*(yB-yA)*K/(d.^2));
%[norm(p0 - p1), rA, norm(p2 - p1), rB]
%fprintf('------------------------------\n');
end