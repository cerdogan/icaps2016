% @name solveL1Distance.m
% @author Can Erdogan
% @date 2015-10-01
% @brief Given two points p0, p2 and distances d01 and d12, find p1 such 
% that |p0p1|^1 = d01 and |p1p2|^1 = d12. Distances in L1 norm.
% We want to solve |x1-c0x| + |y1-c0y| = d01 and |x1-c2x| + |y1-c2y| = d12.
% Two options: 
% x1 same, y1 different.
% y1 same, x1 different.

function p1 = solveL1Distance (p0, p2, d01, d12)
%[p0, p2, d01, d12]
    c0x = p0(1); c0y = p0(2);
    c2x = p2(1); c2y = p2(2);
    coeffs = [[1 1]; [1 -1]; [-1 1]; [-1 -1]];
    for i = 1 : 4
        a0 = coeffs(i,1); b0 = coeffs(i,2);
        for j = 1 : 4
            
            a2 = coeffs(j,1); b2 = coeffs(j,2);
 %           [a0,b0,a2,b2]
            temp = a2 * b0 - a0 * b2;
            if(temp == 0), continue; end;
            rhs = (a2 * d01 - a0 * d12) - a0 * a2 * (c2x - c0x) - a0 * b2 * c2y + a2 * b0 * c0y;
            y1 = rhs / temp;
            %x1 = (d01 - b0 * (y1 - c0y) + a0 * c0x) / a0;
            x1 = (b2*d01 - b0*d12 + a0*b2*c0x - a2*b0*c2x + b0*b2*c0y - b0*b2*c2y)/(a0*b2 - a2*b0);
            y1 = -(a2*d01 - a0*d12 + a0*a2*c0x - a0*a2*c2x + a2*b0*c0y - a0*b2*c2y)/(a0*b2 - a2*b0);
            p1 = [x1, y1];
          %  [p1, norm(p1-p0,1), norm(p1-p2,1)]
            if((norm(p1 - p0, 1) - d01) > 1e-3), continue; end;
            if((norm(p1 - p2, 1) - d12) > 1e-3), continue; end;
  %          fprintf('-------------------------------\n');
   %         norm(p1 - p2, 1)
            return;
        end;
    end;
%     c0x = p0(1); c0y = p0(2);
%     c2x = p2(1); c2y = p2(2);
%     [p0, p2, d01, d12]
%     % Case 1: (x1-c0x) + (y1-c0y) = d01 && (x1-c2x) + (c2y-y1) = d12
%     % (c2x - c0x) + 2y1 - (c0y+c2y) = d01 - d12
%     y1 = ((d01-d12) - (c2x - c0x) + (c0y+c2y)) / 2;
%     x1 = d12 - (y1 - c0y) + c0x;
%     p1 = [x1,y1]
%     norm(p0-p1,1)
end