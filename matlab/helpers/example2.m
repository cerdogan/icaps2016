%%
% @file example2.m
% @author Can Erdogan
% @date 2015-10-10
% @brief Visualizes the second set of inequalities for a challenging ex.

% ------------------------------------------------------------------------
function [ineqs, lims] = example2 (visualize)

  clf
  lims = [-30 30 -30 30];
  
  % Setup the variables
  x = lims(1):1:lims(2);
  syms x1 x2
  global xs ineqFs
  xs = [x1,x2];
  
  % Generate the inequalities
  ineq1 = (-(-(0.2*(-x1-5)).^4 + (0.4 * -x1).^3 + 1.5*(-x1).^2 +10) - x2);
  ineq2 = -(-(0.25*(x1)+3).^2+20 - x2);
  ineq3 = -(-(0.2*x1+0.5).^5 + 2*(0.2*x1+0.5).^4 + 3*(0.2*x1+0.5).^2 - x2);
  ineq4 = (-(0.25*x1+0.5).^3-20 - x2);
  ineqs = [ineq1; ineq2; ineq3; ineq4];
  if(nargin < 1 || not(visualize)), return; end;
  
  % Setup the functions for faster evaluation
  ineqFs = {};
  for i = 1 : numel(ineqs), ineqFs{i} = matlabFunction(ineqs(i)); end;
  
  % Plot the inequalities
  vals1 = coeffs(ineq1,x2);
  plot(x,subs(-vals1(1)/vals1(2),x), '-c'); hold on;
  vals2 = coeffs(ineq2,x2);
  plot(x,subs(-vals2(1)/vals2(2),x), '-g'); hold on;
  vals3 = coeffs(ineq3,x2);
  plot(x,subs(-vals3(1)/vals3(2),x), '-r'); hold on;
  vals4 = coeffs(ineq4,x2);
  plot(x,subs(-vals4(1)/vals4(2),x), '-'); hold on;
  axis(lims); axis equal;
  
  % Perform nonlinear programming to minimize distance to a random point
  % and stay within the feasible spaces of the inequalities
  maxNum = 1000000;
  num = 0;
  for i = 1 : maxNum
    sourcePoint = [(lims(2) - lims(1))*rand + lims(1), (lims(4) - lims(3))*rand + lims(3)];
    [res, ~] = nonlcons(sourcePoint);
    if(sum(res < 0) < 4), continue; end;
    num = num + 1;
    if(mod(i,25) == 1),pause(0.01); end;
  end
  
  V = (lims(2) - lims(1)) * (lims(4) - lims(3)) * (num / maxNum);
  fprintf('Sampled volume: %f\n', V);
  
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