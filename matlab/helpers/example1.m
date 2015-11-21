%%
% @file example1.m
% @author Can Erdogan
% @date 2015-10-10
% @brief Visualizes the first set of inequalities for a simple example.

% ------------------------------------------------------------------------
function [ineqs, lims] = example1 (visualize)

  clf
  lims = [-10 10 0 40];
  
  % Setup the variables
  syms x1 x2
  global xs ineqFs
  xs = [x1,x2];

  % Generate the inequalities
  ineq1 = x1^2 - 40 + x2; % <= 0
  ineq2 = x1^2 - x2; % <= 0
  ineq3 = -4*x1 - 30 + x2; % <= 0 
  ineqs = [ineq1; ineq2; ineq3];
  if(nargin < 1 || not(visualize)), return; end;

  % Setup the functions for faster evaluation
  ineqFs = {};
  for i = 1 : numel(ineqs), ineqFs{i} = matlabFunction(ineqs(i)); end;
  
  % Visualize the inequalities
  x = lims(1):1:lims(2);
  vals1 = coeffs(ineq1,x2);
  plot(x,subs(-vals1(1)/vals1(2),x), '-c'); hold on;
  vals2 = coeffs(ineq2,x2);
  plot(x,subs(-vals2(1)/vals2(2),x), '-g'); hold on;
  vals3 = coeffs(ineq3,x2);
  plot(x,subs(-vals3(1)/vals3(2),x), '-r'); hold on;
  axis(lims); axis equal
 
  % Perform nonlinear programming to minimize distance to a random point
  % and stay within the feasible spaces of the inequalities
  options=optimset('Display','none', 'Algorithm', 'active-set');
  for i = 1 : 10000
  
    % Generate the random source point
    sourcePoint = [(lims(2) - lims(1))*rand + lims(1), (lims(4) - lims(3))*rand + lims(3)];
    fun = @(x) norm(sourcePoint - x);
  
    % Local minima of the distance function is sufficient as long as a
    % feasible point is reached
    while(true)
      initPoint = [(lims(2) - lims(1))*rand + lims(1), (lims(4) - lims(3))*rand + lims(3)];
      [solution, ~, exitflag] = fmincon(fun, initPoint, [], [], [], [], [lims(1), lims(3)], [lims(2), lims(4)], @nonlcons, options);
      if(exitflag ~= -2), break; end;
    end
    
    % Plot the feasible sample
    plot(solution(1), solution(2), 'ko');
    if(mod(i,25) == 1),pause(0.01); end;
  end
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