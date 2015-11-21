% Plot the constraints
clear all;
syms x1 x2;
mainEq = x1^2 - 40 + x2; % <= 0
solution = [
     -1.04960381014068e-05
          2.49999999907264];
    sample = [      -9.79819216474075
          61.6945556831228
]';
solution = [-8;40];
sample = [4;-20];

% y = mx + b; (y1 - mx1) = (y2 - mx2) -> (y2 - y1) = mx2 - mx1 => m = (y2 - y1) / (x2 - x1).
    m = (sample(2) - solution(2)) / (sample(1) - solution(1));
    b = sample(2) - m * sample(1);
    lineEq = m * x1 - x2 + b;

clf
plot(-10:0.1:10, -(-10:0.1:10).^2 +40); hold on;  % y <= -x^2 + 40
plot([sample(1); solution(1)], [sample(2); solution(2)], '-ro'); hold on;  % y <= 4x + 30
axis([-10 10 -60 60]);

% Substitute the line
val = coeffs(lineEq, x2)
temp1 = subs(mainEq, x2, -val(1)/val(2))
sols_x2 = eval(solve(temp1))
sols_y2 = subs(-val(1)/val(2), x1, sols_x2);
ps = [sols_x2, sols_y2]' 
      plot(ps(1,:), ps(2,:), '-go'); hold on;
      distsSample = sqrt(sum(abs(ps - repmat(sample',size(ps,1),1)).^2,2));
v = (solution - sample) / norm(solution - sample);