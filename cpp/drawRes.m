load res;
for i = 1 : 4 : (size(res,1))
	if(abs(res(i+3) +12) < 1e-1)
		p = [res(i), res(i+1)];
		r = res(i+2);
		if(r < 0)
			%drawCircle(p,-r,'r');
		else 
			drawCircle(p,r);
		end
	else
		p1 = [res(i), res(i+1)];
		p2 = [res(i+2), res(i+3)];
		ps = [p1;p2];
		plot(ps(:,1), ps(:,2), '-', 'LineWidth', 3); hold on;
	end
end
%load res;
%for i = 1 : 3 : (size(res,1))
%	p = [res(i), res(i+1)];
%	r = res(i+2);
%	if(r < 0)
%		drawCircle(p,-r,'r');
%	else 
%		drawCircle(p,r);
%	end
%end
