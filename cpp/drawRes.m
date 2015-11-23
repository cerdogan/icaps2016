load res;
for i = 1 : 4 : (size(res,1))
	p = [res(i), res(i+1)];
	if(p(1) < 0)
		p(1) = -p(1); 
		len = abs(res(i+3));
		p0 = [res(i-4), res(i-3)];
		if(res(i-4) < 0), p0(1) = -p0(1); end;
		alpha = atan2(p(2) - p0(2), p(1) - p0(1));
		if(res(i+3) > 0)
			alpha2 = alpha + asin(len/norm(p0-p));
		else 
			alpha2 = alpha - asin(len/norm(p0-p));
		end
		p2 = p0 + [cos(alpha2), sin(alpha2)] * res(i-2);
		p2b = p + [cos(alpha2), sin(alpha2)] * res(i+2);
		plot([p2(1); p2b(1)], [p2(2); p2b(2)], 'r-o'); hold on;
	end
	
	r = res(i+2);
	if(r < 0)
		drawCircle(p,-r,'r');
	else 
		drawCircle(p,r);
	end
end
