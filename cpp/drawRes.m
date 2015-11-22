load res;
for i = 1 : 3 : (size(res,1))
	p = [res(i), res(i+1)];
	r = res(i+2);
	if(r < 0)
		drawCircle(p,-r,'r');
	else 
		drawCircle(p,r);
	end
end
