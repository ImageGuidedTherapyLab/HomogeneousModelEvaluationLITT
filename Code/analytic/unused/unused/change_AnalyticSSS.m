function [temperature]= sammodel1D(u0,ua,k,w,P,r,mua,mus,R1,R2,anfact)
mutr  = mua+mus*(1.0-anfact);
mueff = sqrt(3.0*mua*mutr);
cblood = 3617.0 + 274.15; %value pulled from online, in Kelvin

	h = sqrt(w*cblood/k);	
	f1 = @(r_var) exp(r_var*h)/r_var;
	f2 = @(r_var) exp(-r_var*h)/r_var;
	
	A = f(R1);
	B = f(R2);
	C = exp(h*R2)(h*R2-1)/(R**2);
	D = -exp(-h*R2)(h*R2+1)/(R**2);
	
	E = u0 - (mueff**2)*P*exp(-mueff*r)/(4*pi*(w*cblood-k*(mueff**2)));
	F = (mueff**3)*P*exp(-mueff*r)/(4*pi*(w*cblood-k*(mueff**2)));

	C1 = E/(A*D-B*C);
	C2 = (F-C*C1)/D;
	
	up = (mueff**2)*P*exp(-mueff*r)/(4*pi*(w*cblood-k*(mueff**2)));
	
	temperature = C1*f1(r) + C2*f2(r) + up + ua;


end

