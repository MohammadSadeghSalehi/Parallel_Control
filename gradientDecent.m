function [control,Jc]= gradientDecent(Tl_1,Tl,y0,yref,alpha,c)
N	=41;
t	=linspace(Tl_1,Tl,N);
dt =t(2)-t(1);


A	=[1 0;0 2];
B	=[0 1;1 0];

%alpha 	=1e1;
%y0	=[1;0];
%yref	=[0;1];

maxiter	= 1;
	%c	= zeros(1,N-1);
	rho	= 1e0;

	for iter=1:maxiter
		[J,GradJ]=cost(N,y0,yref,alpha,A,B,c,dt);
		%fprintf(1,'Iter=%i|J=%f|Grad=%f\n',iter,J,norm(GradJ))	
		c	= c-rho*GradJ;	

  end

  Jc=J;
  control=c;
%%cost function%%
function [J,GradJ,x,p]=cost(N,y0,yref,alpha,A,B,c,dt)

x=zeros(2,N);
x(:,1)=y0;
for n=2:N
	x(:,n)=expm(i*A*dt)*expm(i*B*c(n-1)*dt)*x(:,n-1);
end

J=-real(x(:,end)'*yref)+.5*alpha*dt*c*c';

p=zeros(2,N);
%%%%
p(:,N)=-yref;
for n=N:-1:2
	p(:,n-1)=expm(-i*B*c(n-1)*dt)*expm(-i*A*dt)*p(:,n);
end

GradJ=alpha*dt*c;
for n=1:N-1
	GradJ(n)= GradJ(n) + real(p(:,n+1)' ...
				* expm(i*A*dt)*(i*B*dt)*expm(i*B*c(n)*dt)...
				* x(:,n)) ;	
end


end
endfunction
