function parallelmod
%%some initialization%%
time=cputime();
N=6;
NN=201;
T=10;
A	=[1 0;0 2];
B	=[0 1;1 0];

alpha 	=1e1;
tl=linspace(0,T,N);
t=linspace(0,T,NN);
dt=t(2)-t(1);
y0	=[1;0];
yref	=[0;1];
betal=T/(tl(2)-tl(1));
alphal=alpha/betal;
c	= zeros(1,NN-1);
%%calculate p and y on whole interval
y=zeros(2,NN);
y(:,1)=y0;
p=zeros(2,NN);
p(:,NN)=-yref;
for n=1:NN-1
		y(:,n+1)=expm(i*A*dt)*expm(i*B*c(n)*dt)*y(:,n);
	end

for n=NN:-1:2
	p(:,n-1)=expm(-i*B*c(n-1)*dt)*expm(-i*A*dt)*p(:,n);
end
%%%%%%%%%%%%%%%%
for m=1:5
  
%%Calculate lambda^c

lambda=zeros(2,N);
lambda(:,1)=y(:,1);
for j=2:N
  lambda(:,j)= (1-(tl(j)/T))*y(:,(j-1)*40+1)-(tl(j)/T)*p(:,(j-1)*40+1);
endfor
%%%%%%%%%%%%%%%%
J=0;
cl=zeros(1,40);
MY=zeros(2,10);
MP=zeros(2,10);
for j=1:N-1
  for k=1:40
     cl(k)=c((j-1)*40+k);
   endfor
   
      [ctemp,MYtemp,MPtemp]=monotonicmod(tl(j),tl(j+1),lambda(:,j),lambda(:,j+1),alphal,cl,A,B,((NN-1)/(N-1))+1);
   for k=1:40
     c((j-1)*40+k)=ctemp(k);
   endfor
   MY(:,2*j-1)=MYtemp(:,1);
   MY(:,2*j)=MYtemp(:,2);
   MP(:,2*j-1)=MPtemp(:,1);
   MP(:,2*j)=MPtemp(:,2);
endfor

%%%%update y on whole interval
MYtemp=zeros(2,2);
MPtemp=zeros(2,2);

for n=1:N-1
    MYtemp(:,1)=MY(:,2*n-1);
    MYtemp(:,2)=MY(:,2*n);
		y(:,40*n+1)=MYtemp*y(:,40*(n-1)+1);
end

for n=N:-1:2
    MPtemp(:,1)=MP(:,2*(n-1)-1);
    MPtemp(:,2)=MP(:,2*(n-1));
    p(:,40*(n-2)+1)=MPtemp*p(:,40*(n-1)+1);
end

J=-real(y(:,end)'*(yref))+.5*alpha*dt*c*c';
		fprintf(2,'Iter=%i|J=%f\n',m,J);
end
time=cputime()-time
J
endfunction
