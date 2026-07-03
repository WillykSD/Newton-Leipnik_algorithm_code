%%%% Deterministic PZT model (discontinuous case)
%%%% Amplitude response 
%%%% Direct Numerical Method

clear all
clc
global nu1 nu2 mue muP omega1 lambda1 lambda2 gamma0 gamma1 gamma2 beta0 beta1 beta2 delta eta f0 omega omegae omegaP epsilon theta;



%% System parameters %%%%
nu1 = 0.3;  nu2 = 0.05;  mue = 1.0;  muP = 16.0;   omega1 = 0.50 ;                 lambda1 = 0.0045 ;
lambda2 = 0.00002; gamma0 = -.125;  gamma1 = 0.002;  gamma2 = 0.020;  beta0 = 9.996; 
beta1 = -0.16  ;  beta2 = -1.6;  delta = 0.001;   eta = -64.0; f0=.02;  omega=0.95;   

omegae = 50.0; omegaP = 282.843 ;    
epsilon = .99; theta = .99;

h=0.01; N=10000 ;

ff0 = (0.0:.05:0.1);
nv=length(ff0);
Amp=zeros(1,nv);  Ampb=zeros(1,nv); Ampc=zeros(1,nv);
for ii=1:nv

f0=ff0(ii);

%% SOLVER

r1 = omegae^(epsilon-1);   r2 = omegaP^2/omegae^(1+theta);


%% Functions for Newton-Leipnik algorithm
G=@(x,u,v,z,t) (  -( nu1 + nu2*sign(u) )*u - omega1^2*x - lambda1*x^2 - lambda2*x^3 - delta*z - ( gamma0 + gamma1*x + gamma2*x^2 )*v + f0*cos( omega*t ) ) ;
HH=@(x,u,y,v) ( 1/r1*(-mue*v-y-(beta0+beta1*x+beta2*x^2)*u) );

%% Initial conditions
x=0.0000000;  u=0;  y=0; v=0;  z=0; w=0;  t=0;  c0Epsi=1; c0Theta=1;     %%%%%%%%%%%   h=0.001; N=100000 ;

%% Newton-Leipnik algorithm

TABc0Epsi=zeros(1,N); TABV=zeros(1,N);  TABc0Theta=zeros(1,N);   TABZ=zeros(1,N); 

fid = fopen('table1.txt', 'w');

for i=1:N
    
    x0=x; u0=u; v0=v; z0=z; t0=t;
    
    x=x+h*u;                                                              
    u=u+h*G(x,u,v,z,t);                                                  
    y=y+h*v;                                                             
    
    %%%
    c0Epsi=(1-(1+epsilon)/i)*c0Epsi ;    
    TABc0Epsi(1,i)=c0Epsi;
    
    TABV(1,i)=v;
    
    SUMV=0;
    for j=1:i
        SUMV=SUMV+TABc0Epsi(j)*TABV(i+1-j);
    end
    %%
    
    v=HH(x,u,y,v)*h^epsilon-SUMV ;                                        
    
                                           
    %%
    c0Theta=(1-(2-theta)/i)*c0Theta ;    
    TABc0Theta(1,i)=c0Theta;
    
    TABZ(1,i)=z;
     
    SUMZ=0;
    for j=1:i
        SUMZ=SUMZ+TABc0Theta(j)*TABZ(i+1-j);
    end    
    %%%
    
    psi=(z+h*w+SUMZ)*h^(theta-1);                                        
    w=w+h*(-muP*w-r2*psi-eta*G(x0,u0,v0,z0,t0));                          
    z=z+h*w;                                                              
    t=t+h; 
    
    fprintf(fid, '%f %f %f %f %f  %f %f\n',t,x,u,y,v,z,w); 
end

fclose(fid);

clear('TABc0Epsi','TABV','TABc0Theta','TABZ')
clc

fid = fopen('table1.txt','r');
[A,COUNT] = fscanf(fid,'%f',[7,inf]);



%% Amplitude x
 xx0=A(2,N);

 xx1 = xx0;
 xx2 = xx1;

jj=N-1;

while jj<=N && jj>6000
   
   xx2 = xx1;
   xx1 = xx0;
   xx0 = A(2,jj);
    
   if (xx2<=xx1 && xx1>=xx0)
       xm=xx1;
   end
   
   if (xx2>=xx1 && xx1<=xx0)
       xmi=xx1;
   end 

   jj=jj-1;
end

Amp(ii)=(xm-xmi)/2;



end

figure(9);hold on
plot(ff0,Amp,'ko','Markersize', 6);

% figure(2);hold on
% plot(ff0,Ampb,'ko','Markersize', 6);
% 
% figure(3);hold on
% plot(ff0,Ampc,'ko','Markersize', 6);

