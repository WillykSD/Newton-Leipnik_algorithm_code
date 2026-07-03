%%%% Bifurcation diagram 
%%%% Control parameter = theta or epsilon of f0     


clear all
clc
global nu1 nu2 mue muP omega1 omega3 lambda1 lambda2 gamma0 gamma1 gamma2 beta0 beta1 beta2 delta eta f0 omega;

%% System parameters %%%%
nu1 = 0.3;  nu2 = 0.05;  mue = 1.0;  muP = 16.0;   omega1 = 0.50 ;  omega3 = 5.66; lambda1 = 0.0045 ;
lambda2 = 0.00002; gamma0 = -.125;  gamma1 = 0.002;  gamma2 = 0.020;  beta0 = 9.996; 
beta1 = -0.16  ;  beta2 = -1.6;  delta = 0.001;   eta = -64.0; f0=3;  omega=0.5;   

omegae = 50.0; omegaP = 282.843 ;    



epsilon = .8; ttheta = 0.4:.005:0.995;

ne=length(ttheta);

mi=50000;


fid = fopen('bifur11.txt', 'w');

for jj=1:ne
theta = ttheta(jj);
    
    
r1 = omegae^(epsilon-1);   r2 = omegaP^2/omegae^(1+theta);



%% Functions 
G=@(x,u,v,z,t) (  -( nu1 + nu2*sign(u) )*u - omega1^2*x - lambda1*x^2 - lambda2*x^3 - delta*z - ( gamma0 + gamma1*x + gamma2*x^2 )*v + f0*cos( omega*t ) ) ;
HH=@(x,u,y,v) ( 1/r1*(-mue*v-y-(beta0+beta1*x+beta2*x^2)*u) );

%% Initial conditions
x=0;  u=0;  y=0; v=0;  z=0; w=0;  t=0; h=0.01; N=100000 ; c0Epsi=1; c0Theta=1;
xmi=x; xm=x; xmax=x;

%%
TABc0Epsi=zeros(1,N); TABV=zeros(1,N);  TABc0Theta=zeros(1,N);   TABZ=zeros(1,N); 



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
    
                                           
    %%%
    c0Theta=(1-(2-theta)/i)*c0Theta ;    
    TABc0Theta(1,i)=c0Theta;
    
    TABZ(1,i)=z;
     
    SUMZ=0;
    for j=1:i
        SUMZ=SUMZ+TABc0Theta(j)*TABZ(i+1-j);
    end    
    %%
    
    psi=(z+h*w+SUMZ)*h^(theta-1);                                         
    w=w+h*(-muP*w-r2*psi-eta*G(x0,u0,v0,z0,t0));                          
    z=z+h*w;                                                              
    t=t+h; 
    
    xmi=xmax; xmax=xm; xm=x  ;
 
    if (i>=mi && (xmax-xmi)>0.0 && (xm-xmax)<0.0)
        fprintf(fid, '%f  %f\n',theta,xmax);
    end
 
end
display(jj)
end

fclose(fid);


fid = fopen('bifur11.txt','r');
[AA,COUNT] = fscanf(fid,'%f',[2,inf]);
plot(AA(1,:),AA(2,:),'.','Markersize', 6); ST = fclose(fid); 




