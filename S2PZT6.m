%%%% Deterministic PZT model (discontinuous case)
%%%% Numerical solution using Newton-Leipnik algorithm
%%%%       S2PZT7.m    and    S2PZT8.m  are best programs for the purpose
 

clear all
clc
global nu1 nu2 mue muP omega1 omega3 lambda1 lambda2 gamma0 gamma1 gamma2 beta0 beta1 beta2 delta eta f0 omega;

%% System parameters %%%%
nu1 = 0.3;  nu2 = 0.05;  mue = 1.0;  muP = 16.0;   omega1 = 0.50 ;  omega3 = 5.66; lambda1 = 0.0045 ;
lambda2 = 0.00002; gamma0 = -.125;  gamma1 = 0.002;  gamma2 = 0.020;  beta0 = 9.996; 
beta1 = -0.16  ;  beta2 = -1.6;  delta = 0.001;   eta = -64.0; f0=.02;  omega=0.95;   

omegae = 50.0; omegaP = 282.843 ;    
epsilon = .65; theta = .95;   r1 = omegae^(epsilon-1);   r2 = omegaP^2/omegae^(1+theta);



%% Functions for Newton-Leipnik algorithm
G=@(x,u,v,z,t) (  -( nu1 + nu2*sign(u) )*u - omega1^2*x - lambda1*x^2 - lambda2*x^3 - delta*z - ( gamma0 + gamma1*x + gamma2*x^2 )*v + f0*cos( omega*t ) ) ;
HH=@(x,u,y,v) ( 1/r1*(-mue*v-y-(beta0+beta1*x+beta2*x^2)*u) );

%% Initial conditions
x=0;  u=0;  y=0; v=0;  z=0; w=0;  t=0; h=0.01; N=10000 ; c0Epsi=1; c0Theta=1;

%% Newton-Leipnik algorithm


fid1 = fopen('tabc0Epsi.txt', 'w'); fid2 = fopen('tabV.txt', 'w');  fid3= fopen('tabc0Theta.txt', 'w');   fid4 = fopen('tabZ.txt', 'w'); 

fid = fopen('table1.txt', 'w');

for i=1:N
    
    x0=x; u0=u; v0=v; z0=z; t0=t;
    
    x=x+h*u;                                                              %1
    u=u+h*G(x,u,v,z,t);                                                   %2
    y=y+h*v;                                                              %3
    
    %%
    c0Epsi=(1-(1+epsilon)/i)*c0Epsi ;    
    fprintf(fid1, '%f\n',c0Epsi);
    
    fprintf(fid2, '%f\n',v);
    %fclose(fid1); fclose(fid2);
    
    fop1 = fopen('tabc0Epsi.txt','r');
    [TABc0Epsi,COUNT] = fscanf(fop1,'%f',[1,inf]);
    
    fop2 = fopen('tabV.txt','r');
    [TABV,COUNT] = fscanf(fop2,'%f',[1,inf]);
    
    SUMV=0;
    for j=1:i
        SUMV=SUMV+TABc0Epsi(j)*TABV(i+1-j);
    end
    %%
    
    v=HH(x,u,y,v)*h^epsilon-SUMV ;                                        %4
    
                                           
    %%
    c0Theta=(1-(2-theta)/i)*c0Theta ;    
    fprintf(fid3, '%f\n',c0Theta);
    
    fprintf(fid4, '%f\n',z);
    %fclose(fid1); fclose(fid2);
    
    fop3 = fopen('tabc0Theta.txt','r');
    [TABc0Theta,COUNT] = fscanf(fop3,'%f',[1,inf]);
    
    fop4 = fopen('tabZ.txt','r');
    [TABZ,COUNT] = fscanf(fop4,'%f',[1,inf]);
    
    SUMZ=0;
    for j=1:i
        SUMZ=SUMZ+TABc0Theta(j)*TABZ(i+1-j);
    end    
    %%
    
    psi=(z+h*w+SUMZ)*h^(theta-1);                                         %5
    w=w+h*(-muP*w-r2*psi-eta*G(x0,u0,v0,z0,t0));                          %6
    z=z+h*w;                                                              %7
    t=t+h;
    
    fprintf(fid, '%f %f %f %f %f  %f %f\n',t,x,u,y,v,z,w); 
end
fclose(fid);


fid = fopen('table1.txt','r');
[A,COUNT] = fscanf(fid,'%f',[7,inf]);

figure(1);hold on
plot(A(1,:),A(2,:),'k');