%%%% Deterministic PZT model (discontinuous case): Approximate analytical solution 
%%%% (harmonic balance method)
%%%% Amplitude response Newton-Raphson method.    Ampl=f(F0)


clear all
clc
global nu1 nu2 mue muP omega1 omega3 lambda1 lambda2 gamma0 gamma1 gamma2 beta0 beta1 beta2 delta eta f0 omega;

%% System parameters %%%%
nu1 = 0.3;  nu2 = 0.05;  mue = 1.0;  muP = 16.0;   omega1 = 0.50 ;  omega3 = 5.66; lambda1 = 0.0045 ;
lambda2 = 0.00002; gamma0 = -.125;  gamma1 = 0.002;  gamma2 = 0.020;  beta0 = 9.996; 
beta1 = -0.16  ;  beta2 = -1.6;  delta = 0.001;   eta = -64.0; omega=4;  % f0=.02;   

omegae = 50.0; omegaP = 282.843 ;    
epsilon = 1; theta = 1;   r1 = omegae^(epsilon-1);   r2 = omegaP^2/omegae^(1+theta);

%% Other functions definition

varsigma1 = ( omega^(1+epsilon)*cos((1/2)*pi*epsilon) ) ;     
varsigma2 = ( omega^(1+epsilon)*sin((1/2)*pi*epsilon) ) ;

varrho1 = ( omega^(1-theta)*cos((1/2)*pi*theta) ) ;    
varrho2 = ( omega^(1-theta)*sin((1/2)*pi*theta) ) ;

R1 = ( (1-r1*varsigma2)/((1-r1*varsigma2)^2+(mue*omega+r1*varsigma1)^2) ) ;
R2 = ( 1/4*(3*lambda2-2*lambda1*gamma2/gamma1) ) ;
R3 = ( omega1^2-omega^2+delta*eta*omega^2*(r2*varrho2-omega^2)/((r2*varrho2-omega^2)^2+(muP*omega+r2*varrho1)^2) ) ;
R4p = ( -delta*eta*omega^2*(muP*omega+r2*varrho1)/((r2*varrho2-omega^2)^2+(muP*omega+r2*varrho1)^2)+(nu1+nu2)*omega ) ; %%%%
R4m = ( -delta*eta*omega^2*(muP*omega+r2*varrho1)/((r2*varrho2-omega^2)^2+(muP*omega+r2*varrho1)^2)+(nu1-nu2)*omega ) ; %%%%
R5 = ( -(mue*omega+r1*varsigma1)/((1-r1*varsigma2)^2+(mue*omega+r1*varsigma1)^2) ) ;
R6 = ( omega^2*(beta2*gamma0+gamma2*beta0) ) ;

%% Amplitude funtions definition
AFunP=@(a0,f0) ( a0^2*((((1/4)*omega*gamma2*a0^2+gamma0*omega)*((1/4)*omega*beta2*a0^2+beta0*omega)*R1+R2*a0^2+R3)^2+(R4p+((1/4)*omega*gamma2*a0^2+gamma0*omega)*((1/4)*omega*beta2*a0^2+beta0*omega)*R5)^2)-f0^2 ) ;
AFunM=@(a0,f0) ( a0^2*((((1/4)*omega*gamma2*a0^2+gamma0*omega)*((1/4)*omega*beta2*a0^2+beta0*omega)*R1+R2*a0^2+R3)^2+(R4m+((1/4)*omega*gamma2*a0^2+gamma0*omega)*((1/4)*omega*beta2*a0^2+beta0*omega)*R5)^2)-f0^2 ) ;



%% Newton-Raphson method
% a0 = 0.0:0.00001:1; 
% 
% %  f01 = 0.0:0.01:0.5;
% %  f02 = 0.51:0.01:1;
% omeg3a = 0.285:0.002:0.315;
% omeg4a = 0.32:0.005:0.37;
% omeg5a = 0.38:0.01:0.5;
% omeg6a = 0.55:0.05:5.0;   omega=[omeg1a omeg2a omeg3a omeg4a omeg5a omeg6a];

% omega = 0.0:0.01:5.0;   % a0 = 0.0:0.00001:.05; omega = 0.0:0.005:5.0;
       a0 = 0.0:0.0001:6;   f0 = 0.0:0.05:1;                 %%% a0 = 0.0:0.00001:6; f0=[f01 f02];
n=length(f0); m=length(a0);
fid = fopen('table1.txt', 'w');

for i=1:n
    for j=1:m-1
        e1=AFunP(a0(j),f0(i));
        e2=AFunP(a0(j+1),f0(i));
        if (e1*e2<=0)
            x=f0(i); y=( a0(j)+a0(j+1) )/2;
            b=abs(beta0+(1/4)*beta2*y^2)*omega*y/sqrt((r1*varsigma2-1)^2+(mue*omega+r1*varsigma1)^2);
            c=abs(eta)*omega^2*y/sqrt((r2*varrho2-omega^2)^2+(muP*omega+r2*varrho1)^2);
            fprintf(fid, '%f  %f %f %f\n',x,y,b,c); 
        end
    end
    display(i)
end
fclose(fid);

fid = fopen('table1.txt','r');
[AP,COUNT] = fscanf(fid,'%f',[4,inf]);


% figure(1);hold on
% plot(A(1,:),A(2,:),'k.','Markersize', 10);
% 
% figure(2);hold on
% plot(A(1,:),A(3,:),'k.','Markersize', 10);
% 
% figure(3);hold on
% plot(A(1,:),A(4,:),'k.','Markersize', 10);


fid = fopen('table2.txt', 'w');

for i=1:n
    for j=1:m-1
        e1=AFunM(a0(j),f0(i));
        e2=AFunM(a0(j+1),f0(i));
        if (e1*e2<=0)
            x=f0(i); y=( a0(j)+a0(j+1) )/2;
            b=abs(beta0+(1/4)*beta2*y^2)*omega*y/sqrt((r1*varsigma2-1)^2+(mue*omega+r1*varsigma1)^2);
            c=abs(eta)*omega^2*y/sqrt((r2*varrho2-omega^2)^2+(muP*omega+r2*varrho1)^2);
            fprintf(fid, '%f  %f %f %f\n',x,y,b,c); 
        end
    end
    display(i)
end
fclose(fid);

fid = fopen('table2.txt','r');
[AM,COUNT] = fscanf(fid,'%f',[4,inf]);

A=(AP+AM)./2;

figure(4);hold on
plot(A(1,:),A(2,:),'r.','Markersize', 10);

figure(5);hold on
plot(A(1,:),A(3,:),'r.','Markersize', 10);

figure(6);hold on
plot(A(1,:),A(4,:),'r.','Markersize', 10);
