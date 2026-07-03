%%%% Deterministic PZT model (discontinuous case): Approximate analytical solution 
%%%% (harmonic balance method)
%%%% Amplitude response Newton-Raphson method.   
%%%% efficiency = f(epsilon,theta)

clear all
clc
global nu1 nu2 mue muP omega1 omega3 lambda1 lambda2 gamma0 gamma1 gamma2 beta0 beta1 beta2 delta eta f0 omega;

%% System parameters %%%%
nu1 = 0.3;  nu2 = 0.05;  mue = 1.0;  muP = 16.0;   omega1 = 0.50 ;  omega3 = 5.66; lambda1 = 0.0045 ;
lambda2 = 0.00002; gamma0 = -.125;  gamma1 = 0.002;  gamma2 = 0.020;  beta0 = 9.996; 
beta1 = -0.16  ;  beta2 = -1.6;  delta = 0.001;   eta = -64.0; f0=.02;  omega=10;   

omegae = 50.0; omegaP = 282.843 ; 

%% Initial array and varying parameters
epsilon=0.02:.02:1;   ttheta = 0.02:.02:1;

ne=length(ttheta);



a0 = 0.0:0.000001:0.005; % a0 = 0.13:0.00001:0.17; %a0 = 0.0:0.00001:0.04; %a0 = 0.040:0.00001:0.18; 
n=length(epsilon); m=length(a0);

P1maxP=zeros(ne,n);  P2maxP=zeros(ne,n);  AmpaP=zeros(ne,n);
P1maxM=zeros(ne,n);  P2maxM=zeros(ne,n);  AmpaM=zeros(ne,n);
P1max=zeros(ne,n);   P2max=zeros(ne,n);   Ampa=zeros(ne,n);


for jj=1:ne
    
theta = ttheta(jj);    

r1 =@(epsilon) (omegae^(epsilon-1));   r2 = omegaP^2/omegae^(1+theta);

%% Other functions definition

varsigma1 =@(epsilon) ( omega^(1+epsilon)*cos((1/2)*pi*epsilon) ) ;     
varsigma2 =@(epsilon) ( omega^(1+epsilon)*sin((1/2)*pi*epsilon) ) ;

varrho1 = ( omega^(1-theta)*cos((1/2)*pi*theta) ) ;    
varrho2 = ( omega^(1-theta)*sin((1/2)*pi*theta) ) ;

R1 =@(epsilon) ( (1-r1(epsilon)*varsigma2(epsilon))/((1-r1(epsilon)*varsigma2(epsilon))^2+(mue*omega+r1(epsilon)*varsigma1(epsilon))^2) ) ;
R2 = ( 1/4*(3*lambda2-2*lambda1*gamma2/gamma1) ) ;
R3 = ( omega1^2-omega^2+delta*eta*omega^2*(r2*varrho2-omega^2)/((r2*varrho2-omega^2)^2+(muP*omega+r2*varrho1)^2) ) ;
R4p = ( -delta*eta*omega^2*(muP*omega+r2*varrho1)/((r2*varrho2-omega^2)^2+(muP*omega+r2*varrho1)^2)+(nu1+nu2)*omega ) ; %%%%
R4m = ( -delta*eta*omega^2*(muP*omega+r2*varrho1)/((r2*varrho2-omega^2)^2+(muP*omega+r2*varrho1)^2)+(nu1-nu2)*omega ) ; %%%%
R5 =@(epsilon) ( -(mue*omega+r1(epsilon)*varsigma1(epsilon))/((1-r1(epsilon)*varsigma2(epsilon))^2+(mue*omega+r1(epsilon)*varsigma1(epsilon))^2) ) ;
R6 = ( omega^2*(beta2*gamma0+gamma2*beta0) ) ;

%% Amplitude funtions definition
AFunP=@(a0,epsilon) ( a0^2*((((1/4)*omega*gamma2*a0^2+gamma0*omega)*((1/4)*omega*beta2*a0^2+beta0*omega)*R1(epsilon)+R2*a0^2+R3)^2+(R4p+((1/4)*omega*gamma2*a0^2+gamma0*omega)*((1/4)*omega*beta2*a0^2+beta0*omega)*R5(epsilon))^2)-f0^2 ) ;
AFunM=@(a0,epsilon) ( a0^2*((((1/4)*omega*gamma2*a0^2+gamma0*omega)*((1/4)*omega*beta2*a0^2+beta0*omega)*R1(epsilon)+R2*a0^2+R3)^2+(R4m+((1/4)*omega*gamma2*a0^2+gamma0*omega)*((1/4)*omega*beta2*a0^2+beta0*omega)*R5(epsilon))^2)-f0^2 ) ;


%% Newton-Raphson method
for i=1:n
    for j=1:m-1
        %%
        e1=AFunP(a0(j),epsilon(i));
        e2=AFunP(a0(j+1),epsilon(i));
        if (e1*e2<=0)
            x=epsilon(i); y=( a0(j)+a0(j+1) )/2;
            b=abs(beta0+(1/4)*beta2*y^2)*omega*y/sqrt((r1(x)*varsigma2(x)-1)^2+(mue*omega+r1(x)*varsigma1(x))^2);
            c=abs(eta)*omega^2*y/sqrt((r2*varrho2-omega^2)^2+(muP*omega+r2*varrho1)^2);
            P1maxP(jj,i)=b; P2maxP(jj,i)=c; AmpaP(jj,i)=y;
        end
        %
        %
        e1=AFunM(a0(j),epsilon(i));
        e2=AFunM(a0(j+1),epsilon(i));
        if (e1*e2<=0)
            x=epsilon(i); y=( a0(j)+a0(j+1) )/2;
            b=abs(beta0+(1/4)*beta2*y^2)*omega*y/sqrt((r1(x)*varsigma2(x)-1)^2+(mue*omega+r1(x)*varsigma1(x))^2);
            c=abs(eta)*omega^2*y/sqrt((r2*varrho2-omega^2)^2+(muP*omega+r2*varrho1)^2);
            P1maxM(jj,i)=b; P2maxM(jj,i)=c; AmpaM(jj,i)=y;
        end 
    end
    
    P1max(jj,i)=(P1maxP(jj,i)+P1maxM(jj,i))/2; P2max(jj,i)=(P2maxP(jj,i)+P2maxM(jj,i))/2; Ampa(jj,i)=(AmpaP(jj,i)+AmpaM(jj,i))/2;   %amplitudes
    P1max(jj,i)=P1max(jj,i)^2*omega^2; P2max(jj,i)=P2max(jj,i)^2;                                                                   %output powers
    P1max(jj,i)=P1max(jj,i)/( f0*Ampa(jj,i)*omega ) ;    P2max(jj,i)=P2max(jj,i)/( f0*Ampa(jj,i)*omega ) ;                          %efficiencies
    
end

display(jj)

end



[XX,YY] = meshgrid(epsilon,ttheta);


P1max = 0.2500000000e-3*P1max;  P2max = 0.2500000000e-3*P2max;
P1max=100*P1max;   P2max=100*P2max;


figure(1);hold on
mesh(XX,YY,P1max); xlabel('\epsilon'); ylabel('\theta'); zlabel('Efficiency');

figure(2);hold on
mesh(XX,YY,P2max); xlabel('\epsilon'); ylabel('\theta'); zlabel('Efficiency');




