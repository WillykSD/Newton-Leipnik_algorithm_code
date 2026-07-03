%%%% Deterministic PZT model (discontinuous case): Approximate analytical solution 
%%%% (harmonic balance method)
%%%% Amplitude response Newton-Raphson method.   
%%%% Ouput power = f(omega,epsilon)

clear all
clc
global nu1 nu2 mue muP omega1 omega3 lambda1 lambda2 gamma0 gamma1 gamma2 beta0 beta1 beta2 delta eta f0 omega;

%% System parameters %%%%
nu1 = 0.3;  nu2 = 0.05;  mue = 1.0;  muP = 16.0;   omega1 = 0.50 ;  omega3 = 5.66; lambda1 = 0.0045 ;
lambda2 = 0.00002; gamma0 = -.125;  gamma1 = 0.002;  gamma2 = 0.020;  beta0 = 9.996; 
beta1 = -0.16  ;  beta2 = -1.6;  delta = 0.001;   eta = -64.0; f0=.02;  %omega=0.95;   

omegae = 50.0; omegaP = 282.843 ; 

%% Initial array and varying parameters
eepsilon=0.1:.1:1;   theta = 0.9;

ne=length(eepsilon);

 omeg1a = 0.0:0.005:0.5;
 omeg2a = 0.55:0.05:5.0;

% omega = 0.0:0.01:5.0;   % a0 = 0.0:0.00001:.05; omega = 0.0:0.005:5.0;
a0 = 0.0:0.00002:0.2; omega=[omeg1a omeg2a];
n=length(omega); m=length(a0);

P1maxP=zeros(ne,n);  P2maxP=zeros(ne,n); 
P1maxM=zeros(ne,n);  P2maxM=zeros(ne,n);
P1max=zeros(ne,n);   P2max=zeros(ne,n);


for jj=1:ne
    
epsilon = eepsilon(jj);    

r1 = omegae^(epsilon-1);   r2 = omegaP^2/omegae^(1+theta);

%% Other functions definition

varsigma1 =@(omega) ( omega^(1+epsilon)*cos((1/2)*pi*epsilon) ) ;     
varsigma2 =@(omega) ( omega^(1+epsilon)*sin((1/2)*pi*epsilon) ) ;

varrho1 =@(omega) ( omega^(1-theta)*cos((1/2)*pi*theta) ) ;    
varrho2 =@(omega) ( omega^(1-theta)*sin((1/2)*pi*theta) ) ;

R1 =@(omega) ( (1-r1*varsigma2(omega))/((1-r1*varsigma2(omega))^2+(mue*omega+r1*varsigma1(omega))^2) ) ;
R2 =@(omega) ( 1/4*(3*lambda2-2*lambda1*gamma2/gamma1) ) ;
R3 =@(omega) ( omega1^2-omega^2+delta*eta*omega^2*(r2*varrho2(omega)-omega^2)/((r2*varrho2(omega)-omega^2)^2+(muP*omega+r2*varrho1(omega))^2) ) ;
R4p =@(omega) ( -delta*eta*omega^2*(muP*omega+r2*varrho1(omega))/((r2*varrho2(omega)-omega^2)^2+(muP*omega+r2*varrho1(omega))^2)+(nu1+nu2)*omega ) ; %%%%
R4m =@(omega) ( -delta*eta*omega^2*(muP*omega+r2*varrho1(omega))/((r2*varrho2(omega)-omega^2)^2+(muP*omega+r2*varrho1(omega))^2)+(nu1-nu2)*omega ) ; %%%%
R5 =@(omega) ( -(mue*omega+r1*varsigma1(omega))/((1-r1*varsigma2(omega))^2+(mue*omega+r1*varsigma1(omega))^2) ) ;
R6 =@(omega) ( omega^2*(beta2*gamma0+gamma2*beta0) ) ;

%% Amplitude funtions definition
AFunP=@(a0,omega) ( a0^2*((((1/4)*omega*gamma2*a0^2+gamma0*omega)*((1/4)*omega*beta2*a0^2+beta0*omega)*R1(omega)+R2(omega)*a0^2+R3(omega))^2+(R4p(omega)+((1/4)*omega*gamma2*a0^2+gamma0*omega)*((1/4)*omega*beta2*a0^2+beta0*omega)*R5(omega))^2)-f0^2 ) ;
AFunM=@(a0,omega) ( a0^2*((((1/4)*omega*gamma2*a0^2+gamma0*omega)*((1/4)*omega*beta2*a0^2+beta0*omega)*R1(omega)+R2(omega)*a0^2+R3(omega))^2+(R4m(omega)+((1/4)*omega*gamma2*a0^2+gamma0*omega)*((1/4)*omega*beta2*a0^2+beta0*omega)*R5(omega))^2)-f0^2 ) ;


%% Newton-Raphson method
for i=1:n
    for j=1:m-1
        %%
        e1=AFunP(a0(j),omega(i));
        e2=AFunP(a0(j+1),omega(i));
        if (e1*e2<=0)
            x=omega(i); y=( a0(j)+a0(j+1) )/2;
            b=abs(beta0+(1/4)*beta2*y^2)*x*y/sqrt((r1*varsigma2(x)-1)^2+(mue*x+r1*varsigma1(x))^2);
            c=abs(eta)*x^2*y/sqrt((r2*varrho2(x)-x^2)^2+(muP*x+r2*varrho1(x))^2);
            P1maxP(jj,i)=b^2*x^2; P2maxP(jj,i)=c^2; 
        end
        %
        %
        e1=AFunM(a0(j),omega(i));
        e2=AFunM(a0(j+1),omega(i));
        if (e1*e2<=0)
            x=omega(i); y=( a0(j)+a0(j+1) )/2;
            b=abs(beta0+(1/4)*beta2*y^2)*x*y/sqrt((r1*varsigma2(x)-1)^2+(mue*x+r1*varsigma1(x))^2);
            c=abs(eta)*x^2*y/sqrt((r2*varrho2(x)-x^2)^2+(muP*x+r2*varrho1(x))^2);
            P1maxM(jj,i)=b^2*x^2; P2maxM(jj,i)=c^2; 
        end 
    end
    
    P1max(jj,i)=(P1maxP(jj,i)+P1maxM(jj,i))/2; P2max(jj,i)=(P2maxP(jj,i)+P2maxM(jj,i))/2;
    
end

display(jj)

end



[XX,YY] = meshgrid(omega,eepsilon);


figure(1);hold on
mesh(XX,YY,P1max); xlabel('\omega'); ylabel('\epsilon'); zlabel('P_{1max}');

figure(2);hold on
mesh(XX,YY,P2max); xlabel('\omega'); ylabel('\epsilon'); zlabel('P_{2max}');

