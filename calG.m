function [g] = calG(x,t,e,alpha,beta,rou)
%% x,t,e,alpha,beta,rou
%km km/s g/cm3
n0=100; %parameter of integral; the bigger the more accurate
if x(3)~=0
   g=zeros(3,3);
   warning('x3~=0 not supported');return;
end
x1=x(1)-e(1);x2=x(2)-e(2);e3=e(3);
x1=x1*1000;
x2=x2*1000;
e3=e3*1000;
alpha=alpha*1000;
beta=beta*1000;
rou=rou*1000;
t0=t;
[pha,theta,r]=cart2sph(x1,x2,e3);
theta=pi/2-theta;
miu=beta^2*rou;

gp=zeros(3,3);
gs1=zeros(3,3);
gs2=zeros(3,3);
gs3=zeros(3,3);



t=t0;
if t>r/alpha;
    pL0=(t^2/r^2-1/alpha^2)^0.5;
    pL1=[0:n0]/n0*pL0/2;
    pL1R=(pL0^2-pL1.^2).^(0.5);
    pL2=[n0:2*n0]/n0*pL0/2;
    pL2R=-(pL0^2-pL2.^2).^0.5;
    p0=pL1(1);
    M0=calM(alpha,beta, r,theta,pha,t,p0)/pL1R(1);
    for i=2:length(pL1)
        p1=pL1(i);
        M1=calM(alpha,beta, r,theta,pha,t,p1)/pL1R(i);
        gp=gp+real((M1+M0)/2)*(p1-p0);
        p0=p1;M0=M1;
    end
    p0=pL2(1);pR0=pL2R(1);
    M0=calM(alpha,beta, r,theta,pha,t,p0)/p0;
    for i=2:length(pL2)
        p1=pL2(i);pR1=pL2R(i);
        M1=calM(alpha,beta, r,theta,pha,t,p1)/p1;
        gp=gp+real((M1+M0)/2)*(pR1-pR0);
        p0=p1;M0=M1;pR0=pR1;
    end
end

if t>r/beta;
    pL0=(t^2/r^2-1/beta^2)^0.5;
    pL1=[0:n0]/n0*pL0/2;
    pL1R=(pL0^2-pL1.^2).^(0.5);
    pL2=[n0:2*n0]/n0*pL0/2;
    pL2R=-(pL0^2-pL2.^2).^0.5;
    p0=pL1(1);
    N0=calN(alpha,beta, r,theta,pha,t,p0)/pL1R(1);
    for i=2:length(pL1)
        p1=pL1(i);
        N1=calN(alpha,beta, r,theta,pha,t,p1)/pL1R(i);
        gs1=gs1+real((N1+N0)/2)*(p1-p0);
        p0=p1;N0=N1;
    end
    p0=pL2(1);pR0=pL2R(1);
    N0=calN(alpha,beta, r,theta,pha,t,p0)/p0;
    for i=2:length(pL2)
        p1=pL2(i);pR1=pL2R(i);
        N1=calN(alpha,beta, r,theta,pha,t,p1)/p1;
        gs1=gs1+real((N1+N0)/2)*(pR1-pR0);
        p0=p1;N0=N1;pR0=pR1;
    end
end

if sin(theta)>beta/alpha
    t2=r/alpha*sin(theta)+r*cos(theta)*(1/beta^2-1/alpha^2)^0.5;
    p2=(((t/r-(1/beta^2-1/alpha^2)^0.5*cos(theta))/sin(theta))^2-1/alpha^2)^0.5;
    if t>t2 && t<r/beta
       pL0=(1/beta^2-t^2/r^2)^0.5;
       pL1=[0:n0]/n0*(p2);
       pL1R=(pL0^2+pL1.^2).^0.5;
       p0=pL1(1);
       N0=-calN(alpha,beta, r,theta,pha,t,p0)/pL1R(1);
       for i=2:length(pL1)
           p1=pL1(i);
           N1=-calN(alpha,beta,r,theta,pha,t,p1)/pL1R(i);
           gs2=gs2+imag((N1+N0)/2)*(p1-p0);
           p0=p1;N0=N1;
       end

    end
    if t>=r/beta
         pL0=(-1/beta^2+t^2/r^2)^0.5;
         D=p2-pL0;
         pL2=[0:n0]/n0*D/2+pL0;
         pL2R=(pL2.^2-pL0^2).^0.5;
         pL1=[0:n0]/n0*D/2+D/2+pL0;
         pL1R=(pL1.^2-pL0^2).^0.5;
         p0=pL1(1);
         N0=-calN(alpha,beta, r,theta,pha,t,p0)/pL1R(1);
         for i=2:length(pL1)
             p1=pL1(i);
             N1=-calN(alpha,beta, r,theta,pha,t,p1)/pL1R(i);
             gs3=gs3+imag((N1+N0)/2*(p1-p0));
             p0=p1;N0=N1;
         end
            p0=pL2(1);pR0=pL2R(1);
            N0=-calN(alpha,beta, r,theta,pha,t,p0)/p0;
         for i=2:length(pL2)
              p1=pL2(i);pR1=pL2R(i);
              N1=-calN(alpha,beta, r,theta,pha,t,p1)/p1;
              gs3=gs3+imag((N1+N0)/2*(pR1-pR0));
              p0=p1;N0=N1;pR0=pR1;
         end 
    end 
end


gs=gs1+gs2+gs3;
g=(gp+gs1+gs2+gs3)/(pi^2*r*miu);