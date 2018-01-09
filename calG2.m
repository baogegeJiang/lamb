function [g] = calG(x1,x2,t,x3,alpha,beta,rou)
%% x1,x2,t,x3.alpha,beta,rou
x1=x1*1000;
x2=x2*1000;
x3=x3*1000;
alpha=alpha*1000;
beta=beta*1000;
rou=rou*1000;
t0=t;
[pha,theta,r]=cart2sph(x1,x2,x3);
theta=pi/2-theta;
miu=beta^2*rou;
%gp
gp=zeros(3,3);
gs1=zeros(3,3);
gs2=zeros(3,3);
gs3=zeros(3,3);
n0=400;


t=t0;
if t>r/alpha;
    pL0=(t^2/r^2-1/alpha^2)^0.5;
    pL1=[0:n0]/n0*pL0/2;
    pL1R=(pL0^2-pL1.^2).^(0.5);
    pL2=[n0:2*n0]/n0*pL0/2;
    pL2R=-(pL0^2-pL2.^2).^0.5;
    % pL1
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


if sin(theta)>beta/alpha
    t2=r/alpha*sin(theta)+r*cos(theta)*(1/beta^2-1/alpha^2)^0.5;
    p2=(((t/r-(1/beta^2-1/alpha^2)^0.5*cos(theta))/sin(theta))^2-1/alpha^2)^0.5;
else
    t2=r/beta;
    p2=((t/r)^2-1/beta^2)^0.5;
end
%fprintf('%f\n',t2)
if t>real(t2)
   if sin(theta)<beta/alpha
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
          gs1=gs1+real((N1+N0)/2*(p1-p0));
          p0=p1;N0=N1;
      end
      p0=pL2(1);pR0=pL2R(1);
      N0=calN(alpha,beta, r,theta,pha,t,p0)/p0;
      for i=2:length(pL2)
          p1=pL2(i);pR1=pL2R(i);
          N1=calN(alpha,beta, r,theta,pha,t,p1)/p1;
          gs1=gs1+real((N1+N0)/2*(pR1-pR0));
          p0=p1;N0=N1;pR0=pR1;
      end
   else
       if t<=r/beta
          
          pL0=(t^2/r^2-1/beta^2)^0.5;
          pL1=[0:n0]/n0*p2/2;
          pL1R=(pL0^2-pL1.^2).^(0.5);
          pL2=[n0:2*n0]/n0*p2/2;
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
              gs1=gs1+real((N1+N0)/2*(pR1-pR0));
      
              p0=p1;N0=N1;pR0=pR1;
          end 
       else
           
          pL0=(t^2/r^2-1/beta^2)^0.5;
          p3=(t^2/r^2-1/beta^2)^0.5;
         % fprintf('%f %f\n',p3,p2)
          pL1=[0:n0]/n0*p3/2;
          pL1R=(pL0^2-pL1.^2).^(0.5);
          pL2=[0:n0]/n0*(p2-p3/2)+p3/2;
          pL2R=-(pL0^2-pL2.^2).^0.5;
          p0=pL1(1);
          N0=calN(alpha,beta, r,theta,pha,t,p0)/pL1R(1);
          for i=2:length(pL1)
              p1=pL1(i);
              N1=calN(alpha,beta, r,theta,pha,t,p1)/pL1R(i);
              gs1=gs1+real((N1+N0)/2*(p1-p0));
              p0=p1;N0=N1;
          end
          p0=pL2(1);pR0=pL2R(1);
          N0=calN(alpha,beta, r,theta,pha,t,p0)/p0;
          for i=2:length(pL2)
              p1=pL2(i);pR1=pL2R(i);
              N1=calN(alpha,beta, r,theta,pha,t,p1)/p1;
              gs1=gs1+real((N1+N0)/2*(pR1-pR0));
              p0=p1;N0=N1;pR0=pR1;
          end 
       end
   end
end
    


gs=gs1+gs2+gs3;
g=(gp+gs1+gs2+gs3)/(pi^2*r*miu);
%g=(gp+gs)/(pi^2*r*miu);