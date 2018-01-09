function[N]=calN(alpha,beta,r,theta,pha,t,p)
%% alpha,beta,r,theta,pha,t,p
N=zeros(3,3);ii=(-1)^0.5;
%cal q
tauBeta0=r*(1/beta^2+p^2)^0.5;
if t<=tauBeta0
    q=1/r*(-t*sin(theta)+(r^2/beta^2+r^2*p^2-t^2)^0.5*cos(theta));
else
    q=1/r*(-t*sin(theta)+ii*(-(r^2/beta^2+r^2*p^2)+t^2)^0.5*cos(theta));
    
end


etaAlpha=(alpha^(-2)+p^2-q^2)^0.5;
etaBeta=(beta^(-2)+p^2-q^2)^0.5;
gamma=etaBeta^2+p^2-q^2;
sigma=gamma^2+4*etaAlpha*etaBeta*(q^2-p^2);

N(1,1)=1/etaBeta*(etaBeta^2*gamma-(gamma-4*etaAlpha*etaBeta)*...
    ((q^2+p^2)*sin(pha)^2-p^2));
N(1,2)=1/etaBeta*(q^2+p^2)*(gamma-4*etaAlpha*etaBeta)*sin(pha)...
    *cos(pha);
N(1,3)=-q*gamma*cos(pha);
N(2,1)=N(1,2);
N(2,2)=1/etaBeta*(etaBeta^2*gamma-(gamma-4*etaAlpha*etaBeta)*...
    ((p^2+q^2)*cos(pha)^2-p^2));
N(2,3)=-q*gamma*sin(pha);
N(3,1)=-2*q*etaAlpha*etaBeta*cos(pha);
N(3,2)=-2*q*etaAlpha*etaBeta*sin(pha);
N(3,3)=2*etaAlpha*(q^2-p^2);

%involve  eta sigma

N=etaBeta*N/sigma;