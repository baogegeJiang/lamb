function[M]=calM(alpha,beta,r,theta,pha,t,p)
%% alpha,beta,r,theta,pha,t,p
M=zeros(3,3);ii=(-1)^0.5;
%cal q
tauAlpha0=r*(1/alpha^2+p^2)^0.5;
if t<tauAlpha0
    q=1/r*(-t*sin(theta)+(r^2/alpha^2+r^2*p^2-t^2)^0.5*cos(theta));
else
    q=1/r*(-t*sin(theta)+ii*(-(r^2/alpha^2+r^2*p^2)+t^2)^0.5*cos(theta));
end

etaAlpha=(alpha^(-2)+p^2-q^2)^0.5;
etaBeta=(beta^(-2)+p^2-q^2)^0.5;
gamma=etaBeta^2+p^2-q^2;
sigma=gamma^2+4*etaAlpha*etaBeta*(q^2-p^2);

M(1,1)=2*etaBeta*((q^2+p^2)*cos(pha)^2-p^2);
M(1,2)=2*etaBeta*(q^2+p^2)*sin(pha)*cos(pha);
M(1,3)=2*q*etaBeta*etaAlpha*cos(pha);
M(2,1)=M(1,2);
M(2,2)=2*etaBeta*((q^2+p^2)*sin(pha)^2-p^2);
M(2,3)=2*q*etaAlpha*etaBeta*sin(pha);
M(3,1)=q*gamma*cos(pha);
M(3,2)=q*gamma*sin(pha);
M(3,3)=etaAlpha*gamma;

%involve  eta sigma

M=etaAlpha*M/sigma;
