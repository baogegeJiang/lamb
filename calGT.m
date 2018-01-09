x1=10;x2=0;x3=0;t=[0:0.0025:1]*4;e1=0;e2=0;e3=0.2;fn=5;
alpha=8.000;beta=4.620;rou=3.300;
gP=[];g11=[];g12=[];g13=[];g21=[];g22=[];g23=[];g31=[];g32=[];g33=[];

tic
figure(1)
for i=1:length(t)
    g=calG([x1 x2 x3],t(i),[e1 e2 e3],alpha,beta,rou);
    gP(:,:,i)=g;
    g11(i)=g(1,1);
    g12(i)=g(1,2);
    g13(i)=g(1,3);
    g21(i)=g(2,1);
    g22(i)=g(2,2);
    g23(i)=g(2,3);
    g31(i)=g(3,1);
    g32(i)=g(3,2);
    g33(i)=g(3,3);
end
toc

r=((x1-e1)^2+(x2-e2)^2+(x3-e3)^2)^0.5;
ta=r/alpha;
tb=r/beta;
clf
hold on
A=1.5*max(max(max(max(abs([ g11 g31 g22 g13 g33])))));
plot(t,g11/A+5);
plot(t,g31/A+4);
plot(t,g22/A+3);
plot(t,g13/A+2);
plot(t,g33/A+1);
yL=[0:0.2:6.1];
ha=plot(ta+yL*0,yL,'.b');
hb=plot(tb+yL*0,yL,'.r');
legend([ha hb],{'P','S'});
xlim(t([1 end]));
ylim([0 6]);
set(gca,'yTick',[1:5],'yTickLabel',{'g^H33','g^H13','g^H22','g^H31','g^H11'});
title(sprintf('G(%d,%d,0,t;0,0,%3.1f,0) A=%s Jiang',x1,x2,e3,num2str(A)));
xlabel('t/s');
print(gcf,'-djpeg','-r300',sprintf('%d_%3.1fJiang.jpg',x1,e3));

figure(fn)
clf
hold on
A0=1.5*max(max(max(max(abs([ g11 g31 g22 g13 g33])))));
hJ=plot(t,g11/A0+5,'b');
plot(t,g31/A0+4,'b');
plot(t,g22/A0+3,'b');
plot(t,g13/A0+2,'b');
plot(t,g33/A0+1,'b');
yL=[0:0.2:6.1];
ha=plot(ta+yL*0,yL,'.b');
hb=plot(tb+yL*0,yL,'.r');
xlim(t([1 end]));
ylim([0 6]);
set(gca,'yTick',[1:5],'yTickLabel',{'g^H33','g^H13','g^H22','g^H31','g^H11'});
title(sprintf('G(%d,%d,0,t;0,0,%3.1f,0) A=%s J & Z',x1,x2,e3,num2str(A)));
xlabel('t/s');

calGZhang;