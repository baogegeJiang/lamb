%x1=2;x2=0;t=[0:0.0025:1]*4;x3=10;
%alpha=8;beta=4.62;rou=3.300;
gP=[];
g11=[];g12=[];g13=[];g21=[];g22=[];g23=[];g31=[];g32=[];g33=[];
fclose all;
zDir=pwd;
inputFile=[zDir,'/Debug/input.dat'];
outputFile=[zDir,'/Debug/Gij.dat'];
exeDir=[zDir,'/Debug/'];
cd(exeDir);
system('main.exe >> null ');
cd(zDir);
input=readdata(inputFile);
input{1,1}=num2str(x1);
input{2,1}=num2str(x2);
input{3,1}=num2str(x3);
input{4,1}=num2str(e1);
input{5,1}=num2str(e2);
input{6,1}=num2str(e3);
input{7,1}=num2str(alpha);
input{8,1}=num2str(beta);
input{9,1}=num2str(rou);
input{10,1}=num2str(t(end));
input{11,1}=num2str(length(t));
wfile(input,inputFile);
fclose all;

figure(2);
output=readdata(outputFile);
[m,n]=size(output);
for i=1:m
    a=str2num(output{i,1})*0+1;
    g11(i)=str2num(output{i,2})*a;
    g12(i)=str2num(output{i,3})*a;
    g13(i)=str2num(output{i,4})*a;
    g21(i)=str2num(output{i,5})*a;
    g22(i)=str2num(output{i,6})*a;
    g23(i)=str2num(output{i,7})*a;
    g31(i)=str2num(output{i,8})*a;
    g32(i)=str2num(output{i,9})*a;
    g33(i)=str2num(output{i,10})*a;
end
t=t(1:m);
r=((x1-e1)^2+(x2-e2)^2+(x3-e3)^2)^0.5;
ta=r/alpha;
tb=r/beta;
clf
hold on
A=1.5*max(max(abs([ g11 g31 g22 g13 g33])));
plot(t,g11/A+5);
plot(t,g31/A+4);
plot(t,g22/A+3);
plot(t,g13/A+2);
plot(t,g33/A+1);
yL=[0:0.2:5.5];
ha=plot(ta+yL*0,yL,'.b');
hb=plot(tb+yL*0,yL,'.r');
legend([ha hb],{'P','S'});
%xlim(t([1 end]));
ylim([0 6.1]);
set(gca,'yTick',[1:5],'yTickLabel',{'g^H33','g^H13','g^H22','g^H31','g^H11'});
title(sprintf('G(%d,%d,0,t;0,0,%3.1f,0) A=%s Zhang',x1,x2,e3,num2str(A)));
xlabel('t/s');
print(gcf,'-djpeg','-r300',sprintf('%d_%3.1fZhang.jpg',x1,e3));
fclose all;

figure(fn)
%A0=1.5*max(max(max(max(abs([ g11 g31 g22 g13 g33])))));
hZ=plot(t,g11/A0+5,'r');
plot(t,g31/A0+4,'r');
plot(t,g22/A0+3,'r');
plot(t,g13/A0+2,'r');
plot(t,g33/A0+1,'r');
yL=[0:0.2:6.1];
%ha=plot(ta+yL*0,yL,'.b');
%hb=plot(tb+yL*0,yL,'.r');
legend([hJ hZ],{'Jiang','Zhang(ref)'});
%xlim(t([1 end]));
%ylim([0.5 6]);
%set(gca,'yTick',[1:5],'yTickLabel',{'g^H33','g^H13','g^H22','g^H31','g^H11'});
%title(sprintf('G(%d,%d,0,t;0,0,%d,0) A=%s',x1,x2,x3,num2str(A)));
%xlabel('t/s');
print(gcf,'-djpeg','-r300',sprintf('%d_%3.1fJiang_Zhang.jpg',x1,e3));

