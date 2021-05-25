close all;
clear;

syms t;
k=10.*10.^3;
m=331;
c=50;
v=(65.*1.60934.*1000)./3600;
Y=.05;
zeta=c./(2.*sqrt(k.*m));
wb=v.*2.*pi;
wn=sqrt(k./m);
r=wb./wn;
X=Y.*sqrt((1+(2.*zeta.*r).^2)/((1-r.^2).^2+(2.*zeta.*r).^2))
theta1=atan((2.*zeta.*wn.*wb)./(wn.^2-wb.^2))
theta2=atan((wn.^2)./(wn.^2-wb.^2))
xp=X.*cos(wb.*t-theta1-theta2);

figure(1)
fplot(xp,[0,5])