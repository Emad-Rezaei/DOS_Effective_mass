%Developed by Emad
% This code reads TransDOS, differential conductivity, and calculates
% effective masses of electrons and holes and prints out the maximum
% effective mass ratio between two quantities. 
clear all
clc
disp('This is to approximate effective mass around Fermi level through DOS')
cod=input('what code is employed to calculate transport properties? type 1 or 2 or 3 for BoltzTrap1, BoltzTrap2 and Boltzwann, respectively')
qe=1.6*1e-16;
vol=270.2565;
conv=(1/1.8897)^3*1e-30;
vol=vol*conv;
%%%%Opening DOS%%%%
fid=fopen('GaAs_boltzdos.dat','r');
if fid == -1
    disp ('file name is wrong')
else
    if cod==1
        ff=textscan(fid,'%f %f %f'); %BoltzTrap1
    elseif cod==2
        ff=textscan(fid,'%f %f','Headerlines',5); %Boltzdos
    elseif cod==3
        ff=textscan(fid,'%f %f %f %f %f %f %f','Headerlines',5); %Bolttzwann
    else
        ff=textscan(fid,'%f %f %f','Headerlines',1); %QE DOS
    end
end
fclose(fid)
c=size(ff{1},1);
f=zeros(c-1,2);
if cod>1
  Ef=input('please enter the intrinsic Fermi level');  
end
for x=2:c
    if cod==3
    f(x-1,1)=ff{1}(x,1)-Ef; %Boltwann
    f(x-1,2)=ff{2}(x,1);
    elseif cod==1
    %Ef=0;
    %f(x-1,1)=13.605698066*ff{1}(x,1)-Ef; %BTZP1
    f(x-1,1)=13.605698066*ff{1}(x,1)-Ef;
    f(x-1,2)=ff{3}(x,1);
    else 
    %Ef=0;
    f(x-1,1)=ff{1}(x,1)-Ef;
    f(x-1,2)=ff{2}(x,1);
    end
end
%%%%% Defining effective mass estimation range on DOS%%%%%%% 
a=input('how far is final point in valence band? with a negative sign'); %negative sign has to be typed!
b=input('how far is final point in conduction band? with a positive sign'); %positive sign has to be typed!
aa=input('how far is initial point in valence band? with a negative sign');
bb=input('how far is initial point in conduction band? with a positive sign');
dd=zeros(c,2);
ee=zeros(c,2);
g=0;
h=0;
hbar=6.63*10^-34;
me0=9.1094*10^-31;
cc=size(f,1);
for i=1:cc;
    if a<=f(i,1) & f(i,1)<=b
       if bb<=f(i,1)
           g=g+1;
           dd(g,2)=f(i,2);
           dd(g,1)=(abs(f(i,1)))^(1/2);
       elseif f(i,1)<=aa
           h=h+1;
           ee(h,2)=f(i,2);
           ee(h,1)=(abs(f(i,1)))^(1/2);
       end
    end
end
d=zeros(g,2);
e=zeros(h,2);
for v=1:g
    for w=1:2
        for z=1:h
        d(v,w)=dd(v,w);
        e(z,w)=ee(z,w);
        end
    end
end
k=zeros(2,2);
l=0;
m=0;
mc=0;
ymc=0;
xymc=0;
s=zeros(2,1);
for j=1:g;
    l=l+d(j,1);
    m=m+d(j,1).*d(j,1);
    ymc=ymc+d(j,2);
    xymc=xymc+d(j,1).*d(j,2);
end
k(1,1)=g;
k(1,2)=l;
k(2,1)=l;
k(2,2)=m;
s(1,1)=ymc;
s(2,1)=xymc;
r=inv(k)*s;
mc=((hbar^3*abs(r(2,1))*pi^2)/sqrt(2))^(2/3);
mc=mc/(vol^(2/3)*qe);
mc=mc/me0;
disp('electron effective mass is:')
disp(mc)
o=zeros(2,2);
n=0;
q=0;
mv=0;
ymv=0;
xymv=0;
t=zeros(2,1);
for p=1:h;
    n=n+e(p,1);
    q=q+e(p,1).*e(p,1);
    ymv=ymv+e(p,2);
    xymv=xymv+e(p,1).*e(p,2);
end
o(1,1)=h;
o(1,2)=n;
o(2,1)=n;
o(2,2)=q;
t(1,1)=ymv;
t(2,1)=xymv;
u=inv(o)*t;
mv=((hbar^3*abs(u(2,1))*pi^2)/sqrt(2))^(2/3);
mv=mv/(vol^(2/3)*qe);
mv=mv/me0;
disp('hole effective mass is:')
disp(mv)
rmass= abs(r(2,1))/abs(u(2,1));
if rmass<1
    rmass=inv(rmass);
end
disp('mass ratio is')
disp(rmass)
disp(abs(mv/mc))
subplot(3,1,1);
plot(d(:,1),d(:,2),'*',d(:,1),r(2,1).*d(:,1)+r(1,1),':')
title('Electrons')
xlabel('E^{3/2}')
ylabel('DOS(1/ev)')
subplot(3,1,2);
plot(e(:,1),e(:,2),'*',e(:,1),u(2,1).*e(:,1)+u(1,1),':')
title('Holes')
xlabel('E^{3/2}')
ylabel('DOS(1/ev)')
subplot(3,1,3)
plot(f(:,1),f(:,2))
title('NormalDOS')
xlabel('E-E_{f}(ev)')
ylabel('DOS(1/ev)')
xlim([-2 2]);
