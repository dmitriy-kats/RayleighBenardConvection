%Solves Rayleigh Bernard Convection in domain of H by L with dimensionless number Ra
%Nx and Ny specifies the number of grid points
clc;clear;close all;

H=0.5;
L=1;
Ra=5E4;
Nx=80;
Ny=80;


Nu=0;
NuT=[0 Nu];
dx=L/Nx;
dy=H/(Ny-2);
dt=1/Ra;
ds=min(1/4*H*Ra^(-1/2),1/16*H*Ra^(-1/3));


unl=zeros(Nx*Ny,1)*rand(1)^5;
vnl=zeros(Nx*(Ny-1),1)*rand(1)^2;
Tnl=ones(Nx*Ny,1);
for y=1:Ny-1
    Tnl(Nx*y+1:Nx*(y+1))=(H-dy/2-dy*(y-1))*ones(Nx,1);
end

[X,Yd] = meshgrid([dx/2:dx:L-dx/2],[-dy/2:dy:L+dy/2]);
Y=flipud(Yd);


for tt=1:50000
A1=spalloc(Nx*Ny,Nx*Ny,Nx*(Ny-2)*5+Nx*4); 
b1=zeros(Nx*Ny,1);
count=1;
for j=Ny-1:-1:2
for i=1:Nx
  P=count+Nx;
  w=P;
  if i==Nx-1
  e=P+1;
  ee=P+2-Nx;
  elseif i==Nx
  e=P+1-Nx;
  ee=e+1;
  else
  e=P+1;
  ee=e+1;
  end
  en=e+Nx;
  es=e-Nx;
  A1(count,w)=-unl(e)/2/dx-1/dx^2;
  A1(count,e)=dx*dy/dt+2/dx^2+2/dy^2;
  A1(count,ee)=unl(e)/2/dx-1/dx^2;
  vavgtemp=1/4*(vnl(w)+vnl(e)+vnl(w-Nx)+vnl(e-Nx));
  A1(count,en)=vavgtemp/(2*dy)-1/dy^2;
  A1(count,es)=-vavgtemp/(2*dy)-1/dy^2;
  b1(count)=dx*dy*unl(e)/dt;
  count=count+1;  
end
end

for ii=1:Nx
    A1(count,ii)=1;
    A1(count,ii+Nx)=1;
    count=count+1;
    A1(count,Nx*Ny-ii+1)=1;
    A1(count,Nx*Ny-ii+1-Nx)=1;
    count=count+1;
end

usl=A1\b1;
clearvars P w e ee en es count

A2=spalloc(Nx*Ny-Nx,Nx*Ny-Nx,5*(Ny-3)*Nx+Nx);
b2=zeros(Nx*Ny-Nx,1);
count=1;

for j=Ny-1:-1:3
for i=1:Nx
  P=count+Nx;
  n=P;
  nn=P+Nx;
  s=P-Nx;
  ne=P+1;
  nw=P-1;
 if i==Nx
  ne=P+1-Nx;
 elseif i==1
  nw=P+Nx-1;
 end
  e=ne;
  w=P;
  A2(count,s)=-1/dy^2-vnl(n)/(2*dy);
  A2(count,n)=dx*dy/dt+2/dx^2+2/dy^2;
  if nn<Nx*Ny-Nx+1
  A2(count,nn)=vnl(n)/(2*dy)-1/dy^2;
  end
  uavgtemp=1/4*(unl(w)+unl(e)+unl(w+Nx)+unl(e+Nx));
  A2(count,ne)=uavgtemp/(2*dx)-1/dx^2;
  A2(count,nw)=-uavgtemp/(2*dx)-1/dx^2;
  b2(count)=dx*dy*vnl(n)/dt+Ra/2*(Tnl(n)+Tnl(nn));
  count=count+1;
end
end

for ii=1:Nx
    A2(count,ii)=1;
    count=count+1;
    A2(count,Nx*Ny-ii+1-Nx)=1;
    count=count+1;
end
vsl=A2\b2;
clearvars P n nn s ne nw count

A3=spalloc(Nx*Ny,Nx*Ny,5*Nx*(Ny-2)+4*Nx-1);
b3=zeros(Nx*Ny,1);
count=1;

for j=Ny-1:-1:2
for i=1:Nx
  P=count+Nx;
  N=P+Nx;
  S=P-Nx;
  E=P+1;
  W=P-1;
  if i==1
      W=P+Nx-1;
  elseif i==Nx
      E=P+1-Nx;
  end
  w=P;
  e=E;
  n=P;
  s=S;
  if P==round(Ny/1.5*Nx)
  A3(count,P)=1;
  b3(count)=0;
  else
  A3(count,E)=1/dx^2;
  A3(count,P)=-2/dx^2-2/dy^2;
  A3(count,W)=1/dx^2;
  A3(count,N)=1/dy^2;
  A3(count,S)=1/dy^2;
  b3(count)=1/dt*((usl(e)-usl(w))/dx+(vsl(n)-vsl(s))/dy); 
  end
  count=count+1;
end
end

for ii=1:Nx
    A3(count,ii)=1;
    A3(count,ii+Nx)=-1;
    count=count+1;
    A3(count,Nx*Ny-ii+1)=1;
    A3(count,Nx*Ny-ii+1-Nx)=-1;
    count=count+1;
end

Pnpo=A3\b3;
clearvars P E W N S e w n s count


PR=reshape(1:Nx*Ny,[Ny,Nx]);
PRt=PR(1,:);
PR=[PR(2:end,:); PRt];

rtemp=reshape(PR,[1,Nx*Ny]);

for i=1:length(Pnpo)
    gradpx=(Pnpo(rtemp(i))-Pnpo(i))/dx;
    unl(i)=usl(i)-dt*gradpx;
    if i<=Nx*Ny-Nx
    gradpy=(Pnpo(i+Nx)-Pnpo(i))/dy;
    vnl(i)=vsl(i)-dt*gradpy;
    end
end


A4=spalloc(Nx*Ny,Nx*Ny,5*Nx*(Ny-2)+4*Nx);
b4=zeros(Nx*Ny,1);
count=1;

for j=Ny-1:-1:2
for i=1:Nx
  P=count+Nx;
  N=P+Nx;
  S=P-Nx;
  E=P+1;
  W=P-1;
  if i==1
      W=P+Nx-1;
  elseif i==Nx
      E=P+1-Nx;
  end
  w=P;
  e=E;
  n=P;
  s=S;
  uavg=1/2*(unl(w)+unl(e));
  vavg=1/2*(vnl(n)+vnl(s));
  A4(count,E)=uavg/(2*dx)-1/dx^2;
  A4(count,P)=1/dt+2/dx^2+2/dy^2;
  A4(count,W)=-uavg/(2*dx)-1/dx^2;
  A4(count,N)=vavg/(2*dy)-1/dy^2;
  A4(count,S)=-vavg/(2*dy)-1/dy^2;
  b4(count)=1/dt*Tnl(P); 
  count=count+1;
end
end

for ii=1:Nx
    A4(count,ii)=1;
    A4(count,ii+Nx)=1;
    b4(count)=2;
    count=count+1;
    A4(count,Nx*Ny-ii+1)=1;
    A4(count,Nx*Ny-ii+1-Nx)=1;
    count=count+1;
end

Tnl=A4\b4;
clearvars P E W N S e w n s count
Nu=0;
count=1;
for j=Ny-1:-1:2
for i=1:Nx
  P=count+Nx;
  n=P;
  s=P-Nx;
  N=P+Nx;
  S=s;
  VV=(vnl(n)+vnl(s))/2;
  Nu=Nu+dx*dy*(VV*Tnl(P)-(Tnl(N)-Tnl(S))/2/dy);
end
end
Nu=Nu/(L*H)
NuT=[NuT; tt*dt Nu];
% Tsurf=surf(X,Y,flipud(reshape(Tnl,[Nx,Ny])'));
% xlabel('x')
% ylabel('y')
% zlabel('T')
% direction = [0 1 0];
% rotate(Tsurf,direction,150);
if mod(tt,5)==0
drawnow
%imagesc(X(1,:),Y(:,1),flipud(reshape(Tnl,[Nx,Ny])'));
subplot(2,1,1)
pcolor([dx/2:dx:L-dx/2],H+dy/2:-dy:-dy/2,flipud(reshape(Tnl,[Nx,Ny])'));
shading interp
title(['Ra ' num2str(Ra,'%10.5e\n') '  Time is: ' num2str(tt*dt)])
xlabel('X')
ylabel('Y')
colorbar
subplot(2,1,2)
plot(NuT(:,1),NuT(:,2))
xlabel('Time')
ylabel('Nu')
end
end


