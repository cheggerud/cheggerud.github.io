
% things I want to edit
%make its so alphas get finer near the bifurc (save run time:).
clear all;
% Parameters
tic
L=100;
n=100;xstep=L/(n+1);%space
T=10000;m=1000;tstep=T/(m+1);%time
%diffusion of u
% global Unorms Vnorms
% here we make all the same except alphaU. 
D1=1;
D2=1;
alpha=0.1;
advUvec=linspace(-1,1,51);
% advUvec=0.1;
advVvec=linspace(-1,1,51);


d1=0.001; %death rate of u
d2=0.001;%death rate of v


i0=1;

i01=200;
i02=200;
i03=200;
KB=0.001;
k1=0.001;
k2=0.001;
ku1=0.001;
ku2=0.00001;
ku3=0.001;
kv1=0.00001;
kv2=0.001;
kv3=0.001;

m1=1;
a1=10;
m2=1;
a2=10;

numparamsteps=length(advVvec);
xaxis=advVvec;% update x axis for plotting. 

numparamsteps2=length(advUvec);

STOREU=zeros(n,m,floor(numparamsteps/2));
STOREV=zeros(n,m,floor(numparamsteps/2));
STOREUS=zeros(n,m,floor(numparamsteps/2));
STOREVS=zeros(n,m,floor(numparamsteps/2));

parfor kk=1:numparamsteps2

advU=advUvec(kk);
    
    
UU=zeros(n,numparamsteps);
VV=zeros(n,numparamsteps);
UUU=zeros(1,numparamsteps);
VVV=zeros(1,numparamsteps);
UUS=zeros(n,numparamsteps);
VVS=zeros(n,numparamsteps);
UUUS=zeros(1,numparamsteps);
VVVS=zeros(1,numparamsteps);



for ii=1:numparamsteps %alpha param loop:
    
    
    
    %update the parameter (here we update alpha1. )
    
    advV=advVvec(ii);
    
    
    %store every other solution: 
  
  
  
%     if floor(ii/2)==ii/2
%        
%         STOREU(:,:,ii/2)=Theta(1:n,:); STOREV(:,:,ii/2)=Theta(n+1:2*n,:);
%             STOREUS(:,:,ii/2)=ThetaS(1:n,:); STOREVS(:,:,ii/2)=ThetaS(n+1:2*n,:);
%     end
    
% advU=Alphamin+ii*alphastep;%advectionrate of u
% advV=Alphamin+ii*alphastep;%advectionrate of v
u=arrayfun(@Iniu,(xstep:xstep:L-xstep)./L)';%Initial cond of u % arrayfun makes an array evaled at each x in the array
us=arrayfun(@Iniu,(xstep:xstep:L-xstep)./L)';%Initial cond of single u
v=arrayfun(@Iniv,(xstep:xstep:L-xstep)./L)';%Initial cond of v
vs=arrayfun(@Iniv,(xstep:xstep:L-xstep)./L)';%Initial cond of single v

% U=zeros(n,m);%u的积分
% US=zeros(n,m);%single u的积分
% V=zeros(n,m);%v的积分
% VS=zeros(n,m);% single v的积分
Theta=zeros(2*n,m);%解 of (u, v)
ThetaS=zeros(2*n,m);%解 of single (u,v)
A1=Lap(n,L,advU,D1)./xstep^2;%-Laplace/h^2 of u
B1=ADV(n,L,advU,D1)./(2*xstep);%advection of u
A2=Lap(n,L,advV,D2)./xstep^2;%-Laplace/h^2 of v
B2=ADV(n,L,advV,D2)./(2*xstep);%advection of v

for i=1:m % through time!
    
G1=zeros(1,n); G2=zeros(1,n);%Initialize the G1 and G2 function vectors.
G1S=zeros(1,n);G2S=zeros(1,n);
for j=1:n %hrough space
   %parameters of nonlocal


I1=i01.*exp(-KB.*j*xstep-ku1.*sum(u(1:j)).*xstep-kv1.*sum(v(1:j)).*xstep);
I2=i02.*exp(-KB.*j*xstep-ku2.*sum(u(1:j)).*xstep-kv2.*sum(v(1:j)).*xstep);
I3=i03.*exp(-KB.*j*xstep-ku3.*sum(u(1:j)).*xstep-kv3.*sum(v(1:j)).*xstep);

I1us=i01.*exp(-KB.*j*xstep-ku1.*sum(us(1:j)).*xstep);
I1vs=i01.*exp(-KB.*j*xstep-kv1.*sum(vs(1:j)).*xstep);


I2us=i02.*exp(-KB.*j*xstep-ku2.*sum(us(1:j)).*xstep);
I2vs=i02.*exp(-KB.*j*xstep-kv2.*sum(vs(1:j)).*xstep);

I3us=i03.*exp(-KB.*j*xstep-ku3.*sum(us(1:j)).*xstep);
I3vs=i03.*exp(-KB.*j*xstep-kv3.*sum(vs(1:j)).*xstep);
 Iusplit=ku1*I1+ku2*I2+ku3*I3;
 Ivsplit=kv1*I1+kv2*I2+kv3*I3;
 Iussplit=ku1*I1us+ku2*I2us+ku3*I3us;
 Ivssplit=kv1*I1vs+kv2*I2vs+kv3*I3vs;
 
 

%  Iu=i0.*exp(-k0.*j*xstep-k1.*sum(u(1:j)).*xstep-k2.*sum(v(1:j)).*xstep);
%  Ius=i0.*exp(-k0.*j*xstep-k1.*sum(us(1:j)).*xstep);
%  Ivs=i0.*exp(-k0.*j*xstep-k2.*sum(vs(1:j)).*xstep);
%  
 
 %checking why the light split blows up
 
%  IUdiffer(j,i)=Iu-Iusplit;
%  IVdiffer(j,i)=Iu-Ivsplit;
%   IUSdiffer(j,i)=Ius-Iussplit;
%  IVSdiffer(j,i)=Ivs-Ivssplit;
%U(j,i)=sum(u(1:j)).*xstep;
%V(j,i)=sum(v(1:j)).*xstep;%计算u,v积分
%US(j,i)=sum(us(1:j)).*xstep;
%VS(j,i)=sum(vs(1:j)).*xstep;%计算u,v积分
%I=i0.*exp(-k0.*j*h-theta(j));


  Lit1=m1.*Iusplit./(Iusplit+a1);Lit2=m2.*Ivsplit./(Ivsplit+a2);%g1(I) 和g2(I)
 Lit1S=m1.*Iussplit./(Iussplit+a1);Lit2S=m2.*Ivssplit./(Ivssplit+a2);

%  Lit1=m1.*Iu./(Iu+a1);Lit2=m2.*Iu./(Iu+a2);%g1(I) 和g2(I)
%  Lit1S=m1.*Ius./(Ius+a1);Lit2S=m2.*Ivs./(Ivs+a2);%g1(I) 和g2(I)


G1(j)=Lit1;G2(j)=Lit2;
G1S(j)=Lit1S;G2S(j)=Lit2S;
end
u=(eye(n)+tstep.*D1.*A1+tstep.*advU.*B1-tstep.*diag(G1)+tstep.*d1.*eye(n))\u;%implicit finite diff schemes u_{t+1}
us=(eye(n)+tstep.*D1.*A1+tstep.*advU.*B1-tstep.*diag(G1S)+tstep.*d1.*eye(n))\us;% single u
v=(eye(n)+tstep.*D2.*A2+tstep.*advV.*B2-tstep.*diag(G2)+tstep.*d2.*eye(n))\v;%v_{t+1}
vs=(eye(n)+tstep.*D2.*A2+tstep.*advV.*B2-tstep.*diag(G2S)+tstep.*d2.*eye(n))\vs;%single v
Theta(:,i)=[u;v];
ThetaS(:,i)=[us;vs];
end
UU(:,ii)=Theta(1:n,m);VV(:,ii)=Theta(n+1:2*n,m);
UUS(:,ii)=ThetaS(1:n,m);VVS(:,ii)=ThetaS(n+1:2*n,m);
UUU(ii)=sum(UU(:,ii)).*xstep;VVV(ii)=sum(VV(:,ii)).*xstep;
UUUS(ii)=sum(UUS(:,ii)).*xstep;VVVS(ii)=sum(VVS(:,ii)).*xstep;


Unorms(ii,kk)=UUU(ii);
Vnorms(ii,kk)=VVV(ii);


end%alpha param loop:



end
toc
% save  Theta 
% save  ThetaS 
% save  UUU 
% save  VVV 
% save  UUUS 
% save  VVVS
save Unorms
save Vnorms
save advVvec 
save advUvec

% 
% 
% coextinct=zeros(size(Vnorms));
% coexist=zeros(size(Vnorms));
% justU=zeros(size(Vnorms));
% justV=zeros(size(Vnorms));
% comat=zeros(size(Vnorms));
% [n,m]=size(Vnorms);
% for ii=1:m
%     
%     for kk=1:n
%        
%         if Vnorms(kk,ii)<0.000001 & Unorms(kk,ii)<0.000001
%         coextinct(kk,ii)=1;
%         comat(kk,ii)=0;
%         
%         end
%         if Vnorms(kk,ii)>0.000001 & Unorms(kk,ii)>0.000001
%             coexist(kk,ii)=1;
%             comat(kk,ii)=1;
%         end
%          if Vnorms(kk,ii)>0.000001 & Unorms(kk,ii)<0.000001
%             justV(kk,ii)=1;
%             comat(kk,ii)=2;
%          end
%            if Vnorms(kk,ii)<0.000001 & Unorms(kk,ii)>0.000001
%             justU(kk,ii)=1;
%             comat(kk,ii)=3;
%         end
%     end
%     
%     
% end
% figure
% s=pcolor(advVvec,advUvec,comat)
% 
% 
% 
% %Theta(1:n, :)代表的是u(x,t);Theta(n+1:2*n,:)代表的是v(x,t)
% %  figure
% % % 
% % % % surf(VV)
% % % hold on
% % % % alphaa=alphastep:alphastep:Alphamax;
% % % % plot(alphaa,UUU./UUUS,alphaa,VVV./VVVS,'LineWidth',2)
% % 
% % 
% %  plot(xaxis,UUtwoU./UUUS,'b',xaxis,VVV./VVVS,'r','LineWidth',2)
% %  %plot(alphavec,UUU/max(UUU),alphavec,VVV/max(VVV),'LineWidth',2)
% % %set(gcf,'DefaultTextInterpreter','tex')
% % % legend('$\|u^*\| /\|\tilde{u}\|$','$\|v^*\|/ \|\tilde{v}\|$')
% 
% legend({'$\|u^*\| /\|\tilde{u}\|$','$\|v^*\|/ \|\tilde{v}\|$'},'Interpreter','latex')
%  xlabel('advection rate $\alpha$','Interpreter','latex')  
% 
% 

%figure(2)
%alphaa=dalpha:dalpha:Alphamax;
%plot(alphaa,UUU,alphaa,VVV,alphaa,UUUS,alphaa,VVVS,'LineWidth',2)
%legend('\itu', '\itv','\itu^*', '\itv^*')
%xlabel('\alpha')  
%ylabel('L_1 norm')
% 
% figure(2)
% hold on
% plot(alphaa,VVV,'b')
% 
% tplot=tstep:tstep:T-tstep;
% xplot=xstep:xstep:L-xstep;
% figure
% surf(tstep:tstep:T-tstep,xstep:xstep:L-xstep,Theta(1:n,:),'LineStyle','none')%u(x,t)  
% xlabel('Time t')
% ylabel('Depth x')
% title('u(x,t)')
% figure
% surf(tstep:tstep:T-tstep,xstep:xstep:L-xstep,Theta(n+1:2*n,:),'LineStyle','none')%v(x,t)
% 
% 
%  
% xlabel('Time t')
% ylabel('Depth x')
% title('v(x,t)')
%figure(3)
%plot(h:h:L-h,Theta(1:n,500),'r',h:h:L-h,Theta(n+1:2*n,500),'b')%u(x,88),v(x,88)

%figure(4)
%plot(h:h:L-h,Theta(1:n,m),'r',h:h:L-h,Theta(n+1:2*n,m),'b','LineWidth',2)%u(x,88),v(x,88)
%legend('\itu', '\itv')
%xlabel('Depth x')  
%ylabel('Population density')

%caxis([0, max(Theta(1:n,m))])



%(u,v)初值1
function u=Iniu(x)
u=cos(pi*x)+2;
end

%(u,v)初值1
function v=Iniv(x)
v=cos(pi*x)+2;

end




%function v=Iniv(x)
%if x<=1/2
%v=sin(pi.*x);
%else
%v=1;
%end





function advs=ADV(n,L,q,D)
h=L/(n+1);
k=2*h*q/D;
A=[k,0;0,k];
for i=3:n
    D=A;
    C=zeros(i);
    C(1:i-1,1:i-1)=D;
    C(i-1,i-2:i)=[-1,0,1];
    C(i,i-1:i)=[0,k];
    A=C;
end
advs=A;

end

function lap=Lap(n,L,q,D)
h=L/(n+1); %xstep;
k=2*h*q/D;
A=[2+k,-2;-2,2-k]; %this is for BC
for i=3:n
    D=A;
    C=zeros(i);
    C(1:i-1,1:i-1)=D;
    C(i-1,i-2:i)=[-1,2,-1];
    C(i,i-1:i)=[-2,2-k];
    A=C;
end
lap=A;
end



