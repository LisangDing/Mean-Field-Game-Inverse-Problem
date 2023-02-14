clear all
% We consider the 2-D circumstance
% g11=x+4
% g12=x+2
% g22=2x+1
load('forward_data_2D.mat','M','N','rho','vx','vy','cbase','l','t')

dx=l/M;
dt=t/N;

% The feasible observations
rhohat=rho;
vxhat=vx;
vyhat=vy;

% The ground truth
chat=cbase;

rho=ones(M,M,N+1)*rho(1,1,1);
vx=ones(M,M,N)*vx(1,1,1);
vy=ones(M,M,N)*vy(1,1,1);
cbase(2:M,:)=cbase(1,1)*ones(M-1,M);

% Dual variables
psix=zeros(M,M,N);
psiy=zeros(M,M,N);
chi=zeros(M,M,N);
phi=zeros(M,M,N);
thetax=zeros(1,M,N);
thetay=zeros(M,1,N);

% functions: how cbase determines g11, g12, g22
ftildepm=@(x) 1;
ftilde=@(x) x;

f11=@(x) x+4;
f12=@(x) x+2;
f22=@(x) 2*x+1;
f11pm=@(x) 1;
f12pm=@(x) 1;
f22pm=@(x) 2;

% The iteration time and step size setting
n=500000;
tau=0.000001;
sigma=0.000001;

lambda=1000;

% Set the parameters, it might should be adjusted when dealing with
% different datas.
alpha=1./sum(sum(sum(rhohat.^2/M/M/(N+1))))*lambda;
beta=1./sum(sum(sum(vxhat.^2+vyhat.^2)/M/M/N));
gamma=10^(-2)/dt;


p=2;

% Use accelerating parameter (acpm) to change the step size for updating 
% cbase.
acpm=100;

alpha0=0;

% Statistical quantities
rhocri=[];
vxcri=[];
vycri=[];
cbasecri=[];


for i=1:n
    % reserve the variable from last iteration
    rhoold=rho;
    vxold=vx;
    vyold=vy;
    cbaseold=cbase;
    
    phix1=[phi(2:M,:,:);phi(1,:,:)];
    phiy1=[phi(:,2:M,:),phi(:,1,:)];
    
    deltaphix=(phix1-phi)/dx;
    deltaphiy=(phiy1-phi)/dx;
    
    psixx2=[psix(M,:,:);psix(1:M-1,:,:)];
    psiyy2=[psiy(:,M,:),psiy(:,1:M-1,:)];
    
    chix2=[chi(M,:,:);chi(1:M-1,:,:)];
    chiy2=[chi(:,M,:),chi(:,1:M-1,:)];
    
    deltapsixx=(psix-psixx2)/dx;
    deltapsiyy=(psiy-psiyy2)/dx;
    
    g11=f11(cbase);
    g12=f12(cbase);
    g22=f22(cbase);
    
    wx=repmat(g11,[1,1,N]).*vx+repmat(g12,[1,1,N]).*vy;
    wy=repmat(g12,[1,1,N]).*vx+repmat(g22,[1,1,N]).*vy;
    
    derivrho(:,:,1:N)=-(deltaphix.*vx+deltaphiy.*vy);
    derivrho(:,:,2:N)=derivrho(:,:,2:N)+alpha*(rho(:,:,2:N)-rhohat(:,:,2:N))+...
        (phi(:,:,1:N-1)-phi(:,:,2:N))/dt+(deltapsixx(:,:,2:N)+deltapsiyy(:,:,2:N)).*ftildepm(rho(:,:,2:N));
    derivrho(:,:,1)=derivrho(:,:,1)+(alpha+alpha0)*(rho(:,:,1)-rhohat(:,:,1))-phi(:,:,1)/dt;
    derivrho(:,:,N+1)=(alpha+alpha0)*(rho(:,:,N+1)-rhohat(:,:,N+1))+phi(:,:,N)/dt;
    
    gpsix=repmat(g11,[1,1,N]).*psix+repmat(g12,[1,1,N]).*psiy;
    gpsiy=repmat(g12,[1,1,N]).*psix+repmat(g22,[1,1,N]).*psiy;
    
    derivvx=beta*(vx-vxhat)-rho(:,:,1:N).*deltaphix+...
        (chiy2-chi)/dx.*repmat(g11,[1,1,N])-(chix2-chi)/dx.*repmat(g12,[1,1,N])+...
        repmat(thetax,[M,1,1])/dx.*repmat(g11,[1,1,N])+...
        repmat(thetay,[1,M,1])/dx.*repmat(g12,[1,1,N]);
    derivvx(:,:,2:N)=derivvx(:,:,2:N)-(deltapsixx(:,:,2:N)+deltapsiyy(:,:,2:N)).*wx(:,:,2:N);
    derivvx(:,:,2:N-1)=derivvx(:,:,2:N-1)+(gpsix(:,:,2:N-1)-gpsix(:,:,3:N))/dt;
    derivvx(:,:,1)=derivvx(:,:,1)-gpsix(:,:,2)/dt;
    derivvx(:,:,N)=derivvx(:,:,N)+gpsix(:,:,N)/dt;
    
    
    derivvy=beta*(vy-vyhat)-rho(:,:,1:N).*deltaphiy+...
        (chiy2-chi)/dx.*repmat(g12,[1,1,N])-(chix2-chi)/dx.*repmat(g22,[1,1,N])+...
        repmat(thetax,[M,1,1])/dx.*repmat(g12,[1,1,N])+...
        repmat(thetay,[1,M,1])/dx.*repmat(g22,[1,1,N]);
    derivvy(:,:,2:N)=derivvy(:,:,2:N)-(deltapsixx(:,:,2:N)+deltapsiyy(:,:,2:N)).*wy(:,:,2:N);
    derivvy(:,:,2:N-1)=derivvy(:,:,2:N-1)+(gpsiy(:,:,2:N-1)-gpsiy(:,:,3:N))/dt;
    derivvy(:,:,1)=derivvy(:,:,1)-gpsiy(:,:,2)/dt;
    derivvy(:,:,N)=derivvy(:,:,N)+gpsiy(:,:,N)/dt;
    
    
    
    pmA11=-0.5*(deltapsixx+deltapsiyy).*vx.^2;
    pmA12=-0.5*(deltapsixx+deltapsiyy).*vx.*vy*2;
    pmA22=-0.5*(deltapsixx+deltapsiyy).*vy.^2;
    pmB11(:,:,2:N-1)=(psix(:,:,2:N-1)-psix(:,:,3:N))/dt.*vx(:,:,2:N-1);
    pmB12(:,:,2:N-1)=(psix(:,:,2:N-1)-psix(:,:,3:N))/dt.*vy(:,:,2:N-1)+...
        (psiy(:,:,2:N-1)-psiy(:,:,3:N))/dt.*vx(:,:,2:N-1);
    pmB13(:,:,2:N-1)=(psiy(:,:,2:N-1)-psiy(:,:,3:N))/dt.*vy(:,:,2:N-1);
    pmB11(:,:,N)=psix(:,:,N)/dt.*vx(:,:,N);
    pmB12(:,:,N)=psix(:,:,N)/dt.*vy(:,:,N)+psiy(:,:,N)/dt.*vx(:,:,N);
    pmB22(:,:,N)=psiy(:,:,N)/dt.*vy(:,:,N);
    pmB11(:,:,1)=-psix(:,:,2)/dt.*vx(:,:,1);
    pmB12(:,:,1)=-psix(:,:,2)/dt.*vy(:,:,1)-psiy(:,:,2)/dt.*vx(:,:,1);
    pmB22(:,:,1)=-psiy(:,:,2)/dt.*vy(:,:,1);
    pmC11=((chiy2-chi)/dx+repmat(thetax,[M,1,1])/dx).*vx;
    pmC12=((chiy2-chi)/dx+repmat(thetax,[M,1,1])/dx).*vy+...
        ((chi-chix2)/dx+repmat(thetay,[1,M,1])/dx).*vx;
    pmC22=((chi-chix2)/dx+repmat(thetay,[1,M,1])/dx).*vy;
    
    
    cbasex1=[cbase(2:M,:,:);cbase(1,:,:)];
    cbasex2=[cbase(M,:,:);cbase(1:M-1,:,:)];
    cbasey1=[cbase(:,2:M,:),cbase(:,1,:)];
    cbasey2=[cbase(:,M,:),cbase(:,1:M-1,:)];
    
    
    derivcbase=gamma/2/dx^p*p*...
        (abs(cbase-cbasex1).^(p-1).*sign(cbase-cbasex1)+...
        abs(cbase-cbasex2).^(p-1).*sign(cbase-cbasex2)+...
        abs(cbase-cbasey1).^(p-1).*sign(cbase-cbasey1)+...
        abs(cbase-cbasey2).^(p-1).*sign(cbase-cbasey2))+...
        (sum(pmA11(:,:,2:N),3)+sum(pmB11,3)+sum(pmC11,3)).*f11pm(cbase)+...
        (sum(pmA12(:,:,2:N),3)+sum(pmB12,3)+sum(pmC12,3)).*f12pm(cbase)+...
        (sum(pmA22(:,:,2:N),3)+sum(pmB22,3)+sum(pmC22,3)).*f22pm(cbase);
    
    rho=rho-tau*derivrho;
    vx=vx-tau*derivvx;
    vy=vy-tau*derivvy;
    
    cbase(2:M,:)=cbase(2:M,:)-dt*tau*derivcbase(2:M,:)*acpm;
    
    
    % dual step
    rhostar=2*rho-rhoold;
    vxstar=2*vx-vxold;
    vystar=2*vy-vyold;
    cbasestar=2*cbase-cbaseold;
    g11star=f11(cbasestar);
    g12star=f12(cbasestar);
    g22star=f22(cbasestar);
    mxstar=rhostar(:,:,1:N).*vxstar;
    mystar=rhostar(:,:,1:N).*vystar;
    wxstar=repmat(g11star,[1,1,N]).*vxstar+repmat(g12star,[1,1,N]).*vystar;
    wystar=repmat(g12star,[1,1,N]).*vxstar+repmat(g22star,[1,1,N]).*vystar;
    mxstarx2=[mxstar(M,:,:);mxstar(1:M-1,:,:)];
    mystary2=[mystar(:,M,:),mystar(:,1:M-1,:)];
    wxstary1=[wxstar(:,2:M,:),wxstar(:,1,:)];
    wystarx1=[wxstar(2:M,:,:);wxstar(1,:,:)];
    xistar=-ftilde(rhostar(:,:,1:N))+0.5*(...
        repmat(g11star,[1,1,N]).*vxstar.^2+...
        repmat(g12star,[1,1,N]).*vxstar.*vystar*2+...
        repmat(g22star,[1,1,N]).*vystar.^2);
    xistarx1=[xistar(2:M,:,:);xistar(1,:,:)];
    xistary1=[xistar(:,2:M,:),xistar(:,1,:)];
    
    
    phi=phi+sigma*((rhostar(:,:,2:N+1)-rhostar(:,:,1:N))/dt+...
        (mxstar-mxstarx2+mystar-mystary2)/dx);
    psix(:,:,2:N)=psix(:,:,2:N)+sigma*(...
        (xistarx1(:,:,2:N)-xistar(:,:,2:N))/dx+...
        (wxstar(:,:,2:N)-wxstar(:,:,1:N-1))/dt );
    psiy(:,:,2:N)=psiy(:,:,2:N)+sigma*(...
        (xistary1(:,:,2:N)-xistar(:,:,2:N))/dx+...
        (wystar(:,:,2:N)-wystar(:,:,1:N-1))/dt );
    chi=chi+sigma*((wxstary1-wxstar)/dx-(wystarx1-wystar)/dx);
    thetax=thetax+sigma*sum(wxstar,1);
    thetay=thetay+sigma*sum(wystar,2);
    
    % Check whether the iteration has blowed up, if so, stop.
    if mod(i,100)==0
        if isnan(sum(sum(sum(rho))))
            break
        end
    end
    if mod(i,100)==1
        cbasecri=[cbasecri; sqrt( sum(sum((cbase-chat).^2))/M^2 ) ];
    end
end


% The final result is 'cbase', compare it with the ground truth as
% follows.
figure
for k=1:N+1
    subplot(3,3,1)
    surf(rho(:,:,k))
    zlim([min(min(min(rho))) max(max(max(rho)))])
    subplot(3,3,2)
    plot(rhocri)
    subplot(3,3,3)
    surf(vx(:,:,min(k,N)))
    zlim([min(min(min(vx))) max(max(max(vx)))])
    subplot(3,3,4)
    plot(vxcri)
    subplot(3,3,5)
    surf(vy(:,:,min(k,N)))
    zlim([min(min(min(vy))) max(max(max(vy)))])
    pause(0.05)
    subplot(3,3,6)
    plot(vycri)
    subplot(3,3,7)
    mesh(cbase)
    subplot(3,3,8)
    mesh(chat)
    subplot(3,3,9)
    plot(cbasecri)
end
