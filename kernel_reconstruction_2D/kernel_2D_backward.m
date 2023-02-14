clear all
% We consider the 2-D kernel circumstance
load('kernel_forward_data_2D.mat','M','N','rho','vx','vy','g11','g12','g22','A','epsilon','dt','dx')


L=M/2;
rhohat=rho;
vxhat=vx;
vyhat=vy;

xdist=((1:L+1)'*ones(1,L+1)-1)*dx;
ydist=(ones(L+1,1)*(1:L+1)-1)*dx;
Kinit=exp(-(repmat(A(1,1),[L+1,L+1]).*xdist.^2+...
    repmat(A(1,2)+A(2,1),[L+1,L+1]).*xdist.*ydist+...
    repmat(A(2,2),[L+1,L+1]).*ydist.^2)/epsilon);

K=ones(L+1,L+1)*Kinit(1,1);

rho=ones(M,M,N+1)*rho(1,1,1);
vx=ones(M,M,N)*vx(1,1,1);
vy=ones(M,M,N)*vy(1,1,1);


psix=zeros(M,M,N);
psiy=zeros(M,M,N);
chi=zeros(M,M,N);
phi=zeros(M,M,N);
thetax=zeros(1,M,N);
thetay=zeros(M,1,N);



n=20000000;
tau=0.000001;
sigma=0.000001;

lambda=100;

alpha=1./sum(sum(sum(rhohat.^2/M/M/(N+1))))*lambda;
beta=1./sum(sum(sum(vxhat.^2+vyhat.^2)/M/M/N));
gamma=10^(-2)/dt;

p=2;

acpm=10^3;

alpha0=0;

rhocri=[];
vxcri=[];
vycri=[];
Kcri=[];


tic
for i=1:n
    % reserve the variable from last iteration
    rhoold=rho;
    vxold=vx;
    vyold=vy;
    Kold=K;
   
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
    
    wx=repmat(g11,[1,1,N]).*vx+repmat(g12,[1,1,N]).*vy;
    wy=repmat(g12,[1,1,N]).*vx+repmat(g22,[1,1,N]).*vy;
   
   
    derivrho(:,:,1:N)=-(deltaphix.*vx+deltaphiy.*vy);
    derivrho(:,:,2:N)=derivrho(:,:,2:N)+alpha*(rho(:,:,2:N)-rhohat(:,:,2:N))+...
        (phi(:,:,1:N-1)-phi(:,:,2:N))/dt+...
        conv_backward(K,deltapsixx+deltapsiyy,2,N,M,L,N,dx);
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
    
    for q1=1:L+1
        for q2=1:L+1
            if q1<=L && q1>=2 && q2<=L && q2>=2
                pmK(q1,q2,:)=sum(sum((deltapsixx+deltapsiyy).*...
                    (rhotrans(rho(:,:,1:N),q1-1,q2-1,M)+...
                    rhotrans(rho(:,:,1:N),-(q1-1),q2-1,M)+...
                    rhotrans(rho(:,:,1:N),q1-1,-(q2-1),M)+...
                    rhotrans(rho(:,:,1:N),-(q1-1),-(q2-1),M)),1),2)*dx^2;
            elseif mod(q1-1,L)==0 && mod(q2-1,L)==0
                pmK(q1,q2,:)=sum(sum((deltapsixx+deltapsiyy).*...
                    rhotrans(rho(:,:,1:N),q1-1,q2-1,M),1),2)*dx^2;
            elseif mod(q2-1,L)==0
                pmK(q1,q2,:)=sum(sum((deltapsixx+deltapsiyy).*...
                    (rhotrans(rho(:,:,1:N),q1-1,q2-1,M)+...
                    rhotrans(rho(:,:,1:N),-(q1-1),q2-1,M)),1),2)*dx^2;
            else
                pmK(q1,q2,:)=sum(sum((deltapsixx+deltapsiyy).*...
                    (rhotrans(rho(:,:,1:N),q1-1,q2-1,M)+...
                    rhotrans(rho(:,:,1:N),q1-1,-(q2-1),M)),1),2)*dx^2;
            end
        end
    end
    
    derivK=sum(pmK(:,:,2:N),3);
    
    derivK(1:L,:)=derivK(1:L,:)+gamma/2/dx^p*p*(...
        abs(K(2:L+1,:)-K(1:L,:)).^(p-1).*sign(K(1:L,:)-K(2:L+1,:)));
    derivK(:,1:L)=derivK(:,1:L)+gamma/2/dx^p*p*(...
        abs(K(:,2:L+1)-K(:,1:L)).^(p-1).*sign(K(:,1:L)-K(:,2:L+1)));
    derivK(2:L+1,:)=derivK(2:L+1,:)+gamma/2/dx^p*p*(...
        abs(K(2:L+1,:)-K(1:L,:)).^(p-1).*sign(K(2:L+1,:)-K(1:L,:)));
    derivK(:,2:L+1)=derivK(:,2:L+1)+gamma/2/dx^p*p*(...
        abs(K(:,2:L+1)-K(:,1:L)).^(p-1).*sign(K(:,2:L+1)-K(:,1:L)));
    
    rho=rho-tau*derivrho;
    vx=vx-tau*derivvx;
    vy=vy-tau*derivvy;
    K=K-tau*dt*derivK*acpm;
    K(1,1)=Kinit(1,1);
    K=max(K,0);
   
    
    % dual step
    rhostar=2*rho-rhoold;
    vxstar=2*vx-vxold;
    vystar=2*vy-vyold;
    Kstar=2*K-Kold;
    mxstar=rhostar(:,:,1:N).*vxstar;
    mystar=rhostar(:,:,1:N).*vystar;
    wxstar=repmat(g11,[1,1,N]).*vxstar+repmat(g12,[1,1,N]).*vystar;
    wystar=repmat(g12,[1,1,N]).*vxstar+repmat(g22,[1,1,N]).*vystar;
    mxstarx2=[mxstar(M,:,:);mxstar(1:M-1,:,:)];
    mystary2=[mystar(:,M,:),mystar(:,1:M-1,:)];
    wxstary1=[wxstar(:,2:M,:),wxstar(:,1,:)];
    wystarx1=[wxstar(2:M,:,:);wxstar(1,:,:)];
    xistar=-conv_backward(Kstar,rhostar,1,N,M,L,N,dx)+0.5*(...
        repmat(g11,[1,1,N]).*vxstar.^2+...
        repmat(g12,[1,1,N]).*vxstar.*vystar*2+...
        repmat(g22,[1,1,N]).*vystar.^2);
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
   
   
   
    if mod(i,15000)==0
        num=num2str(i);
        file_name=strcat('L2_gamma4',num,'.mat');
        save(file_name)
    end
   
   
    if mod(i,100)==0
        if isnan(sum(sum(sum(rho))))
            break
        end
    end
    if mod(i,1000)==0
        rhocri=[rhocri;sum(sum(sum((rho-rhoold).^2)))];
        vxcri=[vxcri;sum(sum(sum((vx-vxold).^2)))];
        vycri=[vycri;sum(sum(sum((vy-vyold).^2)))];
        Kcri=[Kcri;sum(sum((K-Kold).^2))];
    end
end
toc


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
    % c11
    subplot(3,3,7)
    mesh(K)
    subplot(3,3,8)
    mesh(Kinit)
    subplot(3,3,9)
    plot(Kcri)
end