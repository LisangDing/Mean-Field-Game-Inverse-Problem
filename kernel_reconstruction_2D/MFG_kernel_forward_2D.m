clear all
% We consider the 2-D circumstance
l=1; % total length
t=1; % total time
M=24; % the division of length
N=30; % the division of time
n=2000000; % interation time
sigma=0.0000002; % iteration parameter for CP
tau=0.0000001; % iteration parameter for CP
dx=l/M;
dt=t/N;

acpm=1000;

mx=zeros(M,M,N);
my=zeros(M,M,N);
phi=zeros(M,M,N);

% the initailization of rho
rho0=1+sin(2*pi*((0:M-1)+0.5)/M).^2*0.5;
rhoT=1+sin(pi*((0:M-1)+0.5)/M).^2*0.5;
rho0=rho0/sum(rho0);
rhoT=rhoT/sum(rhoT);
rho=ones(M,M,N+1)*rho0(1)/dx;
rho0=rho0';
rhoT=rhoT';
rho(:,:,1)=rho0*rho0'/dx/dx;
rho(:,:,N+1)=rhoT*rhoT'/dx/dx;

% The ground metric is known
c=(1-sin(pi*(0:M-1)/M).^2*0.6)';
cbase=c*c';
g11=cbase+4;
g12=cbase+2;
g22=cbase*2+1;

% adaptive matrix
A=[3,1;1,2];
epsilon=0.5;

% Statistical quantities
rhocri=[];
mxcri=[];
mycri=[];

tic
for i=1:n
    % resersve the value of rho and m
    rhoold=rho;
    mxold=mx;
    myold=my;
   
   
    phix1=[phi(2:M,:,:);phi(1,:,:)];
    phiy1=[phi(:,2:M,:),phi(:,1,:)];
   
    mxderiv=1./rho(:,:,1:N).*(repmat(g11,[1,1,N]).*mx+...
        repmat(g12,[1,1,N]).*my)+(phi-phix1)/dx;
    myderiv=1./rho(:,:,1:N).*(repmat(g12,[1,1,N]).*mx+...
        repmat(g22,[1,1,N]).*my)+(phi-phiy1)/dx;
    rhoderiv(:,:,2:N)=-0.5*(repmat(g11,[1,1,N-1]).*mx(:,:,2:N).*mx(:,:,2:N)+...
        2*repmat(g12,[1,1,N-1]).*mx(:,:,2:N).*my(:,:,2:N)+...
        repmat(g22,[1,1,N-1]).*my(:,:,2:N).*my(:,:,2:N))./...
        rho(:,:,2:N).^2+conv_forward(A,epsilon,rho,2,N,M,N,dx)+...
        (phi(:,:,1:N-1)-phi(:,:,2:N))/dt;
   
    mx=mx-acpm*tau*mxderiv;
    my=my-acpm*tau*myderiv;
    rho(:,:,2:N)=rho(:,:,2:N)-acpm*tau*rhoderiv(:,:,2:N);
   
    % update phi
    mxstar=2*mx-mxold;
    mystar=2*my-myold;
    rhostar=2*rho-rhoold;
    mxstar2=[mxstar(M,:,:);mxstar(1:M-1,:,:)];
    mystar2=[mystar(:,M,:),mystar(:,1:M-1,:)];
    phi=phi+sigma*((rhostar(:,:,2:N+1)-rhostar(:,:,1:N))/dt+...
        (mxstar-mxstar2)/dx+(mystar-mystar2)/dx);

    if mod(i,10000)==0
        rhocri=[rhocri;sum(sum(sum((rho-rhoold).^2)))];
        mxcri=[mxcri;sum(sum(sum((mx-mxold).^2)))];
        mycri=[mycri;sum(sum(sum((my-myold).^2)))];
        if isnan(sum(rhocri)) || isnan(sum(mxcri))  || isnan(sum(mycri))
            break
        end
    end
end
toc

phix1=[phi(2:M,:,:);phi(1,:,:)];
phiy1=[phi(:,2:M,:),phi(:,1,:)];
wx=(phix1-phi)/dx;
wy=(phiy1-phi)/dx;
vx=mx./rho(:,:,1:N);
vy=my./rho(:,:,1:N);
wx_trans=repmat(g11,[1,1,N]).*vx+repmat(g12,[1,1,N]).*vy;
wy_trans=repmat(g12,[1,1,N]).*vx+repmat(g22,[1,1,N]).*vy;

% plot the forward result
figure
for k=1:N+1
    subplot(3,3,1)
    surf(rho(:,:,k))
    zlim([min(min(min(rho))) max(max(max(rho)))])
    pause(0.05)
    subplot(3,3,2)
    plot(rhocri)
    subplot(3,3,3)
    surf(mx(:,:,min(k,N)))
    zlim([min(min(min(mx))) max(max(max(mx)))])
    pause(0.05)
    subplot(3,3,4)
    plot(mxcri)
    subplot(3,3,5)
    surf(my(:,:,min(k,N)))
    zlim([min(min(min(my))) max(max(max(my)))])
    pause(0.05)
    subplot(3,3,6)
    plot(mycri)
    subplot(3,3,7)
    surf(wx(:,:,min(k,N)))
    zlim([min(min(min(wx))) max(max(max(wx)))])
    pause(0.05)
    subplot(3,3,8)
    surf(wy(:,:,min(k,N)))
    zlim([min(min(min(wy))) max(max(max(wy)))])
    pause(0.05)
end

% Save the datas
save('kernel_forward_data_2D.mat','M','N','rho','vx','vy','g11','g12',...
    'g22','A','epsilon','dt','dx')
