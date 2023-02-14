clear all
% We consider the 1-D circumstance
l=1; % total length
t=1; % total time
M=50; % the division of length
N=30; % the division of time
n=100000; % interation time
sigma=0.002; % iteration parameter for CP
tau=0.001; % iteration parameter for CP

epsilon=0.1;

% double grid
M=2*M;
N=2*N;

dx=l/M;
dt=t/N;

m=zeros(M,N);
phi=zeros(M,N);

% the initail of rho
rho0=1+sin(2*pi*(0:M-1)/M).^2*0.5;
rhoT=1+sin(pi*(0:M-1)/M).^2*0.5;
rho0=rho0/sum(rho0);
rhoT=rhoT/sum(rhoT);
rho=ones(M,N+1)*rho0(1);
rho0=rho0'/dx;
rhoT=rhoT'/dx;
rho(:,1)=rho0;
rho(:,N+1)=rhoT;


g=(1-sin(pi*(0:M-1)/M).^2*0.6)';

K=abs(ones(M,1)*(1:M)-(ones(M,1)*(1:M))')*dx;
K=min(K,1-K);
K=exp(-K.*K/epsilon);

ftilde=@(x)x;
ftildeinv=@(x)x;

rhocri=[];
mcri=[];

tic
for i=1:n
    % resersve the value of rho and m
    rhoold=rho;
    mold=m;
    
    phi1=[phi(2:M,:);phi(1,:)];
    
    mderiv=1./rho(:,1:N).*(g*ones(1,N)).*m+(phi-phi1)/dx;
    rhoderiv(:,2:N)=-0.5*(g*ones(1,N-1)).*m(:,2:N).^2./(rho(:,2:N).^2)+...
        K*rho(:,2:N)*dx+(phi(:,1:N-1)-phi(:,2:N))/dt;
    
    m=m-tau*mderiv;
    rho(:,2:N)=rho(:,2:N)-tau*rhoderiv(:,2:N);
    
    % update phi
    mstar=2*m-mold;
    rhostar=2*rho-rhoold;
    mstar2=[mstar(M,:);mstar(1:M-1,:)];
    phi=phi+sigma*((rhostar(:,2:N+1)-rhostar(:,1:N))/dt+...
        (mstar-mstar2)/dx);

    rhocri=[rhocri;sum(sum((rho-rhoold).^2))];
    mcri=[mcri;sum(sum((m-mold).^2))];
end
toc
phi1=[phi(2:M,:);phi(1,:)];
v=m./rho(:,1:N);



figure
subplot(2,2,1)
mesh(rho)
subplot(2,2,2)
plot(rhocri)
subplot(2,2,3)
mesh(m)
subplot(2,2,4)
plot(mcri)


% figure
% subplot(1,2,1)
% mesh(rho_trans)
% subplot(1,2,2)
% mesh(w_trans)

rho=rho(1:2:M,1:2:N+1);
v=v(1:2:M,1:2:N);
g=g(1:2:M);

M=M/2;
N=N/2;

dx=l/M;
dt=t/N;

K=abs(ones(M,1)*(1:M)-(ones(M,1)*(1:M))')*dx;
K=min(K,1-K);
K=exp(-K.*K/epsilon);

