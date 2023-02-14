clear all

lambda_opt=1;
gamma_opt=[10^(-6);10^(-5);10^(-4);10^(-3);10^(-2);10^(-1)];
lambda_index=1;
gamma_index=1;
save('kernel_control_file.mat')

while true
    clear all
    % We consider the 1-D circumstance
    load('kernel_forward_2x.mat','M','N','rho','v','g','K')
    load('kernel_control_file.mat')
    
    dx=1/M;
    dt=1/N;
    lambda_opt
    lambda=lambda_opt(lambda_index);
    gamma=gamma_opt(gamma_index)/dt;
    
    rhohat=rho;
    vhat=v;

    L=M/2;

    Kinit=K(1,1:L+1)';
    K=ones(L+1,1)*K(1,1);

    % epsilon=1;
    % 
    % rhohat=rhohat+sqrt(sum(sum(rhohat.^2))/M/N)*(rand(M,N+1)-0.5)*epsilon;
    % vhat=vhat+sqrt(sum(sum(vhat.^2))/M/N)*(rand(M,N)-0.5)*epsilon;

    n=3000000; % interation time
    tau=0.001; % iteration parameter for CP
    sigma=0.001; % iteration parameter for CP


    alpha=1/( sum(sum(rhohat.^2))/M/N )*lambda;
    alpha0=alpha/dt;
    beta=1/( sum(sum(vhat.^2))/M/N );



    rho=ones(M,N+1)*rho(1,1);


    rhocri=[];
    vcri=[];
    Kcri=[];


    phi=zeros(M,N);
    psi=zeros(M,N);
    chi=zeros(1,N);

    Kmatrix=zeros(M,M);
    for q=1:L+1
        if q==1
            Kmatrix=Kmatrix+diag(ones(M,1)*K(1));
        elseif q==L+1
            Kmatrix=Kmatrix+diag(ones(M-L,1)*K(q),q-1)+...
                diag(ones(M-L,1)*K(q),1-q);
        else
            Kmatrix=Kmatrix+diag(ones(M-q+1,1)*K(q),q-1)+...
                diag(ones(M-q+1,1)*K(q),1-q);
            Kmatrix=Kmatrix+diag(ones(q-1,1)*K(q),M-q+1)+...
                diag(ones(q-1,1)*K(q),-M+q-1);
        end
    end

    tic
    for i=1:n
        phiold=phi;
        psiold=psi;
        chiold=chi;
        Kold=K;

        rhoold=rho;
        vold=v;

        Kmatrixold=Kmatrix;

        phi1=[phi(2:M,:);phi(1,:)];
        psi2=[psi(M,:);psi(1:M-1,:)];

        derivrho(:,2:N)=alpha*(rho(:,2:N)-rhohat(:,2:N))+...
            (phi(:,1:N-1)-phi(:,2:N))/dt+(phi(:,2:N)-phi1(:,2:N))/dx.*v(:,2:N)+...
            Kmatrixold*(psi(:,2:N)-psi2(:,2:N));
        derivrho(:,1)=(alpha+alpha0)*(rho(:,1)-rhohat(:,1))-...
            phi(:,1)/dt+(phi(:,1)-phi1(:,1))/dx.*v(:,1);
        derivrho(:,N+1)=(alpha+alpha0)*(rho(:,N+1)-rhohat(:,N+1))+phi(:,N)/dt;

        derivv=beta*(v-vhat)+(phi-phi1)/dx.*rho(:,1:N)+...
            g*chi/dx;
        derivv(:,2:N-1)=derivv(:,2:N-1)+(psi2(:,2:N-1)-psi(:,2:N-1))/dx.*v(:,2:N-1).*repmat(g,[1,N-2])+...
            (psi(:,2:N-1)-psi(:,3:N))/dt.*repmat(g,[1,N-2]);
        derivv(:,1)=derivv(:,1)-psi(:,2)/dt.*g;
        derivv(:,N)=derivv(:,N)+(psi2(:,N)-psi(:,N))/dx.*g.*v(:,N)+...
            psi(:,N)/dt.*g;

        for q=2:L+1
            pmb(q)=sum(sum((psi(:,2:N)-psi2(:,2:N)).*...
                ([rho(q:M,2:N);rho(1:q-1,2:N)]+...
                [rho(M+2-q:M,2:N);rho(1:M+1-q,2:N)])));
        end
        pmb(L+1)=pmb(L+1)/2;
        K1=[K(2:L+1);K(1)];
        K2=[K(L+1);K(1:L)];
        derivK=gamma/dx/dx*(2*K-K1-K2);
        derivK(L+1)=gamma/dx/dx*(K(L+1)-K(L));
        derivK=derivK+pmb';

        rho=rho-tau*derivrho;
        v=v-tau*derivv;
        K(2:L+1)=K(2:L+1)-dt*tau*derivK(2:L+1);


        Kmatrix=zeros(M,M);
        for q=1:L+1
            if q==1
                Kmatrix=Kmatrix+diag(ones(M,1)*K(1));
            elseif q==L+1
                Kmatrix=Kmatrix+diag(ones(M-L,1)*K(q),q-1)+...
                    diag(ones(M-L,1)*K(q),1-q);
            else
                Kmatrix=Kmatrix+diag(ones(M-q+1,1)*K(q),q-1)+...
                    diag(ones(M-q+1,1)*K(q),1-q);
                Kmatrix=Kmatrix+diag(ones(q-1,1)*K(q),M-q+1)+...
                    diag(ones(q-1,1)*K(q),-M+q-1);
            end
        end


        % dual step
        rhostar=2*rho-rhoold;
        vstar=2*v-vold;
        Kstar=2*Kmatrix-Kmatrixold;
        mstar=rhostar(:,1:N).*vstar;
        wstar=repmat(g,[1,N]).*vstar;
        xistar=-Kstar*rhostar(:,1:N)*dx+0.5*repmat(g,[1,N]).*vstar.^2;
        mstar2=[mstar(M,:);mstar(1:M-1,:)];
        xistar1=[xistar(2:M,:);xistar(1,:)];

        phi=phi+sigma*((rhostar(:,2:N+1)-rhostar(:,1:N))/dt+(mstar-mstar2)/dx);
        psi(:,2:N)=psi(:,2:N)+sigma*((xistar1(:,2:N)-xistar(:,2:N))/dx+(wstar(:,2:N)-wstar(:,1:N-1))/dt);
        chi=chi+sigma*sum(wstar,1);



        if mod(i,10000)==0
            if isnan(sum(rhocri)) || isnan(sum(vcri)) || isnan(sum(Kcri))
                break
            end
            rhocri=[rhocri;sum(sum((rho-rhoold).^2))];
            vcri=[vcri;sum(sum((v-vold).^2))];
            Kcri=[Kcri;sum((K-Kold).^2)];
        end
    end
    toc

    Kinv=fliplr(K')';
    Kinitinv=fliplr(Kinit')';
    Kplot=[K;Kinv(2:L)];
    Kinitplot=[Kinit;Kinitinv(2:L)];
    
    lambda_text=num2str(lambda_index);
    gamma_text=num2str(gamma_index);
    fig_name=strcat('kernel_','lambda',lambda_text,'gamma',gamma_text,'.jpg');
    file_name=strcat('kernel_','lambda',lambda_text,'gamma',gamma_text,'.mat');
    
    save(file_name)
    

    figure
    subplot(2,3,1)
    plot(Kinitplot)
    hold on
    plot(Kplot)
    hold off
    subplot(2,3,2)
    mesh(rho)
    subplot(2,3,3)
    mesh(v)
    subplot(2,3,4);
    plot(rhocri)
    subplot(2,3,5);
    plot(vcri)
    subplot(2,3,6);
    plot(Kcri)

    saveas(gcf,fig_name)
    
    if lambda_index==length(lambda_opt) && gamma_index==length(gamma_opt)
        break
    end
    if gamma_index==length(gamma_opt)
        gamma_index=1;
        lambda_index=lambda_index+1;
    else
        gamma_index=gamma_index+1;
    end
    save('kernel_control_file.mat','lambda_opt','lambda_index','gamma_opt','gamma_index')
end