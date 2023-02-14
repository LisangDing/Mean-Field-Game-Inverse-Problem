clear all

% try different lambda and gamma as the parameters
lambda_opt=1;
gamma_opt=[10^(-8);10^(-7);10^(-6);10^(-5);10^(-4);10^(-3)];
lambda_index=1;
gamma_index=1;
save('L2_control_file.mat')

% Use the while loop try out all the feasible lambda and gamma set before
while true
    clear all
    % We consider the 1-D circumstance
    load('forward_data_1D.mat','M','N','rho','v','g','dt','dx')
    load('L2_control_file.mat')
    
    lambda=lambda_opt(lambda_index);
    gamma=gamma_opt(gamma_index)/dt;
    
    
    rhohat=rho;
    vhat=v;
    ghat=g;


    n=60000; % interation time
    tau=0.002; % iteration parameter for CP
    sigma=0.001; % iteration parameter for CP

    alpha=1/( sum(sum(rhohat.^2))/M/N )*lambda;
    alpha0=alpha/dt;
    beta=1/( sum(sum(vhat.^2))/M/N );


    g=ones(M,1)*ghat(1);

    rho=ones(M,N+1)*rho(1,1);


    rhocri=[];
    vcri=[];
    gcri=[];


    phi=zeros(M,N);
    psi=zeros(M,N);
    chi=zeros(1,N);


    ftildepm=@(x) 1;
    ftilde=@(x)x;

    tic
    for i=1:n
        phiold=phi;
        psiold=psi;
        chiold=chi;

        rhoold=rho;
        vold=v;
        gold=g;

        phi1=[phi(2:M,:);phi(1,:)];
        psi2=[psi(M,:);psi(1:M-1,:)];

        derivrho(:,2:N)=alpha*(rho(:,2:N)-rhohat(:,2:N))+...
            (phi(:,1:N-1)-phi(:,2:N))/dt+(phi(:,2:N)-phi1(:,2:N))/dx.*v(:,2:N)+...
            (psi2(:,2:N)-psi(:,2:N))/dx.*(-ftildepm(rho(:,2:N)));
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

        g1=[g(2:M);g(1)];
        g2=[g(M);g(1:M-1)];

        pm1=0.5*(psi2-psi)/dx.*v.^2;
        pm2(:,2:N-1)=(psi(:,2:N-1)-psi(:,3:N))/dt.*v(:,2:N-1);
        pm2(:,1)=-psi(:,2)/dt.*v(:,1);
        pm2(:,N)=psi(:,N)/dt.*v(:,N);

        derivg=gamma/dx/dx*(2*g-g1-g2)+sum(pm1(:,2:N),2)+...
            sum(pm2,2)+sum(repmat(chi,[M,1])/dx.*v,2);

        rho=rho-tau*derivrho;
        v=v-tau*derivv;
        g(2:M)=g(2:M)-tau*dt*derivg(2:M);

        % dual step
        rhostar=2*rho-rhoold;
        vstar=2*v-vold;
        gstar=2*g-gold;
        mstar=rhostar(:,1:N).*vstar;
        wstar=repmat(gstar,[1,N]).*vstar;
        xistar=-ftilde(rhostar(:,1:N))+0.5*repmat(gstar,[1,N]).*vstar.^2;
        mstar2=[mstar(M,:);mstar(1:M-1,:)];
        xistar1=[xistar(2:M,:);xistar(1,:)];

        phi=phi+sigma*((rhostar(:,2:N+1)-rhostar(:,1:N))/dt+(mstar-mstar2)/dx);
        psi(:,2:N)=psi(:,2:N)+sigma*((xistar1(:,2:N)-xistar(:,2:N))/dx+(wstar(:,2:N)-wstar(:,1:N-1))/dt);
        chi=chi+sigma*sum(wstar,1);



        if mod(i,1000)==0
            if isnan(sum(rhocri)) || isnan(sum(vcri)) || isnan(sum(gcri))
                break
            end
            rhocri=[rhocri;sum(sum((rho-rhoold).^2))];
            vcri=[vcri;sum(sum((v-vold).^2))];
            gcri=[gcri;sum((g-gold).^2)];
        end
    end
    toc
    
    lambda_text=num2str(lambda_index);
    gamma_text=num2str(gamma_index);
    fig_name=strcat('L2_','lambda',lambda_text,'gamma',gamma_text,'.jpg');
    file_name=strcat('L2_','lambda',lambda_text,'gamma',gamma_text,'.mat');

    save(file_name)
    
    

    figure
    subplot(2,3,1)
    plot(ghat)
    hold on
    plot(g)
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
    plot(gcri)
    
    
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
    save('L2_control_file.mat','lambda_opt','lambda_index','gamma_opt','gamma_index')
    
end
