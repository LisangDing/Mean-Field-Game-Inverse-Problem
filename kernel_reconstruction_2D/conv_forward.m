function rhoconv=conv_forward(A,epsilon,rho,init,term,M,N,dx)
rhoconv=zeros(M,M,N+1);
for i1=1:M
    for i2=1:M
        x=abs(i1-(1:M)'*ones(1,M))*dx;
        y=abs(i2-ones(M,1)*(1:M))*dx;
        xdist=min(x,1-x);
        ydist=min(y,1-y);
        kernel=exp(-(repmat(A(1,1),[M,M]).*xdist.^2+...
            repmat(A(1,2)+A(2,1),[M,M]).*xdist.*ydist+...
            repmat(A(2,2),[M,M]).*ydist.^2)/epsilon);
        rhoconv(i1,i2,init:term)=sum(sum(repmat(kernel,[1,1,term-init+1]).*...
            rho(:,:,init:term)*dx^2,1),2);
    end
end
rhoconv=rhoconv(:,:,init:term);
end