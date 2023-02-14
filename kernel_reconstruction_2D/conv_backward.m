function rhoconv=conv_backward(K,rho,init,term,M,L,N,dx)
revK=zeros(M,M);
revK(1:L+1,1:L+1)=K;
revK(1:L+1,L+2:M)=fliplr(revK(1:L+1,2:L));
revK(L+2:M,:)=flipud(revK(2:L,:));
rhoconv=zeros(M,M,N+1);
for i1=1:M
    for i2=1:M
        rhoconv(i1,i2,init:term)=sum(sum(...
            repmat(revK,[1,1,term-init+1]).*rho(:,:,init:term),1),2)*dx^2;
        revK=[revK(:,M),revK(:,1:M-1)];
    end
    revK=[revK(M,:);revK(1:M-1,:)];
end
rhoconv=rhoconv(:,:,init:term);
end