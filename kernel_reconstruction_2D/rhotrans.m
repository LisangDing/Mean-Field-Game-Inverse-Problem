function rho_trans=rhotrans(rho,x,y,M)
rho_trans=rho;
% x=q1-1, y=q2-1
if x<0
    x=M+x;
end
if y<0
    y=M+y;
end
rho1=[rho(1+x:M,:,:);rho(1:x,:,:)];
rho_trans=[rho1(:,1+y:M,:),rho1(:,1:y,:)];
end