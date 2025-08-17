function v=ascore(beta,W,Y,X)
[p,n]=size(Y);
[~,~,~,K]=size(W);

bX=sparse(kron(X,eye(p)));
Omega=sparse(p*n,p*n);

for i=1:n
    WI(:,:,:)=W(:,:,i,:);
    
    Sigma=Ximat(beta,WI);
    Omega((i-1)*p+(1:p),(i-1)*p+(1:p))=inv(Sigma);
    
end

tildeM=Omega-Omega*(bX')*inv(bX*Omega*(bX'))*bX*Omega;

tildeU=tildeM*reshape(Y,[p*n,1]);

for k=1:K
    for i=1:n
        quar(i)= tildeU((i-1)*p+(1:p))'*W(:,:,i,k)*tildeU((i-1)*p+(1:p))-trace(W(:,:,i,k)*tildeM((i-1)*p+(1:p),(i-1)*p+(1:p)));
    end
    v(k)=sum(quar);
end

