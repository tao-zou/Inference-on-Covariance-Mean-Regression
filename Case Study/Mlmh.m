function AA=Mlmh(beta,W,Y,X)
[p,n]=size(Y);
sizeW=size(W);
K=sizeW(4);



bX=sparse(kron(X,eye(p)));
Omega=sparse(p*n,p*n);


for i=1:n
    WI(:,:,:)=W(:,:,i,:);
    
    Sigma=Ximat(beta,WI);
    Omega((i-1)*p+(1:p),(i-1)*p+(1:p))=inv(Sigma);
end

tildeM=Omega-Omega*(bX')*inv(bX*Omega*(bX'))*bX*Omega;




aa=zeros(K,K);
for k1=1:K
    for k2=1:k1
    Wk1=sparse(p*n,p*n);
    Wk2=sparse(p*n,p*n);
        for i=1:n
            Wk1(((i-1)*p+1):(i*p),((i-1)*p+1):(i*p))=W(:,:,i,k1);
            Wk2(((i-1)*p+1):(i*p),((i-1)*p+1):(i*p))=W(:,:,i,k2);
        end
    aa(k1,k2)=trace(tildeM*Wk1*tildeM*Wk2);
    end
end
    


AA=aa+tril(aa,-1).';



AA=AA/p/n;



