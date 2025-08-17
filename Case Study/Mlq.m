function AA=Mlq(l,W,Y,Sigma,Q)
[p,n]=size(Y);
sizeW=size(W);
K=sizeW(4);

Qb=kron(Q,speye(p));
Sigmab=sparse(p*n,p*n);


for i=1:n
    Sigmad(:,:,i)=Sigma(:,:,i)^(2*l);
    Sigmab(((i-1)*p+1):(i*p),((i-1)*p+1):(i*p))=Sigmad(:,:,i);
end

QSQ=Qb*Sigmab*Qb;


aa=zeros(K,K);
for k1=1:K
    for k2=1:k1
%    aaii=zeros(n,n,n,n);
%    for i=1:n
%        for j1=1:n
%            for j2=1:n
%                for j3=1:n
%                    aaii(i,j1,j2,j3)=trace(Q(i,j1)*Sigmad(:,:,j1)*Q(j1,j2)*W(:,:,j2,k1)*Q(j2,j3)*Sigmad(:,:,j3)*Q(j3,i)*W(:,:,i,k2));
%                end
%            end
%        end
%    end
%    aa(k1,k2)=sum(sum(sum(sum(aaii))));
    Wk1=sparse(p*n,p*n);
    Wk2=sparse(p*n,p*n);
        for i=1:n
            Wk1(((i-1)*p+1):(i*p),((i-1)*p+1):(i*p))=W(:,:,i,k1);
            Wk2(((i-1)*p+1):(i*p),((i-1)*p+1):(i*p))=W(:,:,i,k2);
        end
    aa(k1,k2)=trace(QSQ*Wk1*QSQ*Wk2);
    end
end
    


AA=aa+tril(aa,-1).';



AA=AA/p/n;



