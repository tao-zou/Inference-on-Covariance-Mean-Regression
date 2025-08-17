function [tildeV,AA,BB,AA0,Q2,Q,hatbeta]=OLS(W,Y,X)
[p,n]=size(Y);
sizeW=size(W);
K=sizeW(4);

Q=eye(n)-X'*inv(X*X')*X;
Q2=Q.^2;
tildeV=zeros(p,n);
QY=zeros(p,n);

for i=1:n
    for j=1:n
        QY(:,j)=Q(i,j)*Y(:,j);
    end
    tildeV(:,i)=sum(QY,2);
end


aa=zeros(K,K);
for k1=1:K
    for k2=1:k1
        aaii=zeros(n,n);
        for i1=1:n
            for i2=1:n
                aaii(i1,i2)=Q2(i1,i2)*trace(W(:,:,i1,k1)*W(:,:,i2,k2));
            end
        end
        aa(k1,k2)=sum(sum(aaii));
    end
end

AA=aa+tril(aa,-1).';

    
    
for i=1:n
    aa=zeros(K,K);
    for k1=1:K
        for k2=1:k1
            aa(k1,k2)=trace(W(:,:,i,k1)*W(:,:,i,k2));
        end
    end
    aa=aa+tril(aa,-1).';
    aaa(:,:,i)=aa;
end

AA0=sum(aaa,3);





BB=zeros(K,1);
for k=1:K
    for i=1:n
        bb(i)=tildeV(:,i)'*W(:,:,i,k)*tildeV(:,i);
    end
    BB(k)=sum(bb);
end

hatbeta=AA\BB;
