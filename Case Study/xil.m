function AA=xil(l,W,Y,Sigma)
[p,n]=size(Y);
sizeW=size(W);
K=sizeW(4);


for i=1:n
    Sigmad=Sigma(:,:,i)^(2*l);
    aa=zeros(K,1);
    for k1=1:K
        aa(k1)=trace(Sigmad*W(:,:,i,k1));
    end
    aaa(:,i)=aa;
end

AA=mean(aaa,2)/p;



