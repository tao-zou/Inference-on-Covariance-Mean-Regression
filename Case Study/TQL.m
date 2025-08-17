function [T_QL,mu_QL,mu_QL1,sigma2_QL,sigma2_QL1]=TQL(tildeV,hatSigma,dx,Q2,Q,AA0,AA,W,Y)
[p,n]=size(tildeV);
c=p/n;

I_p=eye(p);
invhatSigma=hatSigma;

for i=1:n
    invhatSigma(:,:,i)=inv(hatSigma(:,:,i));
    tt(i)=trace((tildeV(:,i)*tildeV(:,i)'*invhatSigma(:,:,i)-I_p)^2);
end

T_QL=mean(tt)/p;

mumu=zeros(n,n);

for i=1:n
    for l=1:n
        if l~=i
            mumu(i,l)=Q2(i,l)*trace(invhatSigma(:,:,i)*(hatSigma(:,:,l)-hatSigma(:,:,i)));
        end
    end
end

mu_QL=p+1-2*c*dx+2*sum(sum(mumu))/n;

mu_QL1=p+1-c*(n-sum(diag(Q2)))+2*sum(sum(mumu))/n;

M0=AA0/n/p;
invM0=inv(M0);
M0q=AA/n/p;
invM0q=inv(M0q);

%Mhq=Mlq(0,W,Y,hatSigma,Q);
Mh=Ml(0.5,W,Y,hatSigma);
%Mhq=Mlq(0.5,W,Y,hatSigma,Q);

xih=xil(0.5,W,Y,hatSigma);
ximh=xil(-0.5,W,Y,hatSigma);

sigma2_QL=8*c*(1+ximh'*invM0*Mh*invM0*ximh-2*ximh'*invM0*xih);

sigma2_QL1=8*c*(1+ximh'*invM0q*Mh*invM0q*ximh-2*ximh'*invM0q*xih);







