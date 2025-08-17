function [hatbeta,hatSigma,iter]=OLSp(W,Y,AA,BB,AA0,beta0,Lambda0,mu,epsilon,xi,iterMAX)
[p,n]=size(Y);
sizeW=size(W);
K=sizeW(4);

hatbetaOLS=beta0;


for i=1:n
    WI(:,:,:)=W(:,:,i,:);
    hatSigmaOLS(:,:,i)=Ximat(hatbetaOLS,WI);
    [~,r] = chol(hatSigmaOLS(:,:,i));
    test(i)= ((r == 0) && (rank(hatSigmaOLS(:,:,i)) == p));
end





iter=0;

if prod(test)==1
    hatbeta=hatbetaOLS;
    hatSigma=hatSigmaOLS;
else
    AA=inv(2*mu*AA+AA0);
    
    iter=iter+1;
    
    %Theta Step
    Sigma0=hatSigmaOLS;
    Theta1=Sigma0+mu*Lambda0;
    for i=1:n
        [Gamma,D]=eig(Theta1(:,:,i));
        diagD=diag(D);
        diagD(diagD<epsilon)=epsilon;
        Theta1(:,:,i)=Gamma*diag(diagD)*Gamma';
        if issymmetric(Theta1(:,:,i))==0
            Theta1(:,:,i)=(Theta1(:,:,i)+Theta1(:,:,i)')/2;
        end
    end
    
    %Beta Step
    ThetaL1=Theta1-mu*Lambda0;
    BB1=zeros(K,1);
    for k=1:K
        for i=1:n
            bb(i)=trace(W(:,:,i,k)*ThetaL1(:,:,i));
        end
        BB1(k)=2*mu*BB(k)+sum(bb);
    end
    beta1=AA*BB1;
    
    %Lambda Step
    for i=1:n
        WI(:,:,:)=W(:,:,i,:);
        Sigma1(:,:,i)=Ximat(beta1,WI);
    end
    Lambda1=Lambda0-(Theta1-Sigma1)/mu;
    
    while ((norm(beta1-beta0)>xi) && (iter<iterMAX))
        iter=iter+1;
        beta0=beta1;
        Lambda0=Lambda1;
        Sigma0=Sigma1;
        
        %Theta Step
        Theta1=Sigma0+mu*Lambda0;
        for i=1:n
            [Gamma,D]=eig(Theta1(:,:,i));
            diagD=diag(D);
            diagD(diagD<epsilon)=epsilon;
            Theta1(:,:,i)=Gamma*diag(diagD)*Gamma';
            if issymmetric(Theta1(:,:,i))==0
                Theta1(:,:,i)=(Theta1(:,:,i)+Theta1(:,:,i)')/2;
            end
            %issymmetric(Theta1(:,:,i))
        end 
        

        
        %Beta Step
        ThetaL1=Theta1-mu*Lambda0;
        BB1=zeros(K,1);
        for k=1:K
            for i=1:n
                bb(i)=trace(W(:,:,i,k)*ThetaL1(:,:,i));
            end
            BB1(k)=2*mu*BB(k)+sum(bb);
        end
        beta1=AA*BB1;
 

    
    
    
    
    
        %Lambda Step
        for i=1:n
            WI(:,:,:)=W(:,:,i,:);
            Sigma1(:,:,i)=Ximat(beta1,WI);
        end
        Lambda1=Lambda0-(Theta1-Sigma1)/mu;
        

        
    end
    hatbeta=beta1;
    hatSigma=Sigma1;
end