%% Aug2025
clear all;


%% Read DATA
% Mean Covariates
Xdata=readtable('X.csv');
Xnum=1;%Only consider one covariate (mktrf, based on CAPM) + intercept
dx=1+1;
Xdata=Xdata(:,(1:(dx-1)));

% Response and Covariance Regression Covariates
dataname='dataSTD.csv';
data=readtable(dataname);

% Covariance Regression Covariates
Zindex=3:7;
Zdata=data(:,Zindex);

% Number of time periods
sizeX=size(Xdata);
n=sizeX(1);

% Number of similarity matrices in total
sizeZ=size(Zdata);
K=sizeZ(2);

% Number of stocks
sizedata=size(data);
p=sizedata(1)/n;
sqrtp=sqrt(p);

% Response innitialization
Y=zeros(p,n);

% Covariance regression covariates innitialization
Z=zeros(p,n,K);

for j=1:p
    Y(j,:)=data.ab0((1+(j-1)*n):(j*n))';
    
    for k=1:K
        Z(j,:,k)=table2array(Zdata((1+(j-1)*n):(j*n),k))';
    end
end



%% Pre-setting of tunning parameters
% OLS+
epsilon=1e-5;
mu=0.05;
xi=1e-8;
Lambda0=zeros(p,p,n);
iterMAX=1000;


%% Obtaining Similarity Matrices
W1=speye(p); 
for i=1:n
    W(:,:,i,1)=full(W1);
        for j1=1:p
            for j2=1:p
                W2(j1,j2)=(Z(j1,i,1)==Z(j2,i,1))*1;
            end
        end
        W2=W2-diag(diag(W2));
        W(:,:,i,2)=W2;
        
        for k=2:K
            for j1=1:p
                for j2=1:p
                    W2(j1,j2)=exp(-(Z(j1,i,k)-Z(j2,i,k))^2);
                end
            end
            W2=W2-diag(diag(W2));
            W(:,:,i,(k+1))=W2;                
        end
end
    
 
%% Mean Regression Covariate Matrix
X=ones(dx,n);
X(2:dx,:)=table2array(Xdata)';

%% Mean Regression Coefficient Matrix Innitialization 
B=rand(p,dx);


%% Based on Table 5, the final model only includes two similarity matrices
% I_p and the similarity matrix constructed by IND
% So the similarity matrices used for the final model fitting are given below:
for i=1:n
    Wr(:,:,i,1)=W(:,:,i,1);
    Wr(:,:,i,2)=W(:,:,i,2);
end
    
%% Estimation

%OLS
[tildeV,AA,BB,AA0,Q2,Q,hatbetaOLS]=OLS(Wr,Y,X);
%OLS+
[hatbetaOLSp,hatSigmaOLSp,iterOLSp]=OLSp(Wr,Y,AA,BB,AA0,hatbetaOLS,Lambda0,mu,epsilon,xi,iterMAX);
% Is OLS equal to OLS+?
hatbetaOLS==hatbetaOLSp

    
%% Testing for Covariance Structure 
[T_QL,mu_QL,mu_QL1,sigma2_QL,sigma2_QL1]=TQL(tildeV,hatSigmaOLSp,dx,Q2,Q,AA0,AA,Wr,Y);

TS1=(T_QL-mu_QL1)/sqrt(sigma2_QL);
p_value1=2*(1-cdf('Normal',abs(TS1),0,1));

TS2=(T_QL-mu_QL1)/sqrt(sigma2_QL1);
p_value2=2*(1-cdf('Normal',abs(TS2),0,1));

% Replication of Row 1 in Table 5
[p_value1, p_value2]

% The other rows of Table 5 can be obtained by replacing Wr(:,:,:,:) with
% the corresponding similarity matrices

%% Inference for Covairance Regression Coefficients
% SE of OLS+
M0Q=AA/n/p;
invM0Q=inv(M0Q);
MhQ=Mlq(0.5,Wr,Y,hatSigmaOLSp,Q);

SE=sqrt(diag(2*invM0Q*MhQ*invM0Q/n/p));

% Replication of OLS+ in Table 6
[hatbetaOLSp*1000,SE*1000,hatbetaOLSp./SE]
   

% Adjusted Profile Score Estimation (Z-Estimation) and Inference
fun = @(param)(ascore(param,Wr,Y,X));
[hatbetaS,fvalS,exitflagS]=fsolve(fun,hatbetaOLSp);

MmhM=Mlmh(hatbetaS,Wr,Y,X);

SEz=sqrt(diag(2*inv(MmhM)/n/p));

% Replication of Z-Estimate in Table 6
[hatbetaS*1000,SEz*1000,hatbetaS./SEz]