function [ A1Complete,A1,A2,A3,A4,A5,ZIGMA,U,delta ] = New4covTransformerDecomposer(ZIGMAtemp,m,m1,m2,m3,m4,m5,switch1)
% (this function is used for the 4-partiton upper bound).  With the randomly generated covariance matrix \Sigma, we replace its ith
%largest eigenvalue by three different kinds of generating functions: the first one is constant


%2) This function does eigenvalue decomposition

[V,HH,W] = eig(ZIGMAtemp);
[h,ind] = sort(diag(HH),'descend');
H = HH(ind,ind);
Vs = V(:,ind);

%                 %UPDATE EIGENVALUES

%*************************************************************************
                %Replace the eigenvalues with constant (max):
if switch1 == 1
    
    maxi=max(diag(H));
    H=maxi*eye(m,m);
    %*************************************************************************
    %*************************************************************************
    %                 Replace the eigenvalues with LINEAR:
    
elseif switch1 == 2
    y = linspace(1,.5,m);
    H=diag(y);
    %*************************************************************************
    %                 Replace the eigen values with EXPONENTIAL:
    
elseif switch1 == 3
    
    slope= 0.1;
    
    for i = 1:m
        y(i) = 1-(exp(((-i+m+1)/m)*slope)-1)/(exp(slope)-1);
    end
    H=diag(y);
    
    %*************************************************************************

     %                 Replace the eigen values with EXPONENTIAL (slope 1):
    
elseif switch1 == 4
    
    slope= 1;
    
    for i = 1:m
        y(i) = 1-(exp(((-i+m+1)/m)*slope)-1)/(exp(slope)-1);
    end
    H=diag(y);
    
    %*************************************************************************
%                 Replace the eigen values with EXPONENTIAL (slope 3):
    
elseif switch1 == 5
    
    slope= 5;
    
    for i = 1:m
        y(i) = 1-(exp(((-i+m+1)/m)*slope)-1)/(exp(slope)-1);
    end
    H=diag(y);
    
    %*************************************************************************
%                 Replace the eigen values with EXPONENTIAL (slope 3):
    
elseif switch1 == 6
    
    slope= 15;
    
    for i = 1:m
        y(i) = 1-(exp(((-i+m+1)/m)*slope)-1)/(exp(slope)-1);
    end
    H=diag(y);
    
    %*************************************************************************
  

end

if switch1 == 0   
   ZIGMA = ZIGMAtemp;
else
    ZIGMA= Vs*H*Vs';
end

U = Vs;
delta = H;

A1Complete = Vs*H.^(1/2); %for the original reformulation
A1=Vs(:,1:m1)*H(1:m1,1:m1).^(1/2);
A2 = Vs(:,m1 + 1:m1 + m2)*H(m1 + 1:m1 + m2,m1 + 1:m1 + m2).^(1/2); % This matrix is used only for 2-partition Upper Bound problem
A3 = Vs(:,m1 + m2 + 1:m1 + m2 + m3)*H(m1 + m2 + 1:m1 + m2 + m3,m1 + m2 + 1:m1 + m2 + m3).^(1/2); % This matrix is used only for 3-partition Upper Bound problem
A4 = Vs(:,m1 + m2 + m3 + 1:m1 + m2 + m3 + m4)*H(m1 + m2 + m3 + 1:m1 + m2 + m3 + m4,m1 + m2 + m3 + 1:m1 + m2 + m3 + m4).^(1/2); % This matrix is used only for 4-partition Upper Bound problem
A5 = Vs(:,m1 + m2 + m3 + m4 + 1:m1 + m2 + m3 + m4 + m5)*H(m1 + m2 + m3 + m4 + 1:m1 + m2 + m3 + m4 + m5,m1 + m2 + m3 + m4 + 1:m1 + m2 + m3 + m4 + m5).^(1/2); % This matrix is used only for 4-partition Upper Bound problem

end

