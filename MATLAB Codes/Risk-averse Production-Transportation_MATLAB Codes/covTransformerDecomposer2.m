function [ A1Complete,A1,A1Prime,ZIGMA,U,delta ] = covTransformerDecomposer2(ZIGMAtemp,m,m1,switch1)
%1) % This function is created in response to IJOC Reviewer 2.
% Unlike, the first version of this function, it does not sort 
% the principal components (eigenvectors) based on their variance
% (eigenvalues) in decending orders. This is because we want to see
% how bad the results is when we select the first m_1 principal components
% arbitrary 


%2) This function does eigenvalue decomposition but it does not
% sort the eigenvalue matrix and the eigenvector matri in the 
% decending order of the eigenvalues (variances)

[V,HH,W] = eig(ZIGMAtemp);
%[h,ind] = sort(diag(HH),'descend');
%H = HH(ind,ind);
H = HH;
%Vs = V(:,ind);
Vs = V;
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
    %                 Replace the eigen values with EXPONENTIAL (slope 1):
    
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
A1Prime = Vs(:,m1 + 1:m)*H(m1 + 1:m,m1 + 1:m).^(1/2); % This matrix is used only for Upper Bound problem


end

