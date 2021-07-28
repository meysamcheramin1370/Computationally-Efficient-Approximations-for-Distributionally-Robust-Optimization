function [ f_opt,X_opt,z_opt,CPUTime] = OriginalCombinedProduction( ZIGMAtemp,mu,N,switch1,R0,gamma2,cosiPrime,c,d,alpha_k,beta_k,m,n)
%This function finds the original OPT value (Combined ambiguity
%set)for the Risk-averse Production-Transportation problem.


size=length(mu);

m1 = size; % This parameter is not used. its just for fullfilling m1 in covTransformerDecomposer.

[A1Complete,A1,A1Prime,ZIGMA,U,delta ] = covTransformerDecomposer(ZIGMAtemp,size,m1,switch1);

tic

cvx_begin sdp

variable x(m) nonnegative;
variable z1(size) nonnegative;
variable z2(size) nonnegative;
variable z3(size) nonnegative;
variable z4(size) nonnegative;
variable z5(size) nonnegative;
variables zeta(size,N) y(N);
variable Q(size,size) symmetric;
variable Landa  nonnegative;


%%
minimize( c'*x + Landa*R0 + gamma2*sum(sum(eye(size).*Q)) + (1/N) * sum(y) );
%%
subject to

for i = 1 : N
    
    [Q     1/2*((- alpha_k(1)*z1' + zeta(:,i)')*A1Complete)' ; 1/2*(- alpha_k(1)*z1' + zeta(:,i)')*A1Complete  y(i)-alpha_k(1)*z1'*mu - beta_k(1) + zeta(:,i)'*(mu - cosiPrime(:,i))]>=0 ;
    [Q     1/2*((- alpha_k(2)*z2' + zeta(:,i)')*A1Complete)' ; 1/2*(- alpha_k(2)*z2' + zeta(:,i)')*A1Complete  y(i)-alpha_k(2)*z2'*mu - beta_k(2) + zeta(:,i)'*(mu - cosiPrime(:,i))]>=0 ;
    [Q     1/2*((- alpha_k(3)*z3' + zeta(:,i)')*A1Complete)' ; 1/2*(- alpha_k(3)*z3' + zeta(:,i)')*A1Complete  y(i)-alpha_k(3)*z3'*mu - beta_k(3) + zeta(:,i)'*(mu - cosiPrime(:,i))]>=0 ;
    [Q     1/2*((- alpha_k(4)*z4' + zeta(:,i)')*A1Complete)' ; 1/2*(- alpha_k(4)*z4' + zeta(:,i)')*A1Complete  y(i)-alpha_k(4)*z4'*mu - beta_k(4) + zeta(:,i)'*(mu - cosiPrime(:,i))]>=0 ;
    [Q     1/2*((- alpha_k(5)*z5' + zeta(:,i)')*A1Complete)' ; 1/2*(- alpha_k(5)*z5' + zeta(:,i)')*A1Complete  y(i)-alpha_k(5)*z5'*mu - beta_k(5) + zeta(:,i)'*(mu - cosiPrime(:,i))]>=0 ;
    
end



for i = 1 : N
    for j = 1 : size
        abs(zeta(j,i)) <= Landa;
    end
end


for j = 1 : n
    
    z12 = reshape(z1,n,m)';
    
    sum (z12(:,j)) == d(j);
    
end

for j = 1 : n
    
    z22 = reshape(z2,n,m)';
    
    sum (z22(:,j)) == d(j);
    
end

for j = 1 : n
    
    z32 = reshape(z3,n,m)';
    
    sum (z32(:,j)) == d(j);
    
end

for j = 1 : n
    
    z42 = reshape(z4,n,m)';
    
    sum (z42(:,j)) == d(j);
    
end

for j = 1 : n
    
    z52 = reshape(z5,n,m)';
    
    sum (z52(:,j)) == d(j);
    
end



for i = 1 : m
    
    z13 = reshape(z1,n,m)';
    
    sum (z13(i,:)) == x(i);
    
end

for i = 1 : m
    
    z23 = reshape(z2,n,m)';
    
    sum (z23(i,:)) == x(i);
    
end

for i = 1 : m
    
    z33 = reshape(z3,n,m)';
    
    sum (z33(i,:)) == x(i);
    
end

for i = 1 : m
    
    z43 = reshape(z4,n,m)';
    
    sum (z43(i,:)) == x(i);
    
end

for i = 1 : m
    
    z53 = reshape(z5,n,m)';
    
    sum (z53(i,:)) == x(i);
    
end

x <= 1;




cvx_end

CPUTime = toc;

disp(['Problem is ' cvx_status])
if ~strfind(cvx_status,'Solved')
    return
end
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(['optimal value of cvx:',num2str(cvx_optval)]);
f_opt=cvx_optval;
X_opt=x;
z_opt=[z1,z2,z3,z4,z5];
zetaOpt = zeta;


end

