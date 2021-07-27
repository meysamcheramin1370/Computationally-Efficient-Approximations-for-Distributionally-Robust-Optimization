function [ f_opt,X_opt,z_opt,CPUTime] = WrittenExamNewUB2CombinedProductionExp23( ZIGMAtemp,mu,m1,m2,m3,m4,N,switch1,R0,gamma1,gamma2,cosiPrime,c,d,alpha_k,beta_k,m,n)
%This function finds the new 4-partition upper bound OPT value (Combined ambiguity
%set)for the Risk-averse Production-Transportation problem.


size=length(mu);





[ A1Complete,A1,A2,A3,A4,ZIGMA,U,delta ] = New3covTransformerDecomposer(ZIGMAtemp,size,m1,m2,m3,m4,switch1);

tic

cvx_begin sdp

warning off;

variable x(m) nonnegative;
variable z1(size) nonnegative;
variable z2(size) nonnegative;
variable z3(size) nonnegative;
variable z4(size) nonnegative;
variable z5(size) nonnegative;
variables q(1,size) zeta(size,N) y(N) s11(N) s21(N) s31(N) s41(N) s12(N) s22(N) s32(N) s42(N) s13(N) s23(N) s33(N) s43(N) s14(N) s24(N) s34(N) s44(N) s15(N) s25(N) s35(N) s45(N);
variable Q1(m1,m1) symmetric;
variable Q2(m2,m2) symmetric;
variable Q3(m3,m3) symmetric;
variable Q4(m4,m4) symmetric;
variable Landa  nonnegative;



%%
minimize( c'*x + Landa*R0 + gamma2*(sum(sum(eye(m1).*Q1)) + sum(sum(eye(m2).*Q2)) + sum(sum(eye(m3).*Q3)) + sum(sum(eye(m4).*Q4))) + sqrt(gamma1)*norm(q,2) + (1/N) * sum(y) );
%%
subject to

for i = 1 : N
    
    [Q1     1/2*(q(1:m1)+(- alpha_k(1)*z1' + zeta(:,i)')*A1)' ; 1/2*(q(1:m1)+(- alpha_k(1)*z1' + zeta(:,i)')*A1) s11(i)]>=0 ;
    [Q2     1/2*(q(m1+1:m1+m2)+(- alpha_k(1)*z1' + zeta(:,i)')*A2)' ; 1/2*(q(m1+1:m1+m2)+(- alpha_k(1)*z1' + zeta(:,i)')*A2) s21(i)]>=0 ;
    [Q3     1/2*(q(m1+m2+1:m1+m2+m3)+(- alpha_k(1)*z1' + zeta(:,i)')*A3)' ; 1/2*(q(m1+m2+1:m1+m2+m3)+(- alpha_k(1)*z1' + zeta(:,i)')*A3) s31(i)]>=0 ;
    [Q4     1/2*(q(m1+m2+m3+1:m1+m2+m3+m4)+(- alpha_k(1)*z1' + zeta(:,i)')*A4)' ; 1/2*(q(m1+m2+m3+1:m1+m2+m3+m4)+(- alpha_k(1)*z1' + zeta(:,i)')*A4) s41(i)]>=0 ;
    
    [Q1     1/2*(q(1:m1)+(- alpha_k(2)*z2' + zeta(:,i)')*A1)' ; 1/2*(q(1:m1)+(- alpha_k(2)*z2' + zeta(:,i)')*A1) s12(i)]>=0 ;
    [Q2     1/2*(q(m1+1:m1+m2)+(- alpha_k(2)*z2' + zeta(:,i)')*A2)' ; 1/2*(q(m1+1:m1+m2)+(- alpha_k(2)*z2' + zeta(:,i)')*A2) s22(i)]>=0 ;
    [Q3     1/2*(q(m1+m2+1:m1+m2+m3)+(- alpha_k(2)*z2' + zeta(:,i)')*A3)' ; 1/2*(q(m1+m2+1:m1+m2+m3)+(- alpha_k(2)*z2' + zeta(:,i)')*A3) s32(i)]>=0 ;
    [Q4     1/2*(q(m1+m2+m3+1:m1+m2+m3+m4)+(- alpha_k(2)*z2' + zeta(:,i)')*A4)' ; 1/2*(q(m1+m2+m3+1:m1+m2+m3+m4)+(- alpha_k(2)*z2' + zeta(:,i)')*A4) s42(i)]>=0 ;
    
    [Q1     1/2*(q(1:m1)+(- alpha_k(3)*z3' + zeta(:,i)')*A1)' ; 1/2*(q(1:m1)+(- alpha_k(3)*z3' + zeta(:,i)')*A1) s13(i)]>=0 ;
    [Q2     1/2*(q(m1+1:m1+m2)+(- alpha_k(3)*z3' + zeta(:,i)')*A2)' ; 1/2*(q(m1+1:m1+m2)+(- alpha_k(3)*z3' + zeta(:,i)')*A2) s23(i)]>=0 ;
    [Q3     1/2*(q(m1+m2+1:m1+m2+m3)+(- alpha_k(3)*z3' + zeta(:,i)')*A3)' ; 1/2*(q(m1+m2+1:m1+m2+m3)+(- alpha_k(3)*z3' + zeta(:,i)')*A3) s33(i)]>=0 ;
    [Q4     1/2*(q(m1+m2+m3+1:m1+m2+m3+m4)+(- alpha_k(3)*z3' + zeta(:,i)')*A4)' ; 1/2*(q(m1+m2+m3+1:m1+m2+m3+m4)+(- alpha_k(3)*z3' + zeta(:,i)')*A4) s43(i)]>=0 ;
  
    [Q1     1/2*(q(1:m1)+(- alpha_k(4)*z4' + zeta(:,i)')*A1)' ; 1/2*(q(1:m1)+(- alpha_k(4)*z4' + zeta(:,i)')*A1) s14(i)]>=0 ;
    [Q2     1/2*(q(m1+1:m1+m2)+(- alpha_k(4)*z4' + zeta(:,i)')*A2)' ; 1/2*(q(m1+1:m1+m2)+(- alpha_k(4)*z4' + zeta(:,i)')*A2) s24(i)]>=0 ;
    [Q3     1/2*(q(m1+m2+1:m1+m2+m3)+(- alpha_k(4)*z4' + zeta(:,i)')*A3)' ; 1/2*(q(m1+m2+1:m1+m2+m3)+(- alpha_k(4)*z4' + zeta(:,i)')*A3) s34(i)]>=0 ;
    [Q4     1/2*(q(m1+m2+m3+1:m1+m2+m3+m4)+(- alpha_k(4)*z4' + zeta(:,i)')*A4)' ; 1/2*(q(m1+m2+m3+1:m1+m2+m3+m4)+(- alpha_k(4)*z4' + zeta(:,i)')*A4) s44(i)]>=0 ;
    
    
    [Q1     1/2*(q(1:m1)+(- alpha_k(5)*z5' + zeta(:,i)')*A1)' ; 1/2*(q(1:m1)+(- alpha_k(5)*z5' + zeta(:,i)')*A1) s15(i)]>=0 ;
    [Q2     1/2*(q(m1+1:m1+m2)+(- alpha_k(5)*z5' + zeta(:,i)')*A2)' ; 1/2*(q(m1+1:m1+m2)+(- alpha_k(5)*z5' + zeta(:,i)')*A2) s25(i)]>=0 ;
    [Q3     1/2*(q(m1+m2+1:m1+m2+m3)+(- alpha_k(5)*z5' + zeta(:,i)')*A3)' ; 1/2*(q(m1+m2+1:m1+m2+m3)+(- alpha_k(5)*z5' + zeta(:,i)')*A3) s35(i)]>=0 ;
    [Q4     1/2*(q(m1+m2+m3+1:m1+m2+m3+m4)+(- alpha_k(5)*z5' + zeta(:,i)')*A4)' ; 1/2*(q(m1+m2+m3+1:m1+m2+m3+m4)+(- alpha_k(5)*z5' + zeta(:,i)')*A4) s45(i)]>=0 ;

end

for i = 1 : N
    
    s11(i) + s21(i) + s31(i) + s41(i) == y(i) - alpha_k(1)*z1'*mu - beta_k(1) + zeta(:,i)'*(mu - cosiPrime(:,i));
    
    s12(i) + s22(i) + s32(i) + s42(i) == y(i) - alpha_k(2)*z2'*mu - beta_k(2) + zeta(:,i)'*(mu - cosiPrime(:,i));
    
    s13(i) + s23(i) + s33(i) + s43(i) == y(i) - alpha_k(3)*z3'*mu - beta_k(3) + zeta(:,i)'*(mu - cosiPrime(:,i));
    
    s14(i) + s24(i) + s34(i) + s44(i) == y(i) - alpha_k(4)*z4'*mu - beta_k(4) + zeta(:,i)'*(mu - cosiPrime(:,i));
    
    s15(i) + s25(i) + s35(i) + s45(i) == y(i) - alpha_k(5)*z5'*mu - beta_k(5) + zeta(:,i)'*(mu - cosiPrime(:,i));
    
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
end

