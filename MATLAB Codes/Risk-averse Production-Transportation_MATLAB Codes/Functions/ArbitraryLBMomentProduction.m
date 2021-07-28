function [ f_opt,X_opt,z_opt,CPUTime,Landa1Opt,Landa2Opt,Landa3Opt,Landa4Opt,Landa5Opt ] = ArbitraryLBMomentProduction( ZIGMAtemp,mu,m1,A,b,gamma1,gamma2,switch1,c,d,alpha_k,beta_k,m,n)
%This function finds the lower bound value (General moment-based ambiguity
%set)for Risk-averse Production-Transportation problem.

% This function is created in response to IJOC reviewer's comment
% about how bad the result is if we choose the first m1 components
% arbitrary Not in the decending order of eigenvalues (variances)
size=length(mu);


[ A1Complete,A1,A1Prime,ZIGMA,U,delta ] = covTransformerDecomposer2(ZIGMAtemp,size,m1,switch1);

tic 

cvx_begin sdp
 
                variable x(m) nonnegative;
                variable z1(size) nonnegative;
                variable z2(size) nonnegative;
                variable z3(size) nonnegative;
                variable z4(size) nonnegative;
                variable z5(size) nonnegative;
                variables qrr(m1) s;
                variable Qr(m1,m1) symmetric;
                variable Landa1(2*size) nonnegative;
                variable Landa2(2*size) nonnegative;
                variable Landa3(2*size) nonnegative;
                variable Landa4(2*size) nonnegative;
                variable Landa5(2*size) nonnegative;


                

                                
                %%
                minimize( c'*x + s + sum(sum((gamma2*eye(m1)).*Qr)) + sqrt(gamma1)*norm(qrr,2) );
                %%
                subject to
                
                

                [s-beta_k(1)-Landa1'*(b-A*mu) - alpha_k(1)*z1'*mu     1/2*(qrr + A1'*(A'*Landa1-alpha_k(1)*z1))' ; 1/2*(qrr + A1'*(A'*Landa1-alpha_k(1)*z1)) Qr] >= 0 ;
                [s-beta_k(2)-Landa2'*(b-A*mu) - alpha_k(2)*z2'*mu     1/2*(qrr + A1'*(A'*Landa2-alpha_k(2)*z2))' ; 1/2*(qrr + A1'*(A'*Landa2-alpha_k(2)*z2)) Qr] >= 0 ;
                [s-beta_k(3)-Landa3'*(b-A*mu) - alpha_k(3)*z3'*mu     1/2*(qrr + A1'*(A'*Landa3-alpha_k(3)*z3))' ; 1/2*(qrr + A1'*(A'*Landa3-alpha_k(3)*z3)) Qr] >= 0 ;
                [s-beta_k(4)-Landa4'*(b-A*mu) - alpha_k(4)*z4'*mu     1/2*(qrr + A1'*(A'*Landa4-alpha_k(4)*z4))' ; 1/2*(qrr + A1'*(A'*Landa4-alpha_k(4)*z4)) Qr] >= 0 ;
                [s-beta_k(5)-Landa5'*(b-A*mu) - alpha_k(5)*z5'*mu     1/2*(qrr + A1'*(A'*Landa5-alpha_k(5)*z5))' ; 1/2*(qrr + A1'*(A'*Landa5-alpha_k(5)*z5)) Qr] >= 0 ;

                                                             
               
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
    
    Landa1Opt = Landa1; 
    Landa2Opt = Landa2;
    Landa3Opt = Landa3; 
    Landa4Opt = Landa4;
    Landa5Opt = Landa5; 
    
end

