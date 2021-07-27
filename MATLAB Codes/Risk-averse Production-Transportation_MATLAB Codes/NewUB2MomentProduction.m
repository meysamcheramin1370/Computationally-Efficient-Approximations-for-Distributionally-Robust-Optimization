function [ f_opt,X_opt,z_opt,CPUTime ] = NewUB2MomentProduction( ZIGMAtemp,mu,A,b,m1,gamma1,gamma2,switch1,c,d,alpha_k,beta_k,m,n)
%This function finds the new reformulation of 2-partition upper bound (General moment-based ambiguity
%set)for Risk-averse Production-Transportation problem.

size=length(mu);

m2 = size - m1;

[ A1Complete,A1,A2,ZIGMA,U,delta ] = New1covTransformerDecomposer(ZIGMAtemp,size,m1,m2,switch1);

tic 

cvx_begin sdp

warning off;
 
                variable x(m) nonnegative;
                variable z1(size) nonnegative;
                variable z2(size) nonnegative;
                variable z3(size) nonnegative;
                variable z4(size) nonnegative;
                variable z5(size) nonnegative;
                variables q(size) s s11 s21 s12 s22 s13 s23 s14 s24 s15 s25;
                variable Q1(m1,m1) symmetric;
                variable Q2(m2,m2);
                variable Landa1(2*size) nonnegative;
                variable Landa2(2*size) nonnegative;
                variable Landa3(2*size) nonnegative;
                variable Landa4(2*size) nonnegative;
                variable Landa5(2*size) nonnegative;

                                
                %%
                minimize( c'*x + s + sum(sum((gamma2*eye(m1)).*Q1))+ sum(sum((gamma2*eye(m2)).*Q2)) + sqrt(gamma1)*norm(q,2) );
                %%
                subject to

                [s11     1/2*(q(1:m1)+ A1'*(A'*Landa1 - alpha_k(1)*z1))' ; 1/2*(q(1:m1) + A1'*(A'*Landa1 - alpha_k(1)*z1)) Q1]>=0 ;
                [s21     1/2*(q(m1+1:m1+m2)+ A2'*(A'*Landa1 - alpha_k(1)*z1))' ; 1/2*(q(m1+1:m1+m2) + A2'*(A'*Landa1 - alpha_k(1)*z1)) Q2]>=0 ;
              
                [s12     1/2*(q(1:m1)+ A1'*(A'*Landa2 - alpha_k(2)*z2))' ; 1/2*(q(1:m1) + A1'*(A'*Landa2 - alpha_k(2)*z2)) Q1]>=0 ;
                [s22     1/2*(q(m1+1:m1+m2)+ A2'*(A'*Landa2 - alpha_k(2)*z2))' ; 1/2*(q(m1+1:m1+m2) + A2'*(A'*Landa2 - alpha_k(2)*z2)) Q2]>=0 ;
                
                [s13     1/2*(q(1:m1)+ A1'*(A'*Landa3 - alpha_k(3)*z3))' ; 1/2*(q(1:m1) + A1'*(A'*Landa3 - alpha_k(3)*z3)) Q1]>=0 ;
                [s23     1/2*(q(m1+1:m1+m2)+ A2'*(A'*Landa3 - alpha_k(3)*z3))' ; 1/2*(q(m1+1:m1+m2) + A2'*(A'*Landa3 - alpha_k(3)*z3)) Q2]>=0 ;

                [s14     1/2*(q(1:m1)+ A1'*(A'*Landa4 - alpha_k(4)*z4))' ; 1/2*(q(1:m1) + A1'*(A'*Landa4 - alpha_k(4)*z4)) Q1]>=0 ;
                [s24     1/2*(q(m1+1:m1+m2)+ A2'*(A'*Landa4 - alpha_k(4)*z4))' ; 1/2*(q(m1+1:m1+m2) + A2'*(A'*Landa4 - alpha_k(4)*z4)) Q2]>=0 ;

                [s15     1/2*(q(1:m1)+ A1'*(A'*Landa5 - alpha_k(5)*z5))' ; 1/2*(q(1:m1) + A1'*(A'*Landa5 - alpha_k(5)*z5)) Q1]>=0 ;
                [s25     1/2*(q(m1+1:m1+m2)+ A2'*(A'*Landa5 - alpha_k(5)*z5))' ; 1/2*(q(m1+1:m1+m2) + A2'*(A'*Landa5 - alpha_k(5)*z5)) Q2]>=0 ;

                
                
                s11 + s21 == s- beta_k(1) - alpha_k(1)*z1'*mu - Landa1'*(b-A*mu);
                s12 + s22 == s- beta_k(2) - alpha_k(2)*z2'*mu - Landa2'*(b-A*mu);
                s13 + s23 == s- beta_k(3) - alpha_k(3)*z3'*mu - Landa3'*(b-A*mu);
                s14 + s24 == s- beta_k(4) - alpha_k(4)*z4'*mu - Landa4'*(b-A*mu);
                s15 + s25 == s- beta_k(5) - alpha_k(5)*z5'*mu - Landa5'*(b-A*mu);
                
                
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

