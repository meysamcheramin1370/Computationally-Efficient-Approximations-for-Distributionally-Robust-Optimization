function [ f_opt,X_opt,CPUTime] = GeneralNewUB2CombinedNewsvendorPaperExp22( ZIGMAtemp,mu,m1,m2,m3,N,switch1,R0,cosiPrime,c,v,g,gamma2)
%This function finds the new 3-partition upper bound OPT value (Combined ambiguity
%set)for Multiproduct Newsvendor problem.


m=length(mu);
n = m;




[ A1Complete,A1,A2,A3,ZIGMA,U,delta ] = New2covTransformerDecomposer(ZIGMAtemp,m,m1,m2,m3,switch1);

tic 

cvx_begin sdp
 
                variable X(n) nonnegative;
                variables zeta(m,N) y(N) s11(N) s21(N) s31(N) s12(N) s22(N) s32(N);
                variable Q1(m1,m1) symmetric;
                variable Q2(m2,m2) symmetric;
                variable Q3(m3,m3) symmetric;
                variable Landa  nonnegative;
               
                               
                %%
                minimize( Landa*R0 + gamma2 * (sum(sum(eye(m1).*Q1)) + sum(sum(eye(m2).*Q2))+ sum(sum(eye(m3).*Q3))) + (1/N) * sum(y) );
                %%
                subject to
                
                for i = 1 : N

                    [Q1     1/2*(zeta(:,i)'*A1)' ; 1/2*(zeta(:,i)'*A1)  s11(i)]>=0 ;
                    [Q2     1/2*(zeta(:,i)'*A2)' ; 1/2*(zeta(:,i)'*A2)  s21(i)]>=0 ;
                    [Q3     1/2*(zeta(:,i)'*A3)' ; 1/2*(zeta(:,i)'*A3)  s31(i)]>=0 ;

                    
                    
                    [Q1     1/2*(((v-g)' + zeta(:,i)')*A1)' ; 1/2*((v-g)' + zeta(:,i)')*A1 s12(i)]>=0 ;
                    [Q2     1/2*(((v-g)' + zeta(:,i)')*A2)' ; 1/2*((v-g)' + zeta(:,i)')*A2 s22(i)]>=0 ;
                    [Q3     1/2*(((v-g)' + zeta(:,i)')*A3)' ; 1/2*((v-g)' + zeta(:,i)')*A3 s32(i)]>=0 ;


                end
                
                for i = 1 : N
                    
                   s11(i) + s21(i) + s31(i) == y(i)-(c-v)'*X+zeta(:,i)'*(mu - cosiPrime(:,i));
                   
                   s12(i) + s22(i) + s32(i) == y(i)+(v-g)'*mu - ((c-v)'+(v-g)')*X + zeta(:,i)'*(mu - cosiPrime(:,i));
                
                end
                
                
                                
                for i = 1 : N
                    for j = 1 : m
                        abs(zeta(j,i)) <= Landa;
                    end
                end
                
                
                                      
                      

 cvx_end
    
CPUTime = toc; 
 
    disp(['Problem is ' cvx_status])
    if ~strfind(cvx_status,'Solved')
      return
    end
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp(['optimal value of cvx:',num2str(cvx_optval)]);
    f_opt=cvx_optval;
    X_opt=X;
   

end

