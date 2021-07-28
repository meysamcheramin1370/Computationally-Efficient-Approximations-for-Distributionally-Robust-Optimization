function [ f_opt,X_opt,CPUTime] = NewUB2CombinedNewsvendor( ZIGMAtemp,mu,m1,N,switch1,R0,cosiPrime,c,v,g)
%This function finds the new 2-partition upper bound OPT value (Combined ambiguity
%set)for Multiproduct Newsvendor problem.


m=length(mu);
n = m;

m2 = m - m1;


[ A1Complete,A1,A2,ZIGMA,U,delta ] = New1covTransformerDecomposer(ZIGMAtemp,m,m1,m2,switch1);

tic 

cvx_begin sdp
 
                variable X(n) nonnegative;
                variables zeta(m,N) y(N) s11(N) s21(N) s12(N) s22(N);
                variable Q1(m1,m1) symmetric;
                variable Q2(m2,m2);
                variable Landa  nonnegative;
               
                               
                %%
                minimize( Landa*R0 + sum(sum(eye(m1).*Q1)) + sum(sum(eye(m2).*Q2)) + (1/N) * sum(y) );
                %%
                subject to
                
                for i = 1 : N

                    [Q1     1/2*(zeta(:,i)'*A1)' ; 1/2*(zeta(:,i)'*A1)  s11(i)]>=0 ;
                    [Q2     1/2*(zeta(:,i)'*A2)' ; 1/2*(zeta(:,i)'*A2)  s21(i)]>=0 ;

                    
                    
                    [Q1     1/2*(((v-g)' + zeta(:,i)')*A1)' ; 1/2*((v-g)' + zeta(:,i)')*A1 s12(i)]>=0 ;
                    [Q2     1/2*(((v-g)' + zeta(:,i)')*A2)' ; 1/2*((v-g)' + zeta(:,i)')*A2 s22(i)]>=0 ;

                end
                
                for i = 1 : N
                    
                   s11(i) + s21(i) == y(i)-(c-v)'*X+zeta(:,i)'*(mu - cosiPrime(:,i));
                   
                   s12(i) + s22(i) == y(i)+(v-g)'*mu - ((c-v)'+(v-g)')*X + zeta(:,i)'*(mu - cosiPrime(:,i));
                
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

