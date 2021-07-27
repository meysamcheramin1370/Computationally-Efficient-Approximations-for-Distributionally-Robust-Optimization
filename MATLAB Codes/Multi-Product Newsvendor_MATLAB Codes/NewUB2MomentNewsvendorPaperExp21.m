function [ f_opt,X_opt,CPUTime ] = NewUB2MomentNewsvendorPaperExp21( ZIGMAtemp,mu,A,b,m1,m2,gamma1,gamma2,switch1,c,v,g)
%This function finds the new reformulation of 2-partition upper bound (General moment-based ambiguity
%set)for Multiproduct Newsvendor problem.

m=length(mu);
n = m;

[ A1Complete,A1,A2,ZIGMA,U,delta ] = New1covTransformerDecomposer(ZIGMAtemp,m,m1,m2,switch1);

tic 

cvx_begin sdp
 
                variable X(n) nonnegative;
                variables q(m) s s11 s21 s12 s22;
                variable Q1(m1,m1) symmetric;
                variable Q2(m2,m2) symmetric;
                variable Landa1(2*m) nonnegative;
                variable Landa2(2*m) nonnegative;

                                
                %%
                minimize( s + sum(sum((gamma2*eye(m1)).*Q1))+ sum(sum((gamma2*eye(m2)).*Q2)) + sqrt(gamma1)*norm(q,2) );
                %%
                subject to

                [s11     1/2*(q(1:m1)+ A1'*(A'*Landa1))' ; 1/2*(q(1:m1) + A1'*(A'*Landa1)) Q1]>=0 ;
                [s21     1/2*(q(m1+1:m1+m2)+ A2'*(A'*Landa1))' ; 1/2*(q(m1+1:m1+m2) + A2'*(A'*Landa1)) Q2]>=0 ;

                
                
                [s12    1/2*(q(1:m1) +A1'*(A'*Landa2+(v-g)))' ; 1/2*(q(1:m1) +A1'*(A'*Landa2+(v-g))) Q1]>=0;
                [s22    1/2*(q(m1+1:m1+m2) +A2'*(A'*Landa2+(v-g)))' ; 1/2*(q(m1+1:m1+m2) +A2'*(A'*Landa2+(v-g))) Q2]>=0;

                
                s11 + s21 == s-(c-v)'*X-(Landa1')*(b-A*mu);
                s12 + s22 == s-((c-v)'+(v-g)')*X-Landa2'*(b-A*mu)+(v-g)'*mu;              
               
                             
                
                

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

