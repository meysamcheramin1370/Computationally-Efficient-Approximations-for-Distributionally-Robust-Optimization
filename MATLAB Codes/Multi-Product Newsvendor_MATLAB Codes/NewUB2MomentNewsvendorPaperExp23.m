function [ f_opt,X_opt,CPUTime ] = NewUB2MomentNewsvendorPaperExp23( ZIGMAtemp,mu,A,b,m1,m2,m3,m4,gamma1,gamma2,switch1,c,v,g)
%This function finds the new reformulation of 4-partition upper bound (General moment-based ambiguity
%set)for Multiproduct Newsvendor problem.

m=length(mu);
n = m;

[ A1Complete,A1,A2,A3,A4,ZIGMA,U,delta ] = New3covTransformerDecomposer(ZIGMAtemp,m,m1,m2,m3,m4,switch1);

tic 

cvx_begin sdp
 
                variable X(n) nonnegative;
                variables q(m) s s11 s21 s31 s41 s12 s22 s32 s42;
                variable Q1(m1,m1) symmetric;
                variable Q2(m2,m2) symmetric;
                variable Q3(m3,m3) symmetric;
                variable Q4(m4,m4) symmetric;
                variable Landa1(2*m) nonnegative;
                variable Landa2(2*m) nonnegative;

                                
                %%
                minimize( s + sum(sum((gamma2*eye(m1)).*Q1))+ sum(sum((gamma2*eye(m2)).*Q2))+ sum(sum((gamma2*eye(m3)).*Q3))+ sum(sum((gamma2*eye(m4)).*Q4)) + sqrt(gamma1)*norm(q,2) );
                %%
                subject to

                [s11     1/2*(q(1:m1)+ A1'*(A'*Landa1))' ; 1/2*(q(1:m1) + A1'*(A'*Landa1)) Q1]>=0 ;
                [s21     1/2*(q(m1+1:m1+m2)+ A2'*(A'*Landa1))' ; 1/2*(q(m1+1:m1+m2) + A2'*(A'*Landa1)) Q2]>=0 ;
                [s31     1/2*(q(m1+m2+1:m1+m2+m3)+ A3'*(A'*Landa1))' ; 1/2*(q(m1+m2+1:m1+m2+m3)+ A3'*(A'*Landa1)) Q3]>=0 ;
                [s41     1/2*(q(m1+m2+m3+1:m1+m2+m3+m4)+ A4'*(A'*Landa1))' ; 1/2*(q(m1+m2+m3+1:m1+m2+m3+m4)+ A4'*(A'*Landa1)) Q4]>=0 ;

                
                
                [s12    1/2*(q(1:m1) +A1'*(A'*Landa2+(v-g)))' ; 1/2*(q(1:m1) +A1'*(A'*Landa2+(v-g))) Q1]>=0;
                [s22    1/2*(q(m1+1:m1+m2) +A2'*(A'*Landa2+(v-g)))' ; 1/2*(q(m1+1:m1+m2) +A2'*(A'*Landa2+(v-g))) Q2]>=0;
                [s32    1/2*(q(m1+m2+1:m1+m2+m3) +A3'*(A'*Landa2+(v-g)))' ; 1/2*(q(m1+m2+1:m1+m2+m3) +A3'*(A'*Landa2+(v-g))) Q3]>=0;
                [s42    1/2*(q(m1+m2+m3+1:m1+m2+m3+m4) +A4'*(A'*Landa2+(v-g)))' ; 1/2*(q(m1+m2+m3+1:m1+m2+m3+m4) +A4'*(A'*Landa2+(v-g))) Q4]>=0;

                
                s11 + s21 + s31 + s41 == s-(c-v)'*X-(Landa1')*(b-A*mu);
                s12 + s22 + s32 + s42 == s-((c-v)'+(v-g)')*X-Landa2'*(b-A*mu)+(v-g)'*mu;              
               
                             
                
                

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

