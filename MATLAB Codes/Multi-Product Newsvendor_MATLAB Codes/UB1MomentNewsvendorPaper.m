function [ f_opt,X_opt,CPUTime ] = UB1MomentNewsvendorPaper( ZIGMAtemp,mu,m1,A,b,gamma1,gamma2,switch1,c,v,g)
%This function finds the first upper bound OPT value (General moment-based ambiguity
%set)for Multiproduct Newsvendor problem.

m=length(mu);
n = m;

[ A1Complete,A1,A1Prime,ZIGMA,U,delta ] = covTransformerDecomposer(ZIGMAtemp,m,m1,switch1);

tic 

cvx_begin sdp
 
                variable X(n) nonnegative;
                variables q(m) s;
                variable Qr(m1,m1) symmetric;
                variable Landa1(2*m) nonnegative;
                variable Landa2(2*m) nonnegative;

                                
                %%
                minimize( s + sum(sum((gamma2*eye(m1)).*Qr)) + sqrt(gamma1)*norm(q,2) );
                %%
                subject to

                [s-(c-v)'*X-(Landa1')*(b-A*mu)     1/2*(q(1:m1)+ A1'*(A'*Landa1))' ; 1/2*(q(1:m1) + A1'*(A'*Landa1)) Qr]>=0 ;
                [s-((c-v)'+(v-g)')*X-Landa2'*(b-A*mu)+(v-g)'*mu    1/2*(q(1:m1) +A1'*(A'*Landa2+(v-g)))' ; 1/2*(q(1:m1)+A1'*(A'*Landa2+(v-g))) Qr]>=0;
                
                q(m1 + 1:m) + A1Prime'*(A'*Landa1) == 0;
                q(m1 + 1:m) + A1Prime'*(A'*Landa2 + (v-g)) == 0;
                         

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
    Landa1Opt = Landa1; 
    Landa2Opt = Landa2;

end

