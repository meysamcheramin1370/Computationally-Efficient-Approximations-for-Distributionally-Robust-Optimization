function [ f_opt,X_opt,CPUTime] = GeneralOriginalCombinedNewsvendorPaper( ZIGMAtemp,mu,N,switch1,R0,cosiPrime,c,v,g,gamma2)
%This function finds the original OPT value (Combined ambiguity
%set)for Multiproduct Newsvendor problem.


m=length(mu);
n = m;
m1 = m; % This parameter is not used. its just for fullfilling m1 in covTransformerDecomposer. 

[A1Complete,A1,A1Prime,ZIGMA,U,delta ] = covTransformerDecomposer(ZIGMAtemp,m,m1,switch1);

tic 

cvx_begin sdp
 
                variable X(n) nonnegative;
                variables zeta(m,N) y(N);
                variable Q(m,m) symmetric;
                variable Landa  nonnegative;
               
                               
                %%
                minimize( Landa*R0 + gamma2 * sum(sum(eye(m).*Q)) + (1/N) * sum(y) );
                %%
                subject to
                
                for i = 1 : N

                    [Q     1/2*(zeta(:,i)'*A1Complete)' ; 1/2*(zeta(:,i)'*A1Complete)  y(i)-(c-v)'*X+zeta(:,i)'*(mu - cosiPrime(:,i))]>=0 ;
                    [Q     1/2*(((v-g)' + zeta(:,i)')*A1Complete)' ; 1/2*((v-g)' + zeta(:,i)')*A1Complete y(i)+(v-g)'*mu - ((c-v)'+(v-g)')*X + zeta(:,i)'*(mu - cosiPrime(:,i))]>=0 ;

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
    zetaOpt = zeta; 
    

end

