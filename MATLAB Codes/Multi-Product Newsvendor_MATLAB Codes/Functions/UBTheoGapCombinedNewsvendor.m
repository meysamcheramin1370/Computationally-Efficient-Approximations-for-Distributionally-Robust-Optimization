function [ TheoreticalGap ] = UBTheoGapCombinedNewsvendor(zetaOpt,N,ZIGMAtemp,m1,m,switch1,f_opt2,v,g)
%This function calculates the Relative theoretical upper bound of the GAP between
% Combined LB and original DRO
%   XOpt,Landa1OPT,Landa2OPT = Optimal solutions of the Combined LP problem
% f_opt2 = The optimal solution value of the original DRO

[A1Complete,A1,A1Prime,ZIGMA,U,delta] = covTransformerDecomposer(ZIGMAtemp,m,m1,switch1);

summation = 0;

for i = 1 : N
    summation = summation + sqrt((zetaOpt(:,i)'*A1Prime)*(zetaOpt(:,i)'*A1Prime)') + sqrt((((v-g)'+zetaOpt(:,i)')*A1Prime)*(((v-g)'+zetaOpt(:,i)')*A1Prime)');
end


TheoreticalGap = abs((1/sqrt(N))*(summation/f_opt2))*100;

end

