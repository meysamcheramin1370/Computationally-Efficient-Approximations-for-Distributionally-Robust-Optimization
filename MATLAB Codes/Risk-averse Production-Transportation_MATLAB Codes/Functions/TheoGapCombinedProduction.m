function [ TheoreticalGap ] = TheoGapCombinedProduction(z_opt2,zetaOpt,N,ZIGMAtemp,m1,size,switch1,f_opt1,gamma2,alpha_k)
%This function calculates the Relative theoretical upper bound of the GAP between
% Combined LB and original DRO for the Risk-averse Production_transportation problem
%   XOpt,Landa1OPT,Landa2OPT = Optimal solutions of the Combined LP problem
% f_opt1 = The optimal solution value of the original DRO

[A1Complete,A1,A1Prime,ZIGMA,U,delta] = covTransformerDecomposer(ZIGMAtemp,size,m1,switch1);

summation = 0;

for i = 1 : N
    summation = summation + sqrt(((-alpha_k(1)*z_opt2(:,1)'+zetaOpt(:,i)')*A1Prime)*(((-alpha_k(1)*z_opt2(:,1)'+zetaOpt(:,i)')*A1Prime))')+ sqrt(((-alpha_k(2)*z_opt2(:,2)'+zetaOpt(:,i)')*A1Prime)*(((-alpha_k(2)*z_opt2(:,2)'+zetaOpt(:,i)')*A1Prime))')+ sqrt(((-alpha_k(3)*z_opt2(:,3)'+zetaOpt(:,i)')*A1Prime)*(((-alpha_k(3)*z_opt2(:,3)'+zetaOpt(:,i)')*A1Prime))')+ sqrt(((-alpha_k(4)*z_opt2(:,4)'+zetaOpt(:,i)')*A1Prime)*(((-alpha_k(4)*z_opt2(:,4)'+zetaOpt(:,i)')*A1Prime))')+ sqrt(((-alpha_k(5)*z_opt2(:,5)'+zetaOpt(:,i)')*A1Prime)*(((-alpha_k(5)*z_opt2(:,5)'+zetaOpt(:,i)')*A1Prime))');
end


TheoreticalGap = abs(sqrt(gamma2/N)*(summation/f_opt1))*100;

end

