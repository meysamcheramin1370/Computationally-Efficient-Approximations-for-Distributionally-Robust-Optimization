function [ TheoreticalGap ] = TheoGapMomentNewsvendor( Landa1OPT,Landa2OPT,ZIGMAtemp,m1,m,switch1,f_opt1,A,gamma2,v,g)
%This function calculates the Relative theoretical upper bound of the GAP between
% moment-based LB and original DRO for the Multi-product Newsvendor problem
%   XOpt,Landa1OPT,Landa2OPT = Optimal solutions of LP problem
% f_opt1 = The optimal solution value of the original DRO

[A1Complete,A1,A1Prime,ZIGMA,U,delta] = covTransformerDecomposer(ZIGMAtemp,m,m1,switch1);

summation1 = 0;

for i = m1+1 : m
    summation1 = summation1 + delta(i,i)*((A'*Landa1OPT)'*U(:,i))^2;
end

summation2 = 0;

for i = m1+1 : m
    summation2 = summation2 + delta(i,i)*((A'*Landa2OPT+(v-g))'*U(:,i))^2;
end



TheoreticalGap = abs((sqrt(gamma2)*(sqrt(summation1) + sqrt(summation2)))/f_opt1)*100;
end

