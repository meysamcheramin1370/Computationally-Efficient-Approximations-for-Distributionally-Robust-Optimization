function [ TheoreticalGap ] = ArbitraryTheoGapMomentProduction( z_opt2,Landa1OPT,Landa2OPT,Landa3OPT,Landa4OPT,Landa5OPT,ZIGMAtemp,m1,size,switch1,f_opt1,A,gamma2,alpha_k)
%This function calculates the Relative theoretical upper bound of the GAP between
% moment-based LB and original DRO for the Risk-averse Production_transportation problem
%   XOpt,Landa1OPT,Landa2OPT = Optimal solutions of LP problem
% f_opt1 = The optimal solution value of the original DRO

[A1Complete,A1,A1Prime,ZIGMA,U,delta] = covTransformerDecomposer2(ZIGMAtemp,size,m1,switch1);

z1 = z_opt2(:,1);
z2 = z_opt2(:,2);
z3 = z_opt2(:,3);
z4 = z_opt2(:,4);
z5 = z_opt2(:,5);

summation1 = 0;

for i = m1+1 : size
    summation1 = summation1 + delta(i,i)*((A'*Landa1OPT - alpha_k(1)*z1 )'*U(:,i))^2;
end

summation2 = 0;

for i = m1+1 : size
    summation2 = summation2 + delta(i,i)*((A'*Landa2OPT - alpha_k(2)*z2 )'*U(:,i))^2;
end

summation3 = 0;

for i = m1+1 : size
    summation3 = summation3 + delta(i,i)*((A'*Landa3OPT - alpha_k(3)*z3 )'*U(:,i))^2;
end

summation4 = 0;

for i = m1+1 : size
    summation4 = summation4 + delta(i,i)*((A'*Landa4OPT - alpha_k(4)*z4 )'*U(:,i))^2;
end

summation5 = 0;

for i = m1+1 : size
    summation5 = summation5 + delta(i,i)*((A'*Landa5OPT - alpha_k(5)*z5 )'*U(:,i))^2;
end





TheoreticalGap = abs((sqrt(gamma2)*(sqrt(summation1) + sqrt(summation2) + sqrt(summation3) + sqrt(summation4)+ sqrt(summation5)))/f_opt1)*100;
end

