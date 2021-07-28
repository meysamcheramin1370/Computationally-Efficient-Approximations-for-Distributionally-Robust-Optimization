% **************************************************************************
%      Generating parametes (problem Data) for test problems
% **************************************************************************
clear;
clc;

switch1 = 0; %covariance matrix type (0 = No transformation, 1 = constant, 2 = linear, 3 = exponential)
switch2 = 1; %support type (1= 2sigma, 2 = 3sigma, 3 = 4sigma)
switch3 = 1; %Number of principal components (1 = %10,2=%25,3=%50,4=%75,5= %100)

gamma1 = 1; % The parameter of moment-based ambiguity set (Zhang's paper)
gamma2 = 2; % The parameter of moment-based ambiguity set (Zhang's paper)

m=8; % The number of suppliers
n=25; % The number of customers

size = m*n; % the size of the problem

m1_list = [.1,.25,.5,.75,1]*size;% The number of principal components




%%%%%%%%%%%%%%%%%%%%%% To solve general moment-based LB of CVaR %%%%%%%
Number = 10; % The number of randomly generated test problems for each case (We consider the avarage performance of all of them)

record = zeros (500,12);% 1=ID 2=Support type  3=m1  4=case number 5= original Value  6=original CPU  7= UB1 Value 8= UB1 CPU  9= UB1 Relative gap  10= UB2 Value   11= UB2 CPU  12= UB2 Relative gap



ID = 1;%Test Problem Counter

for switch2 = 1 : 3
    
    for iterationNumber = 1 : Number
        
              
%%%%%%%%%% Certain parameters of the Risk-avers production-transportation problem %%%%%%%%%%
        
        Faci=rand(2,m); %Meysam: It has two rows because it shows the cordinations of supplier locations
        %(X,Y)
        Dema=rand(2,n);%Meysam: It has two rows because it shows the cordinations of demand locations
        %(X,Y)
        
        
        for i=1:m
            
            for j=1:n
                
                dd((i-1)*n+j)=norm(Faci(:,i)-Dema(:,j)); % It calculates the distance between
                % any pair of demand and
                % facility locations
                %the components of dd are \overhead{\xi_{ij}} for all i and j
            end
            
        end
        
        Samp=zeros(10000,n*m);
        
        for t=1:10000
            
            Samp(t,:)=rand(1,n*m).*dd+0.5*dd; % Generating 10000 samples (Scenarios)
            % of random variable \xi
            
        end
        
        mu = mean(Samp)'; % \mu of randome vectore \xi (It creates the mean vector (\mu) from 10000 sample \xi )
        
        ZIGMAtemp = cov(Samp); % It creates covariance matrix (\Sigma) from 10000 sample \xi
        
        %%% To approximate the disutility function as a piece-wise linear function%%
        
        xx=0:0.2:1; % X = [0,.2,.4,.6,.8,1]
        
        yy=0.25*(exp(2*xx)-1); % Disutility function
        
        alpha_k=(yy(2:6)-yy(1:5))/0.2; % Slopes of disutility functiuon (a_1, a_2, ... , a_5)
        
        beta_k=yy(1:5)-alpha_k.*xx(1:5);  % Intercepts of disutility functiuon (b_1, b_2, ... , b_5)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%% To randomly generate support S %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        SDcosi = std(Samp)'; % (MATLAB) S = std(A): If A is a matrix whose columns are random variables and whose rows are
        % observations, then S is a row vector containing the standard deviations corresponding to each
        % column.
        
        I= eye(size);
        A= [I ; -1*I];
        b_list={};
        b_list{1} = [2*SDcosi+mu; 2*SDcosi-mu];
        b_list{2} = [3*SDcosi+mu; 3*SDcosi-mu];
        b_list{3} = [4*SDcosi+mu; 4*SDcosi-mu];
        b = b_list{switch2};
        
        %%%%%%%%%% Certain parameters of Newsvendor problem %%%%%%%%%%%%%%%
        
        c=rand(m,1)*mean(mu)+0.5*mean(mu); % Production cost
        
        d=(rand(n,1)*0.5+0.5)*m/n; % the amount of demand
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% Solving original and LB problems and saving their results %%%%%%
        [ f_opt1,X_opt1,z_opt1,CPUTime1 ] = OriginalMomentProduction( ZIGMAtemp,mu,A,b,gamma1,gamma2,switch1,c,d,alpha_k,beta_k,m,n);
        record(ID,4) = iterationNumber;
        record(ID,5) = f_opt1;
        record(ID,6) = CPUTime1;
        
        for switch3 = 1 : 5
            
            m1 = m1_list(switch3);        
            record(ID,1) = ID;
            record(ID,2) = switch2;
            record(ID,3) = m1;
            
            [ f_opt2,X_opt2,z_opt2,CPUTime2] = UB1MomentProduction( ZIGMAtemp,mu,m1,A,b,gamma1,gamma2,switch1,c,d,alpha_k,beta_k,m,n);           
            record(ID,7) = f_opt2;
            record(ID,8) = CPUTime2;
            record(ID,9) = abs(((f_opt2 - f_opt1 )/ f_opt1) * 100); % Relative Gap
            
            
            [ f_opt3,X_opt3,z_opt3,CPUTime3 ] = NewUB2MomentProduction( ZIGMAtemp,mu,A,b,m1,gamma1,gamma2,switch1,c,d,alpha_k,beta_k,m,n);
            record(ID,10) = f_opt3;
            record(ID,11) = CPUTime3;
            record(ID,12) = abs(((f_opt3 - f_opt1 )/ f_opt1) * 100); % Relative Gap
            ID = ID + 1;
         end
    end
end



