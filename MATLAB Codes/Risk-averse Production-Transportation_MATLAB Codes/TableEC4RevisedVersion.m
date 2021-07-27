%**************************************************************************
%      Generating parametes (problem Data) for test problems
%**************************************************************************
clear;
clc;

switch1 = 0; %covariance matrix type (0 = No transformation, 1 = constant, 2 = linear, 3 = exponential)
m=8; % The number of suppliers
n=25; % The number of customersgamma1 = 1
size = m*n; % the size of the problem
gamma1 = 1; % The parameter of moment-based ambiguity set (Zhang's paper)
gamma2 = 2; % The parameter of moment-based ambiguity set (Zhang's paper)
Number = 1; % The number of randomly generated Original DRO problems for each case (We consider the avarage performance of all of them)
R0=30; % distance in Waserestein metric
N = 10; % The number of data points for the random vector of nominal distribution ({\xi^{i}}^{'}) (Gao's paper)
cosiPrime = zeros(size,N); % to save sample \xi^i

record = zeros (100,12);% 1=ID 2= original Value  3=original CPU  4= UB1 Value 5= UB1 CPU  6= UB1 Relative gap  7= UB2 Value   8= UB2 CPU  9= UB2 Relative gap 10= UB3 Value 11= UB3 CPU  12= UB3 Relative gap

ID = 1;%Test Problem Counter

for iterationNumber = 1 : Number
        
        record(ID,1) = ID;

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

        %%%%%%%%%%%%%% To generate N data points (Scenarios) for the random vector of nominal distribution ({\xi^{i}}^{'})%%%%%%%%%%%%
        for i=1:N
            cosiPrime(:,i) = (rand(1,n*m).*dd+0.5*dd)'; % Generating N data points (Scenarios) of random variable \xi
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        mu = mean(Samp)'; % \mu of randome vectore \xi (It creates the mean vector (\mu) from 10000 sample \xi )

        ZIGMAtemp = cov(Samp); % It creates covariance matrix (\Sigma) from 10000 sample \xi

        %%% To approximate the disutility function as a piece-wise linear function%%

        xx=0:0.2:1; % X = [0,.2,.4,.6,.8,1]

        yy=0.25*(exp(2*xx)-1); % Disutility function

        alpha_k=(yy(2:6)-yy(1:5))/0.2; % Slopes of disutility functiuon (a_1, a_2, ... , a_5)

        beta_k=yy(1:5)-alpha_k.*xx(1:5);  % Intercepts of disutility functiuon (b_1, b_2, ... , b_5)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        %%%%%%%%%% Certain parameters of Production problem %%%%%%%%%%%%%%%%%%%

        c=rand(m,1)*mean(mu)+0.5*mean(mu); % Production cost

        d=(rand(n,1)*0.5+0.5)*m/n; % the amount of demand

    
        %%%%%%%%% Solving original and LB problems and saving their results %%%%%%

        record(ID,2) = iterationNumber;

        
        [ f_OPT1,X_OPT1,z_OPT1,CPUTIME1 ] = OriginalCombinedProduction( ZIGMAtemp,mu,N,switch1,R0,gamma2,cosiPrime,c,d,alpha_k,beta_k,m,n);
        record(ID,3) = f_OPT1;
        
        [ f_opt1,X_opt1,z_opt1,CPUTime1 ] = WrittenExamOriginalCombinedProduction( ZIGMAtemp,mu,N,switch1,R0,gamma1,gamma2,cosiPrime,c,d,alpha_k,beta_k,m,n);
        record(ID,4) = f_opt1;
                
        
        m1 =  size/2;
        m2 =  size/2;
        
        [ f_OPT2,X_OPT2,z_OPT2,CPUTIME2] = NewUB2CombinedProductionExp21( ZIGMAtemp,mu,m1,m2,N,switch1,R0,gamma2,cosiPrime,c,d,alpha_k,beta_k,m,n);        
        record(ID,5) = f_OPT2;
        
        [ f_opt2,X_opt2,z_opt2,CPUTime2] = WrittenExamNewUB2CombinedProductionExp21( ZIGMAtemp,mu,m1,m2,N,switch1,R0,gamma1,gamma2,cosiPrime,c,d,alpha_k,beta_k,m,n);        
        record(ID,6) = f_opt2;
        
        m1 =  size/4;
        m2 =  size/4;
        m3 =  size/4;
        m4 =  size/4;

        
        [ f_OPT3,X_OPT3,z_OPT3,CPUTIME3] = NewUB2CombinedProductionExp23( ZIGMAtemp,mu,m1,m2,m3,m4,N,switch1,R0,gamma2,cosiPrime,c,d,alpha_k,beta_k,m,n);
        record(ID,7) = f_OPT3;
        
        [ f_opt3,X_opt3,z_opt3,CPUTime3] = WrittenExamNewUB2CombinedProductionExp23( ZIGMAtemp,mu,m1,m2,m3,m4,N,switch1,R0,gamma1,gamma2,cosiPrime,c,d,alpha_k,beta_k,m,n);
        record(ID,8) = f_opt3;
        
        m1 =  size/5;
        m2 =  size/5;
        m3 =  size/5;
        m4 =  size/5;
        m5 =  size/5;

       
        [ f_OPT4,X_OPT4,z_OPT4,CPUTIME4] = NewUB2CombinedProductionExp24( ZIGMAtemp,mu,m1,m2,m3,m4,m5,N,switch1,R0,gamma2,cosiPrime,c,d,alpha_k,beta_k,m,n);
        record(ID,9) = f_OPT4;
        
        [ f_opt4,X_opt4,z_opt4,CPUTime4] = WrittenExamNewUB2CombinedProductionExp24( ZIGMAtemp,mu,m1,m2,m3,m4,m5,N,switch1,R0,gamma1,gamma2,cosiPrime,c,d,alpha_k,beta_k,m,n);
        record(ID,10) = f_opt4;
                
        ID = ID + 1;

end



            


