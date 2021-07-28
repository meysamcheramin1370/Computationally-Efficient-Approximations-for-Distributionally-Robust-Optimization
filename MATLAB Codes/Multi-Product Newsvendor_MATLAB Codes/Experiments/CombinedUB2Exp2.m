%**************************************************************************
%      Generating parametes (problem Data) for test problems
%**************************************************************************
clear;
clc;

switch1 = 0; %covariance matrix type (0 = No transformation, 1 = constant, 2 = linear, 3 = exponential)
n = 120; %The size of decision variable x
m = 120; %The size of random vector \xi
N = 10; % The number of data points for the random vector of nominal distribution ({\xi^{i}}^{'}) (Gao's paper)
R0=700; % distance in Waserestein metric
cosiPrime = zeros(m,N); % to save sample \xi^i
gamma1 = 1; % The parameter of moment-based ambiguity set (Zhang's paper)
gamma2 = 2; % The parameter of moment-based ambiguity set (Zhang's paper)
Number = 10; % The number of randomly generated Original DRO problems for each case (We consider the avarage performance of all of them)

record = zeros (100,13);% 1=ID 2= original Value  3=original CPU  4= UB1 Value 5= UB1 CPU  6= UB1 Relative gap  7= UB2 Value   8= UB2 CPU  9= UB2 Relative gap 10= UB3 Value 11= UB3 CPU  12= UB3 Relative gap

ID = 1;%Test Problem Counter

for iterationNumber = 1 : Number

        record(ID,1) = ID;
        mu = 10 * rand(m,1); % \mu of randome vectore \xi (\mu \in Uniform[0,10])
        SDcosi = 2 * rand(m,1); % The standard deviation of random vector \xi (uniformly from [0,2])
        %%%%%%%%%%%%%%%To generate covariance matrix \Sigma randomly%%%%%%%%%%%%%%
        MULTI = (SDcosi)*(SDcosi');
        ZIGMA1 = gallery('randcorr',n); % Correlation matrix
        ZIGMAtemp = ZIGMA1.*MULTI; % To change the correlation matrix into covarinace matrix (Refer to the definitions of corelation and covariance)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% Certain parameters of Newsvendor problem %%%%%%%%%%
        c = zeros(n,1); % to save purchase (wholesale) prices
        v = zeros(n,1); % to save selling (retail) prices
        g = zeros(n,1); % to save salvage prices

        for i = 1 : n
            c(i,1) = 0.1*(5+i-1);
            v(i,1) = 0.15*(5+i-1);
            g(i,1) = 0.05*(5+i-1);
        end
        %The following part is only for combined DRO problem.
        %%%%%%%%%%%%%% To generate N data points for the random vector of nominal distribution ({\xi^{i}}^{'})%%%%%%%%%%%%
        for i=1:N
            cosiPrime(:,i)= 10 * rand(m,1) + 2 * rand(m,1);%cosiPrime(:,i)= mu + SDcosi
        end            
               
        %%%%%%%%% Solving original and LB problems and saving their results %%%%%%

        [ f_opt1,X_opt1,CPUTime1 ] = GeneralOriginalCombinedNewsvendorPaper( ZIGMAtemp,mu,N,switch1,R0,cosiPrime,c,v,g,gamma2);
        record(ID,2) = iterationNumber;
        record(ID,3) = f_opt1;
        record(ID,4) = CPUTime1;
                
        m1 =  m/2;
        m2 =  m/2;
        [ f_opt2,X_opt2,CPUTime2] = GeneralNewUB2CombinedNewsvendorPaperExp21( ZIGMAtemp,mu,m1,m2,N,switch1,R0,cosiPrime,c,v,g,gamma2);        
        record(ID,5) = f_opt2;
        record(ID,6) = CPUTime2;
        record(ID,7) = abs(((f_opt2 - f_opt1 )/ f_opt1) * 100); % Relative Gap

        m1 =  m/3;
        m2 =  m/3;
        m3 =  m/3;
        [ f_opt3,X_opt3,CPUTime3] = GeneralNewUB2CombinedNewsvendorPaperExp22( ZIGMAtemp,mu,m1,m2,m3,N,switch1,R0,cosiPrime,c,v,g,gamma2);
        record(ID,8) = f_opt3;
        record(ID,9) = CPUTime3;
        record(ID,10) = abs(((f_opt3 - f_opt1 )/ f_opt1) * 100); % Relative Gap

        m1 =  m/4;
        m2 =  m/4;
        m3 =  m/4;
        m4 =  m/4;
        [ f_opt4,X_opt4,CPUTime4] = GeneralNewUB2CombinedNewsvendorPaperExp23( ZIGMAtemp,mu,m1,m2,m3,m4,N,switch1,R0,cosiPrime,c,v,g,gamma2);
        record(ID,11) = f_opt4;
        record(ID,12) = CPUTime4;
        record(ID,13) = abs(((f_opt4 - f_opt1 )/ f_opt1) * 100); % Relative Gap

        ID = ID + 1;
end



            


