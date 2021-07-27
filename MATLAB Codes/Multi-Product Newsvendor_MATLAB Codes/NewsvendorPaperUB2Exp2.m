%**************************************************************************
%      Generating parametes (problem Data) for test problems
%**************************************************************************
clear;
clc;

switch1 = 0; %covariance matrix type (0 = No transformation, 1 = constant, 2 = linear, 3 = exponential)
n = 120; %The size of decision variable x
m = 120; %The size of random vector \xi

gamma1 = 1; % The parameter of moment-based ambiguity set (Zhang's paper)
gamma2 = 2; % The parameter of moment-based ambiguity set (Zhang's paper)
Number = 10; % The number of randomly generated Original DRO problems for each case (We consider the avarage performance of all of them)

record = zeros (100,12);% 1=ID 2= original Value  3=original CPU  4= UB1 Value 5= UB1 CPU  6= UB1 Relative gap  7= UB2 Value   8= UB2 CPU  9= UB2 Relative gap 10= UB3 Value 11= UB3 CPU  12= UB3 Relative gap

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
        %%%%%%%%%%%%% To randomly generate support S %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I= eye(m);
        A= [I ; -1*I];
        b = [3*SDcosi+mu; 3*SDcosi-mu];
        
        %%%%%%%%%% Certain parameters of Newsvendor problem %%%%%%%%%%
        c = zeros(n,1); % to save purchase (wholesale) prices
        v = zeros(n,1); % to save selling (retail) prices
        g = zeros(n,1); % to save salvage prices

        for i = 1 : n
            c(i,1) = 0.1*(5+i-1);
            v(i,1) = 0.15*(5+i-1);
            g(i,1) = 0.05*(5+i-1);
        end

        %%%%%%%%% Solving original and p-partition UB problems and saving their results %%%%%%

        [ f_opt1,X_opt1,CPUTime1 ] = OriginalMomentNewsvendorPaper( ZIGMAtemp,mu,A,b,gamma1,gamma2,switch1,c,v,g);
        record(ID,2) = f_opt1;
        record(ID,3) = CPUTime1;
                
        m1 =  m/2;
        m2 =  m/2;
        [ f_opt2,X_opt2,CPUTime2] = NewUB2MomentNewsvendorPaperExp21( ZIGMAtemp,mu,A,b,m1,m2,gamma1,gamma2,switch1,c,v,g);
        record(ID,4) = f_opt2;
        record(ID,5) = CPUTime2;
        record(ID,6) = abs(((f_opt2 - f_opt1 )/ f_opt1) * 100); % Relative Gap

        m1 =  m/3;
        m2 =  m/3;
        m3 =  m/3;
        [ f_opt3,X_opt3,CPUTime3] = NewUB2MomentNewsvendorPaperExp22( ZIGMAtemp,mu,A,b,m1,m2,m3,gamma1,gamma2,switch1,c,v,g);
        record(ID,7) = f_opt3;
        record(ID,8) = CPUTime3;
        record(ID,9) = abs(((f_opt3 - f_opt1 )/ f_opt1) * 100); % Relative Gap

        m1 =  m/4;
        m2 =  m/4;
        m3 =  m/4;
        m4 =  m/4;
        [ f_opt4,X_opt4,CPUTime4] = NewUB2MomentNewsvendorPaperExp23( ZIGMAtemp,mu,A,b,m1,m2,m3,m4,gamma1,gamma2,switch1,c,v,g);
        record(ID,10) = f_opt4;
        record(ID,11) = CPUTime4;
        record(ID,12) = abs(((f_opt4 - f_opt1 )/ f_opt1) * 100); % Relative Gap
        ID = ID + 1;

end



            


