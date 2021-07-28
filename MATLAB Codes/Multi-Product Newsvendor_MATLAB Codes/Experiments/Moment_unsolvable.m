%**************************************************************************
%      Generating parametes (problem Data) for test problems
%**************************************************************************
clear;
clc;

switch1 = 0; %covariance matrix type (0 = No transformation, 1 = constant, 2 = linear, 3 = exponential)
switch2 = 1; %The size of random vector \xi and decision variable x  (1= 700, 2 = 800, 3 = 900, 4 = 1000)
n_list = [300,400,500,600,700]; %The size of decision variable x
m_list = [300,400,500,600,700]; %The size of random vector \xi
gamma1 = 1; % The parameter of moment-based ambiguity set (Zhang's paper)
gamma2 = 2; % The parameter of moment-based ambiguity set (Zhang's paper)

%%%%%%%%%%%%%%%%%%%%%% To solve general moment-based LB of CVaR %%%%%%%
Number = 10; % The number of randomly generated Original DRO problems for each case (We consider the avarage performance of all of them)

record = zeros (500,12);% 1=ID  2=m  3=Main problem number  4 = LB valu 5= LB CPU  6=UB of Theoretical Gap  7= UB1 Value 8= UB1 CPU  9= 4-partition value  10= 4-partition CPU  11= LB-UB1 Gap  12=LB- 4partition UB Gap

ID = 1;%Test Problem Counter



for switch2 = 3 : 5
    
    m = m_list(switch2);
    n = n_list(switch2);
    m1 = 0.25*m; % The number of principal components
    
    record(ID,1) = ID;
    record(ID,2) = m;
    
    
    for iterationNumber = 1 : Number
        
        record(ID,3) = iterationNumber;
        
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
        
        %%%%%%%%% Solving original and LB problems and saving their results %%%%%%
        
        
        [ f_opt2,X_opt2,CPUTime2,Landa1OPT,Landa2OPT] = LBMomentNewsvendorPaper( ZIGMAtemp,mu,m1,A,b,gamma1,gamma2,switch1,c,v,g);
        record(ID,4) = f_opt2;
        record(ID,5) = CPUTime2;
        
        UBTheoreticalGap = UBTheoGapMomentNewsvendor( Landa1OPT,Landa2OPT,ZIGMAtemp,m1,m,switch1,f_opt2,A,gamma2,v,g);% Upper bound of Relative Theoretical Gap
        record(ID,6) = UBTheoreticalGap;
        
        [ f_opt3,X_opt3,CPUTime3] = UB1MomentNewsvendorPaper( ZIGMAtemp,mu,m1,A,b,gamma1,gamma2,switch1,c,v,g);
        record(ID,7) = f_opt3;
        record(ID,8) = CPUTime3;
        
        record(ID,9) = abs(((f_opt2 - f_opt3)/ f_opt2) * 100); % Relative Gap between LB and the first UB
        
        
        if m <= 700
            
            m1 =  m/4;
            m2 =  m/4;
            m3 =  m/4;
            m4 =  m/4;
            [ f_opt4,X_opt4,CPUTime4] = NewUB2MomentNewsvendorPaperExp23( ZIGMAtemp,mu,A,b,m1,m2,m3,m4,gamma1,gamma2,switch1,c,v,g);
            record(ID,10) = f_opt4;
            record(ID,11) = CPUTime4;
            
            
            record(ID,12) = abs(((f_opt2 - f_opt4)/ f_opt2) * 100); % Relative Gap between LB and 4-partition UB
            
        end
        
        
        ID = ID + 1;
        
    end
end




