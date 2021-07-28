%**************************************************************************
%      Generating parametes (problem Data) for test problems
%**************************************************************************
clear;
clc;

switch1 = 0; %covariance matrix type (0 = No transformation, 1 = constant, 2 = linear, 3 = exponential)
switch2 = 1; %support type (1= 2sigma, 2 = 3sigma, 3 = 4sigma)
switch3 = 1; %Number of principal components (1 = %10,2=%25,3=%50,4=%75,5= %100)
switch4 = 1; %gamma2 value
gamma1 = 1; % The parameter of moment-based ambiguity set (Zhang's paper)
gamma2_list = [1,2,3,4,5]; % The parameter of moment-based ambiguity set (Zhang's paper)

n = 200; %The size of decision variable x
m = 200; %The size of random vector \xi

m1_list = [.1,.25,.5,.75,1]*m;% The number of principal components

%%%%%%%%%%%%%%%%%%%%%% To solve general moment-based LB of CVaR %%%%%%%
Number = 5; % The number of randomly generated test problems for each case (We consider the avarage performance of all of them)

record = zeros (500,10);% 1=ID 2=Support type  3=m1  4=case number 5= original Value  6=original CPU  7= LB Value 8= LB CPU  9= Relative gap  10= Relative Theoretical gap

ID = 1;%Test Problem Counter

for switch2 = 2 : 2
    
    for iterationNumber = 1 : Number
        
             
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
        b_list={};
        b_list{1} = [2*SDcosi+mu; 2*SDcosi-mu];
        b_list{2} = [3*SDcosi+mu; 3*SDcosi-mu];
        b_list{3} = [4*SDcosi+mu; 4*SDcosi-mu];
        b = b_list{switch2};
        
        %%%%%%%%%% Certain parameters of Newsvendor problem %%%%%%%%%%
        c = zeros(n,1); % to save purchase (wholesale) prices
        v = zeros(n,1); % to save selling (retail) prices
        g = zeros(n,1); % to save salvage prices
        
        for i = 1 : n
            c(i,1) = 0.1*(5+i-1);
            v(i,1) = 0.15*(5+i-1);
            g(i,1) = 0.05*(5+i-1);
        end
        
        for switch4 = 1 : 5
            
            gamma2 = gamma2_list(switch4);
            record(ID,3) = gamma2;
            %%%%%%%%% Solving original and LB problems and saving their results %%%%%%
            [ f_opt1,X_opt1,CPUTime1 ] = OriginalMomentNewsvendorPaper( ZIGMAtemp,mu,A,b,gamma1,gamma2,switch1,c,v,g);

            record(ID,5) = f_opt1;
            record(ID,6) = CPUTime1;

            for switch3 = 1 : 5

                m1 = m1_list(switch3);       
                record(ID,1) = ID;
                record(ID,2) = switch2;
                record(ID,4) = m1;

                [ f_opt2,X_opt2,CPUTime2,Landa1OPT,Landa2OPT] = LBMomentNewsvendorPaper( ZIGMAtemp,mu,m1,A,b,gamma1,gamma2,switch1,c,v,g);
                record(ID,7) = f_opt2;
                record(ID,8) = CPUTime2;

                record(ID,9) = abs(((f_opt1 - f_opt2)/ f_opt1) * 100); % Relative Gap

                TheoreticalGap = TheoGapMomentNewsvendor( Landa1OPT,Landa2OPT,ZIGMAtemp,m1,m,switch1,f_opt1,A,gamma2,v,g);% Relative Theoretical Gap
                record(ID,10) = TheoreticalGap;            
                ID = ID + 1;
            end
        end
    end
end



