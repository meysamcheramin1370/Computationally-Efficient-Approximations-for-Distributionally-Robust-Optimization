%**************************************************************************
%      Generating parametes (problem Data) for test problems
%**************************************************************************
clear;
clc;

switch1 = 0; %covariance matrix type (0 = No transformation, 1 = constant, 2 = linear, 3 = exponential)
switch3 = 1; %Number of principal components (1 = %10,2=%25,3=%50,4=%75,5= %100)
n = 200; %The size of decision variable x
m = 200; %The size of random vector \xi
m1_list = [.1,.25,.5,.75,1]*m;% The number of principal components
N = 5; % The number of data points for the random vector of nominal distribution ({\xi^{i}}^{'}) (Gao's paper)
R0_list = [700,800,900,1000]; % distance in Waserestein metric (Gao's paper - Example7)
cosiPrime = zeros(m,N); % to save sample \xi^i
Number = 2; % The number of randomly generated Original DRO problems for each case (We consider the avarage performance of all of them)
record = zeros (200,10);% 1=ID   2=Main problem number  3=R_0 4 = m1 5= original Value  6=original CPU  7= LB Value 8= LB CPU  9= Relative gap  10= Relative Theoretical gap
ID = 1;%Test Problem Counter


for iterationNumber = 1 : Number
    
    
    mu = 10 * rand(m,1); % \mu of randome vectore \xi (\mu \in Uniform[0,10])
    SDcosi = 2 * rand(m,1); % The standard deviation of random vector \xi (uniformly from [0,2])
    %%%%%%%%%%%%%%%To generate covariance matrix \Sigma randomly%%%%%%%%%%%%%%
    MULTI = (SDcosi)*(SDcosi');
    ZIGMA1 = gallery('randcorr',n); % Correlation matrix
    ZIGMAtemp = ZIGMA1.*MULTI; % To change the correlation matrix into covarinace matrix (Refer to the definitions of corelation and covariance)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %The following part is only for combined DRO problem.
    %%%%%%%%%%%%%% To generate N data points for the random vector of nominal distribution ({\xi^{i}}^{'})%%%%%%%%%%%%
    for i=1:N
        cosiPrime(:,i)= 10 * rand(m,1) + 2 * rand(m,1); %cosiPrime(:,i)= mu + SDcosi;
    end
    
    %%%%%%%%%% Certain parameters of Newsvendor problem %%%%%%%%%%
    c = zeros(n,1); % to save purchase (wholesale) prices
    v = zeros(n,1); % to save selling (retail) prices
    g = zeros(n,1); % to save salvage prices
    
    for i = 1 : n
        c(i,1) = 0.1*(5+i-1);
        v(i,1) = 0.15*(5+i-1);
        g(i,1) = 0.05*(5+i-1);
    end
    
    
    for switch2 = 1 : 4
        
        R0 = R0_list(switch2);
        record(ID,3) = R0;
        
        %%%%%%%%% Solving original and LB problems and saving their results %%%%%%
        [ f_opt1,X_opt1,CPUTime1 ] = OriginalCombinedNewsvendorPaper( ZIGMAtemp,mu,N,switch1,R0,cosiPrime,c,v,g);
        record(ID,2) = iterationNumber;
        record(ID,5) = f_opt1;
        record(ID,6) = CPUTime1;
        
        for switch3 = 1 : 5
            
            m1 = m1_list(switch3);
            record(ID,1) = ID;
            record(ID,4) = m1;
            
            
            [ f_opt2,X_opt2,CPUTime2,zetaOpt] = LBCombinedNewsvendorPaper( ZIGMAtemp,mu,m1,N,switch1,R0,cosiPrime,c,v,g);
            record(ID,7) = f_opt2;
            record(ID,8) = CPUTime2;
            record(ID,9) = abs(((f_opt1 - f_opt2)/ f_opt1) * 100); % Relative Gap
            TheoreticalGap = TheoGapCombinedNewsvendor(zetaOpt,N,ZIGMAtemp,m1,m,switch1,f_opt1,v,g);% Relative Theoretical Gap
            record(ID,10) = TheoreticalGap;
            
            ID = ID + 1;
            
        end
        
    end
end




