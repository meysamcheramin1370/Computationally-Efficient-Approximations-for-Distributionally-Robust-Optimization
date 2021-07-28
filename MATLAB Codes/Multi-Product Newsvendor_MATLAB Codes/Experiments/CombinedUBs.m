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
N = 10; % The number of data points for the random vector of nominal distribution ({\xi^{i}}^{'}) (Gao's paper)
gamma2 = 2; 
R0= 700; % distance in Waserestein metric (Gao's paper - Example7)
Number = 8; % The number of randomly generated Original DRO problems for each case (We consider the avarage performance of all of them)
cosiPrime = zeros(m,N); % to save sample \xi^i
record = zeros (100,8);% 1=ID   2=m1  3=case number 4= original Value  5=original CPU  6= UB2 Value 7= UB2 CPU  8= Relative gap

ID = 1;%Test Problem Counter


for iterationNumber = 1 : Number
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %The following part is only for combined DRO problem.
    %%%%%%%%%%%%%% To generate N data points for the random vector of nominal distribution ({\xi^{i}}^{'})%%%%%%%%%%%%
    for i=1:N
        cosiPrime(:,i)= 10 * rand(m,1) + 2 * rand(m,1);%cosiPrime(:,i)= mu + SDcosi
    end
    
    
    %%%%%%%%% Solving original and UB problems and saving their results %%%%%%
    [ f_opt1,X_opt1,CPUTime1] = GeneralOriginalCombinedNewsvendorPaper( ZIGMAtemp,mu,N,switch1,R0,cosiPrime,c,v,g,gamma2);
    
    record(ID,3) = iterationNumber;
    record(ID,4) = f_opt1;
    record(ID,5) = CPUTime1;
    
    for switch3 = 1 : 5
        
        m1 = m1_list(switch3);        
        record(ID,1) = ID;
        record(ID,2) = m1;
        
        [ f_opt3,X_opt3,CPUTime3] = GeneralNewUB2CombinedNewsvendor( ZIGMAtemp,mu,m1,N,switch1,R0,cosiPrime,c,v,g,gamma2);
        record(ID,6) = f_opt3;
        record(ID,7) = CPUTime3;
        record(ID,8) = abs(((f_opt3 - f_opt1 )/ f_opt1) * 100); % Relative Gap
        ID = ID + 1;
        
    end
end




