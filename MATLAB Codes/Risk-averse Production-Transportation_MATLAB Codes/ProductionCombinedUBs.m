%**************************************************************************
%      Generating parametes (problem Data) for test problems
%**************************************************************************
clear;
clc;

switch1 = 0; %covariance matrix type (0 = No transformation, 1 = constant, 2 = linear, 3 = exponential)
switch3 = 1; %Number of principal components (1 = %10,2=%25,3=%50,4=%75,5= %100)
m=8; % The number of suppliers
n=25; % The number of customers
size = m*n; % the size of the problem
m1_list = [.1,.25,.5,.75,1]*size;% The number of principal components
gamma2 = 2; 
R0=30; % distance in Waserestein metric (Gao's paper - Example7)
N = 10; % The number of data points for the random vector of nominal distribution ({\xi^{i}}^{'}) (Gao's paper)
cosiPrime = zeros(size,N); % to save sample \xi^i
Number = 1; % The number of randomly generated Original DRO problems for each case (We consider the avarage performance of all of them)cosiPrime = zeros(m,N); % to save sample \xi^i
record = zeros (100,8);% 1=ID   2=m1  3=case number 4= original Value  5=original CPU  6= UB2 Value 7= UB2 CPU  8= Relative gap

ID = 1;%Test Problem Counter


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
    
    
    %%%%%%%%% Solving original and UB problems and saving their results %%%%%%
    [ f_opt1,X_opt1,z_opt1,CPUTime1 ] = OriginalCombinedProduction( ZIGMAtemp,mu,N,switch1,R0,gamma2,cosiPrime,c,d,alpha_k,beta_k,m,n);
    
    record(ID,3) = iterationNumber;
    record(ID,4) = f_opt1;
    record(ID,5) = CPUTime1;
    
    for switch3 = 1 : 5
        
        m1 = m1_list(switch3);      
        record(ID,1) = ID;
        record(ID,2) = m1;
        
        [ f_opt3,X_opt3,z_opt3,CPUTime3] = NewUB2CombinedProduction( ZIGMAtemp,mu,m1,N,switch1,R0,gamma2,cosiPrime,c,d,alpha_k,beta_k,m,n);
        record(ID,6) = f_opt3;
        record(ID,7) = CPUTime3;
        record(ID,8) = abs(((f_opt3 - f_opt1 )/ f_opt1) * 100); % Relative Gap
        ID = ID + 1;
        
    end
end



