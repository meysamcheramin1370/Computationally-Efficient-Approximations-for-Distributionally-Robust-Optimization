%**************************************************************************
%      Generating parametes (problem Data) for test problems
%**************************************************************************
clear;
clc;

switch1 = 0; %covariance matrix type (0 = No transformation, 1 = constant, 2 = linear, 3 = exponential)
switch2 = 1; %The size of random vector \xi and decision variable x  (1= 240, 2 = 320, 3 = 400, 4 = 440, 5 = 480)
m_list = [8,8,8,11,12]; % The number of suppliers
n_list = [30,40,50,40,40]; % The number of customers
N = 10;% The number of principal components; % The number of data points for the random vector of nominal distribution ({\xi^{i}}^{'}) (Gao's paper)
R0=100; % distance in Waserestein metric (Gao's paper - Example7)
Number = 1; % The number of randomly generated Original DRO problems for each case (We consider the avarage performance of all of them)
record = zeros (500,9);% 1=ID  2=m  3=Main problem number  4= LB valu 5= LB CPU  6=UB of Theoretical Gap   7= 4-partition value  8= 4-partition CPU  9=LB- 4partition UB Gap
gamma1 = 1; % The parameter of moment-based ambiguity set (Zhang's paper)
gamma2 = 2; % The parameter of moment-based ambiguity set (Zhang's paper)

ID = 1;%Test Problem Counter

for switch2 = 1 : 5
    
    m = m_list(switch2);
    n = n_list(switch2);
    size = m*n; % the size of the problem
    m1 = 0.25*size; % The number of principal components
    
    cosiPrime = zeros(size,N); % to save sample \xi^i
    
    
    for iterationNumber = 1 : Number
        
        record(ID,1) = ID;
        record(ID,2) = size;
        record(ID,3) = iterationNumber;
        
        
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
        %%%%%%%%%% Certain parameters of Production-Transportation problem %%%%
        
        c=rand(m,1)*mean(mu)+0.5*mean(mu); % Production cost
        
        d=(rand(n,1)*0.5+0.5)*m/n; % the amount of demand
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% Solving original and LB problems and saving their results %%%
        
        
        
        [ f_opt2,X_opt2,z_opt2,CPUTime2,zetaOpt2] = LBCombinedProduction( ZIGMAtemp,mu,m1,N,switch1,R0,gamma2,cosiPrime,c,d,alpha_k,beta_k,m,n);
        
        record(ID,4) = f_opt2;
        record(ID,5) = CPUTime2;
        
        UBTheoreticalGap = UBTheoGapCombinedProduction(z_opt2,zetaOpt2,N,ZIGMAtemp,m1,size,switch1,f_opt2,gamma2,alpha_k);% Upper bound of Relative Theoretical Gap
        record(ID,6) = UBTheoreticalGap;
        
        m1 =  size/4;
        m2 =  size/4;
        m3 =  size/4;
        m4 =  size/4;
        
        [ f_opt4,X_opt4,z_opt4,CPUTime4] = NewUB2CombinedProductionExp23( ZIGMAtemp,mu,m1,m2,m3,m4,N,switch1,R0,gamma2,cosiPrime,c,d,alpha_k,beta_k,m,n);
        record(ID,7) = f_opt4;
        record(ID,8) = CPUTime4;
        
        
        record(ID,9) = abs(((f_opt2 - f_opt4)/ f_opt2) * 100); % Relative Gap between LB and 4-partition UB
        
        ID = ID + 1;
        
    end
end




