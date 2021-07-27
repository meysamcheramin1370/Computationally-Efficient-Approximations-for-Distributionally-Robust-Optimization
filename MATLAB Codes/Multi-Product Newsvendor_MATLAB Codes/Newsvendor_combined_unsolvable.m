%**************************************************************************
%      Generating parametes (problem Data) for test problems
%**************************************************************************
clear;
clc;

switch1 = 0; %covariance matrix type (0 = No transformation, 1 = constant, 2 = linear, 3 = exponential)
switch2 = 1; %The size of random vector \xi and decision variable x  (1= 240, 2 = 320, 3 = 400, 4 = 440, 5 = 480)
n_list = [240,320,400,440,480]; %The size of decision variable x
m_list = [240,320,400,440,480]; %The size of random vector \xi
N = 10;% The number of principal components; % The number of data points for the random vector of nominal distribution ({\xi^{i}}^{'}) (Gao's paper)
R0=1700; % distance in Waserestein metric (Gao's paper - Example7)
Number = 9; % The number of randomly generated Original DRO problems for each case (We consider the avarage performance of all of them)
record = zeros (500,9);% 1=ID  2=m  3=Main problem number  4= LB valu 5= LB CPU  6=UB of Theoretical Gap   7= 4-partition value  8= 4-partition CPU  9=LB- 4partition UB Gap
gamma1 = 1; % The parameter of moment-based ambiguity set (Zhang's paper)
gamma2 = 2; % The parameter of moment-based ambiguity set (Zhang's paper)

ID = 1;%Test Problem Counter

for switch2 = 1 : 4
    
    m = m_list(switch2);
    n = n_list(switch2);
    m1 = 0.25*m; % The number of principal components
    
    cosiPrime = zeros(m,N); % to save sample \xi^i
    
    
    for iterationNumber = 1 : Number
        
        record(ID,1) = ID;
        record(ID,2) = m;        
        record(ID,3) = iterationNumber;

        
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
        
        
                
        %%%%%%%%% Solving original and LB problems and saving their results %%%%%%
              
                 
                        
            [ f_opt2,X_opt2,CPUTime2,zetaOpt] = GeneralLBCombinedNewsvendorPaper(ZIGMAtemp,mu,m1,N,switch1,R0,cosiPrime,c,v,g,gamma2);%This function finds the lower bound (Combined ambiguity;
            
            record(ID,4) = f_opt2;
            record(ID,5) = CPUTime2;
            
            UBTheoreticalGap = UBTheoGapCombinedNewsvendor(zetaOpt,N,ZIGMAtemp,m1,m,switch1,f_opt2,v,g);% Upper bound of Relative Theoretical Gap
            record(ID,6) = UBTheoreticalGap;
            
            m1 =  m/4;
            m2 =  m/4;
            m3 =  m/4;
            m4 =  m/4;
            
            [ f_opt4,X_opt4,CPUTime4] = GeneralNewUB2CombinedNewsvendorPaperExp23( ZIGMAtemp,mu,m1,m2,m3,m4,N,switch1,R0,cosiPrime,c,v,g,gamma2);
            record(ID,7) = f_opt4;
            record(ID,8) = CPUTime4;
        
        
            record(ID,9) = abs(((f_opt2 - f_opt4)/ f_opt2) * 100); % Relative Gap between LB and 4-partition UB
            
            ID = ID + 1;
        
    end
end




