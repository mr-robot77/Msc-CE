classdef Macro < handle
    properties
        M;       % number of mobile nodes in Macrocell
        N;       % number of Subchannels in Macrocell
        X;       % cell width
        Y;       % cell width
        mnodes;  % mobile nodes in Macrocell
        bs;
        noise;
        Ai;      % Path loss Parameter
        Bi;      % Path loss Parameter
        Ci;      % Path loss Parameter
        fc;      % fequency
        WL;      % Path loss Parameter
        nWL;     % Path loss Parameter
        tsinr;  
        Ru;      % Min Target Rate for each user
        Sm;      % Macrocell constellation size
        nSC;     % number of subchannel that each user needs
        nASC;    % number of assigned subchannel to each user         
        A;       % Subchannel assignmet matrix
        P;       % powre assignmet matrix
        h;       % Path gain matrix
        I;       % Interference matrix
        I_th;       % Interference thereshold in each subchannel
        P_th;       % Max Power thereshold  in each subchannel for each Cell
        sinr;    % Sinr matrix
        TR;      %Throughput
        Cell_TR;    % Cell Throughput
        T;       % Number of assigned SC
        It;      % Iteration counter
        Scel;
        nScel;
        nSU;
        Pmax;
%         Alfa;
%         Teta;
%         Delta;
        Beta;
%         nistar;
%         mistar;
        avg_sinr;
        sum_sinr;  % SINR of each user per Iteration
        cost;
        W;
        Miyu;        %Lagranzh Parameter
        Beta3;         %for computing  Lagranzh Parameter
        
    end
    methods
       
        function Mac=Macro(m,n,x,y,ns,nSU)
            Mac.M=m;
            Mac.N=n;
            Mac.X=x;
            Mac.Y=y;
            Mac.bs=complex((Mac.X/2),(Mac.Y/2));
            Mac.noise=10^(-13);
            Mac.Ai=36;
            Mac.Bi=40;
            Mac.Ci=20;
            Mac.fc=2.5 * (10^9);
            Mac.WL=5;
            Mac.nWL=1;
            Mac.It=1;
            Mac.T=floor(Mac.N/Mac.M);
         
            for i=1 :Mac.M
%                 Mac.tsinr(i)=9.55;
                  Mac.tsinr(i)= 2^(1000000/180000) -1;
            end
            Mac.Sm=4;
            for i=1 :Mac.M
                Mac.nSC(i)=Mac.T;  % Compute nSC based on Ru of each user and Constellation size
            end
            
            for i=1 :Mac.M
                Mac.Pmax(i)= 0.05;
            end            
            
%             Mac.Alfa= 1;
%             Mac.Teta=1;
            Mac.mnodes=test(Mac.M,Mac.X,Mac.Y);
            Mac.nScel=ns;
            Mac.nSU=nSU;
            Mac=CreateSmall(Mac);
            Mac.Beta3(1)=10;
            for n=1 : Mac.N
                Mac.Miyu(n,1)= 1;
            end 
        end
       
        function Mac=CreateSmall(Mac)
          if Mac.nScel>0  
            for i=1 : Mac.nScel
                ttt(i) = Small(Mac.nSU,Mac.N,Mac.X,Mac.Y ,Mac ,i+1);
            end
            Mac.Scel=ttt;
          end
        end
        function Mac=Interference(Mac,t)
            %t= Mac.It +1;
            for i=1 : Mac.M
                for n= 1: Mac.N
                    Mac.I(i,n,t)=Mac.noise;
                    % Interference IntraCell
                    for j=1 : Mac.M
                        if j~= i
                            Mac.I(i,n,t)=Mac.I(i,n,t)+ (Mac.A(j,n,t)* Mac.P(j,n,t)* Mac.h(j,n,1));
                        end
                    end
                    % Interference InterCell
                    for k=1 : Mac.nScel
                        for j=1 : Mac.Scel(k).M
                            Mac.I(i,n,t)=Mac.I(i,n,t)+ Mac.Scel(k).A(j,n,t)* Mac.Scel(k).P(j,n,t)* Mac.Scel(k).h(j,n,1);
                        end  
   
                    end 
                end
            end
        end
        function Mac=SINR(Mac,t)
            Mac.Cell_TR(t)=0;
            for i=1 : Mac.M
                Mac.TR(i,t)=0;
                Mac.sum_sinr(i,t)=0;
                c=0;
                for n=1:Mac.N
                    if Mac.A(i,n,t)==1
                      if t>1   
                        Mac.sinr(i,n,t)= Mac.A(i,n,t)* Mac.P(i,n,t)* Mac.h(i,n,1)/Mac.I(i,n,t-1);
                      else
                        Mac.sinr(i,n,t)= Mac.A(i,n,t)* Mac.P(i,n,t)* Mac.h(i,n,1)/Mac.noise;
                      end
                      Mac.sum_sinr(i,t) = Mac.sum_sinr(t) + Mac.sinr(i,n,t);
                      %if Mac.sinr(i,n,t)>= Mac.tsinr(i)
                          Mac.TR(i,t)=Mac.TR(i,t)+ log2( 1 + Mac.sinr(i,n,t));
                          Mac.Cell_TR(t)=Mac.Cell_TR(t)+ log2( 1 + Mac.sinr(i,n,t));
                      %end
                      c=c+1;
                    end
                end
                Mac.nASC(i)=c;
                Mac.avg_sinr(i,t)= Mac.sum_sinr(i,t)/Mac.nASC(i);
            end
         end
        
        function Mac=PathGain(Mac)
          %Path gain to Macro BS
       
 
            for i = 1 : Mac.M
                u = rand(Mac.N, 1); % generating uniform variates
                sigma = 1; % the parameter
                x = sigma * sqrt(-2 * log(u)); % generating Rayleigh-distributed variates
                for n = 1 : Mac.N
%                    db= (Mac.Ai * log10(abs(Mac.mnodes(i)-Mac.bs ))) + Mac.Bi + (Mac.Ci * log10(Mac.fc/5))+ (Mac.nWL*Mac.WL);
%                    Mac.h(i,n,1)= 10^(-db/10);
                     Mac.h(i,n,1) = x(n) * (abs(Mac.mnodes(i)-Mac.bs )^(-3));
                end
            end
            %Path gain to Femto BSs
            for i = 1 : Mac.M
               for k= 1 : Mac.nScel 
                   u = rand(Mac.N, 1); % generating uniform variates
                   sigma = 1; % the parameter
                   x = sigma * sqrt(-2 * log(u)); % generating Rayleigh-distributed variates
                   for n = 1 : Mac.N
%                         db= (Mac.Ai * log10(abs(Mac.mnodes(i)-Mac.bs ))) + Mac.Bi + (Mac.Ci * log10(Mac.fc/5))+ (Mac.nWL*Mac.WL);
%                         Mac.h(i,n,k+1)=10^(-db/10);
                          Mac.h(i,n,k+1) = x(n) * (abs(Mac.mnodes(i)-Mac.Scel(k).bs )^(-3));
                    end
                end
            end   
            for i=1 : Mac.nScel
               Mac.Scel(i) = PathGain (Mac.Scel(i));
            end
        end
        
        function Mac=Update(Mac,t)  % Start from here
            if t==1
                Mac=Init(Mac);
            else
                Mac=UpdatePower(Mac);
            end
            
            Mac=Interference(Mac,Mac.It) ;
            Mac=SINR(Mac,Mac.It);
            for i=1 : Mac.nScel
                Mac.Scel(i) = Update (Mac.Scel(i));
            end
            for i=1 : Mac.nScel
              Mac.Scel(i)=Update_Parameter(Mac.Scel(i));
            end
            Mac=INT_Notify(Mac);
            Mac.It=Mac.It + 1;
            for i=1 : Mac.nScel
              Mac.Scel(i)=It_Increment(Mac.Scel(i));
            end
    
        end
        
        function Mac=UpdatePower(Mac) % Next step after Update function
            t=Mac.It ;
            for i=1 : Mac.M
                for n=1:Mac.N
                    Mac.cost(i,n,t)=(Mac.tsinr(i)*Mac.I(i,n,t-1))/Mac.h(i,n,1);
                end
            end
            Mac=ScAlloc(Mac,t);
            for i=1 : Mac.M
                for n=1:Mac.N
                    Mac.P(i,n,t)=Mac.cost(i,n,t)*Mac.A(i,n,t);
                end
            end 
            
            for i=1 : Mac.M
                Mac.Beta(i)=0;
                for n=1:Mac.N
                    Mac.P(i,n,t)=(Mac.tsinr(i)*Mac.I(i,n,t-1)*Mac.A(i,n,t))/Mac.h(i,n,1);
                    Mac.Beta(i)= Mac.Beta(i) + (Mac.P(i,n,t)*Mac.A(i,n,t));
                end
                Mac.Beta(i) = Mac.Beta(i)/Mac.Pmax(i);
            end
            
            for i=1 : Mac.M
                if Mac.Beta(i)>1
                    %Mac= SC_FU_maxPower(Mac,i,t);
                    for n=1:Mac.N
                        Mac.P(i,n,t) = Mac.P(i,n,t)/Mac.Beta(i);
                    end
                end
            end
            for i=1 : Mac.nScel

                Mac.Scel(i) = UpdatePower (Mac.Scel(i));
            end
           
        end
        
             
        function Mac=Init(Mac)
            for i=1 : Mac.M
                for n=1:Mac.N
                    Mac.cost(i,n,Mac.It)=0.000001;
                end
            end
            Mac=ScAlloc(Mac,Mac.It);
            for i=1 :Mac.M
                for n=1:Mac.N
                     Mac.P(i,n,Mac.It)=Mac.cost(i,n,Mac.It)*Mac.A(i,n,Mac.It);
                end
            end
             
         
            for i=1 : Mac.nScel
%                Mac.Scel(i) = MaxInterfereceNotify(Mac.Scel(i),Mac.I_th);                
               Mac.Scel(i) = Init (Mac.Scel(i));
            end  
            
        end
        
        function Mac=INT_Notify(Mac)  % Compute Maximmum interference for each subchannel
            t=Mac.It;
            
            for n=1 : Mac.N
                s(n)=0;
                for k=1 : Mac.nScel
                    for i=1 : Mac.Scel(k).M
                      s(n) = s(n) +  (Mac.Scel(k).A(i,n,t)* Mac.Scel(k).P(i,n,t)*Mac.Scel(k).h(i,n,1));
                    end
                end
                
            end
            
            for n= 1: Mac.N
                for i=1 : Mac.M
                    if Mac.A(i,n,t)==1
%                         Mac.I_th(n)= (((Mac.Pmax(i)/Mac.nSC(i))* Mac.h(i,n,1)/Mac.tsinr(i))- Mac.noise);
                          Mac.I_th(n)= ((Mac.P(i,n,t)* Mac.h(i,n,1)/Mac.tsinr(i))- Mac.noise);
                    end
                end

                    
                Mac.Miyu(n,t+1)= Mac.Miyu(n,t)- Mac.Beta3(t)*(Mac.I_th(n)-s(n));
                if Mac.Miyu(n,t+1) <0
                    Mac.Miyu(n,t+1)=0;
                end
            end
            Mac.Beta3(t+1) =10/ sqrt(t+1);
            
%             Mac.Beta3(t+1) =10000/ sqrt(t+1);
            for i=1 : Mac.nScel
                Mac.Scel(i) = MaxInterfereceNotify(Mac.Scel(i),Mac.Miyu);
            end            
            
            
        end

        function Mac=ScAlloc(Mac,t)
            for i=1 : Mac.M
                for j=1 : Mac.N
                    Mac.A(i,j,t)=0;
                end
            end
            
            
            k=0;
            for i=1 : Mac.M
                if i>1 
                  k=Mac.nSC(i-1)+k;
                end
                for j= 1 : Mac.nSC(i)
                    for n=1 : Mac.N
                        cost_matrix((k + j),n)=  Mac.cost(i,n,t);
                    end
                end
            end
            [X,Mac.W(t)]= main(cost_matrix );
             n=1;
            for i=1 : Mac.M
                for j=n :Mac.nSC(i)+n-1
                    Mac.A(i,X(j),t)=1;
                end
                n=n+Mac.nSC(i);
            end
            
       
        end
        
        function disp(Mac)
          t= Mac.It;
          for n= 1: Mac.N
              
              %plot(1:1:t-1, Mac.sinr(1,n,:), 1:1:t-1,Mac.sinr(2,n,:), 1:1:t-1,Mac.sinr(3,n,:),1:1:t-1,Mac.sinr(4,n,:),1:1:t-1,Mac.sinr(5,n,:));
              for i=1 : Mac.M
                  for j=1 : t-1
                      r(i,j)=Mac.sinr(i,n,j);
 
                  end
              end
              for i=1 : Mac.M
                  plot(1:1:t-1, r(i,:));%, 1:1:t-1,r(i,:), 1:1:t-1,r(i,:),1:1:t-1,r(i,:),1:1:t-1,r(i,:));
                  hold;
              end
          end
            plot(1:1:t-1, Mac.sinr(1,1,:), 1:1:t-1,Mac.sinr(1,2,:));
            %, 1:1:t-1,Mac.sinr(3,n,:),1:1:t-1,Mac.sinr(4,n,:),1:1:t-1,Mac.sinr(5,n,:));

            xlabel('Time(Iteration)')
            ylabel('SINR')
            legend('user 1','user 2','user 3','user 4','user 5')
            grid
            
            figure
            plot(1:1:t-1, Mac.P(1,n,:), 1:1:t-1,Mac.P(2,n,:), 1:1:t-1,Mac.P(3,n,:),1:1:t-1,Mac.P(4,n,:),1:1:t-1,Mac.P(5,n,:));
            xlabel('Time(Iteration)')
            ylabel('Power')
            legend('user 1','user 2','user 3','user 4','user 5')
            grid
         % end  
        end
        function calc(Mac)
            k=0;  % index of non supported users
            l=0;    % index of supported users
            for i=1 : Mac.M
                if ((Mac.tsinr(i)-Mac.sinr(i,n,t-1 ))< 0.00001 )
                    l=l+1;
                    sp(l)=i;
                else
                    k=k+1;
                    nsp(k)=i;
                end
            end
            disp(Mac.T(t-1))
            disp(k)
            %disp(nsp)
            disp(l)
            %disp(sp)
            outage= k/(l+k)
            
            %figure

        end
        
        function Mac=Reset_It(Mac)
           Mac.It=1;
            for i=1 : Mac.nScel
                Mac.Scel(i) = Reset_It(Mac.Scel(i));
            end
        end
        function Mac=Set_conf(Mac,Pmx,Tsnr) % Set Maximum power and Target Sinr and other Configurationd for users 
            for i=1 :Mac.M 
                Mac.tsinr(i)=Tsnr;
            end
                        
            for i=1 :Mac.M
                Mac.Pmax(i)= Pmx;
            end       
        end
    end 
end
 
