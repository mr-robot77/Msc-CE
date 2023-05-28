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
        Sm;      % Ma crocell constellation size
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
        Fcel;
        nFcel;
        nFU;
        Pmax;
        QAM_Rate;
        QAM_Sinr;
        MinRate;
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
       
        function Mac=Macro(m,n,x,y,nf,nfu)
            Mac.M=m;
            Mac.N=n;
            Mac.X=x;
            Mac.Y=y;
            Mac.bs=complex((Mac.X/2),(Mac.Y/2));
            Mac.noise=10^(-13);
            Mac.Ai=36;
            Mac.Bi=40;
            Mac.Ci=20;
            Mac.fc=2 * (10^9);
            Mac.WL=5;
            Mac.nWL=1;
            Mac.It=1;
            Mac.T=floor(Mac.N/Mac.M);
             Mac.QAM_Rate(1)=2;    %4QAM
            Mac.QAM_Rate(2)=4;    %16QAM
            Mac.QAM_Rate(3)=6;    %64QAM
            Mac.QAM_Rate(4)=8;    %256QAM
            Mac.QAM_Rate(5)=10;   %1024QAM
            Mac.QAM_Sinr(1)= 9.55;
            Mac.QAM_Sinr(2)= 45.11;
            Mac.QAM_Sinr(3)= 179.85;
            Mac.QAM_Sinr(4)= 694.17;               
            Mac.QAM_Sinr(5)= 2667.32;            
         
            for i=1 :Mac.M
                 Mac.Ru(i)=1024000/180000;
            end
 
         
            Mac.Sm=4;
            for i=1 :Mac.M
                Mac.nSC(i)=  ceil(Mac.Ru(i)/Mac.QAM_Rate(2));   % Compute nSC based on Ru of each user and Constellation size
            end
            for i=1 :Mac.M
                Mac.tsinr(i)=2^(Mac.Ru(i)/Mac.nSC(i))-1;
                
            end  
            for i=1 :Mac.M
                Mac.Pmax(i)= 0.05;
            end            
 
%             Mac.Alfa= 1;
%             Mac.Teta=1;
            Mac.mnodes=test(Mac.M,Mac.X,Mac.Y);
            Mac.nFcel=nf;
            Mac.nFU=nfu;
            Mac=CreateFemto(Mac);
            Mac.Beta3(1)=0.1;
            for n=1 : Mac.N
                Mac.Miyu(n,1)= 0.000001;
            end 
        end
       
        function Mac=CreateFemto(Mac)
          if Mac.nFcel>0  
            for i=1 : Mac.nFcel
                ttt(i) = Femto(Mac.nFU,Mac.N,Mac.X,Mac.Y ,Mac ,i+1);
            end
            Mac.Fcel=ttt;
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
                    for k=1 : Mac.nFcel
                        for j=1 : Mac.Fcel(k).M
                            Mac.I(i,n,t)=Mac.I(i,n,t)+ Mac.Fcel(k).A(j,n,t)* Mac.Fcel(k).P(j,n,t)* Mac.Fcel(k).h(j,n,1);
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
               for k= 1 : Mac.nFcel 
                   u = rand(Mac.N, 1); % generating uniform variates
                   sigma = 1; % the parameter
                   x = sigma * sqrt(-2 * log(u)); % generating Rayleigh-distributed variates
                   for n = 1 : Mac.N
%                         db= (Mac.Ai * log10(abs(Mac.mnodes(i)-Mac.bs ))) + Mac.Bi + (Mac.Ci * log10(Mac.fc/5))+ (Mac.nWL*Mac.WL);
%                         Mac.h(i,n,k+1)=10^(-db/10);
                          Mac.h(i,n,k+1) = x(n) * (abs(Mac.mnodes(i)-Mac.Fcel(k).bs )^(-3));
                    end
                end
            end   
            for i=1 : Mac.nFcel
               Mac.Fcel(i) = PathGain (Mac.Fcel(i));
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
            for i=1 : Mac.nFcel
                Mac.Fcel(i) = Update (Mac.Fcel(i));
            end
            for i=1 : Mac.nFcel
              Mac.Fcel(i)=Update_Parameter(Mac.Fcel(i));
            end
            Mac=INT_Notify(Mac);
            Mac.It=Mac.It + 1;
            for i=1 : Mac.nFcel
              Mac.Fcel(i)=It_Increment(Mac.Fcel(i));
            end
    
        end
        
        function Mac=UpdatePower(Mac) % Next step after Update function
           
            t=Mac.It ;   

            
             Mac= Macro_Update(Mac,t);  %compute power and sub-channel allocation through linearization
            for i=1 : Mac.nFcel

                Mac.Fcel(i) = UpdatePower (Mac.Fcel(i));
            end
           
        end
        
             
        function Mac=Init(Mac)
            
            Mac= Macro_Update(Mac,1);
         
            for i=1 : Mac.nFcel
%                Mac.Fcel(i) = MaxInterfereceNotify(Mac.Fcel(i),Mac.I_th);                
               Mac.Fcel(i) = Init (Mac.Fcel(i));
            end  
            
        end
        
        function Mac=INT_Notify_(Mac,If)  % Compute Maximmum interference for each subchannel
            t=Mac.It;
            for n= 1: Mac.N
              Sc(n,1)=0;
              Sc(n,2)=0;
            end
            for n= 1: Mac.N
                for i=1 : Mac.M
                    if Mac.A(i,n,t)==1
                        Sc(n,1)=1;
                        Sc(n,2)=i;

                    end
                    
                end
            end
            for n= 1: Mac.N
              if(Sc(n,1)==1 )
                Mac.I_th(n)= If(Sc(n,2),n);
              else
                  Mac.I_th(n)=Mac.Pmax(1)*1000000;
              end
            
            end
           end
        
        function Mac=INT_Notify(Mac)  % Compute Maximmum interference for each subchannel
            t=Mac.It;
            
            for n=1 : Mac.N
                s(n)=0;
                for k=1 : Mac.nFcel
                    for i=1 : Mac.Fcel(k).M
                      s(n) = s(n) +  (Mac.Fcel(k).A(i,n,t)* Mac.Fcel(k).P(i,n,t)*Mac.Fcel(k).h(i,n,1));
                    end
                end
                
            end
            
            
            
            for n= 1: Mac.N

%                 for i=1 : Mac.M
%                     if Mac.A(i,n,t)==1
% %                             Mac.I_th(n)= (((Mac.Pmax(i)/Mac.nSC(i))* Mac.h(i,n,1)/Mac.tsinr(i))- Mac.noise);
%                            Mac.I_th(n)= ((Mac.P(i,n,t)* Mac.h(i,n,1)/Mac.tsinr(i))- Mac.noise);
% %                           Mac.I_th(n)= (Mac.P(Sc(n,2),n,t)* Mac.h(Sc(n,2),n,1)/Mac.tsinr(Sc(n,2)))- Mac.noise;                          
%                     end
%                  end
                
                Mac.Miyu(n,t+1)= Mac.Miyu(n,t)- Mac.Beta3(t)*(Mac.I_th(n)-s(n));
                if Mac.Miyu(n,t+1) <0
                    Mac.Miyu(n,t+1)=0;
                end
            end
            Mac.Beta3(t+1) =1/ sqrt(t+1);
            
%             Mac.Beta3(t+1) =10000/ sqrt(t+1);
            for i=1 : Mac.nFcel
                Mac.Fcel(i) = MaxInterfereceNotify(Mac.Fcel(i),Mac.Miyu);
            end            
            
            
        end

        function Mac=Macro_Update(Mac,t)
            %param(1)= 1 ; %number of Macrocells
            param(1)= Mac.M;    %number of Macr users
            param(2)=Mac.N;  %number of subchannels
            
               
                for i=1 : Mac.M
                    for n=1 : Mac.N
                       p(i,n)=0;
                       If(i,n)=0;
                       a(i,n)=0;
%                         eff_I(i,n)= Mac.eff_I(i,n); % effectice interference of each user in each subchannel
                        Mbs_h(i,n)= Mac.h(i,n,1);   % Path gain of each user to Macro Bs
                    end
                end
           
                [If,p,a]=mex_macro(param,Mac.Pmax(1),Mac.noise,Mbs_h,Mac.tsinr,Mac.nSC);
%                 [If,p,a]=mex_macro(param,Mac.Pmax(1),Mac.noise,Mbs_h,Mac.Ru,Mac.nSC);                
%                mex_macro(param,Mac.Pmax(1),Mac.noise,Mbs_h,Mac.tsinr);

         
             for i=1 : Mac.M
               for n=1 : Mac.N
                     Mac.P(i,n,t)= p(i,n);
                     Mac.A(i,n,t)=a(i,n);
               end
             end
          
          
         Mac=INT_Notify_(Mac,If);           
        end
        
        
 
 
        
        function Mac=Reset_It(Mac)
           Mac.It=1;
            for i=1 : Mac.nFcel
                Mac.Fcel(i) = Reset_It(Mac.Fcel(i));
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
        
        function Mac=Set_conf2(Mac,macro_MinRate, qam) % Set Maximum MinRate and other Configurationd for users
            for i=1 :Mac.M
                 Mac.MinRate(i)=macro_MinRate;
            end
         
            for i=1 :Mac.M
                 Mac.Ru(i)=Mac.MinRate(i)/180000;
            end
 
         
            Mac.Sm=4;
            for i=1 :Mac.M
                Mac.nSC(i)=  ceil(Mac.Ru(i)/Mac.QAM_Rate(qam));   % Compute nSC based on Ru of each user and Constellation size
            end
            for i=1 :Mac.M
                Mac.tsinr(i)=2^(Mac.Ru(i)/Mac.nSC(i))-1;
            end             
        end
    end 
end
 
