classdef Femto < handle
    properties
        M;
        N;
        X;
        Y;
        mnodes;
        bs;
        noise;
        Ai;
        Bi;
        Ci;
        fc;
        WL;
        nWL;
        tsinr;
        sf;           % femtocell constellation size
        nSC;          % number of subchannel that each user needs
        nASC;         % number of assigned subchannel to each user
        A;
        P;
        h;
        I;
        sinr;
        TR;
        Cell_TR;      % Cell Throughput
        SP_Eff;
        T;
        It;
        Mcel;
        ID;
        Pmax;
        Pmax_SC;     % Max power for each Subchannel
        I_th;        % Interference thereshold in each subchannel
        w;
        W;
        sum_sinr;    % SINR of each user per Iteration
        avg_sinr;
        Lambda;      %Lagranzh Parameter
        Nou;         %Lagranzh Parameter
        Miyu;        %Lagranzh Parameter
        Beta;         %for computing  Lagranzh Parameter
        DT_U;        % Delay Tolerent Users
        DS_U;        % Delay sensitive Users
        Ru;          % Minimum Rate for delay sensitive users;
        H;           % for Compute SC allocation
        
    end
    methods
        function F=Femto(m,n,x,y,Mac,id)
            F.M=m;
            F.N=n;
            F.X=x;
            F.Y=y;
            F.noise=10^(-13);
            F.Ai=25;
            F.Bi=45;
            F.Ci=20;
            F.fc=2 * (10^9);
            F.WL=5;
            F.nWL=1;
            F.It=1;
            for i=1 :F.M
                F.tsinr(i)=179.85;
            end
            F.sf= 64;
            for i=1 :F.M
                F.Pmax(i)= 0.05;
            end
            
            
            F.T=floor(F.N/F.M);
            F.Mcel=Mac;
            F.ID=id;
            F=Create_BS_UE(F);
            
%             for i=1 : (F.M/2)
%                 F.DS_U(i)=i;
%             end
%             
            for i=1 : (F.M)
                F.DT_U(i)= i; % + (F.M/2);
            end
            for i=1 :F.M
                F.Lambda(i,1)= 0.1;
                
            end
            for l=1 : length(F.DS_U)
                i=F.DS_U(l);
                F.Nou(i,1)= 0.1;
            end
            for l=1 : length(F.DT_U)
                i=F.DT_U(l);
                F.Nou(i,1)= 0;
            end
            for i=1 :F.M
                F.Ru(i)= 9;
            end
            F.Beta(1,1)=10;
            F.Beta(2,1)=10;
            
        end
        
        function F=Create_BS_UE(F)
            bsx=unifrnd(1, F.X);
            bsy=unifrnd(1, F.Y);
            rx = unifrnd(1:F.M, 20);
            ry = unifrnd(1:F.M, 20);
            for j=1 : F.M
                rx(j)= rx(j)+ bsx;
                ry(j)= ry(j)+ bsy;
            end
            for j=1 : F.M
                F.mnodes(j)= complex(rx(j), ry(j));
            end
            F.bs = complex(bsx+10,bsy+10);
            
        end
        function F=PathGain(F)
            %Path gain to Macro BS
            for i = 1 : F.M
                u = rand(F.N, 1); % generating uniform variates
                sigma = 1; % the parameter
                x = sigma * sqrt(-2 * log(u)); % generating Rayleigh-distributed variates
                for n = 1 : F.N
                    %                     db=((F.Mcel.Ai * log10(abs(F.mnodes(i)-F.Mcel.bs ))) + F.Mcel.Bi + (F.Mcel.Ci * log10(F.Mcel.fc/5))+ (F.Mcel.nWL*F.Mcel.WL));
                    %                     F.h(i,n,1)=10^(-db/10);
                    F.h(i,n,1) = x(n) * (abs(F.mnodes(i)-F.Mcel.bs )^(-3));
                end
            end
            %Path gain to Femto BSs
            for i = 1 : F.M
                for k= 1 : F.Mcel.nFcel
                    u = rand(F.N, 1); % generating uniform variates
                    sigma = 1; % the parameter
                    x = sigma * sqrt(-2 * log(u)); % generating Rayleigh-distributed variates
                    for n = 1 : F.N
                        % db=((F.Mcel.Fcel(k).Ai * log10(abs(F.mnodes(i)-F.Mcel.Fcel(k).bs ))) + F.Mcel.Fcel(k).Bi + (F.Mcel.Fcel(k).Ci * log10(F.Mcel.Fcel(k).fc/5))+ (F.Mcel.Fcel(k).nWL*F.Mcel.Fcel(k).WL));
                        %  F.h(i,n,k+1)= 10^(-db/10);
                        F.h(i,n,k+1) = x(n) * (abs(F.mnodes(i)-F.Mcel.Fcel(k).bs )^(-3));
                    end
                end
            end
        end
        
        function F=Interference(F,t)
            %t= Mac.It +1;
            for i=1 : F.M
                for n= 1: F.N
                    F.I(i,n,t)=F.noise;
                    %                     Interference IntraCell
%                     for j=1 : F.M
%                         if j~= i
%                             F.I(i,n,t)=F.I(i,n,t)+ (F.A(j,n,t)* F.P(j,n,t)* F.h(j,n,F.ID));
%                         end
%                     end
                    %                     Interference InterCell from Femtocells
%                     for k=1 : F.Mcel.nFcel
%                         if( k+1 )~= F.ID
%                             for j=1 : F.Mcel.Fcel(k).M
%                                 %if F.Mcel.Fcel(k).It>t
%                                   F.I(i,n,t)=F.I(i,n,t)+ F.Mcel.Fcel(k).A(j,n,t)* F.Mcel.Fcel(k).P(j,n,t)* F.Mcel.Fcel(k).h(j,n,F.ID);
%                                 %end
%                             end
%                         end
%                     end
%                     
                    % Interference InterCell from Macrocell
                    for j=1 : F.Mcel.M
                        F.I(i,n,t)=F.I(i,n,t)+ F.Mcel.A(j,n,t)* F.Mcel.P(j,n,t)* F.Mcel.h(j,n,F.ID);
                    end
                    
                end
            end
        end
        
        function F=SINR(F,t)
            F.Cell_TR(t)=0;
            for i=1 : F.M
                F.TR(i,t)=0;
                F.SP_Eff(i,t)=0;
                c=0;
                F.sum_sinr(i,t)=0;
                for n=1:F.N
                    if F.A(i,n,t)==1
                        if t> 1
                          F.sinr(i,n,t)= F.P(i,n,t)*F.A(i,n,t)* F.h(i,n,F.ID)/F.I(i,n,t-1);
                        else
                          F.sinr(i,n,t)= F.P(i,n,t)*F.A(i,n,t)* F.h(i,n,F.ID)/F.noise;  
                        end
                        %                         if F.sinr(i,n,t)>= F.tsinr(i)
                        F.TR(i,t)=F.TR(i,t)+ (log2( 1 + F.sinr(i,n,t)));
                        F.Cell_TR(t)= F.Cell_TR(t)+log2( 1 + F.sinr(i,n,t));
                        F.SP_Eff(i,t)= F.SP_Eff(i,t) + (log2( F.sf)/F.N);
                        %                        end
                        F.sum_sinr(i,t) = F.sum_sinr(i,t) + F.sinr(i,n,t);
                        c=c+1;
                    end
                end
                F.nASC(i)=c;
                F.avg_sinr(i,t)= F.sum_sinr(i,t)/F.nASC(i);
            end
            
            
            
        end
        
        function F=Update(F)
            F=Interference(F,F.It) ;
            F=SINR(F,F.It);
            %             F.It=F.It + 1;
        end
        
        function F=MaxInterfereceNotify(F,Miyu)
            F.Miyu= Miyu;
        end
        function F=SC_MaxPower(F)
            t= F.It;
            for i=1 : F.M
                for n=1:F.N
                    if F.A(i,n,t)==1
                        F.Pmax_SC(n)=F.I_th(n)/F.h(i,n,1);
                    end
                end
            end
            
        end
        
        function F=UpdatePower(F)
            t=F.It ;
%             F=Interference(F,t);
            
            for l=1 : length(F.DS_U)
                i=F.DS_U(l);
                for n=1:F.N
                    F.P(i,n,t)= (1/log(2)) * ((1+F.Nou(i,t))/(F.Lambda(i,t)+F.Miyu(n)* F.h(i,n,1))) - (F.I(i,n,t-1)/F.h(i,n,F.ID));
                    if F.P(i,n,t)<0
                        F.P(i,n,t)=0;
                    end
                    %                     if F.P(i,n,t)>F.Pmax(i)
                    %                         F.P(i,n,t)=F.Pmax(i);
                    %                     end
                end
            end
            for l=1 : length(F.DT_U)
                i=F.DT_U(l);
                for n=1:F.N
                    F.P(i,n,t)= (1/log(2)) * ( 1 /(F.Lambda(i,t)+F.Miyu(n)* F.h(i,n,1))) - (F.I(i,n,t-1)/F.h(i,n,F.ID));
                    if F.P(i,n,t)<0
                        F.P(i,n,t)=0;
                    end
                    %                     if F.P(i,n,t)>F.Pmax(i)
                    %                         F.P(i,n,t)=F.Pmax(i);
                    %                     end
                    %
                end
            end
            
            F= ScAlloc (F,t);
            
            for i=1 : F.M
                for n=1:F.N
                    F.P(i,n,t)= F.P(i,n,t)* F.A(i,n,t);
                end
            end
        end
        
        function F=Init(F)
            t= F.It;
            
            % Assignn Equal power to all subchannells
            for i=1 : F.M
                for n=1:F.N
                    F.P(i,n,t)=F.Pmax(i)/F.N;
                end
            end
            
            for i=1 : F.M
                for j=1 : F.N
                    F.A(i,j,t)=0;
                end
            end
            %              F=Interference(F,t);
            F=SINR(F,t);
            
            for n=1 : F.N
                Sc(n,1)=n;
                Sc(n,2)=0;   % Flag that show a SC Used or not Used
            end
            
            for l=1 : length(F.DS_U)
                DS(l,1)=F.DS_U(l);
                DS(l,2)=0;   %Set a flag for each user
            end
            while(1==1)
                flag=0;
                %find a DS User
                for l=1 : length(F.DS_U)
                    if DS(l,2)==0
                        L=l;
                        i=DS(l,1);
                        flag=1;
                        break;
                    end
                end
                if flag==0
                    break;
                end
                
                %find a unused subchannel
                flag2=0;
                for n=1 : F.N
                    if Sc(n,2)==0
                        best_SC=n;
                        flag2=1;
                        break;
                    end
                end
                if flag2==0
                    break;  % can not found unused SC
                end
                
                for n=1 : F.N
                    if Sc(n,2)==0  %if SC in not used
                        if (F.h(i,n,F.ID))>(F.h(i,best_SC,F.ID))
                            best_SC=n;
                        end
                    end
                end
                F.A(i,best_SC,t)=1;
                Sc(best_SC,2)=1;
                
                %                      F=Interference(F,t);
                F=SINR(F,t);
                if F.TR(i,t)>= F.Ru(i)
                    DS(L,2)=1;
                end
                
            end
            
            while(1==1)
                flag=0;
                %find a unused subchannel
                for n=1 : F.N
                    if Sc(n,2)==0
                        best_SC=n;
                        flag=1;
                        break;
                    end
                end
                if flag==0
                    break;
                end
                best_user=1;
                for i=2 : F.M
                    if (F.h(i,best_SC,F.ID))>(F.h(best_user,best_SC,F.ID))
                        best_user=i;
                    end
                end
                F.A(best_user,best_SC,t)=1;
                Sc(best_SC,2)=1;
            end
            
            F=SINR(F,t);
            
            
            % count number of assigned Subchannel to each user
            for i=1 : F.M
                F.nASC(i)=0;
                for n=1 : F.N
                    if F.A(i,n,t)==1
                        F.nASC(i)=F.nASC(i)+1;
                    end
                end
            end
        end
        
        
        function F=ScAlloc(F,t)  % Algorithm2 Zhang
            
            % assign eaxh Subchannel to best user on it
            for i=1 : F.M
                for j=1 : F.N
                    F.A(i,j,t)=0;
                end
            end
            %              F=Interference(F,t);
            F=SINR(F,t);
            
            if t>1
                for n=1 : F.N
                    Sc(n,1)=n;
                    Sc(n,2)=0;   % Flag that show a SC Used or not Used
                end
                
                for l=1 : length(F.DS_U)
                    DS(l,1)=F.DS_U(l);
                    DS(l,2)=0;   %Set a flag for each user
                end
                while(1==1)
                    flag=0;
                    %find a DS User
                    for l=1 : length(F.DS_U)
                        if DS(l,2)==0
                            L=l;
                            i=DS(l,1);
                            flag=1;
                            break;
                        end
                    end
                    if flag==0
                        break;
                    end
                    
                    %find a unused subchannel
                    flag2=0;
                    for n=1 : F.N
                        if Sc(n,2)==0
                            best_SC=n;
                            flag2=1;
                            break;
                        end
                    end
                    if flag2==0
                        break;  % can not found unused SC
                    end
                    
                    for n=1 : F.N
                        if Sc(n,2)==0  %if SC in not used
                            if (F.h(i,n,F.ID)/F.I(i,n,t-1))>(F.h(i,best_SC,F.ID)/F.I(i,best_SC,t-1))
                                best_SC=n;
                            end
                        end
                    end
                    F.A(i,best_SC,t)=1;
                    Sc(best_SC,2)=1;
                    
                    %                      F=Interference(F,t);
                    F=SINR(F,t);
                    if F.TR(i,t)>= F.Ru(i)
                        DS(L,2)=1;
                    end
                    
                end
                
                while(1==1)
                    flag=0;
                    %find a unused subchannel
                    for n=1 : F.N
                        if Sc(n,2)==0
                            best_SC=n;
                            flag=1;
                            break;
                        end
                    end
                    if flag==0
                        break;
                    end
                    best_user=1;
                    for i=2 : F.M
                        if (F.h(i,best_SC,F.ID)/F.I(i,best_SC,t-1))>(F.h(best_user,best_SC,F.ID)/F.I(best_user,best_SC,t-1))
                            best_user=i;
                        end
                    end
                    F.A(best_user,best_SC,t)=1;
                    Sc(best_SC,2)=1;
                    %                      F=Interference(F,t);
                    
                end
                
                F=SINR(F,t);
                
            end
            % count number of assigned Subchannel to each user
            for i=1 : F.M
                F.nASC(i)=0;
                for n=1 : F.N
                    if F.A(i,n,t)==1
                        F.nASC(i)=F.nASC(i)+1;
                    end
                end
            end
        end
        
        function calc(F)
            k=0;  % index of non supported users
            l=0;    % index of supported users
            for i=1 : F.M
                if ((F.tsinr(i)-F.sinr(i,n,t-1 ))< 0.00001 )
                    l=l+1;
                    sp(l)=i;
                else
                    k=k+1;
                    nsp(k)=i;
                end
            end
            disp(F.T(t-1))
            disp(k)
            %disp(nsp)
            disp(l)
            %disp(sp)
            outage= k/(l+k)
            
            %figure
            
        end
        
        function disp(F)
            t= F.It;
            for n= 1: F.N
                plot(1:1:t-1, F.sinr(1,:), 1:1:t-1,F.sinr(2,:), 1:1:t-1,F.sinr(3,:),1:1:t-1,F.sinr(4,:),1:1:t-1,F.sinr(5,:),1:t-1,F.sinr(6,:));
                xlabel('Time(Iteration)')
                ylabel('SINR')
                legend('user 1','user 2','user 3','user 4','user 5','user 6')
                grid
                
                figure
                plot(1:1:t-1, F.P(1,:), 1:1:t-1,F.P(2,:), 1:1:t-1,F.P(3,:),1:1:t-1,F.P(4,:),1:1:t-1,F.P(5,:),1:t-1,F.P(6,:));
                xlabel('Time(Iteration)')
                ylabel('Power')
                legend('user 1','user 2','user 3','user 4','user 5','user 6')
                grid
            end
        end
        
        function F=Reset_It(F)
            F.It=1;
        end
        function F=Set_conf(F,Pmx,Tsnr)  % Set Maximum power and Target Sinr and other Configurationd for users
            for i=1 :F.M
                F.tsinr(i)=Tsnr;
            end
            if Tsnr == 9.55
                F.sf= 4;
            end
            
            if Tsnr == 45.11
                F.sf= 16;
            end
            
            if Tsnr == 179.85
                F.sf= 64;
            end
            
            if Tsnr == 694.17
                F.sf= 256;
            end
            if Tsnr == 2667.32
                F.sf= 1024;
            end
            
            
            for i=1 :F.M
                F.Pmax(i)= Pmx;
            end
        end
        function F=Update_Parameter(F)
            %Update Lambda for each user
            t=F.It;
            for i=1 : F.M
                s(i)=0;
                for n=1 : F.N
                    s(i) = s(i) +( F.A(i,n,t) * F.P(i,n,t));
                end
                F.Lambda(i,t+1)= F.Lambda(i,t)- (F.Beta(1,t)*(F.Pmax(i)- s(i))) ;
                if F.Lambda(i,t+1)<0
                    F.Lambda(i,t+1)=0;
                end
            end
            
            %Update Nou for each Delay bsensitive user
            for l=1 : length(F.DS_U)
                i=F.DS_U(l);
                F.Nou(i,t+1)= F.Nou(i,t)-(F.Beta(2,t)*(F.TR(i,t)- F.Ru(i)) );
                if F.Nou(i,t+1)<0
                    F.Nou(i,t+1)=0;
                end
            end
            for l=1 : length(F.DT_U)
                i=F.DT_U(l);
                F.Nou(i,t+1)= 0;
            end
            F.Beta(1,t+1) = 10/ sqrt(t+1);
            F.Beta(2,t+1) =10/ sqrt(t+1);
        end
        
        
        function F=It_Increment(F)
            F.It= F.It +1;
        end
    end
end

