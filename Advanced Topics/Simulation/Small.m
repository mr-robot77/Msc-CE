classdef Small < handle
    properties
        M;       % number of mobile nodes in Small-cell
        N;       % number of Subchannels in Small-cell
        X;       % cell width
        Y;       % cell width
        mnodes;  % mobile nodes in Small-cell
        bs;
        noise;
        Ai;      % Path loss Parameter
        Bi;      % Path loss Parameter
        Ci;      % Path loss Parameter
        fc;      % fequency
        WL;      % Path loss Parameter
        nWL;     % Path loss Parameter
        tsinr;  
        ss;           % femtocell constellation size
        nSC;          % number of subchannel that each user needs
        nASC;         % number of assigned subchannel to each user
        A;       % Subchannel assignmet matrix
        P;       % power assignmet matrix
        h;       % Path gain matrix
        I;       % Interference matrix
        sinr;    % Sinr matrix
        TR;      %Throughput
        Cell_TR;    % Cell Throughput
        SP_Eff;
        T;       % Number of assigned SC
        It;      % Iteration counter
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
        function S=Small(m,n,x,y,Mac,id)
            S.M=m;
            S.N=n;
            S.X=x;
            S.Y=y;
            S.noise=10^(-13);
            S.Ai=25;
            S.Bi=45;
            S.Ci=20;
            S.fc=2 * (10^9);
            S.WL=5;
            S.nWL=1;
            S.It=1;
            for i=1 :S.M
                S.tsinr(i)=179.85;
            end
            S.ss= 64;
            for i=1 :S.M
                S.Pmax(i)= 0.05;
            end
            
            
            S.T=floor(S.N/S.M);
            S.Mcel=Mac;
            S.ID=id;
            S=Create_BS_UE(S);
            
%             for i=1 : (F.M/2)
%                 F.DS_U(i)=i;
%             end
%             
            for i=1 : (S.M)
                S.DT_U(i)= i; % + (F.M/2);
            end
            for i=1 :S.M
                S.Lambda(i,1)= 0.1;
                
            end
            for l=1 : length(S.DS_U)
                i=S.DS_U(l);
                S.Nou(i,1)= 0.1;
            end
            for l=1 : length(S.DT_U)
                i=S.DT_U(l);
                S.Nou(i,1)= 0;
            end
            for i=1 :S.M
                S.Ru(i)= 9;
            end
            S.Beta(1,1)=10;
            S.Beta(2,1)=10;
            
        end
        
        function S=Create_BS_UE(S)
            bsx=unifrnd(1, S.X);
            bsy=unifrnd(1, S.Y);
            rx = unifrnd(1:S.M, 20);
            ry = unifrnd(1:S.M, 20);
            for j=1 : S.M
                rx(j)= rx(j)+ bsx;
                ry(j)= ry(j)+ bsy;
            end
            for j=1 : S.M
                S.mnodes(j)= complex(rx(j), ry(j));
            end
            S.bs = complex(bsx+10,bsy+10);
            
        end
        function S=PathGain(S)
            %Path gain to Macro BS
            for i = 1 : S.M
                u = rand(S.N, 1); % generating uniform variates
                sigma = 1; % the parameter
                x = sigma * sqrt(-2 * log(u)); % generating Rayleigh-distributed variates
                for n = 1 : S.N
                    %                     db=((F.Mcel.Ai * log10(abs(F.mnodes(i)-F.Mcel.bs ))) + F.Mcel.Bi + (F.Mcel.Ci * log10(F.Mcel.fc/5))+ (F.Mcel.nWL*F.Mcel.WL));
                    %                     F.h(i,n,1)=10^(-db/10);
                    S.h(i,n,1) = x(n) * (abs(S.mnodes(i)-S.Mcel.bs )^(-3));
                end
            end
            %Path gain to Femto BSs
            for i = 1 : S.M
                for k= 1 : S.Mcel.nScel
                    u = rand(S.N, 1); % generating uniform variates
                    sigma = 1; % the parameter
                    x = sigma * sqrt(-2 * log(u)); % generating Rayleigh-distributed variates
                    for n = 1 : S.N
                        % db=((F.Mcel.Fcel(k).Ai * log10(abs(F.mnodes(i)-F.Mcel.Fcel(k).bs ))) + F.Mcel.Fcel(k).Bi + (F.Mcel.Fcel(k).Ci * log10(F.Mcel.Fcel(k).fc/5))+ (F.Mcel.Fcel(k).nWL*F.Mcel.Fcel(k).WL));
                        %  F.h(i,n,k+1)= 10^(-db/10);
                        S.h(i,n,k+1) = x(n) * (abs(S.mnodes(i)-S.Mcel.Scel(k).bs )^(-3));
                    end
                end
            end
        end
        
        function S=Interference(S,t)
            %t= Mac.It +1;
            for i=1 : S.M
                for n= 1: S.N
                    S.I(i,n,t)=S.noise;
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
                    for j=1 : S.Mcel.M
                        S.I(i,n,t)=S.I(i,n,t)+ S.Mcel.A(j,n,t)* S.Mcel.P(j,n,t)* S.Mcel.h(j,n,S.ID);
                    end
                    
                end
            end
        end
        
        function S=SINR(S,t)
            S.Cell_TR(t)=0;
            for i=1 : S.M
                S.TR(i,t)=0;
                S.SP_Eff(i,t)=0;
                c=0;
                S.sum_sinr(i,t)=0;
                for n=1:S.N
                    if S.A(i,n,t)==1
                        if t> 1
                          S.sinr(i,n,t)= S.P(i,n,t)*S.A(i,n,t)* S.h(i,n,S.ID)/S.I(i,n,t-1);
                        else
                          S.sinr(i,n,t)= S.P(i,n,t)*S.A(i,n,t)* S.h(i,n,S.ID)/S.noise;  
                        end
                        %                         if F.sinr(i,n,t)>= F.tsinr(i)
                        S.TR(i,t)=S.TR(i,t)+ (log2( 1 + S.sinr(i,n,t)));
                        S.Cell_TR(t)= S.Cell_TR(t)+log2( 1 + S.sinr(i,n,t));
                        S.SP_Eff(i,t)= S.SP_Eff(i,t) + (log2( S.ss)/S.N);
                        %                        end
                        S.sum_sinr(i,t) = S.sum_sinr(i,t) + S.sinr(i,n,t);
                        c=c+1;
                    end
                end
                S.nASC(i)=c;
                S.avg_sinr(i,t)= S.sum_sinr(i,t)/S.nASC(i);
            end
            
            
            
        end
        
        function S=Update(S)
            S=Interference(S,S.It) ;
            S=SINR(S,S.It);
            %             F.It=F.It + 1;
        end
        
        function S=MaxInterfereceNotify(S,Miyu)
            S.Miyu= Miyu;
        end
        function S=SC_MaxPower(S)
            t= S.It;
            for i=1 : S.M
                for n=1:S.N
                    if S.A(i,n,t)==1
                        S.Pmax_SC(n)=S.I_th(n)/S.h(i,n,1);
                    end
                end
            end
            
        end
        
        function S=UpdatePower(S)
            t=S.It ;
%             F=Interference(F,t);
            
            for l=1 : length(S.DS_U)
                i=S.DS_U(l);
                for n=1:S.N
                    S.P(i,n,t)= (1/log(2)) * ((1+S.Nou(i,t))/(S.Lambda(i,t)+S.Miyu(n)* S.h(i,n,1))) - (S.I(i,n,t-1)/S.h(i,n,S.ID));
                    if S.P(i,n,t)<0
                        S.P(i,n,t)=0;
                    end
                    %                     if F.P(i,n,t)>F.Pmax(i)
                    %                         F.P(i,n,t)=F.Pmax(i);
                    %                     end
                end
            end
            for l=1 : length(S.DT_U)
                i=S.DT_U(l);
                for n=1:S.N
                    S.P(i,n,t)= (1/log(2)) * ( 1 /(S.Lambda(i,t)+S.Miyu(n)* S.h(i,n,1))) - (S.I(i,n,t-1)/S.h(i,n,S.ID));
                    if S.P(i,n,t)<0
                        S.P(i,n,t)=0;
                    end
                    %                     if F.P(i,n,t)>F.Pmax(i)
                    %                         F.P(i,n,t)=F.Pmax(i);
                    %                     end
                    %
                end
            end
            
            S= ScAlloc (S,t);
            
            for i=1 : S.M
                for n=1:S.N
                    S.P(i,n,t)= S.P(i,n,t)* S.A(i,n,t);
                end
            end
        end
        
        function S=Init(S)
            t= S.It;
            
            % Assignn Equal power to all subchannells
            for i=1 : S.M
                for n=1:S.N
                    S.P(i,n,t)=S.Pmax(i)/S.N;
                end
            end
            
            for i=1 : S.M
                for j=1 : S.N
                    S.A(i,j,t)=0;
                end
            end
            %              F=Interference(F,t);
            S=SINR(S,t);
            
            for n=1 : S.N
                Sc(n,1)=n;
                Sc(n,2)=0;   % Flag that show a SC Used or not Used
            end
            
            for l=1 : length(S.DS_U)
                DS(l,1)=S.DS_U(l);
                DS(l,2)=0;   %Set a flag for each user
            end
            while(1==1)
                flag=0;
                %find a DS User
                for l=1 : length(S.DS_U)
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
                for n=1 : S.N
                    if Sc(n,2)==0
                        best_SC=n;
                        flag2=1;
                        break;
                    end
                end
                if flag2==0
                    break;  % can not found unused SC
                end
                
                for n=1 : S.N
                    if Sc(n,2)==0  %if SC in not used
                        if (S.h(i,n,S.ID))>(S.h(i,best_SC,S.ID))
                            best_SC=n;
                        end
                    end
                end
                S.A(i,best_SC,t)=1;
                Sc(best_SC,2)=1;
                
                %                      F=Interference(F,t);
                S=SINR(S,t);
                if S.TR(i,t)>= S.Ru(i)
                    DS(L,2)=1;
                end
                
            end
            
            while(1==1)
                flag=0;
                %find a unused subchannel
                for n=1 : S.N
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
                for i=2 : S.M
                    if (S.h(i,best_SC,S.ID))>(S.h(best_user,best_SC,S.ID))
                        best_user=i;
                    end
                end
                S.A(best_user,best_SC,t)=1;
                Sc(best_SC,2)=1;
            end
            
            S=SINR(S,t);
            
            
            % count number of assigned Subchannel to each user
            for i=1 : S.M
                S.nASC(i)=0;
                for n=1 : S.N
                    if S.A(i,n,t)==1
                        S.nASC(i)=S.nASC(i)+1;
                    end
                end
            end
        end
        
        
        function S=ScAlloc(S,t)  % Algorithm2 Zhang
            
            % assign eaxh Subchannel to best user on it
            for i=1 : S.M
                for j=1 : S.N
                    S.A(i,j,t)=0;
                end
            end
            %              F=Interference(F,t);
            S=SINR(S,t);
            
            if t>1
                for n=1 : S.N
                    Sc(n,1)=n;
                    Sc(n,2)=0;   % Flag that show a SC Used or not Used
                end
                
                for l=1 : length(S.DS_U)
                    DS(l,1)=S.DS_U(l);
                    DS(l,2)=0;   %Set a flag for each user
                end
                while(1==1)
                    flag=0;
                    %find a DS User
                    for l=1 : length(S.DS_U)
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
                    for n=1 : S.N
                        if Sc(n,2)==0
                            best_SC=n;
                            flag2=1;
                            break;
                        end
                    end
                    if flag2==0
                        break;  % can not found unused SC
                    end
                    
                    for n=1 : S.N
                        if Sc(n,2)==0  %if SC in not used
                            if (S.h(i,n,S.ID)/S.I(i,n,t-1))>(S.h(i,best_SC,S.ID)/S.I(i,best_SC,t-1))
                                best_SC=n;
                            end
                        end
                    end
                    S.A(i,best_SC,t)=1;
                    Sc(best_SC,2)=1;
                    
                    %                      F=Interference(F,t);
                    S=SINR(S,t);
                    if S.TR(i,t)>= S.Ru(i)
                        DS(L,2)=1;
                    end
                    
                end
                
                while(1==1)
                    flag=0;
                    %find a unused subchannel
                    for n=1 : S.N
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
                    for i=2 : S.M
                        if (S.h(i,best_SC,S.ID)/S.I(i,best_SC,t-1))>(S.h(best_user,best_SC,S.ID)/S.I(best_user,best_SC,t-1))
                            best_user=i;
                        end
                    end
                    S.A(best_user,best_SC,t)=1;
                    Sc(best_SC,2)=1;
                    %                      F=Interference(F,t);
                    
                end
                
                S=SINR(S,t);
                
            end
            % count number of assigned Subchannel to each user
            for i=1 : S.M
                S.nASC(i)=0;
                for n=1 : S.N
                    if S.A(i,n,t)==1
                        S.nASC(i)=S.nASC(i)+1;
                    end
                end
            end
        end
        
        function calc(S)
            k=0;  % index of non supported users
            l=0;    % index of supported users
            for i=1 : S.M
                if ((S.tsinr(i)-S.sinr(i,n,t-1 ))< 0.00001 )
                    l=l+1;
                    sp(l)=i;
                else
                    k=k+1;
                    nsp(k)=i;
                end
            end
            disp(S.T(t-1))
            disp(k)
            %disp(nsp)
            disp(l)
            %disp(sp)
            outage= k/(l+k)
            
            %figure
            
        end
        
        function disp(S)
            t= S.It;
            for n= 1: S.N
                plot(1:1:t-1, S.sinr(1,:), 1:1:t-1,S.sinr(2,:), 1:1:t-1,S.sinr(3,:),1:1:t-1,S.sinr(4,:),1:1:t-1,S.sinr(5,:),1:t-1,S.sinr(6,:));
                xlabel('Time(Iteration)')
                ylabel('SINR')
                legend('user 1','user 2','user 3','user 4','user 5','user 6')
                grid
                
                figure
                plot(1:1:t-1, S.P(1,:), 1:1:t-1,S.P(2,:), 1:1:t-1,S.P(3,:),1:1:t-1,S.P(4,:),1:1:t-1,S.P(5,:),1:t-1,S.P(6,:));
                xlabel('Time(Iteration)')
                ylabel('Power')
                legend('user 1','user 2','user 3','user 4','user 5','user 6')
                grid
            end
        end
        
        function S=Reset_It(S)
            S.It=1;
        end
        function S=Set_conf(S,Pmx,Tsnr)  % Set Maximum power and Target Sinr and other Configurationd for users
            for i=1 :S.M
                S.tsinr(i)=Tsnr;
            end
            if Tsnr == 9.55
                S.ss= 4;
            end
            
            if Tsnr == 45.11
                S.ss= 16;
            end
            
            if Tsnr == 179.85
                S.ss= 64;
            end
            
            if Tsnr == 694.17
                S.ss= 256;
            end
            if Tsnr == 2667.32
                S.ss= 1024;
            end
            
            
            for i=1 :S.M
                S.Pmax(i)= Pmx;
            end
        end
        function S=Update_Parameter(S)
            %Update Lambda for each user
            t=S.It;
            for i=1 : S.M
                s(i)=0;
                for n=1 : S.N
                    s(i) = s(i) +( S.A(i,n,t) * S.P(i,n,t));
                end
                S.Lambda(i,t+1)= S.Lambda(i,t)- (S.Beta(1,t)*(S.Pmax(i)- s(i))) ;
                if S.Lambda(i,t+1)<0
                    S.Lambda(i,t+1)=0;
                end
            end
            
            %Update Nou for each Delay bsensitive user
            for l=1 : length(S.DS_U)
                i=S.DS_U(l);
                S.Nou(i,t+1)= S.Nou(i,t)-(S.Beta(2,t)*(S.TR(i,t)- S.Ru(i)) );
                if S.Nou(i,t+1)<0
                    S.Nou(i,t+1)=0;
                end
            end
            for l=1 : length(S.DT_U)
                i=S.DT_U(l);
                S.Nou(i,t+1)= 0;
            end
            S.Beta(1,t+1) = 10/ sqrt(t+1);
            S.Beta(2,t+1) =10/ sqrt(t+1);
        end
        
        
        function S=It_Increment(S)
            S.It= S.It +1;
        end
    end
end

