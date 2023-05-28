clear;
K=5;   % Totall repeat of simulation
nFemto=8;
nMusr=6;
N=18;
X=1000;
Y=1000;
nFusr=4;
T=10;

for t=1 : T
    for j=1 : nFemto+1
        avg_cell_TR(j,t)=0;
    end
    for j=1 : nMusr
        mue_avg_rate(j,t)=0;
    end
end


for k=1 : K
    macrocell(k) = Macro(nMusr,N,X,Y,nFemto,nFusr);
    macrocell(k)=PathGain(macrocell(k));
    macrocell(k)=Update(macrocell(k),1);
    avg_cell_TR(1,1)=avg_cell_TR(1,1)+ macrocell(k).Cell_TR(1);
    for j=1 : nMusr
        mue_avg_rate(j,1)=mue_avg_rate(j,1)+ macrocell(k).TR(j,1);
    end
    
    for j=2 : nFemto+1 
        avg_cell_TR(j,1)=avg_cell_TR(j,1)+ macrocell(k).Fcel(j-1).Cell_TR(1);
    end
    for t=2 : T
        macrocell(k)=Update(macrocell(k),0);
        avg_cell_TR(1,t)=avg_cell_TR(1,t)+ macrocell(k).Cell_TR(t);
        for j=2 : nFemto+1
            avg_cell_TR(j,t)=avg_cell_TR(j,t)+ macrocell(k).Fcel(j-1).Cell_TR(t);
        end
        for j=1 : nMusr
           mue_avg_rate(j,t)=mue_avg_rate(j,t)+ macrocell(k).TR(j,t);
        end        
    end
    
     
end

for t=1 : T
    femto_av(t)=0;
    for j=1 : nFemto+1
        avg_cell_TR(j,t)=avg_cell_TR(j,t)/K;
        if j>1
           femto_av(t)=femto_av(t)+ avg_cell_TR(j,t);
        end
    end
    femto_av(t)= femto_av(t)/nFemto;
    for j=1 : nMusr
           mue_avg_rate(j,t)=mue_avg_rate(j,t)/K;
    end        
end
linespec = {'-+b','-sb','-ob','-*b','-.b','-db'};
plot(1:1:t, avg_cell_TR(1,:), 1:1:t,avg_cell_TR(2,:), 1:1:t,avg_cell_TR(3,:),1:1:t,avg_cell_TR(4,:),1:1:t,avg_cell_TR(5,:),1:1:t,avg_cell_TR(6,:),1:1:t,avg_cell_TR(7,:));%,1:1:t);,avg_cell_TR(8,:),1:1:t,avg_cell_TR(9,:));
xlabel('Time(Iteration)')
ylabel('Throughput(bps/HZ)')
legend('Macrocell','Femtocell 1','Femtocell 2','Femtocell 3','Femtocell 4','Femtocell 5','Femtocell 6','Femtocell 7','Femtocell 8')
grid 
figure

plot(1:1:t, avg_cell_TR(1,:),linespec{6}, 1:1:t,femto_av(:),linespec{2});
xlabel('Time(Iteration)')
ylabel('Average Throughput(bps/HZ)')
legend('Macrocell','Femtocells');
grid 
figure

plot(1:1:t, mue_avg_rate(1,:),linespec{1}, 1:1:t,mue_avg_rate(2,:),linespec{2}, 1:1:t,mue_avg_rate(3,:),linespec{3},1:1:t,mue_avg_rate(4,:),linespec{4},1:1:t,mue_avg_rate(5,:),linespec{5},1:1:t,mue_avg_rate(6,:),linespec{6});
xlabel('Time(Iteration)')
ylabel('Throughput(bps/HZ)')
legend('MUE l','MUE 2','MUE 3','MUE 4','MUE 5','MUE 6')
grid 
