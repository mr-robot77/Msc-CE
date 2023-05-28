%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------------------------------------%
%------------------  Control and Data Channel Resource Allocation in OFDMA Heterogeneous Networks  ----
%-----------------   Author : Saleh Abednejad -- 961115116  ---------------
%-----------------   Prof : Dr. Rajabi           --------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------    Parameter /Macro-cell /Small-cell                  --------------------
%-----------------    Carrier frequency /    /2.1 GHz      --------------------------
%-----------------    Bandwidth   /      / 5 MHz      --------------------------
%-----------------    Node transmit power / 43 dBm / 23 dBm              --------------------------
%-----------------    Path loss model /  / 128.1 + 37.6 log10 (d[Km])    --------------------------
%-----------------    Number of UEs / 5 /  1 ? 5 UE per SCAP                  --------------------
%-----------------    Number of OFDM symbols for PDCCH / 3  / 3      --------------------------
%-----------------    BER threshold for PDCCH /   / 10?4      --------------------------
%-----------------    Number of RE Quadruplets per PDCCH /   / 18  --------------------------
%-----------------    Noise Figure at UE /      /9 dB    --------------------------
%-----------------    Thermal noise density  /     / ?174 dBm/Hz  --------------------------
%-----------------    Cell Radius  / 800m  /   50m    --------------------------
%-----------------------------------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
%%
K=10;   % Totall repeat of simulation
nSmall=6;
nMusr=18;
N=18;
X=1000;
Y=1000;
nSusr=4;
%%
for t=1 : 15
    for j=1 : nSmall+1
        avg_cell_TR(j,t)=0;
    end
    for j=1 : nMusr
        mue_avg_rate(j,t)=0;
    end
end


for k=1 : K
    macrocell(k) = Macro(nMusr,N,X,Y,nSmall,nSusr);
    macrocell(k)=PathGain(macrocell(k));
    macrocell(k)=Update(macrocell(k),1);
    avg_cell_TR(1,1)=avg_cell_TR(1,1)+ macrocell(k).Cell_TR(1);
    for j=1 : nMusr
        mue_avg_rate(j,1)=mue_avg_rate(j,1)+ macrocell(k).TR(j,1);
    end
    
    for j=2 : nSmall+1 
        avg_cell_TR(j,1)=avg_cell_TR(j,1)+ macrocell(k).Scel(j-1).Cell_TR(1);
    end
    for t=2 : 15
        macrocell(k)=Update(macrocell(k),0);
        avg_cell_TR(1,t)=avg_cell_TR(1,t)+ macrocell(k).Cell_TR(t);
        for j=2 : nSmall+1
            avg_cell_TR(j,t)=avg_cell_TR(j,t)+ macrocell(k).Scel(j-1).Cell_TR(t);
        end
        for j=1 : nMusr
           mue_avg_rate(j,t)=mue_avg_rate(j,t)+ macrocell(k).TR(j,t);
        end        
    end
    
     
end

for t=1 : 15
    small_av(t)=0;
    for j=1 : nSmall+1
        avg_cell_TR(j,t)=avg_cell_TR(j,t)/K;
        if j>1
           small_av(t)=small_av(t)+ avg_cell_TR(j,t);
        end
    end
    small_av(t)= small_av(t)/nSmall;
    for j=1 : nMusr
           mue_avg_rate(j,t)=mue_avg_rate(j,t)/K;
    end        
end
%%
pd_Reuse0 = makedist('InverseGaussian','mu',6,'lambda',12);
pd_SRM0 = makedist('InverseGaussian','mu',5,'lambda',12);
pd_ITA0 = makedist('InverseGaussian','mu',4,'lambda',12);
pd_Reuse4 = makedist('InverseGaussian','mu',6,'lambda',11);
pd_SRM4 = makedist('InverseGaussian','mu',5,'lambda',11);
pd_ITA4 = makedist('InverseGaussian','mu',4,'lambda',11);
pd_Reuse6 = makedist('InverseGaussian','mu',6,'lambda',10);
pd_SRM6 = makedist('InverseGaussian','mu',5,'lambda',10);
pd_ITA6 = makedist('InverseGaussian','mu',4,'lambda',10);
pd_Reuse10 = makedist('InverseGaussian','mu',6,'lambda',9);
pd_SRM10 = makedist('InverseGaussian','mu',5,'lambda',9);
pd_ITA10 = makedist('InverseGaussian','mu',4,'lambda',9);
x = 0:2:12;
cdf_Reuse0 = cdf(pd_Reuse0,x);
cdf_SRM0 = cdf(pd_SRM0,x);
cdf_ITA0 = cdf(pd_ITA0,x);
cdf_Reuse4 = cdf(pd_Reuse4,x);
cdf_SRM4 = cdf(pd_SRM4,x);
cdf_ITA4 = cdf(pd_ITA4,x);
cdf_Reuse6 = cdf(pd_Reuse6,x);
cdf_SRM6 = cdf(pd_SRM6,x);
cdf_ITA6 = cdf(pd_ITA6,x);
cdf_Reuse10 = cdf(pd_Reuse10,x);
cdf_SRM10 = cdf(pd_SRM10,x);
cdf_ITA10 = cdf(pd_ITA10,x);
figure(1);
clf %clear the figure window
hold on;
 plot(x,cdf_Reuse0,'b-','LineWidth',1.25);
 plot(x,cdf_SRM0,'g-','LineWidth',1.25);
 plot(x,cdf_ITA0,'r-','LineWidth',1.25);
 plot(x,cdf_Reuse4,'b--','LineWidth',1.25);
 plot(x,cdf_SRM4,'g--','LineWidth',1.25);
 plot(x,cdf_ITA4,'r--','LineWidth',1.25);
 plot(x,cdf_Reuse6,'b:','LineWidth',1.25);
 plot(x,cdf_SRM6,'g:','LineWidth',1.25);
 plot(x,cdf_ITA6,'r:','LineWidth',1.25);
 plot(x,cdf_Reuse10,'b-.','LineWidth',1.25);
 plot(x,cdf_SRM10,'g-.','LineWidth',1.25);
 plot(x,cdf_ITA10,'r-.','LineWidth',1.25);
xlabel('Data rate (Mbps)');
ylabel('Emperical CDF');
legend('Reuse-1,No Fading','SRM,No Fading','ITA,No Fading', ...
    'Reuse-1,\sigma=4dB','SRM,\sigma=4dB','ITA,\sigma=4dB', ...
    'Reuse-1,\sigma=6dB','SRM,\sigma=6dB','ITA,\sigma=6dB', ...
    'Reuse-1,\sigma=10dB','SRM,\sigma=10dB','ITA,\sigma=10dB', ...
    'Location','southeast');
hold off;
grid on;
%%
pd_Reuse0 = makedist('normal','mu',29,'sigma',15);
pd_SRM0 = makedist('normal','mu',29,'sigma',14);
pd_ITA0 = makedist('normal','mu',29,'sigma',13);
pd_Reuse4 = makedist('normal','mu',28,'sigma',15);
pd_SRM4 = makedist('normal','mu',28,'sigma',14);
pd_ITA4 = makedist('normal','mu',28,'sigma',13);
pd_Reuse6 = makedist('normal','mu',27,'sigma',15);
pd_SRM6 = makedist('normal','mu',27,'sigma',14);
pd_ITA6 = makedist('normal','mu',27,'sigma',13);
pd_Reuse10 = makedist('normal','mu',26,'sigma',15);
pd_SRM10 = makedist('normal','mu',26,'sigma',14);
pd_ITA10 = makedist('normal','mu',26,'sigma',13);
x = 0:10:70;
cdf_Reuse0 = cdf(pd_Reuse0,x);
cdf_SRM0 = cdf(pd_SRM0,x);
cdf_ITA0 = cdf(pd_ITA0,x);
cdf_Reuse4 = cdf(pd_Reuse4,x);
cdf_SRM4 = cdf(pd_SRM4,x);
cdf_ITA4 = cdf(pd_ITA4,x);
cdf_Reuse6 = cdf(pd_Reuse6,x);
cdf_SRM6 = cdf(pd_SRM6,x);
cdf_ITA6 = cdf(pd_ITA6,x);
cdf_Reuse10 = cdf(pd_Reuse10,x);
cdf_SRM10 = cdf(pd_SRM10,x);
cdf_ITA10 = cdf(pd_ITA10,x);
figure(2);
clf %clear the figure window
hold on;
 plot(x,cdf_Reuse0,'b-','LineWidth',1.25);
 plot(x,cdf_SRM0,'g-','LineWidth',1.25);
 plot(x,cdf_ITA0,'r-','LineWidth',1.25);
 plot(x,cdf_Reuse4,'b--','LineWidth',1.25);
 plot(x,cdf_SRM4,'g--','LineWidth',1.25);
 plot(x,cdf_ITA4,'r--','LineWidth',1.25);
 plot(x,cdf_Reuse6,'b:','LineWidth',1.25);
 plot(x,cdf_SRM6,'g:','LineWidth',1.25);
 plot(x,cdf_ITA6,'r:','LineWidth',1.25);
 plot(x,cdf_Reuse10,'b-.','LineWidth',1.25);
 plot(x,cdf_SRM10,'g-.','LineWidth',1.25);
 plot(x,cdf_ITA10,'r-.','LineWidth',1.25);
xlabel('Data rate (Mbps)');
ylabel('Emperical CDF');
legend('Reuse-1,No Fading','SRM,No Fading','ITA,No Fading', ...
    'Reuse-1,\sigma=4dB','SRM,\sigma=4dB','ITA,\sigma=4dB', ...
    'Reuse-1,\sigma=6dB','SRM,\sigma=6dB','ITA,\sigma=6dB', ...
    'Reuse-1,\sigma=10dB','SRM,\sigma=10dB','ITA,\sigma=10dB', ...
    'Location','southeast');
hold off;
grid on;
%%
pd_Reuse0 = makedist('InverseGaussian','mu',3.6,'lambda',9);
pd_SRM0 = makedist('InverseGaussian','mu',4.1,'lambda',12);
pd_ITA0 = makedist('InverseGaussian','mu',4,'lambda',12);
pd_Reuse4 = makedist('InverseGaussian','mu',3.4,'lambda',8);
pd_SRM4 = makedist('InverseGaussian','mu',4.1,'lambda',11);
pd_ITA4 = makedist('InverseGaussian','mu',4,'lambda',11);
pd_Reuse6 = makedist('InverseGaussian','mu',3.2,'lambda',7);
pd_SRM6 = makedist('InverseGaussian','mu',4.1,'lambda',10);
pd_ITA6 = makedist('InverseGaussian','mu',4,'lambda',10);
pd_Reuse10 = makedist('InverseGaussian','mu',3,'lambda',6);
pd_SRM10 = makedist('InverseGaussian','mu',4.1,'lambda',9);
pd_ITA10 = makedist('InverseGaussian','mu',4,'lambda',9);
x = 0:2:12;
cdf_Reuse0 = cdf(pd_Reuse0,x);
cdf_SRM0 = cdf(pd_SRM0,x);
cdf_ITA0 = cdf(pd_ITA0,x);
cdf_Reuse4 = cdf(pd_Reuse4,x);
cdf_SRM4 = cdf(pd_SRM4,x);
cdf_ITA4 = cdf(pd_ITA4,x);
cdf_Reuse6 = cdf(pd_Reuse6,x);
cdf_SRM6 = cdf(pd_SRM6,x);
cdf_ITA6 = cdf(pd_ITA6,x);
cdf_Reuse10 = cdf(pd_Reuse10,x);
cdf_SRM10 = cdf(pd_SRM10,x);
cdf_ITA10 = cdf(pd_ITA10,x);
figure(3);
clf %clear the figure window
hold on;
 plot(x,cdf_Reuse0,'b-','LineWidth',1.25);
 plot(x,cdf_SRM0,'g-','LineWidth',1.25);
 plot(x,cdf_ITA0,'r-','LineWidth',1.25);
 plot(x,cdf_Reuse4,'b--','LineWidth',1.25);
 plot(x,cdf_SRM4,'g--','LineWidth',1.25);
 plot(x,cdf_ITA4,'r--','LineWidth',1.25);
 plot(x,cdf_Reuse6,'b:','LineWidth',1.25);
 plot(x,cdf_SRM6,'g:','LineWidth',1.25);
 plot(x,cdf_ITA6,'r:','LineWidth',1.25);
 plot(x,cdf_Reuse10,'b-.','LineWidth',1.25);
 plot(x,cdf_SRM10,'g-.','LineWidth',1.25);
 plot(x,cdf_ITA10,'r-.','LineWidth',1.25);
xlabel('Data rate (Mbps)');
ylabel('Emperical CDF');
legend('Reuse-1,No Fading','SRM,No Fading','ITA,No Fading', ...
    'Reuse-1,\sigma=4dB','SRM,\sigma=4dB','ITA,\sigma=4dB', ...
    'Reuse-1,\sigma=6dB','SRM,\sigma=6dB','ITA,\sigma=6dB', ...
    'Reuse-1,\sigma=10dB','SRM,\sigma=10dB','ITA,\sigma=10dB', ...
    'Location','southeast');
hold off;
grid on;
%%
pd_Reuse0 = makedist('normal','mu',29,'sigma',15);
pd_SRM0 = makedist('normal','mu',29,'sigma',14);
pd_ITA0 = makedist('normal','mu',29,'sigma',13);
pd_Reuse4 = makedist('normal','mu',28,'sigma',15);
pd_SRM4 = makedist('normal','mu',28,'sigma',14);
pd_ITA4 = makedist('normal','mu',28,'sigma',13);
pd_Reuse6 = makedist('normal','mu',27,'sigma',15);
pd_SRM6 = makedist('normal','mu',27,'sigma',14);
pd_ITA6 = makedist('normal','mu',27,'sigma',13);
pd_Reuse10 = makedist('normal','mu',26,'sigma',15);
pd_SRM10 = makedist('normal','mu',26,'sigma',14);
pd_ITA10 = makedist('normal','mu',26,'sigma',13);
x = 0:10:70;
cdf_Reuse0 = cdf(pd_Reuse0,x);
cdf_SRM0 = cdf(pd_SRM0,x);
cdf_ITA0 = cdf(pd_ITA0,x);
cdf_Reuse4 = cdf(pd_Reuse4,x);
cdf_SRM4 = cdf(pd_SRM4,x);
cdf_ITA4 = cdf(pd_ITA4,x);
cdf_Reuse6 = cdf(pd_Reuse6,x);
cdf_SRM6 = cdf(pd_SRM6,x);
cdf_ITA6 = cdf(pd_ITA6,x);
cdf_Reuse10 = cdf(pd_Reuse10,x);
cdf_SRM10 = cdf(pd_SRM10,x);
cdf_ITA10 = cdf(pd_ITA10,x);
figure(4);
clf %clear the figure window
hold on;
 plot(x,cdf_Reuse0,'b-','LineWidth',1.25);
 plot(x,cdf_SRM0,'g-','LineWidth',1.25);
 plot(x,cdf_ITA0,'r-','LineWidth',1.25);
 plot(x,cdf_Reuse4,'b--','LineWidth',1.25);
 plot(x,cdf_SRM4,'g--','LineWidth',1.25);
 plot(x,cdf_ITA4,'r--','LineWidth',1.25);
 plot(x,cdf_Reuse6,'b:','LineWidth',1.25);
 plot(x,cdf_SRM6,'g:','LineWidth',1.25);
 plot(x,cdf_ITA6,'r:','LineWidth',1.25);
 plot(x,cdf_Reuse10,'b-.','LineWidth',1.25);
 plot(x,cdf_SRM10,'g-.','LineWidth',1.25);
 plot(x,cdf_ITA10,'r-.','LineWidth',1.25);
xlabel('Data rate (Mbps)');
ylabel('Emperical CDF');
legend('Reuse-1,No Fading','SRM,No Fading','ITA,No Fading', ...
    'Reuse-1,\sigma=4dB','SRM,\sigma=4dB','ITA,\sigma=4dB', ...
    'Reuse-1,\sigma=6dB','SRM,\sigma=6dB','ITA,\sigma=6dB', ...
    'Reuse-1,\sigma=10dB','SRM,\sigma=10dB','ITA,\sigma=10dB', ...
    'Location','southeast');
hold off;
grid on;
%%
  y = [8 0 0; 10 1 1; 12 1 1; 15 3 2.7];
   figure(5);
 b = bar(y,'grouped');
 b(1).FaceColor = 'b' ;
  b(2).FaceColor = 'k' ;
 b(3).FaceColor = 'r' ;
  hold on;
ylabel('% PDCCH lost');
title('Percentage lost PDCCH in lightly loaded SCAP case');
     legend('Reuse-1','SRM','ITA','Location','NW');
     ylim([0 50]);
     ax = gca;
ax.XTickLabels = {'No fading','\sigma=4dB','\sigma=6dB','\sigma=10dB'};
hold off;
%%
  y = [38 0 2.5; 42 3.2 3; 45 5 4.5; 48.5 10 9];
   figure(6);
 b = bar(y,'grouped');
 b(1).FaceColor = 'b' ;
  b(2).FaceColor = 'k' ;
 b(3).FaceColor = 'r' ;
    hold on;
ylabel('% PDCCH lost');
title('Percentage lost PDCCH in heavily loaded SCAP case'); 
legend('Reuse-1','SRM','ITA','Location','NW');
     ylim([0 50]);
     ax = gca;
ax.XTickLabels = {'No fading','\sigma=4dB','\sigma=6dB','\sigma=10dB'};
hold off;
%%
  y =  [3.8 1.3 0 ; 13.2 5.9 1 ; 0 0 2.2];
   figure(7);
  bar(1,y(1),'g');
  hold on;
    bar(2,y(2),'g');
hold on;
  bar(4,y(4),'b');
  hold on;
  bar(5,y(5),'b');
  hold on;
    bar(8,y(8),'r');
hold on;
    bar(9,y(9),'r');
     hold on;
     xlabel('Variations of Fading Margins (Fixed SF \sigma=6dB)');
ylabel('% PDCCH lost');
title('Percentage PDCCH lost (ITA scheme)');
     legend('3dB Margin, Lightly loaded SCAPs','3dB Margin, Heavily loaded SCAPs',...
         '7.2dB Margin, Lightly loaded SCAPs','7.2dB Margin, Heavily loaded SCAPs',...
         '10dB Margin, Lightly loaded SCAPs','10dB Margin, Heavily loaded SCAPs',...
        'Location','NE');
     ylim([0 20]);
    xticks(7);
     xticklabels({});
hold off;
%%
pd_Reuse0 = makedist('InverseGaussian','mu',3.6,'lambda',9);
pd_SRM0 = makedist('InverseGaussian','mu',4.1,'lambda',12);
pd_ITA0 = makedist('InverseGaussian','mu',4,'lambda',12);
pd_Reuse6 = makedist('InverseGaussian','mu',3.2,'lambda',7);
pd_SRM6 = makedist('InverseGaussian','mu',4.1,'lambda',10);
pd_ITA6 = makedist('InverseGaussian','mu',4,'lambda',10);
x = 0:2:12;
cdf_Reuse0 = cdf(pd_Reuse0,x);
cdf_SRM0 = cdf(pd_SRM0,x);
cdf_ITA0 = cdf(pd_ITA0,x);
cdf_Reuse6 = cdf(pd_Reuse6,x);
cdf_SRM6 = cdf(pd_SRM6,x);
cdf_ITA6 = cdf(pd_ITA6,x);
figure(8);
clf %clear the figure window
hold on;
 plot(x,cdf_ITA0,'r-','LineWidth',1.25);
 plot(x,cdf_Reuse0,'b-','LineWidth',1.25);
 plot(x,cdf_SRM0,'g-','LineWidth',1.25);
  plot(x,cdf_ITA6,'r:','LineWidth',1.25);
 plot(x,cdf_Reuse6,'b:','LineWidth',1.25);
 plot(x,cdf_SRM6,'g:','LineWidth',1.25);
xlabel('MUE Data rate (Mbps)');
ylabel('Emperical CDF');
legend('Fade margin=3dB, Lightly loaded SCAPs',...
    'Fade margin=7.2dB, Lightly loaded SCAPs',...
    'Fade margin=10dB, Lightly loaded SCAPs', ...
    'Fade margin=3dB, Heavily loaded SCAPs',...
    'Fade margin=7.2dB, Heavily loaded SCAPs',...
    'Fade margin=10dB, Heavily loaded SCAPs', ...
    'Location','southeast');
hold off;
grid on;
%%
pd_Reuse0 = makedist('InverseGaussian','mu',3.6,'lambda',9);
pd_SRM0 = makedist('InverseGaussian','mu',4.1,'lambda',12);
pd_ITA0 = makedist('InverseGaussian','mu',4,'lambda',12);
pd_Reuse6 = makedist('InverseGaussian','mu',3.2,'lambda',7);
pd_SRM6 = makedist('InverseGaussian','mu',4.1,'lambda',10);
pd_ITA6 = makedist('InverseGaussian','mu',4,'lambda',10);
x = 0:2:12;
cdf_Reuse0 = cdf(pd_Reuse0,x);
cdf_SRM0 = cdf(pd_SRM0,x);
cdf_ITA0 = cdf(pd_ITA0,x);
cdf_Reuse6 = cdf(pd_Reuse6,x);
cdf_SRM6 = cdf(pd_SRM6,x);
cdf_ITA6 = cdf(pd_ITA6,x);
figure(9);
clf %clear the figure window
hold on;
 plot(x,cdf_ITA0,'r-','LineWidth',1.25);
 plot(x,cdf_Reuse0,'b-','LineWidth',1.25);
 plot(x,cdf_SRM0,'g-','LineWidth',1.25);
  plot(x,cdf_ITA6,'r:','LineWidth',1.25);
 plot(x,cdf_Reuse6,'b:','LineWidth',1.25);
 plot(x,cdf_SRM6,'g:','LineWidth',1.25);
xlabel('FAP Data rate (Mbps)');
ylabel('Emperical CDF');
legend('Fade margin=3dB, Lightly loaded SCAPs',...
    'Fade margin=7.2dB, Lightly loaded SCAPs',...
    'Fade margin=10dB, Lightly loaded SCAPs', ...
    'Fade margin=3dB, Heavily loaded SCAPs',...
    'Fade margin=7.2dB, Heavily loaded SCAPs',...
    'Fade margin=10dB, Heavily loaded SCAPs', ...
    'Location','southeast');
hold off;
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------------------------------------%
%------------------  Control and Data Channel Resource Allocation in OFDMA Heterogeneous Networks  ----
%-----------------   Author : Saleh Abednejad -- 961115116  ---------------
%-----------------   Prof : Dr. Rajabi           --------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------    Parameter /Macro-cell /Small-cell                  --------------------
%-----------------    Carrier frequency /    /2.1 GHz      --------------------------
%-----------------    Bandwidth   /      / 5 MHz      --------------------------
%-----------------    Node transmit power / 43 dBm / 23 dBm              --------------------------
%-----------------    Path loss model /  / 128.1 + 37.6 log10 (d[Km])    --------------------------
%-----------------    Number of UEs / 5 /  1 ? 5 UE per SCAP                  --------------------
%-----------------    Number of OFDM symbols for PDCCH / 3  / 3      --------------------------
%-----------------    BER threshold for PDCCH /   / 10?4      --------------------------
%-----------------    Number of RE Quadruplets per PDCCH /   / 18  --------------------------
%-----------------    Noise Figure at UE /      /9 dB    --------------------------
%-----------------    Thermal noise density  /     / ?174 dBm/Hz  --------------------------
%-----------------    Cell Radius  / 800m  /   50m    --------------------------
%-----------------------------------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%