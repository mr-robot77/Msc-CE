%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------------------------------------%
%------------------  HW1 : TPC     ----------------------------------------
%-----------------   Author : Saleh Abednejad -- 961115116  ---------------
%-----------------   Prof : Dr. Rajabi           --------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------    Number of users Nu = 5                  --------------------
%-----------------    Background noise power 2 = 10?10      --------------------------
%-----------------    Maximum power of each user Pi = 1mW      --------------------------
%-----------------    Target SINR ^= 0:05              --------------------------
%-----------------    Path gain hi = 0:09d?3           --------------------------
%-----------------------------------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
%%
noOfNodes  = 5;
L=500;
rx = rand(1,noOfNodes)*L; 
ry = rand(1,noOfNodes)*L;
bsx=L/2;
bsy=L/2;
bs=[bsx 0 0 0 0;bsy 0 0 0 0];
figure(1);
clf;
hold on;
grid on;
plot(bsx, bsy, '*r');
title('CDMA Network');
xlim ([0 500]);
ylim ([0 500]);
for i = 1:noOfNodes
    plot(rx(i), ry(i), 'p');
    text(rx(i), ry(i), num2str(i));
end;
%hold off;
for j=1:noOfNodes
     z(j)=complex(rx(j),ry(j)); %User
 end     
for i=1:noOfNodes
    H(i)=0.09*(abs(z(i)-(bs(i)))^(-3)); %Path Gain
end
%hold off;
Gamm=[0.05, 0.05, 0.05,0.05,0.05]; %target SIR 
P=[.001 .001 .001 .001 .001]; %  Transmit Power
N=[1e-10, 1e-10, 1e-10,1e-10,1e-10]; % Noise 
SIR=[H(1)*P(1)/(H(5)*P(5)+H(4)*P(4)+H(3)*P(3)+H(2)*P(2)+N(1)),...
      H(2)*P(2)/(H(5)*P(5)+H(4)*P(4)+H(3)*P(3)+H(1)*P(1)+N(2)),...
       +H(3)*P(3)/(H(5)*P(5)+H(4)*P(4)+H(1)*P(1)+H(2)*P(2)+N(3)),...
         H(4)*P(4)/(H(5)*P(5)+H(1)*P(1)+H(3)*P(3)+H(2)*P(2)+N(4)),...
          +H(5)*P(5)/(H(1)*P(1)+H(4)*P(4)+H(3)*P(3)+H(2)*P(2)+N(5))];% initial SIR at each reciever
%%
%algorithm starts here
iterations=1;
while iterations<20      
    P=(Gamm./SIR(iterations,:)).*P; % New power used by transmitters
    iterations=iterations+1;
    power(iterations,:)=P;
    SIR(iterations,:)=[H(1)*P(1)/(H(5)*P(5)+H(4)*P(4)+H(3)*P(3)+H(2)*P(2)+N(1)),...
                        H(2)*P(2)/(H(5)*P(5)+H(4)*P(4)+H(3)*P(3)+H(1)*P(1)+N(2)),...
                          H(3)*P(3)/(H(5)*P(5)+H(4)*P(4)+H(1)*P(1)+H(2)*P(2)+N(3)),...
                            H(4)*P(4)/(H(5)*P(5)+H(1)*P(1)+H(3)*P(3)+H(2)*P(2)+N(4)),...
                              H(5)*P(5)/(H(1)*P(1)+H(4)*P(4)+H(3)*P(3)+H(2)*P(2)+N(5))];  
end
    power;
%%
%ploting 
figure(2);
hold on;
plot(1:iterations,SIR(:,1),'-.sb');
plot(1:iterations,SIR(:,2),'-*k');
plot(1:iterations,SIR(:,3),':or');
plot(1:iterations,SIR(:,4),'-pg');
plot(1:iterations,SIR(:,5),'-dm');
xlabel('Iterations');
ylabel('SIR');
title('SIR vs number of Iterations');
legend(' SIR of user 1',' SIR of user 2',' SIR of user 3',' SIR of user 4',' SIR of user 5');
grid on;
%%
figure(3);
plot(1:iterations,power(:,1),'-.sb');
hold on;
plot(1:iterations,power(:,2),'-*k');
plot(1:iterations,power(:,3),':or');
plot(1:iterations,power(:,4),'-pg');
plot(1:iterations,power(:,5),'-dm');
xlabel('Iterations');
ylabel('power');
title('power ');
legend(' power of user 1',' power of user 2',' power of user 3',' power of user 4',' power of user 5');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------------------------------------%
%------------------  HW1 : TPC     ----------------------------------------
%-----------------   Author : Saleh Abednejad -- 961115116  ---------------
%-----------------   Prof : Dr. Rajabi           --------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------    Number of users Nu = 5                  --------------------
%-----------------    Background noise power 2 = 10?10      --------------------------
%-----------------    Maximum power of each user Pi = 1mW      --------------------------
%-----------------    Target SINR ^= 0:05              --------------------------
%-----------------    Path gain hi = 0:09d?3           --------------------------
%-----------------------------------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%