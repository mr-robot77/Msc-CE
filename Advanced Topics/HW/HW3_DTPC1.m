%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------------------------------------%
%------------------  HW3 : DTPC - 1                    ------------------------------------
%-----------------   Author : Saleh Abednejad -- 961115116      ---------------
%-----------------   Prof : Dr. Rajabi                      -------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------    Background noise power 2 = 10?10          ---------
%-----------------    OPC constant  = 10?4      --------------------------
%-----------------    Path gain hi = 0:09d?3.      --------------------------
%-----------------------------------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; 
clear;
close all;
%%
noOfNodes  = 6;
L=500;
rx = rand(1,noOfNodes)*L; 
ry = rand(1,noOfNodes)*L;
bsx=L/2;
bsy=L/2;
bs=[bsx 0 0 0 0 0;bsy 0 0 0 0 0];
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
G=[0.05, 0.05, 0.05,0.05,0.05,0.05]; %target SINR
E=[0.0001, 0.0001, 0.0001,0.0001,0.0001,0.0001]; %eta
P=[0.001 ,0.001 ,0.001, 0.001 ,0.001, 0.001]; %  Transmit Power
N=[1e-10, 1e-10, 1e-10,1e-10,1e-10,1e-10]; % Noise 
SIR=[H(1)*P(1)/(H(5)*P(5)+H(4)*P(4)+H(3)*P(3)+H(2)*P(2)+N(1)),...
      H(2)*P(2)/(H(5)*P(5)+H(4)*P(4)+H(3)*P(3)+H(1)*P(1)+N(2)),...
       +H(3)*P(3)/(H(5)*P(5)+H(4)*P(4)+H(1)*P(1)+H(2)*P(2)+N(3)),...
         H(4)*P(4)/(H(5)*P(5)+H(1)*P(1)+H(3)*P(3)+H(2)*P(2)+N(4)),...
          +H(5)*P(5)/(H(1)*P(1)+H(4)*P(4)+H(3)*P(3)+H(2)*P(2)+N(5)),...
            H(6)*P(6)/(H(1)*P(1)+H(4)*P(4)+H(3)*P(3)+H(2)*P(2)+H(5)*P(5)+N(6))];% initial SINR 
%%
%algorithm 
iterations=1;
while iterations<25  
    P1=(G./SIR(iterations,:)).*P; % New power TPC
    P2=(E.*SIR(iterations,:))./P; % New power OPC
    P=max(P1,P2);
    if P<1   
    iterations=iterations+1;
    power(iterations,:)=P;
SIR(iterations,:)=[H(1)*P(1)/(H(5)*P(5)+H(4)*P(4)+H(3)*P(3)+H(2)*P(2)+N(1)),...
                    H(2)*P(2)/(H(5)*P(5)+H(4)*P(4)+H(3)*P(3)+H(1)*P(1)+N(2)),...
                      H(3)*P(3)/(H(5)*P(5)+H(4)*P(4)+H(1)*P(1)+H(2)*P(2)+N(3)),...
                        H(4)*P(4)/(H(5)*P(5)+H(1)*P(1)+H(3)*P(3)+H(2)*P(2)+N(4)),...
                          H(5)*P(5)/(H(1)*P(1)+H(4)*P(4)+H(3)*P(3)+H(2)*P(2)+N(5)),...
                            H(6)*P(6)/(H(1)*P(1)+H(4)*P(4)+H(3)*P(3)+H(2)*P(2)+H(5)*P(5)+N(6))];
    else
    iterations=iterations+1;
    power(iterations,:)=1;
SIR(iterations,:)=[H(1)*P(1)/(H(5)*P(5)+H(4)*P(4)+H(3)*P(3)+H(2)*P(2)+N(1)),...
                    H(2)*P(2)/(H(5)*P(5)+H(4)*P(4)+H(3)*P(3)+H(1)*P(1)+N(2)),...
                     H(3)*P(3)/(H(5)*P(5)+H(4)*P(4)+H(1)*P(1)+H(2)*P(2)+N(3)),...
                      H(4)*P(4)/(H(5)*P(5)+H(1)*P(1)+H(3)*P(3)+H(2)*P(2)+N(4)),...
                       H(5)*P(5)/(H(1)*P(1)+H(4)*P(4)+H(3)*P(3)+H(2)*P(2)+N(5)),...
                        H(6)*P(6)/(H(1)*P(1)+H(4)*P(4)+H(3)*P(3)+H(2)*P(2)+H(5)*P(5)+N(6))];
    end
end
    power;
%%
%ploting 
figure(2);
hold on;
plot(1:iterations,SIR(:,1),'-.sc');
plot(1:iterations,SIR(:,2),'-*b');
plot(1:iterations,SIR(:,3),':ok');
plot(1:iterations,SIR(:,4),'-^r');
plot(1:iterations,SIR(:,5),'-pg');
plot(1:iterations,SIR(:,6),'--dm');
xlabel('Iterations');
ylabel('SIR');
title('SIR vs number of Iterations');
legend(' SIR of user 1',' SIR of user 2',' SIR of user 3',' SIR of user 4',' SIR of user 5',' SIR of user 6');
%%
figure(3);
plot(1:iterations,power(:,1),'-.sc');
hold on;
plot(1:iterations,power(:,2),'-*b');
plot(1:iterations,power(:,3),':ok');
plot(1:iterations,power(:,4),'-^r');
plot(1:iterations,power(:,5),'-pg');
plot(1:iterations,power(:,6),'--dm');
xlabel('Iterations');
ylabel('power');
title('power ');
legend(' power of user 1',' power of user 2',' power of user 3',' power of user 4',' power of user 5',' SINR of user 6');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------------------------------------%
%------------------  HW3 : DTPC - 1                    ------------------------------------
%-----------------   Author : Saleh Abednejad -- 961115116      ---------------
%-----------------   Prof : Dr. Rajabi                      -------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------    Background noise power 2 = 10?10          ---------
%-----------------    OPC constant  = 10?4      --------------------------
%-----------------    Path gain hi = 0:09d?3.      --------------------------
%-----------------------------------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%