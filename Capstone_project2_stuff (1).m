clear
clc

load Stuff2Data_Names.mat
load Stuff2Data_Numbers.mat
%% Initialize Constants
Num = Stuff2Data;
Name = Stuff2Data1;

L = Num(:,6);            %[ft] longest flow path
Sl = Num(:,3);           %[unitless] average slope along the longest flow path
P = Num(:,8);            %[in] rainfall depth (tc value in hrs at 50 yr event form NRCC)
A = Num(:,1)*1000^2;     %[m^2]                   
N = length(L);

CN = 70*0.4 + 85*0.23 + 74*0.2 + 83*0.15 + 98*0.02; 
S = 1000/CN-10;    %[in] watershed storage
Ia = 0.05*S;     %[in] initial abstraction

%% CN Synthetic Triangular Hydrograph Method
tc = zeros(N,1);     %[hrs] Kirpich (1940)
Q_tri = zeros(N,1);  %[m] runoff depth
qp_tri = zeros(N,1); %[m^3/s] 
x_tri = zeros(N,3);
y_tri = zeros(N,3);

for i = 1:N
    tc(i) = (0.0078 * L(i)^(0.77) * Sl(i)^(-0.385)) / 60; 
    Q_tri(i) = (P(i)-Ia)^2 / (P(i)-Ia+S)*0.0254; 
    qp_tri(i) = 2*Q_tri(i)*A(i) / (2.937*tc(i)*3600);
    
end

tp_tri = 1.1*tc;                 %[hrs] peak time, Kirpich
tr_tri = 1.67*tp_tri;            %[hrs] recession time, Kirpich
%style = ['--','b-','r:','r-.','b--','g--','c--','m--','y-','k--',...
%    'b:','g:','m:','c:','k:','m-.','g-.'];

for j = 1:N
    x_tri(j,2) = tp_tri(j);
    x_tri(j,3) = tp_tri(j)+tr_tri(j);
    y_tri(j,2) = qp_tri(j);
    plot(x_tri(j,:), y_tri(j,:))
    hold on
end

hold off

legend(Name(1), Name(2), Name(3), Name(4), Name(5), Name(6), Name(7),...
    Name(8), Name(9), Name(10), Name(11), Name(12), Name(13),Name(14),...
    Name(15), Name(16), Name(17))
xlabel('Time [hrs]')
ylabel('Storm Runoff [m^3/s]')
title('Synthetic Triangular Hydrograph of 17 Different Stream Crossings for a 50 year Return Period')

