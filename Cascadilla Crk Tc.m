%Initialize Constants

L = 48000;                 %[ft] longest flow path
h1 = 380;                   %[ft] lowest elevation point
h2 = 1870;                  %[ft] highest elevation point
H = h2-h1;                  %[ft] elevation difference along longest flow path
S = H/L;                  % average slope along the longest flow path

CN = 70*0.4 + 85*0.23 + 74*0.2 + 83*0.15 + 98*0.02; 
                       %note row crop CN for straight row, 1/4 acre for residential
C = 0.4*(0.02*1.27) + 0.23*(0.47*1.09) + 0.2*(0.02*1.21) + 0.15*0.5 + 0.02*0.95;
                       %residential- single family units and upper limit,
                       %pavement- upper limit of all options
                       
%Find rainfall intensity
t = 15;             %years
p = 0.25;           %exceedence probability
EP = 1-((1-p)^(1/t));
T = 1/EP;           % 50-yr rainfall event with 2-hr duration
intensity = 2.40/2*25.4          %[mm/yr] found using NRCC data

t_c1 = (0.0078 * L^(0.77) * S^(-0.385)) / 60;      %[hrs] Kirpich (1940)
t_c2 = L^1.15 / (7700 * H^0.38);                  %[hrs] Soil Conservation Service (1972)
t_c3 = (10 * L^0.8 * ((1000/CN)-9)^0.7 / (1900*S^0.5)) / 60; 
                                                   %[hrs] SCS Lag Equation (1973)
t_c4 = (1.8 * (1.1-C) * L^0.5 * S^(-0.333)) / 60;  %[hrs] FAA (1970)

%% Rational Method 
A = 3.5*10^7;                       %[m^2] area of watershed 
i_m_sec = (intensity*0.001)/(3.154*10^7);   %Intensity in [m/sec]
qp_rational = C*i_m_sec*A          % Peak runoff in [m^3/sec;] 

%% CN Synthetic Triangular Hydrograph Method
P = 2.4;            %[in] rainfall depth
St = 1000/CN-10;    %[in] watershed storage
Ia_1 = 0.2*St;      %[in] initial abstraction #1
Ia_2 = 0.05*St;     %[in] initial abstraction #2

Q_1 = ((P-Ia_1)^2 / (P-Ia_1+St))*0.0254; %[m] runoff depth #1
Q_2 = (P-Ia_2)^2 / (P-Ia_2+St)*0.0254; %[m] runoff depth #2

qp_1 = 2*Q_1*A / (2.937*t_c1*3600);  %[m^3/s] I_a #1, Kirpich
qp_2 = 2*Q_1*A / (2.937*t_c2*3600);  %[m^3/s] I_a #1, SCS
qp_3 = 2*Q_2*A / (2.937*t_c1*3600);  %[m^3/s] I_a #2, Kirpich
qp_4 = 2*Q_2*A / (2.937*t_c2*3600);  %[m^3/s] I_a #2, SCS

%Just in case here's the graphs
tp_tri_1 = 1.1*t_c1;                 %[hrs] peak time, Kirpich
tp_tri_2 = 1.1*t_c2;                 %[hrs] peak time, SCS
tr_tri_1 = 1.67*tp_tri_1;            %[hrs] recession time, Kirpich
tr_tri_2 = 1.67*tp_tri_2;            %[hrs] recession time, SCS

x_tri_1 = [0, tp_tri_1, tp_tri_1+tr_tri_1];     %Kirpich
x_tri_2 = [0, tp_tri_2, tp_tri_2+tr_tri_2];     %SCS

y_tri_1 = [0, qp_1, 0];
y_tri_2 = [0, qp_2, 0];
y_tri_3 = [0, qp_3, 0];
y_tri_4 = [0, qp_4, 0];

figure(1)
plot(x_tri_1, y_tri_1, 'r')
hold on
plot(x_tri_2, y_tri_2, 'b')
hold on
plot(x_tri_1, y_tri_3, 'm')
hold on
plot(x_tri_2, y_tri_4, 'g')
hold off

legend('Ia = 0.2S, Kirpich', 'Ia = 0.2S, SCS', 'Ia = 0.05S, Kirpich',...
    'Ia = 0.05S, SCS')
xlabel('Time [hrs]')
ylabel('Storm Runoff [m^3/s]')
title('Synthetic Triangular Hydrograph for a 50 year, 2 hour Precipitation Event in Cascadilla Creek')

%% CN Unit Hydrograph Method
tc=2;        %Rounded tc value for both Kirpich and NRCS 1972
delta = 0.133*2;
tp= (delta/2)+0.6*tc;
tr=1.67*tp; 
Qv =(0.001)*A;   %[m^3]
qp_unit = (Qv/((0.5)*(tp+tr)))/3600;     %[m^3/sec] 
slope1 = qp_unit/tp;                     %[m^3/sec/hr]
slope2 = qp_unit/tr;                     %[m^3/sec/hr]
step1= tp/delta;
u_step = qp_unit/step1;
counter1 = [0 1 2 3 4 5];
x1 = delta*counter1;
y1=u_step*counter1;
counter2 = [5 6 7 8 9 10 11 12 13];
x2 = delta*counter2;
b = slope2*delta;
counter3 = [1 2 3 4 5 6 7 8] ;
different = b*counter3;
ones = [1 1 1 1 1 1 1 1]*qp_unit;
y2(1) = qp_unit;
y2(2:9)= ones - different;

figure(2)
plot(x1,y1)
hold on
plot(x2,y2)
hold off

xlabel('Time [hrs]')
ylabel('Unit Runoff Rate [m^3/s/mm]')
title('Unit Hydrograph for Cascadilla Creek')

%% Making the Composite Hydrograph 

del = [ x1 x2 ]' % delta values

U = [ y1 y2]' % u values

index = find (del >2,1) % fnd index for over tc, that way you know where P ends up 
deltaP = P/ index 
 
P_val = [1 2 3 4 5 6 7 8 9 10 ]* deltaP % ascending P values up to the P for tc 
S = 1000/CN-10    %[in] watershed storage
I = 0.2*S      %[in] initial abstraction 

Q = (P_val - I).^2 ./ (P_val- I +S) %[in]

index_P = find (P_val >I,1) % fnd index for over tc, that way you know where P ends up 

Q(1:index_P) = 0 % this demonstrates the vector where we assume no runoff is generated until the cumulative rainfall is ~Ia (i.e., 0.2S)
delta_Q = diff(Q')

index_Q_del = find(delta_Q > 0, 1)
% 
num_of_q = length(delta_Q) - index_Q_del -1 
matrix = zeros(length(delta_Q) +length(U), num_of_q)
k = 1
for i = 1: num_of_q
matrix(index_Q_del : index_Q_del + length(U) -1 ,k) = (U .*  delta_Q(index_Q_del)')
 k = k+1
 index_Q_del = index_Q_del + 1 

end 

q_total = sum(matrix, 2) * 25.4 %[m^3/sec/mm]

T_values = (1:1:length(q_total)) * delta 

figure(3)
plot(T_values, q_total)
hold on

xlabel('Time [hrs]')
ylabel('Storm Runoff [m^3/s]')
title('Composite Unit Hydrograph for a 50 year, 2 hour Precipitation Event in Cascadilla Creek')

%% Composite Hydrograph with Ia = 0.05S

del = [ x1 x2 ]' % delta values

U = [ y1 y2]' % u values

index = find (del >2,1) % fnd index for over tc, that way you know where P ends up 
deltaP = P/ index 
 
P_val = [1 2 3 4 5 6 7 8 9 10 ]* deltaP % ascending P values up to the P for tc 
S = 1000/CN-10    %[in] watershed storage
I = 0.05*S      %[in] initial abstraction 

Q = (P_val - I).^2 ./ (P_val- I +S) 

index_P = find (P_val >I,1) % fnd index for over tc, that way you know where P ends up 

Q(1:index_P) = 0 % this demonstrates the vector where we assume no runoff is generated until the cumulative rainfall is ~Ia (i.e., 0.2S)
delta_Q = diff(Q')

index_Q_del = find(delta_Q > 0, 1)
% 
num_of_q = length(delta_Q) - index_Q_del -1 
matrix = zeros(length(delta_Q) +length(U), num_of_q)
k = 1
for i = 1: num_of_q
matrix(index_Q_del : index_Q_del + length(U) -1 ,k) = (U .*  delta_Q(index_Q_del)')
 k = k+1
 index_Q_del = index_Q_del + 1 

end 

q_total = sum(matrix, 2) * 25.4  %[m^3/sec/mm]

T_values = (1:1:length(q_total)) * delta 

plot(T_values, q_total)
hold off

legend('Ia = 0.2S', 'Ia = 0.05S')
