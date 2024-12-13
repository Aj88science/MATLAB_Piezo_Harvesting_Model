%defining the constants
u0 = 0.05; % Base excitation value
i = sqrt(-1); % Imaginary Number constant
% for 1st and 2nd coupled elements
m1 = 0.03121 ; % Proof mass of oscillator one
m2 = 0.01763 ; % Proof mass of oscillator two
k1 = 409.5 ; % stiffness constant for oscillator one
k2 = 136.4 ; % stiffness constant for oscillator two
c1 = 30.72*10e-5 ; % damping constant of oscillator one
c2 = 30.72*10e-5 ; % damping constant of oscillator two
a1 = 4.16*10e-10 ; % fprce factor of piezo one
a2 = 4.16*10e-10 ; % force factor of piezo two
Cp1 = 1.08*10e-8 ; % clamped capacitance of piezo one
Cp2 = 7.569*10e-9; % clamped capacitance of piezo two
R1 = 1*10e6 ; % load resistance for piezo one
R2 = 1*10e6 ; % load resistance for piezo two
% for 3rd Element
m3 = 0.01896; % Proof mass of the oscillator in kg
k3 = 306.936; % stiffnes constant of the oscillator in N/m
c3 = 30.72*10e-5; % damping constant of the oscillator in Ns/m
a3 = 4.16*10e-10; % force factor of the piezo transducer in N/V
Cp3 = 7.94*10e-9; % clamped capacitance of the piezo transducer in F
R3 = 1*10e6; % Load resistance for piezo output in Ohm
% for 4th Element
m4 = 0.02161; % Proof mass of the oscillator in kg
k4 = 288.521; % stiffnes constant of the oscillator in N/m
c4 = 30.72*10e-5; % damping constant of the oscillator in Ns/m
a4 = 4.16*10e-10; % force factor of the piezo transducer in N/V
Cp4 = 8.654*10e-9; % clamped capacitance of the piezo transducer in F
R4 = 1*10e6; % Load resistance for piezo output in Ohm
% for 5th and 6th coupled elements
m5 = 0.02504 ; % Proof mass of oscillator 5
m6 = 0.01498 ; % Proof mass of oscillator 6
k5 = 348.1; % stiffness constant for oscillator 5
k6 = 158.5 ; % stiffness constant for oscillator 6
c5 = 30.72*10e-5 ; % damping constant of oscillator 5
c6 = 30.72*10e-5 ; % damping constant of oscillator 6
a5 = 4.16*10e-10 ; % fprce factor of piezo 5
a6 = 4.16*10e-10 ; % force factor of piezo 6
Cp5 = 9.48*10e-9 ; % clamped capacitance of piezo 5
Cp6 = 6.755*10e-9 ; % clamped capacitance of piezo 6
R5 = 1*10e6 ; % load resistance for piezo 5
R6 = 1*10e6 ; % load resistance for piezo 6
% defining frequency sweep range
f_range = 0:0.1:30 ;
w_range = 2*pi*f_range ;

% defining arrays for storing relative displacement values
u1 = zeros(size(w_range));
u2 = zeros(size(w_range));
u3 = zeros(size(w_range));
u4 = zeros(size(w_range));
u5 = zeros(size(w_range));
u6 = zeros(size(w_range));
% defining arrays for storing output voltage values
v1 = zeros(size(w_range));
v2 = zeros(size(w_range));
v3 = zeros(size(w_range));
v4 = zeros(size(w_range));
v5 = zeros(size(w_range));
v6 = zeros(size(w_range));
v_total = zeros(size(w_range));
% defining arrays for storing harvested power values
P1 = zeros(size(w_range));
P2 = zeros(size(w_range));
P3 = zeros(size(w_range));
P4 = zeros(size(w_range));
P5 = zeros(size(w_range));
P6 = zeros(size(w_range));
P_total = zeros(size(w_range));
for idx = 1:length(w_range)
w = w_range(idx);
% defining the matrix element equations for matrix M1 and N1
M1A1 = -(m1*w^2) + c1*w*1i + c2*w*1i + k1 + k2 ;
M1A2 = -(c2*w*1i + k2) ;
M1A3 = a1 ;
M1A4 = -a2 ;
M1B1 = -(c2*w*1i + k2) ;
M1B2 = -(m2*w^2) + c2*w*1i + k2 ;
M1B3 = 0 ;
M1B4 = a2 ;
M1C1 = -a1*w*1i ;
M1C2 = 0 ;
M1C3 = Cp1*w*1i + (1/R1) ;
M1C4 = 0 ;
M1D1 = a2*w*1i ;
M1D2 = -a2*w*1i ;
M1D3 = 0 ;
M1D4 = Cp2*w*1i + (1/R2) ;

N1E1 = k1 + c1*w*1i ;
N1E2 = 0 ;
N1E3 = -a1*w*1i ;
N1E4 = 0 ;
% Defining the matrix element equations for Matrix M3 and N3
M3A1 = -(m3*w^2) + k3 + (c3*1i*w) ;
M3A2 = a3 ;
M3B1 = (a3*1i*w);
M3B2 = -((1/R3) + Cp3*1i*w);
N3E1 = k3 + c3*1i*w;
N3E2 = a3*1i*w;
% Defining the matrix element equations for Matrix M4 and N4
M4A1 = -(m4*w^2) + k4 + (c4*1i*w) ;
M4A2 = a4 ;
M4B1 = (a4*1i*w);
M4B2 = -((1/R4) + Cp4*1i*w);
N4E1 = k4 + c4*1i*w;
N4E2 = a4*1i*w;
% defining the matrix element equations for matrix M5 and N5
M5A1 = -(m5*w^2) + c5*w*1i + c6*w*1i + k5 + k6 ;
M5A2 = -(c6*w*1i + k6) ;
M5A3 = a5 ;
M5A4 = -a6 ;
M5B1 = -(c6*w*1i + k6) ;
M5B2 = -(m6*w^2) + c6*w*1i + k6 ;
M5B3 = 0 ;
M5B4 = a2 ;
M5C1 = -a5*w*1i ;
M5C2 = 0 ;
M5C3 = Cp5*w*1i + (1/R5) ;
M5C4 = 0 ;
M5D1 = a6*w*1i ;
M5D2 = -a6*w*1i ;
M5D3 = 0 ;
M5D4 = Cp6*w*1i + (1/R6) ;
N5E1 = k5 + c5*w*1i ;
N5E2 = 0 ;
N5E3 = -a5*w*1i ;
N5E4 = 0 ;
% defining the matrices for coupled element 1 and 2

M1 = [
M1A1, M1A2, M1A3, M1A4;
M1B1, M1B2, M1B3, M1B4;
M1C1, M1C2, M1C3, M1C4;
M1D1, M1D2, M1D3, M1D4
];
N1 = [N1E1; N1E2; N1E3; N1E3];
% defining the matrix for element 3
M3 = [
M3A1, M3A2;
M3B1, M3B2
];
N3 = [N3E1;N3E2];
% defining the matrix for element 4
M4 = [
M4A1, M4A2;
M4B1, M4B2
];
N4 = [N4E1;N4E2];
% defining the matrices for coupled element 5 and 6
M5 = [
M5A1, M5A2, M5A3, M5A4;
M5B1, M5B2, M5B3, M5B4;
M5C1, M5C2, M5C3, M5C4;
M5D1, M5D2, M5D3, M5D4
];
N5 = [N5E1; N5E2; N5E3; N5E3];
% defining the matrix operations
X12 = u0*(M1\N1);
X3 = u0*(M3\N3) ;
X4 = u0*(M4\N4) ;
X56 = u0*(M5\N5) ;
%storing the generaed values in the resective arrays
u1(idx) = X12(1)/100;
u2(idx) = X12(2)/100;

u3(idx) = X3(1)/100;
u4(idx) = X4(1)/100;
u5(idx) = X56(1)/100;
u6(idx) = X56(2)/100;
v1(idx) = X12(3);
v2(idx) = X12(4);
v3(idx) = X3(2);
v4(idx) = X4(2);
v5(idx) = X56(3);
v6(idx) = X56(4);
v_total(idx) = X12(3) + X12(4) + X3(2) + X4(2) + X56(3) + X56(4);
P1(idx) = (abs(X12(3))^2)/R1 ;
P2(idx) = (abs(X12(4))^2)/R2 ;
P3(idx) = (abs(X3(2))^2)/R3 ;
P4(idx) = (abs(X4(2))^2)/R4 ;
P5(idx) = (abs(X56(3))^2)/R5 ;
P6(idx) = (abs(X56(4))^2)/R6 ;
P_total(idx) = (abs(X12(3))^2)/R1 + (abs(X12(4))^2)/R2 + (abs(X3(2))^2)/R3 + (abs(X4
(2))^2)/R4 + (abs(X56(3))^2)/R5 + (abs(X56(4))^2)/R6 ;
end
figure;
subplot(3,2,1);
plot(f_range, abs(u1));
title("frequency vs Rel. displacement 1");
xlabel("frequency");
ylabel("|u1|");
subplot(3,2,2);
plot(f_range, abs(u2));
title("frequency vs Rel. displacement 2");
xlabel("frequency");
ylabel("|u2|");
subplot(3,2,3);
plot(f_range, abs(u3));
title("frequency vs Rel. displacement 3");
xlabel("frequency");
ylabel("|u3|");
subplot(3,2,4);
plot(f_range, abs(u4));
title("frequency vs Rel. displacement 4");
xlabel("frequency");
ylabel("|u4|");

subplot(3,2,5);
plot(f_range, abs(u5));
title("frequency vs Rel. displacement 5");
xlabel("frequency");
ylabel("|u5|");
subplot(3,2,6);
plot(f_range, abs(u6));
title("frequency vs Rel. displacement 6");
xlabel("frequency");
ylabel("|u6|");
hold on;
figure;
subplot(3,2,1);
plot(f_range, abs(v1));
title("frequency vs Output voltage 1");
xlabel("frequency");
ylabel("|v1|");
subplot(3,2,2);
plot(f_range, abs(v2));
title("frequency vs Output voltaget 2");
xlabel("frequency");
ylabel("|v2|");
subplot(3,2,3);
plot(f_range, abs(v3));
title("frequency vs Output voltage 3");
xlabel("frequency");
ylabel("|v3|");
subplot(3,2,4);
plot(f_range, abs(v4));
title("frequency vs Output voltage 4");
xlabel("frequency");
ylabel("|v4|");
subplot(3,2,5);
plot(f_range, abs(v5));
title("frequency vs Output voltage 5");
xlabel("frequency");
ylabel("|v5|");
subplot(3,2,6);
plot(f_range, abs(v6));
title("frequency vs Output voltage 6");
xlabel("frequency");
ylabel("|v6|");

hold on;
figure;
subplot(3,2,1);
plot(f_range, P1);
title("frequency vs Harvested Power 1");
xlabel("frequency");
ylabel("|P1|");
subplot(3,2,2);
plot(f_range, P2);
title("frequency vs Harvested Power 2");
xlabel("frequency");
ylabel("|P2|");
subplot(3,2,3);
plot(f_range, P3);
title("frequency vs Harvested Power 3");
xlabel("frequency");
ylabel("|P3|");
subplot(3,2,4);
plot(f_range, P4);
title("frequency vs Harvested Power 4");
xlabel("frequency");
ylabel("|P4|");
subplot(3,2,5);
plot(f_range, P5);
title("frequency vs Harvested power 5");
xlabel("frequency");
ylabel("|P5|");
subplot(3,2,6);
plot(f_range, P6);
title("frequency vs Harvested Power 6");
xlabel("frequency");
ylabel("|P6|");
hold on
figure;
subplot(2,1,1);
plot(f_range, abs(v_total));
title("frequency vs Total Voltage Output");
xlabel("frequency");
ylabel("|V_total|");

subplot(2,1,2);
plot(f_range, P_total);
title("frequency vs Total Harvested Power");
xlabel("frequency");
ylabel("|P_total|");