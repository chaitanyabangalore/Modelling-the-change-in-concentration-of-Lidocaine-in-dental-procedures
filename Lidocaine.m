%Analytical Methods Project 1

%Modelling the change in concentration of Lidocaine used as local anesthesia in dental procedures using System of Ordinary Differential Equations

%Submitted by Chaitanya Bangalore, Pujitha Vasantha Gopal

clc; clear all;

%initializing constant parameters
k1 = 1.25*(10^-10);            %diffusivity of lidocaine in the tissue
Vmax = 1.45*(10^-7);           %Michelis Menten constant for Lidocaine (M/s^5)
Km = 143*(10^-6);              %Michelis Menten constant for Lidocaine (micrometer)
k2 = Vmax/Km; 

%input of system of differential equations
syms A(t) C(t) M(t)
ode1 = diff(A) == -k1*A;       %Injecting the drug into the gums
ode2 = diff(C) == k1*A - k2*C; %Circulation of the drug in the blood stream
ode3 = diff(M) == k2*C;        %Metabolism of the drug at the target area
odes = [ode1; ode2; ode3];

%solving the differential equations uding the dsolve function that returns
%the solution as a structure
S = dsolve(odes);
ISol(t) = S.A;
TSol(t) = S.C;
MSol(t) = S.M;

%input of the initial conditions to solve the system of equations
cond1 = A(0) == 0.08535;         %Initial dose of Lidocaine administered in mol/l
cond2 = C(0) == 0;               %Initial concentration of Lidocaine in bloodstream
cond3 = M(0) == 0;               %Initial concentration of Lidocaine at target area for drug to metabolize
conds = [cond1; cond2; cond3];
[ASol(t), CSol(t), MSol(t)] = dsolve(odes,conds);

%use of Variable Precision Arithmetic Function to get precise solutions to
%the differential equations
A(t) = vpa(ASol)
C(t) = vpa(CSol)
M(t) = vpa(MSol)

%Plotting the solutions against time
fplot(ASol,[0,3600],'b')        %time period:3600s
grid on
xlabel('time in seconds') 
ylabel('A(t)')
title('Change of Lidocaine concentration at the gums')
figure

fplot(CSol,[0,3600],'m')       %time period:3600s
grid on
xlabel('time in seconds') 
ylabel('C(t)') 
title('Circulation of Lidocaine in the target area')
figure

fplot(MSol,[0,3600],'r')       %time period:3600s
grid on
xlabel('time in seconds') 
ylabel('M(t)') 
title('Metabolism of Lidocaine at the target area')

