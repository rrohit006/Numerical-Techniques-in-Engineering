clc
close all
clear
N=1000;
%Initial values
i_0=0; %A
i_n=132.56*10^(-3) %A
t_0=0; %s
t_n=50*10^(-3); %s
%The parameters are:
Vs=5 %V
R=0.1 %ohm
L=10^(-3) %H
C=10^(-3) %F
h=(t_n-t_0)/(N+1);
%Coefficient Values
a=(L-(R/2)*h);
b=(((h^2)/C)-(2*L));
c=(L+(R/2)*h);
d=0;
A=diag(b*ones(1,N))+diag(c*ones(1,N-1),1)+diag(a*ones(1,N-1),-1);
B_1=repmat(d,N-2,1);
B=[d-a*i_0
 B_1
 d-c*i_n];
T=inv(A)*B;
t=linspace(t_0,t_n,N+2)
i=[i_0
 T
 i_n];
%Exact Solution (03)
alp=R/(2*L);
w=1/(sqrt(L*C));
s_1=-alp+sqrt((alp^2)-(w^2));
s_2=-alp-sqrt((alp^2)-(w^2));
A=-Vs/(L*(s_1-s_2));
i_exact=A*exp(s_1*t)-A*exp(s_2*t);
%Find and Plot (04)
%04(a)
x_L=w*L;
V_L=i.*x_L
figure(1)
plot(t,V_L);
xlabel('t');
ylabel('V_L');
grid on
%04(b)
PR=(i.^2)*R;
figure(2)
plot(t,PR);
xlabel('t');
ylabel('PR');
grid on
%04(c)
PL=(i.^2)*x_L;
figure(3)
plot(t,PL);
xlabel('t');
ylabel('PL');
grid on
%Visualizations and plots
%05(a)
figure(4)
plot(t,i);
xlabel('t');
ylabel('i');
grid on;
figure(5)
plot(t,i); hold on
figure(5)
plot(t,i_exact,'r--');
plot(t,i,'g');
xlabel('t');
ylabel('i & i_exact');
grid on
ei=i_exact'-i;
En=sqrt (sum(ei.^2)/(N+2));
%05(b)
%N=100
figure(6)
plot(t,V_L);
xlabel('t');
ylabel('V_L');
grid on
%05(c)
%N=100
figure(7)
plot(t,PR);
xlabel('t');
ylabel('PR');
grid on
figure(8)
plot(t,PR,'r'); hold on
PR_exact=(i_exact.^2)*R;
figure(8)
plot(t,PR_exact,'r--');
plot(t,PR,'g');
xlabel('t');
ylabel('PR & PR_exact');
grid on
