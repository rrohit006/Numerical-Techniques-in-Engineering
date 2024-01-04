clc
close all
clear
n=[75:25:1000 2000 3500 5000];
for j=1:length(n)
N=n(j)
 
%Initial values
i_0=0; %A
i_n=132.56*10^(-3); %A
t_0=0; %s
t_n=50*10^(-3); %s
%The parameters are:
Vs=5; %V
R=0.1; %ohm
L=10^(-3); %H
C=10^(-3); %F
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
t=linspace(t_0,t_n,N+2);
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
ei=i_exact'-i;
En(j)=sqrt (sum(ei.^2)/(N+2));
end
figure(1)
semilogy(n,En,'o-');
xlabel('n');
ylabel('En');
grid on;
