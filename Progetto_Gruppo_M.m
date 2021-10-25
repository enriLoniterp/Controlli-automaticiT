%%%%%%% PROGETTO GRUPPO M - 1B %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                            %
%   MANFREDA FILIPPO                   %
%   ANDRINI ENRICO                     %
%   CRISTAUDO GIUSEPPE                 %
%   SARNERI ENRICO                     %
%   SANTANDREA PIETRO                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DEFINIZIONE DEL SISTEMA

% x1_dot = x2
% x2_dot = (-MpgLsin(x1)-K(x1-x3)-ro(x2-x4))/Jp
% x3_dot = x4
% x4_dot =(K(x1-x3)+ro(x2-x4)+u)/I
% y = x1

%% SPECIFICHE DEL PROGETTO

K = 3; ro = 0.2; L = 1; W=5;
I = 0.0075; Jp = 0.02; Mp = 0.07; g = 6.67*(10^-11);
w_n = 100; A_n = 0.02; Bn = 30;
h = 5; Tah = 0.5; 
Ta0 = 0.25;
x_1r=pi/2;
x_1g=90;
u_=Mp*g*L*sin(x_1r);
x_2=0;
x_3=(u_/K)+x_1g;
x_4=0;

%% MATRICI LINEARIZZATE

A=[0,1,0,0;
   -K/Jp,-ro/Jp,K/Jp,ro/Jp;
   0,0,0,1;
   K/I,ro/I,-K/I,-ro/I];
B=[0;0;0;1/I];
C=[1,0,0,0];
D=0;

%% DEFINIZIONE FUNZIONE DI TRASFERIMENTO

s = tf('s');
G = ((K+ro*s)/s^2)*(1/(I*(s*(Jp*s+ro)+K)+Jp*(K+ro*s)));

zpk(G)
w_plot_min=10^(-2);
w_plot_max=10^5;

[Mag,phase,w]=bode(G,{w_plot_min,w_plot_max});

%% SPECIFICHE DI CONTROLLO

%Ingresso
wt=W*1;

%Margine id Fase
xi=sqrt(log(0.01)^2/(pi^2+log(0.01)^2));
Mf=xi*100;
Ta5=0.5;

%W min e max
wc_min=3/0.5/xi;
wc_max=w_n;
patch([w_plot_min,wc_min,wc_min,w_plot_min],[-300,-300,0,0],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([wc_max,w_plot_max,w_plot_max,wc_max],[0,0,200,200],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([wc_max,w_plot_max,w_plot_max,wc_max],[-29.5,-29.5,200,200],'y','FaceAlpha',0.3,'EdgeAlpha',0);
hold on;
margin(Mag,phase,w);
grid on;
patch([wc_min,wc_max,wc_max,wc_min],[-180+Mf,-180+Mf,-270,-270],'y','FaceAlpha',0.3,'EdgeAlpha',0);

%% PARTE OPERATIVA

%Regolatore statico
R_s=24;

%Sistema esteso
Ge=G*R_s;

%Regolatore dinamico
R_d =((s*(s^2 + 36.67*s + 550))/((s+41)^3));

%Regolatore
R=R_d*R_s;

%Funzione d'anello
Ls=G*R;
[Mag,phase,w]=bode(Ls,{w_plot_min,w_plot_max});
margin(Ls);
hold on;
grid on;

%Funzione di sensitività complementare
F=Ls/(1+Ls);

hold on;
margin(F);
figure();
step(F);
figure();
rlocus(Ls);
open("Simulink_Gruppo_M.slx");
