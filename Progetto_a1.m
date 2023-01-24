%% Progetto traccia a1
% Gruppo Q
%   Clarissa Bovo
%   Fabio Cangiano
%   Francesca Porzia Fedi
%   Antonio Zara

clear all; close all; clc

%% SPECIFICHE
% -)Errore a regime minore o uguale a 0.01 in risposta a un gradino
%   w(t)=2*1(t) e d(t)=2*1(t)
%
% -)Mf maggiore o uguale a 30°
%
% -)S% minore o uguale a 1%
%
%-)Ta,1=1
%
% -)Attenuazione di almeno 35db per d(t) tra
%       [omega_d_min,omega_d_max] = [0,0.1]
%
% -)Attenuazione di almeno 80db per n(t) tra
%       [omega_n_min,omega_n_max] = [10^3,10^6]

% Ampiezza gradini
WW = 2;
DD = 2;

% Errore a regime
e_star = 0.01;

% Sovraelongazione massima e Tempo di assetamento all'1%
S_100_spec = 0.01;
T_a1_spec = 1;

% Attenuazione disturbo sull'uscita
A_d = 35;
omega_d_MAX = 0.1;

% Attenuazione disturbo di misura
A_n = 80;
omega_n_min = 1e3;
omega_n_MAX = 1e6;

%% Creazione sistema

% Parametri forniti
J=10; % [R] indica il momento di inerzia del drone rispetto all'asse di rotazione passante per il baricentro

beta=0.5; % [R] indica il coefficiente di attrito dovuto alla presenza dell'aria

a=0.01; % [R] indica la semi ampiezza planare del drone

F_v=-5; % [R] indica la forza costante do vuta all'azione del vento

theta_e=pi/3; %Valore di equilibrio per x1

% Matrici per la costruzione del sistema linearizzato
A=[0, 1; ((a*F_v)/(2*J))*cos(theta_e), -beta/J];
B=[0; a/J];
C=[1, 0];
D=0;

% Creazione del sistema linearizzato nel dominio del tempo
SYS = ss(A,B,C,D);

% Passaggio nel dominio di Laplace
s = tf('s'); %Variabile di traferimento

GG = tf(SYS); %Funzione di trasferimento G(s)

%omega_n=sqrt(0.00125)
%psi=0.05/(2*omega_n)
%mu=0.001/(omega_n)^2

%% Diagramma di Bode

figure(1);
bode(GG);
xlim([1e-3, 1e4]);
grid on, zoom on;

% Stop qui per diagramma di bode G(s)
if 0
    return
end

%% Regolatore statico - proporzionale senza poli nell'origine

% Valore minimo prescritto per L(0)
mu_s = (DD+WW)/e_star;

% Guadagno minimo del regolatore ottenuto come L(0)/G(0)
G_0 = abs(evalfr(GG,0));
RR_s = mu_s / G_0; 

% Sistema esteso Ge(s)
GG_e = RR_s*GG;

%% Diagrammi di Bode di Ge con specifiche

figure(2);
hold on;


% Calcolo specifiche S% => Margine di fase
logsq = (log(S_100_spec))^2;
xi = sqrt(logsq/(pi^2+logsq));
Mf_spec = xi*100;
omega_Ta_MAX = 460/(Mf_spec*T_a1_spec); 

% Specifiche su d
omega_d_min = 0.0001; % lower bound per il plot
Bnd_d_x = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
Bnd_d_y = [A_d; A_d; -150; -150];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche su n
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 100; 100];
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche tempo d'assestamento (minima pulsazione di taglio)
omega_Ta_min = 1e-4; % lower bound per il plot
omega_Ta_MAX = 460/(Mf_spec*T_a1_spec); % omega_c >= 460/(Mf*T^*) ~ 4.6
Bnd_Ta_x = [omega_Ta_min; omega_Ta_MAX; omega_Ta_MAX; omega_Ta_min];
Bnd_Ta_y = [0; 0; -150; -150];
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_mag = ["A_d", "A_n", "\omega_{c,min}", "G(j\omega)"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(GG);
xlim([1e-3, 1e6]);
grid on; zoom on;

% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_MAX;
omega_c_MAX = omega_n_min;

phi_spec = Mf_spec - 180;
phi_low = -270; % lower bound per il plot

Bnd_Mf_x = [omega_c_min; omega_c_MAX; omega_c_MAX; omega_c_min];
Bnd_Mf_y = [phi_spec; phi_spec; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;

% Legenda colori
Legend_arg = ["G(j\omega)"; "M_f"];
legend(Legend_arg);

% Stop qui per diagrama di bode Ge(s)
if 0
    return
end

%% Regolatore dinamico - Rete Anticipatrice

% Parametri per la creazione della rete
Mf_star = Mf_spec+5; % Margine di fase per L(s) finale

omega_c_star = omega_c_min+1; % Pulsazione alla quale vorremo |L(j*omega_c_star)|db = 0db

% Formule di inversione per il calcolo di alpha e tau
[mag_omega_c_star, arg_omega_c_star, ~] = bode(GG_e, omega_c_star); % Calcolo ampezza e fase rispetto a omega_c_star di Ge(s)
mag_omega_c_star_db = 20*log10(mag_omega_c_star);

M_star = 10^(-mag_omega_c_star_db/20);
phi_star = Mf_star - 180 - arg_omega_c_star;

alpha_tau = (cos(phi_star*pi/180) - 1/M_star)/(omega_c_star*sin(phi_star*pi/180)); %è passato da gradi a radianti facendo pi/180
tau = (M_star - cos(phi_star*pi/180))/(omega_c_star*sin(phi_star*pi/180));
alpha = alpha_tau / tau;


check_flag = min(tau, alpha_tau);
if check_flag < 0
    disp('Errore: polo/zero positivo');
    return;
end

%Costruzione Rete Anticipatrice e L(s)
R_d = (1 + tau*s)/(1 + alpha_tau*s); % Rete Anticipatrice

% Aggiungiamo un polo a parte reale negativa per tirare giù l'ampiezza e
% mantenere la stabilità del sistema

Polo = 1/(1+(1/80)*s);

LL = R_d*Polo*GG_e; % L(s)

%% Diagrammi di Bode con specifiche includendo regolatore dinamico

figure(3);
hold on;


% Specifiche su ampiezza
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(LL);
xlim([1e-3, 1e6])
grid on; zoom on;

% Specifiche su fase
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;
legend(Legend_arg);

% Stop qui per diagramma di bode di L(s)
if 0
    return 
end

%% Check prestazioni in anello chiuso

% Funzione di sensitività complementare
FF = LL/(1+LL);

% Risposta al gradino
figure(4);

% Gradino w(t)=2*1(t)

T_simulation = 5;
[y_step,t_step] = step(WW*FF, T_simulation);
plot(t_step,y_step,'b');
grid on, zoom on, hold on;

% Vincolo sovraelongazione
patch([0,T_simulation,T_simulation,0],[WW*(1+S_100_spec),WW*(1+S_100_spec),WW+1,WW+1],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
ylim([0, WW+1]);

% Vincolo tempo di assestamento all'1%
LV = abs(evalfr(WW*FF,0)); % valore limite gradino: W*F(0)
patch([T_a1_spec,T_simulation,T_simulation,T_a1_spec],[LV*(1-0.01),LV*(1-0.01),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a1_spec,T_simulation,T_simulation,T_a1_spec],[LV*(1+0.01),LV*(1+0.01),LV+1,LV+1],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

Legend_step = ["Risposta al gradino"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);

%% Diagramma di bode sistema in anello chiuso
figure(5);
bode(FF)
grid on, zoom on

%% Check al gradino

% Risposta al gradino
figure(6);

WW_check = -2; % Gradino w(t)=-2*1(t)

tt = (0:1e-2:5);
T_simulation = 5;
[y_step,t_step] = step(WW_check*FF, T_simulation);
grid on, zoom on, hold on;
plot(t_step,y_step,'b');
plot([0,5],[WW_check,WW_check]);

legend('y_w', 'ww')
%% Check disturbo in uscita

% Funzione di sensitività
SS = 1/(1+LL);
figure(7);

% Testiamo il sistema su d(t)
tt = (0:1e-1:1e2);
dd = 0;
for kk = 1:4   
    dd=0.3*sin(0.025*kk*tt)+dd;
end
y_d = lsim(SS,dd,tt);
hold on, grid on, zoom on
plot(tt,dd,'m')
plot(tt,y_d,'b')
grid on
legend('dd','y_d')

%% Check disturbo di misura

% Funzione di sensitività complementare negata
FF_neg = -FF;
figure(8);

% Testiamo il sistema su n(t)
tt = (0:1e-8:1e-2)';
nn = 0;
for kk = 1:4   
    nn=0.2*sin(1e3*kk*tt)+nn;
end
y_n = lsim(FF_neg,nn,tt);
hold on, grid on, zoom on
plot(tt,nn,'m')
plot(tt,y_n,'b')
grid on
legend('nn','y_n')

%% Utilità per Simulink

x=[1,2,3,4];
DD_ampiezza=0.3;
DD_frequenza=0.025;
NN_ampiezza=0.2;
NN_frequenza=1e3;

%parte sistema lineare
R=R_d*RR_s*Polo;
[nuR,deR]=tfdata(R);
numR=nuR{1};
denR=deR{1};
[nuG,deG]=tfdata(GG);
numG=nuG{1};
denG=deG{1};

%parte sistema non lineare
x01=theta_e;
x02=0;
ye=x01;
u_e=2.17;

%% File Simulink punto 4 del progetto

%open('Progetto_a1_simulink_lineare.slx');

%% File Simulink punto 5 del progetto

%open('Progetto_a1_simulink_non_lineare.slx');