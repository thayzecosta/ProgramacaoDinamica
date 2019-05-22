% Arquivo de entrada para a programacao dinamica

function pars=dados_serra_da_mesa()
disp('Serra da Mesa');

% dados da simulacao
pars.mes_inicial=5;
pars.mes_final=3;
pars.ano_inicial = 2;
pars.ano_final = 63;

% horizonte de planejamento (numero de meses)
anos = pars.ano_final-pars.ano_inicial;
T = anos*12;

% eliminar vertimentos turbinaveis? 0 - nao, 1 - sim
pars.eliminar_vert_turb = 1;

% horizonte do MCP
pars.H = 4*12;

% tipo de previsao para o MCP: 1 - MLT, 2 - PAR1
pars.prev = 1; 

%serie sintetica? 0 - nao, 1 - sim
pars.sintetica = 0;

% Interpolacao de custo futuro: 1-linear; 2-spline
pars.interp=1;

% tolerancia de convergencia da programacao dinamica
pars.tol=0.01;

% Discretizacoes
pars.num_disc_x=100;
pars.num_disc_w=10;

% Custo [$/MW^2]
pars.c = 0.0002; % 106.6;

% quantidade de dias por mes
dias = [31 28 31 30 31 30 31 31 30 31 30 31];
% numero de segundos por mes
N_sec = 86400*dias;
% beta
pars.beta=N_sec/10^6;
%pars.beta=2.628*ones(1,12);

% Demanda
pars.D = 1275*ones(1,T);

% Geracao final
pars.geracao_final=0;

% Polinomio cota volume a0+a1*x+a2*x^2+a3*x^3+a4*x^4 %ok
pars.a0 =  3.914048E+02;
pars.a1 =  2.772160E-03;
pars.a2 = -4.357250E-08;
pars.a3 =  2.903040E-13;
pars.a4 = 0;

% Polinomio cota jusante b0+b1*x+b2*x^2+b3*x^3+b4*x^4
pars.b0 =  3.327979E+02;    %cota media= 359.65; 
pars.b1 =  1.342970E-03;
pars.b2 =  8.819558E-08;
pars.b3 = -1.627669E-11;
pars.b4 = 0;

% Volume maximo e minimo[hm^3]
pars.xmax = 54400;
pars.xmin = 11150;
pars.x_0 = pars.xmax; 

% Turbinagem maxima e minima [m3/s]
pars.qmax = 1074.6; %descontando teif e ip
pars.qmin = 0;

% Constante k [MW/(m3/s).m]
pars.k = 0.009124;

% % Leitura das vazoes
VAZ= xlsread('vazoes_serra_da_mesa.xlsx');
pars.VAZ=VAZ(1:85,:);
pars.W_medio = mean(pars.VAZ);

% funcoes
%cota montante
pars.fm = @(x) pars.a0 + pars.a1*x + pars.a2*x.^2 + pars.a3*x.^3 + pars.a4*x.^4; 
%cota jusante
pars.fj = @(q,v) pars.b0 + pars.b1*(q+v) + pars.b2*(q+v).^2 + pars.b3*(q+v).^3 + pars.b4*(q+v).^4;
%geracao
pars.fg = @(x,q,v) pars.k*(pars.fm(x) - pars.fj(q,v)).*q;
%custo presente
pars.fc = @(x,q,v,mes) ((pars.beta(mes)*10^6)/3600)*pars.c*(pars.D(1) - pars.fg(x,q,v)).^2  ; % Calcula o custo do estagio
