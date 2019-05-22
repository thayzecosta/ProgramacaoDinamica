% Arquivo de entrada para a programacao dinamica

function pars=dados_sobradinho()
disp('Sobradinho');

% dados da simulacao
pars.mes_inicial=5;
pars.mes_final=3;
pars.ano_inicial = 2;
pars.ano_final = 77;

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
pars.D = 1200*ones(1,T);
%pars.D = 1050*ones(1,T);

% Geracao final
pars.geracao_final=0;

% Geracao final
pars.geracao_final=0;

% Polinomio cota volume a0+a1*x+a2*x^2+a3*x^3+a4*x^4 %ok
pars.a0 =  3.741790E+02;
pars.a1 =  1.396690E-03;
pars.a2 = -5.351590E-08;
pars.a3 =  1.155989E-12;
pars.a4 = -9.545989E-18;

% Polinomio cota jusante b0+b1*x+b2*x^2+b3*x^3+b4*x^4
pars.b0 =  3.596538E+02;    %cota media= 359.65; 
pars.b1 =  1.964010E-03;
pars.b2 = -2.968730E-07;
pars.b3 =  2.508280E-11;
pars.b4 = -7.702299E-16;

% Volume maximo e minimo[hm^3]
pars.xmax = 34116;
pars.xmin = 5447;
pars.x_0 = pars.xmax; %pars.Rmin + (pars.Rmax - pars.Rmin)/2;

% Turbinagem maxima e minima [m3/s]
pars.qmax = 3860; %descontando teif e ip
pars.qmin = 0;%102;

% Constante k [MW/(m3/s).m]
pars.k = 0.009023;

% Leitura das vazoes
VAZ= xlsread('vazoes_sobradinho.xlsx');
pars.VAZ= VAZ(1:85,:);
pars.W_medio = mean(pars.VAZ);

% funcoes
%cota montante
pars.fm = @(x) pars.a0 + pars.a1*x + pars.a2*x.^2 + pars.a3*x.^3 + pars.a4*x.^4; 
%cota jusante
pars.fj = @(q,v) pars.b0 + pars.b1*(q+v) + pars.b2*(q+v).^2 + pars.b3*(q+v).^3 + pars.b4*(q+v).^4;
%geracao
pars.fg = @(x,q,v) pars.k*(pars.fm(x) - pars.fj(q,v)).*q;
%custo presente
pars.fc = @(x,q,v,mes) pars.c*((pars.beta(mes)*10^6)/3600)*(pars.D(1) - pars.fg(x,q,v)).^2  ; % Calcula o custo do est?gio
