% Arquivo de entrada para a programacao dinamica

function pars=dados_furnas()
disp('Furnas')

% dados da simulacao
pars.mes_inicial=5;
pars.mes_final=4;
pars.ano_inicial = 2;
pars.ano_final = 79;


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
pars.D = 1312*ones(1,T);

% Geracao final
pars.geracao_final=0;

% Polinomio cota volume a0+a1*x+a2*x^2+a3*x^3+a4*x^4
pars.a0 =  7.352458E+02;
pars.a1 =  3.496580E-03;
pars.a2 = -1.974370E-07;
pars.a3 =  6.917049E-12;
pars.a4 = -9.773650E-17;

% Polinomio cota jusante b0+b1*x+b2*x^2+b3*x^3+b4*x^4
pars.b0 = 672.9; %0;% 6.716328E+02;    %cota media= 672.9; 
pars.b1 =1.017380E-03;
pars.b2 =-1.799719E-07;
pars.b3 =2.513280E-11;
pars.b4 = 0;

% Volume maximo e minimo[hm^3]
pars.xmax = 22950;
pars.xmin = 5733;
pars.x_0 = pars.xmax;

% Turbinagem maxima e minima [m3/s]
pars.qmax = 1231.60; %descontando teif e ip
pars.qmin = 0;%102;

% Constante k [MW/(m3/s).m]
pars.k = 0.008633;

% % Leitura das vazoes
VAZ = xlsread('vazoes_furnas.xlsx');
pars.VAZ = VAZ(1:85,:);
pars.W_medio = mean(pars.VAZ);

% funcoes
%cota montante
pars.fm = @(R) pars.a0 + pars.a1.*R + pars.a2*R.^2 + pars.a3*R.^3 + pars.a4*R.^4; 
%cota jusante
pars.fj = @(x,v) pars.b0 + pars.b1.*(x+v) + pars.b2*(x+v).^2 + pars.b3*(x+v).^3 + pars.b4*(x+v).^4;
%geracao
pars.fg = @(R,x,v) pars.k*(pars.fm(R) - pars.fj(x,v)).*x;
%custo presente
pars.fc = @(R,x,v,mes) pars.c*pars.beta(mes)*((10^6)/3600)*(pars.D(1) - pars.fg(R,x,v)).^2  ; % Calcula o custo do est?gio

