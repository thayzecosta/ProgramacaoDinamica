%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Planejamento da operacao energetica 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% Escolha da politica:
% algo = 1 : Programacao Dinamica Deterministica
% algo = 2 : Programacao Dinamica Estocastica Independente
% algo = 3 : Programacao Dinamica Estocastica Markoviana


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clear all;
tic  %inicia a contagem de tempo
algo =3;

% Os parâmetros do modelo são escolhidos nos seguinte arquivos, conforme a usina escolhida: 
pars=dados_furnas();
%pars = dados_sobradinho();
%pars=dados_emborcacao();
%pars=dados_serra_da_mesa();

%pars=dados_furnas_sint();
%pars=dados_sobradinho_sint();
%pars=dados_emborcacao_sint();
%pars=dados_serra_da_mesa_sint();

sintetica = pars.sintetica;

if algo ==1     %PDD
    regras_decisao_PDD_cont(:,:) = PDD_xcontinuo2(pars);
   
elseif algo ==2 %PDEI 
        regras_decisao_PDEI_cont(:,:) = PDEI_continuo(pars);
        
elseif algo ==3 %PDEM
     regras_decisao(:,:,:) = PDEM_continuo(pars);
     
%elseif algo ==4 %MCP 
    
else
    error('Opcao nao existente!');
end

disp('Simulando...')
if algo==1
    [geracao, CP , x ,v,q_decisao,vert_turb,produtividade]  = simulador_geral(pars, regras_decisao_PDD_cont,algo);
    disp('=====================================================================')

    disp('Simulacao de programacao dinamica deterministica ');
    disp(['Inicio: ano ', int2str(1930+pars.ano_inicial), ' mes: ', int2str(pars.mes_inicial)]);
    disp(['Final: ano ', int2str(1930 + pars.ano_final), ' mes: ', int2str(pars.mes_final)] );

    disp('======================')

elseif algo==2

    [geracao, CP , x ,v,q_decisao,vert_turb,produtividade] = simulador_geral(pars, regras_decisao_PDEI_cont,algo);
    disp('=====================================================================')

    disp('Simulacao de programacao dinamica estocastica independente ');
    disp(['Inicio: ano ', int2str(1930+pars.ano_inicial), ' mes: ', int2str(pars.mes_inicial)]);
    disp(['Final: ano ', int2str(1930 + pars.ano_final), ' mes: ', int2str(pars.mes_final)] );

    disp('======================')

elseif algo==3
    [geracao, CP , x ,v,q_decisao,vert_turb,produtividade] = simulador_geral(pars, regras_decisao,algo);


elseif algo==4
    [geracao, CP , x ,v,q_decisao,vert_turb,produtividade] = simulador_geral(pars,0,algo);
end

%===================================================================
a1 = (pars.ano_inicial-1)*12+pars.mes_inicial;
a2 = (pars.ano_final-1)*12+pars.mes_final;
VAZ=[];
if ~sintetica
    for t=1:length(pars.VAZ(:,1))
        VAZ = [VAZ, pars.VAZ(t,:)];
    end
else
    for t=1:length(pars.VAZ_sint(:,1))
        VAZ = [VAZ, pars.VAZ_sint(t,:)];
    end
end



m1 = mean(geracao(a1:a2));  % geracao media
m5 = std(geracao(a1:a2));
m2 = mean(CP(a1:a2));  % custo medio
m6 = std(CP(a1:a2));
m3 = mean(produtividade(a1:a2));
m4 = mean(v(a1:a2));


disp(['Custo medio: $ ', num2str(m2)]);
disp(['Custo desvio padrao: $ ', num2str(m6)]);
disp(['Geracao media: ', num2str(m1), ' MW medio']);
disp(['Geracao desvio padrao: ', num2str(m5), ' MW medio']);
disp(['Produtividade media: ', num2str(m3), ' MW/m^3/s']);
disp(['Vertimento medio: ', num2str(m4), ' m^3/s']);


figure(3); % trajetoria do reservatorio
subplot(2,1,1)
plot(x(a1:end),'g'); title('Trajetoria de armazenamento'); hold on

plot(pars.xmax*ones(a2-a1+1,1),'k--');
plot(pars.xmin*ones(a2-a1+1,1),'c--');
ylim([0,pars.xmax+1000]);
xlabel ('Estagio t')
ylabel ('Volume [hm^3]');
legend('Armazenamento PDD','Volume maximo','Volume minimo');

subplot(2,1,2)
plot(q_decisao(a1:a2),'b'); hold on
plot(VAZ(a1:a2),'r'); title('Turbinagem, Vazoes e Vertimento'); hold on
plot(v(a1:a2),'g');
plot(pars.qmax*ones(a2-a1+1,1),'k--');
plot(pars.qmin*ones(a2-a1+1,1),'c--');
xlabel ('Estagio t');
ylabel ('Vazao [m^3/s]');
legend('Turbinagem','Vazao afluente','Vertimento','Eng. maximo','Def. minima')

figure(4); % custos e gera??o
subplot(3,1,1)
plot(geracao(a1:a2)); title('Geracao'); hold on
plot((pars.D(1)*ones(1,a2-a1+1)-geracao(a1:a2)),'r'); hold on
plot(pars.D(1)*ones(a2-a1,1),'k--');
%ylim([0,pars.Rmax+1000]);
xlabel('Estagio t');
ylabel('Geracao [MW medio]');
legend('Geracao hidroeletrica','Geracao termoeletrica','Capacidade instalada');

subplot(3,1,2)
plot(CP(a1:a2)); title('Custo mensal')
xlabel ('Estagio t');
ylabel ('Custo [$]');

subplot(3,1,3)
plot(produtividade(a1:a2)); title('Produtividade')
xlabel ('Estagio t');
ylabel ('Produtividade [MW/m^3/s]');

mix_t = mean((pars.D(1)-geracao)/pars.D(1))
mix_h = 1-mix_t


toc  %finaliza a contagem de tempo
