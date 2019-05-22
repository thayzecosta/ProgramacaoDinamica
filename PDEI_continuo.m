%==========================================================================
function regras_decisao2=PDEI_continuo(pars)
%==========================================================================
% Funcao que retorna as tabelas de decisao da PDEI. 
%==========================================================================

for mes=1:12
    P{mes}=get_prob_ind_log_sk(pars.VAZ,mes,pars.num_disc_w );
    P{mes}=P{mes}';
end

%Discretizacao do volume
x=(linspace(pars.xmin,pars.xmax,pars.num_disc_x))';

% Inicializacao das matrizes de custo futuro e presente, e de turbinagem
% futura e presente
geracao_final=0;
V_fut   = geracao_final*ones(pars.num_disc_x,12);
V_atual = geracao_final*ones(pars.num_disc_x,12);
q_fut   = zeros(pars.num_disc_x,12);
q_atual = zeros(pars.num_disc_x,12); 

figure(5); hold;

% programacao dinamica (sentido backward)
t=0;
mes=1;
while true
    t=t+1;
    disp(['iteracao: ', int2str(t)]);
    %mes_fut=mes;
    mes=mes-1;
    if mes==0
        mes=12;
    end   
    % varredura dos volumes
    for i=1:pars.num_disc_x       
        if i==1
            lb=pars.qmin;
        else
            lb=q_atual(i-1,mes);
        end
        ub=pars.qmax;
        x_custo= [x,V_fut];
        [q_atual(i,mes),V_atual(i)] = secao_aurea_est(lb,ub,x(i),pars,x_custo,mes,P{mes});
    end
    V_fut = V_atual;
    
    % teste de convergencia das turbinagens:
    if t>24 && mes==1        
        if max(max(abs(q_fut-q_atual))) < pars.tol
            break; 
        end
    end
    q_fut(:,mes) = q_atual(:,mes);    
    
    shg
    plot(V_atual);drawnow
    title('Custo Futuro')
    xlabel('Armazenamento')
    ylabel('Custo Futuro')
    xlim([1 pars.num_disc_x])
end

figure(5)
for i=1:12
    sb=subplot(3,4,i);
    %plot(regras_decisao_PDD_cont{i}(:,2),'r'); hold on
    plot(q_fut(:,i),'r'); hold on
    title(['Mes: ', int2str(i)]);
    legend(sb,'PDEI')
   % regras_decisao2(:,i)=regras_decisao{i}(:,2);
end

regras_decisao2 = q_fut;

