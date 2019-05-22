% Programacao Dinamica Deterministica

function regras_decisao2 = PDD_xcontinuo2(pars)

%Discretizacao do volume
x=(linspace(pars.xmin,pars.xmax,pars.num_disc_x))';

% Inicializacao das matrizes de custo futuro e presente, e de turbinagem
% futura e presente
geracao_final=0;
V_fut   = geracao_final*ones(pars.num_disc_x,12);
V_atual = geracao_final*ones(pars.num_disc_x,12);
q_fut   = zeros(pars.num_disc_x,12);
q_atual = zeros(pars.num_disc_x,12); 

figure(1);hold;



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
        [q_atual(i,mes),V_atual(i)] = secao_aurea(lb,ub,x(i),pars,x_custo,mes,t);
        CP(t,i)=V_atual(i);
    end
    V_fut = V_atual;
    
    % teste de convergencia das turbinagens:
    if t>24 && mes==1        
        if max(max(abs(q_fut-q_atual))) < pars.tol
            break; 
        end
    end
    q_fut(:,mes) = q_atual(:,mes);  
    custo_medio = mean((CP)./t);
    
    figure(1);
    %shg
    %subplot(1,2,1)
    %plot(custo_medio);  drawnow; 
    %hold on;
    
    %figure(2)
    shg
    %subplot(1,2,2)
    plot(V_atual);drawnow;
    hold on
    title('Custo Futuro')
    xlabel('Armazenamento')
    ylabel('Custo Futuro')
    xlim([1 pars.num_disc_x])
end

figure(10)
for i=1:12
    sb=subplot(3,4,i);
    %plot(regras_decisao_PDD_cont{i}(:,2),'r'); hold on
    plot(q_fut(:,i),'r'); hold on
    title(['Mes: ', int2str(i)]);
    legend(sb,'PDD')
   % regras_decisao2(:,i)=regras_decisao{i}(:,2);
end

regras_decisao2 = q_fut;
