function regras_decisao=PDEM_continuo(pars)

P_matriz = P_dependente(pars);

%Discretizacao do volume
x=(linspace(pars.xmin,pars.xmax,pars.num_disc_x))';

W=zeros(1,pars.num_disc_w,12);
P=zeros(pars.num_disc_w,pars.num_disc_w,12);

figure(5)
regras_decisao = zeros(pars.num_disc_x,pars.num_disc_w,12);
% contem as decisoes otimas para cada discretizacao de W anterior e cada
% mes

for mes=1:12
    W(:,:,mes) = P_matriz(1, 2:pars.num_disc_w+1, mes);
    P(:,:,mes) = P_matriz(2:(pars.num_disc_w+1),(2:pars.num_disc_w+1),mes);
end    

V_fut = pars.geracao_final*ones(pars.num_disc_x,pars.num_disc_w,12);
V_atual = pars.geracao_final*ones(pars.num_disc_x,pars.num_disc_w,12);
q_fut = zeros(pars.num_disc_x,pars.num_disc_w,12);
q_atual = zeros(pars.num_disc_x,pars.num_disc_w,12);   

%parpool
%% DISCRETIZACAO DO TEMPO  
t=0;
while true
    t=t+1;
    disp(['iteracao: ', int2str(t)]);
    mes=mes-1;
    if mes==0
        mes=12;
    end   
    
    V_at = V_atual;
    
    %% DISCRETIZACAO DA VAZAO ATUAL
    for b=1: pars.num_disc_w
        %% DISCRETIZACAO DO VOLUME
        %for a=1:pars.num_disc_x    
        parfor a=1:pars.num_disc_x        
            %% DISCRETIZACAO DA VAZAO FUTURA
            %if a==1
                lb = pars.qmin;
            %else
                %lb = q_atual(a-1,b,mes);
            %end
                ub = pars.qmax;
                [x_min,custo_min] =secao_aurea_PDEM(lb,ub,x(a),pars,mes,P,W, V_at, V_fut,x,b);
                V_atual(a,b,mes) = custo_min;
                q_atual(a,b,mes) = x_min;                    
        end        
    end
    
    if t>24&& mes==1
        disp('Teste')
        if max(max(max(abs(q_fut - q_atual))))<pars.tol
            break;
        end    
        max(max(max(abs(q_fut - q_atual))))
    end
    figure(5)
    if mes==1 && t>12
        V_fut(:,:,:)=V_atual(:,:,:);
        q_fut(:,:,:)=q_atual(:,:,:);
        regras_decisao(:,:,:)= q_atual(:,:,:);
        plot(regras_decisao(:,1,1));
        hold on
    end
    
end
