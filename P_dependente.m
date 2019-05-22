%clear all

function P=P_dependente(pars)
%====================================================================
%% ESCOLHA A USINA:
%====================================================================
%pars = dados_furnas();
%pars= dados_sobradinho();
%pars= dados_serra_da_mesa();
%pars= dados_foz_do_areia();
%pars= dados_emborcacao();

div_log=0; %0: partes iguais
VAZ = pars.VAZ;

%==========================================================================
%% INICIALIZACAO DAS VARIAVEIS
P = zeros(pars.num_disc_w+1,pars.num_disc_w+1,12);
%matriz_lim=zeros(pars.num_disc_W,2,12);
n_anos = length(VAZ);
rho_x = zeros(1,12);
phi=zeros(1,12);
rho_y = zeros(1,12);
desvio_y= zeros(1,12);
%min1=zeros(1,12);
%max1=zeros(1,12);
tau = zeros(1,12);
media_y=zeros(1,12);
desvio_y_cond=zeros(1,12);
n_faixas = pars.num_disc_w;

%% C�LCULOS B�SICOS E GR�FICOS DA S�RIE ORIGINAL

media = mean(VAZ);
desvio = std(VAZ);
sk = skewness(VAZ);
curtose = kurtosis(VAZ);

% figure(11)
% x=(1:1:12);
% e=errorbar(x,media, desvio,'MarkerEdgeColor','red','MarkerSize',6);
% e.Marker = '*'; e.LineWidth = 1.5;
% %title('M�dia e desvio padr�o')
% xlim([0.5 12.5]);
% xticks(1:1:12)
% xticklabels({'Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez'})
% xlabel('M�s')
% ylabel('Vaz�o [m^3/s]')


figure(1)
subplot(2,2,1)
x=(1:1:12);
e=errorbar(x,media, desvio,'MarkerEdgeColor','red','MarkerSize',6);
e.Marker = '*'; e.LineWidth = 1.5;
title('M�dia e desvio padr�o')
xlim([0.5 12.5])

subplot(2,2,2)
bar(sk); hold on; title('Assimetria');
xlim([0.5 12.5])

subplot(2,2,3)
bar(curtose); hold on; title('Curtose');
xlim([0.5 12.5])

% Correla��o:
for mes=1:12
    if mes==1
        A2 = VAZ(1:n_anos-1,12) - media(12)*ones(n_anos-1,1);
        A1 = VAZ(2:n_anos,1) - media(1)*ones(n_anos-1,1);
        c1 = (1/(n_anos-1)) * A1'*A2;
        c0 = desvio(1)*desvio(12);
        rho_x(1) = c1/c0;
    else
        A1 = VAZ(:,mes) - media(mes)*ones(n_anos,1);
        A2 = VAZ(:,mes-1) - media(mes-1)*ones(n_anos,1);
        c1 = (1/(n_anos)) * A1'*A2;
        c0 = desvio(mes)*desvio(mes-1);
        rho_x(mes) = c1/c0;
    end
end

subplot(2,2,4)
bar(rho_x); hold on; title('Correla��o com o m�s anterior');
xlim([0.5 12.5])

%% ESTAT�STICAS DA S�RIE TRANSFORMADA

% dados da s�rie y
for mes=1:12
    a = 1 + sk(mes)^2/2;
    b = sk(mes)^2 + sk(mes)^4/4;
    phi(mes) = (a + sqrt(b))^(1/3) + (a - sqrt(b))^(1/3) -1;
    rho_y(mes) = (log(rho_x(mes)*(phi(mes)-1)+1))/log(phi(mes));
    desvio_y(mes) = sqrt(log(phi(mes)));
    media_y(mes) = 0.5*log(desvio(mes)^2/(phi(mes)^2 - phi(mes)));
    tau(mes) = media(mes) - sqrt(desvio(mes)^2/(phi(mes) -1));
    
    
%     %% PLOT DAS FUN��ES DENSIDADE DE PROBABILIDADE INDEPENDENTE
    figure(2)
    subplot(4,3,mes)

    m1=media_y(mes); s1=desvio_y(mes);

    log_min(mes) = m1 - 3*s1;
    log_max(mes) = m1 + 3*s1;
%     
    novo_min(mes) = exp(log_min(mes)) + tau(mes);
    novo_max(mes) = exp(log_max(mes)) + tau(mes);
%     
%     if mes==7
%         novo_min(mes)=min(VAZ(:,7));
%         tau(7) = novo_min(7)- exp(log_min(7)) ;
%     end


    x = (novo_min(mes):1:novo_max(mes));
    y = lognpdf(x-tau(mes),media_y(mes),desvio_y(mes));
    pl=plot(x,y,'b');    
    xlabel('Vaz�o (W_{t-1})','FontSize',8)
    ylabel('Densidade','FontSize',9)
    pl.LineWidth = 1.5;
    title(['M�s: ', int2str(mes)]);
    hold on;
    yL = get(gca,'YLim');
   % line([media(mes) media(mes)],yL,'Color', 'r', 'linewidth',1.5);
    
 end

for mes=1:12
    %===============================================================================
    %% DIVIS�O DAS FAIXAS
    %===============================================================================
    if div_log==1
        % Divis�o ,log
        
        matriz_lim(1,1,mes) = log_min(mes);
        delta = (log_max(mes)- log_min(mes))/(pars.num_disc_W);
        
        for i=1:pars.num_disc_W
            matriz_lim(i,2,mes) = matriz_lim(i,1,mes) + delta;
            if i<pars.num_disc_W
                matriz_lim(i+1,1,mes) = matriz_lim(i,2,mes);
            end
        end
        matriz_lim(:,:,mes) = exp(matriz_lim(:,:,mes))+tau(mes);
        Wmes(:,mes) = 0.5*(matriz_lim(:,1,mes) + matriz_lim(:,2,mes)); % m�dia dos limites
        
    %============================================================================
    %Divis�o comum (partes iguais)

    else 
        max1(mes) = max(VAZ(:,mes));
        %max1(mes) = max(setdiff(VAZ(:,mes),max3(mes)));
        min1(mes) = min(VAZ(:,mes));
        %min1(mes) = min(setdiff(VAZ(:,mes), min2(mes)));
        
        matriz_lim(1,1,mes) = min1(mes);
        delta = (max1(mes)-min1(mes))/(pars.num_disc_w);
        for i=1:pars.num_disc_w
            matriz_lim(i,2,mes) = matriz_lim(i,1,mes) + delta;
            if i<pars.num_disc_w
                matriz_lim(i+1,1,mes) = matriz_lim(i,2,mes);
            end
        end
        Wmes(:,mes) = 0.5*(matriz_lim(:,1,mes) + matriz_lim(:,2,mes)); % m�dia dos limites
        
    end
    
    
end


    % Preenchimento de P da markoviana
    for mes=1:12
        mes_ant = mes-1;
        if mes_ant==0
            mes_ant=12;
        end
        
        P(1,2:pars.num_disc_w+1,mes) = Wmes(:,mes)';
        P(2:pars.num_disc_w+1,1,mes) = Wmes(:,mes_ant);
        mes_ant = mes-1;
        if mes_ant==0
            mes_ant=12;
        end
        
        for k=1:pars.num_disc_w %n_faixas de discretiza��es do W anterior
            desvio_y_cond(mes) = desvio_y(mes)*sqrt(1-rho_y(mes)^2);
             teste=Wmes(k,mes_ant)-tau(mes_ant);
            if teste<1
                teste=40;
                disp(['teste',num2str(mes)]);
            end
            y_ant(k) = log(teste);
            media_y_cond(k,mes) = media_y(mes) + rho_y(mes)*(desvio_y(mes)/desvio_y(mes_ant))*(y_ant(k)-media_y(mes_ant));
            
            
            %===================================================================================================
%             % plot das fdp's
            figure(2)
            subplot(4,3,mes)
            %x = (log_min(mes):1:log_max(mes))';
            x = (novo_min(mes):10:novo_max(mes))';
            y = lognpdf(x-tau(mes),media_y_cond(k,mes),desvio_y_cond(mes));
            
            plot(x,y, 'Color',[0.5 0.7 0.7])            
            title(['M�s: ', int2str(mes)]);
            hold on;
            
            yL = get(gca,'YLim');
            line([media(mes) media(mes)],yL,'Color','r');
            if k == pars.num_disc_w
                hold off
            end
            z_surf{mes}(k,:)= y;
            
            %===================================================================================================
            
            p = logncdf([-Inf matriz_lim(1,2,mes)-tau(mes)], media_y_cond(k,mes), desvio_y_cond(mes));
            P(k+1,2,mes) = p(2)-p(1);
            
            %p = normcdf([log(matriz_lim(n_faixas-1,2)-tau) Inf], media_y, desvio_y);
            p = logncdf([matriz_lim(pars.num_disc_w,1,mes)-tau(mes) Inf], media_y_cond(k,mes), desvio_y_cond(mes));
            %P(pars.num_disc_W,2,pars.num_disc_W) = p(2)-p(1);
            P(k+1,pars.num_disc_w+1,mes) = p(2)-p(1);
            
            for j=3:(pars.num_disc_w)%-1
                %p = normcdf([log(matriz_lim(j,1)-tau)  log(matriz_lim(j,2)-tau)], media, desvio);
                p = logncdf([matriz_lim(j-1,1,mes)-tau(mes) matriz_lim(j-1,2,mes)-tau(mes)], media_y_cond(k,mes), desvio_y_cond(mes));
                P(k+1,j,mes) = p(2)-p(1);
                %P(j,2,pars.num_disc_W) = p(2)-p(1);
            end
        end
        
% %         
%             figure(3)
%             lg = length(z_surf{mes}(1,:));
%             delta1 = (novo_max(mes) - novo_min(mes))/(lg-1);
%             delta2 = (novo_max(mes_ant) - novo_min(mes_ant))/(lg-1);
%             subplot(3,4,mes)
%             x_surf = (novo_min(mes):delta1:novo_max(mes))';
%             y_surf = Wmes(:,mes);
%             %y_surf = (novo_min(mes_ant):delta2:novo_max(mes_ant))';
%         
%             surf(x_surf,y_surf,z_surf{mes}(:,:));%,P(2:pars.num_disc_W+1,2:pars.num_disc_W+1,mes))
%             title(['FDP cond. no m�s: ', int2str(mes)]);
%             xlim([0 novo_max(mes_ant)]); xlabel('W ant.');
%             ylim([0 novo_max(mes)]); ylabel('W atual')
% %           xlim(P(1,2:pars.num_disc_W+1))
%  %         xlabel ('Aflu�ncia')
%   %        ylabel ('Aflu�ncia anterior');
%         
%         
%         figure(4)
%         for z=1:pars.num_disc_W
%             subplot(3,4,mes)
%             plot(P(1,2:pars.num_disc_W+1,mes), P(z+1,2:pars.num_disc_W+1,mes));%,P(2:pars.num_disc_W+1,2:pars.num_disc_W+1,mes))
%             %  xlim(P(2:pars.num_disc_W+1,1,mes))
%             %  ylim(P(1,2:pars.num_disc_W+1,mes))
%             xlim([novo_min(mes) novo_max(mes)]); xlabel('mes atual')
%             hold on
%         end
%         
%         figure(5)
%         subplot(3,4,mes)
%         surf(P(2:pars.num_disc_W+1,1,mes),P(1,2:pars.num_disc_W+1,mes),P(2:pars.num_disc_W+1,2:pars.num_disc_W+1,mes));%,P(2:pars.num_disc_W+1,2:pars.num_disc_W+1,mes))
%         title(['Prob. cond. no m�s: ', int2str(mes)]);
%         xlim([0 novo_max(mes_ant)]); xlabel('mes ant.');
%         ylim([0 novo_max(mes)]); ylabel('mes atual')
     end
