function [geracao, CP, R ,v,q_decisao,vert_turb,produtividade] = simulador_geral(pars, regras_decisao,algo)

%R_seguinte = pars.Rmax;        %nivel inicial do reservatorio
% horizonte
mes_inicial = pars.mes_inicial;
mes_final=pars.mes_final;

nW=pars.num_disc_w;

if algo==4
    H=pars.H;
    previsao =pars.prev; %1:MLT, 2:PAR1
end
    beta=[];
%============================================================================
serie_sintetica=pars.sintetica; % colocar 1(sim) ou 0 (nao)

%============================================================================

if algo ==4 %(MCP)
    V=[]; 
    for a=1:pars.anos
        if serie_sintetica==0
            V = [V pars.VAZ(a,1:12)];
        else
            V = [V pars.VAZ_sint(a,1:12)];
        end
        beta = [beta pars.beta];
    end
    for a=1:pars.H/12
        beta=[beta pars.beta];
    end
    VAZ=pars.VAZ;    
    media_x   = mean(VAZ);          %media

    if previsao==2
        desvpad_x = std(VAZ);           %desvio padrao
        sk = skewness(VAZ);
        anos=length(VAZ(:,1));
        % correlacao
        for i=1:12
            if i==1
                c1 = (VAZ(1:anos-1,12) - media_x(12)*ones(anos-1,1))'*(VAZ(2:anos,1) - media_x(1)*ones(anos-1,1))*(1/(anos-1));
                c2 = desvpad_x(12)*desvpad_x(1);
                rho_x(1) = c1/c2;
            else
                c1 = (VAZ(:,i-1) - media_x(i-1)*ones(anos,1))'*(VAZ(:,i) - media_x(i)*ones(anos,1))*(1/anos);
                c2 = desvpad_x(i-1)*desvpad_x(i);
                rho_x(i) = c1/c2;
            end
    
             a = 1 + sk(i)^2/2;
             b = sk(i)^2 + sk(i)^4/4;
             phi(i) = (a + sqrt(b))^(1/3) + (a - sqrt(b))^(1/3) -1;
             rho_y(i) = (log(rho_x(i)*(phi(i)-1)+1))/log(phi(i));
             desvio_y(i) = sqrt(log(phi(i)));
             media_y(i) = 0.5*log(desvpad_x(i)^2/(phi(i)^2 - phi(i)));
             tau(i) = media_x(i) - sqrt(desvpad_x(i)^2/(phi(i) -1));
        end
        
    end
    %sd=std(pars.VAZ(1:85,:));
    %m_log = log(  pars.W_medio./(sqrt(1+(sd.^2)./(pars.W_medio.^2)))  ) ;
    %s_log = sqrt(log( 1 + (sd.^2)./(pars.W_medio.^2 )  ));
   % valor_esperado = exp( media_y + (desvio_y.^2)/2)+tau
end




%============================================================================

%anos = pars.ano_final;
%custo = 0;
%d = pars.D(1);
grafico=zeros(1,pars.num_disc_w+1);

R=zeros((pars.ano_final-pars.ano_inicial)*12,1);         % trajetoria do resevatorio
v=zeros((pars.ano_final-pars.ano_inicial)*12,1);         % vertimento
q_decisao=zeros((pars.ano_final-pars.ano_inicial)*12,1); % turbinagem otima
CP=zeros((pars.ano_final-pars.ano_inicial)*12,1);        % custo de geracao
geracao=zeros((pars.ano_final-pars.ano_inicial)*12,1);   % geracao hidreletrica
if algo==3 
    P=P_dependente(pars);
end
%Discretizacao dos estados R
R_disc = linspace(pars.xmin,pars.xmax,pars.num_disc_x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if serie_sintetica==1
    vazao_anterior = pars.VAZ_sint(pars.ano_inicial,pars.mes_inicial-1);
else
    vazao_anterior = pars.VAZ(pars.ano_inicial,pars.mes_inicial-1); 
end    
barra = waitbar(0,'Simulando...');    

for t=pars.ano_inicial:pars.ano_final 
    
    waitbar((t-pars.ano_inicial)/(pars.ano_final));
    for mes=1:12
        if algo==3
            Want = P(2:pars.num_disc_w+1,1,mes);
        end
        if (t==pars.ano_inicial && mes<mes_inicial) || (t==pars.ano_final && mes>mes_final)

        else
            % Define vazao W_{t+1}
            if serie_sintetica==1
                vazao = pars.VAZ_sint(t,mes);
            else
                vazao = pars.VAZ(t,mes);
            end
            % Define estagio t
            estagio = (t-1)*12+mes;
            disp(['estagio: ', num2str(estagio)]);
            % inicializa vertimento turbinavel eliminado
            vert_turb(estagio)=0;
            % Define armazenamento R_t
            if estagio - ((pars.ano_inicial-1)*12+mes_inicial)==0
                R1 = pars.xmax;
                R(estagio-1) = R1;
            else    
                R1 = R(estagio-1);
            end
            
            if algo==1 || algo==2           
                if R1> (pars.xmax-0.1)
                    q_decisao(estagio) = regras_decisao(end,mes);
                    
                elseif R1<=pars.xmin
                    q_decisao(estagio) = regras_decisao(1,mes);
                else
                    q_decisao(estagio) = interp1( R_disc, regras_decisao(:,mes),R1 ,'linear');  
                end
            elseif algo==3
                %===========================
                % para PDEM
                
                for j=1:nW+1
                    if j==1
                        %delta=(Want(2)-Want(1))/2;
                        if (vazao_anterior<=Want(1))
                            pos_Want=1;
                            grafico(1,j)=grafico(1,j)+1;
                            %achou_j = 1;
                            break
                        end
                    else
                        if j==nW+1% && vazao_anterior >= Want(j-1)
                            pos_Want=j;
                            grafico(1,j)=grafico(1,j)+1;
                            break
                        else
                            %delta1 = (Want(j)-Want(j-1))/2; delta2=(Want(j+1)-Want(j))/2;
                            if (Want(j-1) <= vazao_anterior)&&(vazao_anterior <= Want(j))
                                pos_Want=j;
                                grafico(1,j)=grafico(1,j)+1;
                                %achou_j = 1;
                                break
                            end
                        end
                    end
                end
                for k=1:pars.num_disc_x-1
                    achou=0;
                    % Interpolacao para obter a decis?o de turbinagem
                    if (R_disc(k) <= R1)&&(R1 <= R_disc(k+1))
                        achou =1;
                        
                        if pos_Want>1 && pos_Want<nW
                            % Parte 1
                            y1 = regras_decisao(k,pos_Want,mes);
                            y2 = regras_decisao(k+1,pos_Want,mes);
                            x1 = R_disc(k);
                            x2 = R_disc(k+1);
                            %x_decisao = y1 + ((y2 - y1)*(R1- x1))/(x2 -x1);
                            A = (x2-R1)/(x2-x1);
                            x_decisao2= A*y1 + (1-A)*y2;
                            %achou =1;
                            %pos_Want
                            % Parte 2
                            y1 = regras_decisao(k,pos_Want-1,mes);
                            y2 = regras_decisao(k+1,pos_Want-1,mes);
                            %x1 = R_disc(k);
                            %x2 = R_disc(k+1);
                            %x_decisao = y1 + ((y2 - y1)*(R1- x1))/(x2 -x1);
                            A = (x2-R1)/(x2-x1);
                            x_decisao1 = A*y1 + (1-A)*y2;
                            
                            q_decisao(estagio)= interp1([Want(pos_Want-1);Want(pos_Want)],[x_decisao1;x_decisao2],vazao_anterior,'linear');
                        else
                            if pos_Want==1
                                
                                y1 = regras_decisao(k,pos_Want,mes);
                                y2 = regras_decisao(k+1,pos_Want,mes);
                                x1 = R_disc(k);
                                x2 = R_disc(k+1);
                                %x_decisao = y1 + ((y2 - y1)*(R1- x1))/(x2 -x1);
                                A = (x2-R1)/(x2-x1);
                                q_decisao(estagio)= A*y1 + (1-A)*y2;
                            else
                                
                                y1 = regras_decisao(k,pos_Want-1,mes);
                                y2 = regras_decisao(k+1,pos_Want-1,mes);
                                x1 = R_disc(k);
                                x2 = R_disc(k+1);
                                %x_decisao = y1 + ((y2 - y1)*(R1- x1))/(x2 -x1);
                                A = (x2-R1)/(x2-x1);
                                q_decisao(estagio)= A*y1 + (1-A)*y2;
                            end
                            %break
                        end
                        
                        
                        if achou
                            break
                        end
                    end
                end
               
                %==============================
            elseif algo==4
                anos_prev = H/12;
                V_prev=[];
                if previsao==2
                    vazao_ant=V(estagio-1);
                    mesp=mes;
                    for k=1:pars.H
                        mes_ant = mesp-1;
                        if mes_ant==0
                            mes_ant=12;
                        end
                        V_prev_y(k) = media_y(mesp) + rho_y(mesp)*(desvio_y(mesp)/desvio_y(mes_ant))*(log(vazao_ant-tau(mes_ant)) - media_y(mes_ant)) + desvio_y(mesp)^2*(1-rho_y(mesp)^2)/2;
                        V_prev(k) = exp(V_prev_y(k))+tau(mesp);
                        if ~(isreal(V_prev(k)))
                            V_prev(k) = media_x(mesp);
                            complexo = 1;
                        end
                            
                        %V_prev(k) = media_x(mesp) + rho_x(mesp)*(desvpad_x(mesp)/desvpad_x(mes_ant))*(vazao_ant - media_x(mes_ant)) ;
                        mesp=mesp+1;
                        if mesp==13
                            mesp=1;
                        end
                        vazao_ant=V_prev(k);
                    end
%                     if anos_prev>1
%                         for an=2:anos_prev
%                             % disp('aqui')
%                             if mes==1
%                                 V_prev=[V_prev  1.1*pars.W_medio];
%                                 %V_prev=[V_prev  pars.W_esperado];
%                             else
%                                 V_prev=[V_prev pars.W_medio(mes:12) 1.1*pars.W_medio(1:mes-1)];
%                                 %V_prev=[V_prev pars.W_esperado(mes:12) pars.W_esperado(1:mes-1)];
%                             end
%                         end
%                     end
                    
                elseif previsao==1

                    for an=1:anos_prev
                        if mes==1
                            V_prev= [V_prev pars.W_medio];
                        else
                            V_prev= [V_prev pars.W_medio(mes:12) pars.W_medio(1:mes-1)];
                        end
                    end
                end
                
                [ q_decisao(estagio), exitflag(estagio)]= MCP_opt_2019(R1,pars,H,V_prev,beta(estagio:estagio+H-1));
                q_decisao(estagio)=min(pars.qmax,q_decisao(estagio));
            end
            
            % C�lculo do armazenamento no est�gio t+1
            R(estagio) = R1 + pars.beta(mes)*(vazao - q_decisao(estagio));
            
            % C�lculo do vertimento e atualiza��o do armazenamento: 
            if R(estagio)>pars.xmax
                v(estagio) = (R(estagio)-pars.xmax)/pars.beta(mes);
                R(estagio) = pars.xmax;
                
            elseif R(estagio)<=pars.xmin
                q_decisao(estagio) = q_decisao(estagio) - (pars.xmin-R(estagio))/pars.beta(mes);
                R(estagio) = R1 + pars.beta(mes)*(vazao - q_decisao(estagio));
                %R(estagio)=pars.Rmin;
                v(estagio) = 0;
            else
                v(estagio)=0;
            end    
            
            % eliminando vertimento turbinavel
            if pars.eliminar_vert_turb ==1 && v(estagio)>0 && q_decisao(estagio)<pars.qmax
                diferenca = pars.qmax - q_decisao(estagio);
                if diferenca>v(estagio)
                    diferenca = v(estagio);
                    %v(estagio)=0;
                end
                vert_turb(estagio)=diferenca;
                q_decisao(estagio) = q_decisao(estagio) + vert_turb(estagio);
                v(estagio) = v(estagio)-vert_turb(estagio);
                % Balanco hidrico
                R(estagio) = R1 + pars.beta(mes)*(vazao - q_decisao(estagio)-v(estagio));                
            end
             
            %Calculo do custo
            %CP(estagio) = custo_presente((R1+R(estagio))/2,x_decisao(estagio)+v(estagio),pars,d);
            %CP(estagio) = pars.fc((R1+R(estagio))/2,x_decisao(estagio),v(estagio),mes);
            CP(estagio) = pars.fc((R1),q_decisao(estagio),v(estagio),mes);
            
            %CP(estagio) = custo_presente(R1,x_decisao(estagio),v(estagio),pars,d);
            %custo = custo + CP(estagio);

            %Calculo da geracao
            %cm = cota_montante(pars,(R1+R(estagio))/2);
            %cm = cota_montante(pars,R1);
            %cm = pars.fm(R1);
            
            %cj = cota_jusante(pars,x_decisao(estagio) + v(estagio));
            
            %altura_liquida = cm - cj;
            %geracao(estagio) = pars.k*altura_liquida*x_decisao(estagio);        
            %CP(estagio) = pars.c*(d - geracao(estagio))^2;
            geracao(estagio) = pars.fg((R1),q_decisao(estagio),v(estagio));
            %geracao(estagio) = pars.fg(R1,x_decisao(estagio),v(estagio));
            
            vazao_anterior = vazao;
        end
    end

end
close(barra);
produtividade = geracao./q_decisao;
if algo ==3
    figure(11)
    bar(grafico)
end