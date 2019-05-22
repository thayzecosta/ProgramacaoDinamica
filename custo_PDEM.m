function custo = custo_PDEM(pars,R, mes, x, C_atual,C_seguinte,P,W,R_vetor,b1)
%parpool
%logic2=0;
custo1=zeros(1,pars.num_disc_w);


R_seguinte = R + pars.beta(mes)*(W(1,:,mes) - x);
% vert = (max((x_seguinte-pars.xmax)/pars.beta(mes),0))';
% deft = (max((pars.xmin-x_seguinte)/pars.beta(mes),0))';
% x_seguinte = x_seguinte -vert +deft;
% %R_seguinte = min(pars.xmax,max(pars.xmin,R_seguinte));
% 
% 
% %gh_x = pars.fg(x,q,0); 
% %gh_min = pars.fg(x,,0);
% 
% 
% %cp_min = pars.fc(x*ones(pars.num_disc_w,1),q*ones(pars.num_disc_w,1)-deft,vert,mes);
% 
% %cp = pars.c*((pars.beta(mes)*10^6)/3600)*(pars.D(1) - (gh_min-gh_delta))^2 ;
% 
% %cc = ; %penalidade por ficar abaixo do xmin
% 
% %delta_x = ( pars.xmin - x_seguinte(i))/pars.beta(mes);
%         % x que leva o reservatorio no minimo
% %        q_min = q - delta_x;
%         
%         gh_x = pars.fg(x,q,vert); 
%         
%         
%         gh_min = pars.fg(x,q-deft,vert);
%         gh = gh_min - abs(gh_x-gh_min);
%         
%         cp = pars.c*((pars.beta(mes)*10^6)/3600)*(pars.D(1) - gh).^2 ;
%         
% 
% 
% 
% %     hjus_min = pars.b0 + pars.b1*(x_min_pen) + pars.b2*(x_min_pen)^2 + pars.b3*(x_min_pen)^3 + pars.b4*(x_min_pen)^4;%%
% %     gh_min = pars.k*(hmon - hjus_min - pars.xi)*x_min_pen;
% %     gh_delta = abs(gh-gh_min);
% %     cp = pars.c*((pars.beta(mes)*10^6)/3600)*(pars.D(1) - (gh_min-gh_delta))^2 ;
% 
% % for d=1:pars.num_disc_w
% %     if mes~=12
% %         cf=interp1(R_vetor, C_atual(:,d,mes+1), R_seguinte(d),'linear');
% %     else
% %         cf=interp1(R_vetor, C_seguinte(:,d,1), R_seguinte(d),'linear');
% %     end
% %     custo1(d) = P(b1,d,mes)*(cp(d)+cf);                 
% % end 
% 
% for d=1:pars.num_disc_w
%     if mes~=12
%         cf=interp1(R_vetor, C_atual(:,d,mes+1), x_seguinte(d),'linear');
%     else
%         cf=interp1(R_vetor, C_seguinte(:,d,1), x_seguinte(d),'linear');
%     end
%     custo1(d) = P(b1,d,mes)*(cp(d)+cf);                 
% end 
% 
% 
% 
% custo=sum(custo1);  




%% DISCRETIZA��O DA VAZ�O FUTURA
                for d=1: pars.num_disc_w
                    %% CALCULO DO VOLUME SEGUINTE
                    logic1=0;
                    logic2=0;
                    %% VERIFICACAO E REALOCA??O DE VOLUMES, TURBINAGENS E VERTIMENTOS
                    if (R_seguinte(d)<=pars.xmax)&&(R_seguinte(d)>=pars.xmin)
                        vert=0;
                        %cc=0;                   
                    elseif (R_seguinte(d)>pars.xmax)
                        vert=(R_seguinte(d)-pars.xmax)/pars.beta(mes);
                        R_seguinte(d) = pars.xmax;
                        logic1=0;
                        %cc=0;               
                    else 
                        vert=0;
                        delta_x = (pars.xmin - R_seguinte(d))/pars.beta(mes);%%%
                        x_min_pen = x - delta_x; %%%
                        logic2=1;
                        R_seguinte(d) = pars.xmin;
                    end    
                    
                    %% CALCULOS     
                    % VOLUME Rmed(volume tabelado, R_seguinte)
                    Rmed = R;%+R_seguinte)/2;%(R(b)+R_seguinte(d))/2;%
                
                    % ALTURA DE MONTANTE
                    hmon = pars.a0 + pars.a1*Rmed + pars.a2*Rmed^2 + pars.a3*Rmed^3 + pars.a4*Rmed^4;
                     
                    % ALTURA DE JUSANTE
                    hjus = pars.b0 + pars.b1*(x+vert) + pars.b2*(x+vert)^2 + pars.b3*(x+vert)^3 + pars.b4*(x+vert)^4;
                    gh = pars.k*(hmon - hjus)*x;
                    cp = pars.c*((pars.beta(mes)*10^6)/3600)*(pars.D(1) - gh)^2 ;
                    %hjus = pars.b0 + pars.b1*(x(c)) + pars.b2*(x(c))^2 + pars.b3*(x(c))^3 + pars.b4*(x(c))^4;
                    if logic2
                        %x_min_novo = x_min_pen - delta_x;
                        hjus_min = pars.b0 + pars.b1*(x_min_pen) + pars.b2*(x_min_pen)^2 + pars.b3*(x_min_pen)^3 + pars.b4*(x_min_pen)^4;%%
                        gh_min = pars.k*(hmon - hjus_min)*x_min_pen;
                        gh_delta = abs(gh-gh_min);
                        cp = pars.c*((pars.beta(mes)*10^6)/3600)*(pars.D(1) - (gh_min-gh_delta))^2 ;                      
                    end
                    % GERACAO
                    %gh = pars.k*(hmon - hjus - pars.xi)*x_min_novo;
                    %if logic2
                        
                    %    gh_pen2 = pars.k*(hmon - hjus_xmin - pars.xi)*x_min_pen; %%
                   % end                
%                     if gh>pars.potnom
%                         gh=pars.potnom;
%                     end                
                    %% CUSTO(mes, demanda, termicas, hidraulicas)
                   % cp = pars.c*((pars.D(1) - gh)*(pars.beta(mes)*10^6)/3600)^2 ;
                    %if logic2
                    %    cp2 = (pars.c*((pars.D(1) - gh_pen2)*(pars.beta(mes)*10^6)/3600)^2);%%
                       % cc = abs(cp2-cp);
                   %     logic2=0;
                   % end
                    %% CALCULO DOS CUSTOS TOTAIS INTERPOLADOS
                    %if t==pars.T
                    %    cf=0;
                    %else
                        if logic1~=1 && logic2~=1
                            if pars.interp==1
                                if mes~=12
                                    cf=interp1(R_vetor, C_atual(:,d,mes+1), R_seguinte(d),'linear');
                                else
                                    cf=interp1(R_vetor, C_seguinte(:,d,1), R_seguinte(d),'linear');
                                end
                            else
                                if mes~=12
                                    cf=interp1(R_vetor, C_atual(:,d,mes+1), R_seguinte(d),'spline');
                                else
                                    cf=interp1(R_vetor, C_seguinte(:,d,1), R_seguinte(d),'spline');
                                end
                            end
                        elseif logic1==1
                            if mes~=12
                                cf= C_atual(end,d,mes+1);
                            else
                                cf=C_seguinte(end,d,1);
                            end
                        elseif logic2==1
                            if mes~=12
                                cf= C_atual(1,d,mes+1);
                            else
                                cf=C_seguinte(1,d,1);
                            end
                            
                        end
                    %end
                    custo1(d) = P(b1,d,mes)*(cp+cf); %%%%%%
                    
                end
                custo=sum(custo1);     
                
