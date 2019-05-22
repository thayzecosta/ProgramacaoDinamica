% Calculo do custo estoc�stico de acordo com a decisao

function c = custo(pars,x,q,x_custo,mes,P)

x_seguinte = x*ones(pars.num_disc_w,1) + pars.beta(mes)*(P(1,:)' - q*ones(pars.num_disc_w,1));
c = 0;
vert=0;

for i=1:pars.num_disc_w    
    if x_seguinte(i) >= pars.xmax % reservatorio acima do maximo
        vert = (x_seguinte(i)-pars.xmax)/pars.beta(mes);
        cf = x_custo(pars.num_disc_x,2);
        cp = pars.fc(x,q,vert,mes);
        c1 = P(2,i)*(cf+cp);
        
        
    elseif x_seguinte(i) <= pars.xmin % reservatorio abaixo do minimo
        %custo futuro
        cf = x_custo(1,2);
        %custo da correcao
        delta_x = ( pars.xmin - x_seguinte(i))/pars.beta(mes);
        % x que leva o reservatorio no minimo
        q_min = q - delta_x;
        
        gh_x = pars.fg(x,q,0); gh_min = pars.fg(x,q_min,0);
        gh = gh_min - abs(gh_x-gh_min);
        
        cp = pars.c*((pars.beta(mes)*10^6)/3600)*(pars.D(1) - gh)^2 ;
        
        % custo presente
        c1 =  P(2,i)*(cf+cp);                
    else
        cp = pars.fc(x,q,vert,mes);
        for k=1:pars.num_disc_x
            if (x_custo(k,1) <= x_seguinte(i))&&(x_seguinte(i) <= x_custo(k+1,1))
                if pars.interp==2
                    cf=spline(x_custo(k:k+1,1),x_custo(k:k+1,2),x_seguinte(i));
                else
                    cf = x_custo(k,2) +  (x_custo(k+1,2) - x_custo(k,2))*(x_seguinte(i) - x_custo(k,1))/(x_custo(k+1,1)- x_custo(k,1));
                end
                break
            end   
        end
        c1 = P(2,i)*(cf+cp);
    end
    c = c + c1;
end
