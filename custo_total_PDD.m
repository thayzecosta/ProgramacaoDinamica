function c = custo_total_PDD(pars,x,q,x_custo,mes,t)

x_seguinte = x + pars.beta(mes)*(pars.W_medio(mes) - q);
cc = 0;

%vert=0;
%for i=1:pars.num_disc_W
if x_seguinte >= pars.xmax % reservatorio acima do maximo
    vert = (x_seguinte - pars.xmax)/pars.beta(mes);
    cp = pars.fc(x,q,vert,mes);
    %cp = custo_presente(R,x,vert,pars,d);
    cf = x_custo(pars.num_disc_x,2);
 
elseif x_seguinte <= pars.xmin % reservatorio abaixo do minimo
    %custo futuro
    cf = x_custo(1,2);
    
    
    %custo da correcao
    delta_x = ( pars.xmin - x_seguinte)/pars.beta(mes);%N_sec(mes)
    % x que leva o reservatorio no minimo
    x_min = q - delta_x;
     
    cc = abs(pars.fc(x,x_min,0,mes) - pars.fc(x,q,0,mes));
    %cc = abs(custo_presente(R,x_min,0,pars,d)-custo_presente(R,x,0,pars,d));
    %cp = custo_presente(R,x_min,0,pars,d);
    cp = pars.fc(x,x_min,0,mes);
else
    cp = pars.fc(x,q,0,mes);
    %cp = custo_presente(R,x,0,pars,d);
    if pars.interp==2
        cf=spline(x_custo(:,1),x_custo(:,2),x_seguinte);
        
    else
        cf = interp1(x_custo(:,1),x_custo(:,2),x_seguinte);
        
        %cf = R_Custos(k,2) +  (R_Custos(k+1,2) - R_Custos(k,2))*(R_seguinte - R_Custos(k,1))/(R_Custos(k+1,1)- R_Custos(k,1));
    end
end

c = ((t)*cf +(cc + cp))/(t+1);


