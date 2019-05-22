% ------------------------GOLDEN SECTION METHOD----------------------------
function [x_min,cf]=secao_aurea(a,b,x,pars,x_custo,mes,t)

%a=0;                            % start of interval
%b=2;                            % end of interval
epsilon=0.001;                   % accuracy value
iter= 10000;                     % maximum number of iterations
alpha=double((sqrt(5)-1)/2);     % golden proportion coefficient, around 0.618
k=0;                             % number of iterations

lambda=a+(1-alpha)*(b-a);        % computing x values
mu=a+alpha*(b-a);

% Interpolacoes de custo futuro
f_x1 = custo_total_PDD(pars,x,lambda,x_custo,mes,t);
f_x2 = custo_total_PDD(pars,x,mu,x_custo,mes,t);
% R_seguinte1 = R + pars.beta(mes)*(pars.W_medio(mes)- lambda) ;
% f_x1 = custo_presente(R,lambda,pars,d) + custo_futuro(R_Custos,R_seguinte1,pars);   % computing values in x points

% R_seguinte2 = R + pars.beta(mes)*(pars.W_medio(mes)- mu) ;
% f_x2=custo_presente(R,mu,pars,d)+ custo_futuro(R_Custos,R_seguinte2,pars);

while ((abs(b-a)>epsilon) && (k<iter))
    k=k+1;
    if(f_x1<f_x2)
        b=mu;
        mu=lambda;
        lambda=a+(1-alpha)*(b-a);
        f_x1 = custo_total_PDD(pars,x,lambda,x_custo,mes,t);
        f_x2 = custo_total_PDD(pars,x,mu,x_custo,mes,t);
        
%         R_seguinte1 = R + pars.beta(mes)*(pars.W_medio(mes)- lambda) ;
%         f_x1=custo_presente(R,lambda,pars,d)+custo_futuro(R_Custos,R_seguinte1,pars);
%         R_seguinte2 = R + pars.beta(mes)*(pars.W_medio(mes)- mu) ;
%         f_x2=custo_presente(R,mu,pars,d)+custo_futuro(R_Custos,R_seguinte2,pars);
    else
        a=lambda;
        lambda=mu;
        mu=a+alpha*(b-a);
        f_x1 = custo_total_PDD(pars,x,lambda,x_custo,mes,t);
        f_x2 = custo_total_PDD(pars,x,mu,x_custo,mes,t);
%         R_seguinte1 = R + pars.beta(mes)*(pars.W_medio(mes)- lambda) ;
%         f_x1=custo_presente(R,lambda,pars,d)+custo_futuro(R_Custos,R_seguinte1,pars);
%         R_seguinte2 = R + pars.beta(mes)*(pars.W_medio(mes)- mu) ;
%         f_x2=custo_presente(R,mu,pars,d)+custo_futuro(R_Custos,R_seguinte2,pars);
    end    
end

% chooses minimum point
if(f_x1<=f_x2)
    x_min=lambda;
    cf = custo_total_PDD(pars,x,lambda,x_custo,mes,t);
    %cf=custo_presente(R,x_min,pars,d)+custo_futuro(R_Custos,R_seguinte1,pars);
else
    x_min=mu;
    cf = custo_total_PDD(pars,x,mu,x_custo,mes,t);
    %cf=custo_presente(R,x_min,pars,d)+custo_futuro(R_Custos,R_seguinte2,pars);
end
