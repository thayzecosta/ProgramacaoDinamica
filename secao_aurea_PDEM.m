% secao aurea estocastica

function [x_min,cf] = secao_aurea_PDEM(a,b,R,pars,mes,P,W, C_atual, C_seguinte,R_vetor,b1)
 
%a=0;                            % start of interval
%b=2;                            % end of interval
epsilon=0.01;                      % accuracy value
iter= 10000;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);       % golden proportion coefficient, around 0.618
k=0;                             % number of iterations

x1=a+(1-tau)*(b-a);              % computing x values
x2=a+tau*(b-a); 

f_x1 = custo_PDEM(pars,R, mes, x1, C_atual,C_seguinte,P,W,R_vetor,b1);
f_x2 = custo_PDEM(pars,R, mes, x2, C_atual,C_seguinte,P,W,R_vetor,b1);

while ((abs(b-a)>epsilon) && (k<iter))
    k=k+1;
    if(f_x1<=f_x2)
        b=x2;
        x2=x1;
        x1=a+(1-tau)*(b-a);
        f_x1 = custo_PDEM(pars,R, mes, x1, C_atual,C_seguinte,P,W,R_vetor,b1);
        f_x2 = custo_PDEM(pars,R, mes, x2, C_atual,C_seguinte,P,W,R_vetor,b1);   
    else
        a=x1;
        x1=x2;
        x2=a+tau*(b-a);
        f_x1 = custo_PDEM(pars,R, mes, x1, C_atual,C_seguinte,P,W,R_vetor,b1);
        f_x2 = custo_PDEM(pars,R, mes, x2, C_atual,C_seguinte,P,W,R_vetor,b1);         
    end 
end

% chooses minimum point
if(f_x1<=f_x2)
    x_min=x1;
    cf = custo_PDEM(pars,R, mes, x_min, C_atual,C_seguinte,P,W,R_vetor,b1);
else
    x_min=x2;
    cf = custo_PDEM(pars,R, mes, x_min, C_atual,C_seguinte,P,W,R_vetor,b1);
end


