% secao aurea estocastica

function [x_min,cf] = secao_aurea_est(a,b,x,pars,x_custo,mes,P)

%a=0;                            % start of interval
%b=2;                            % end of interval
epsilon=0.01;                      % accuracy value
iter= 5000;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);       % golden proportion coefficient, around 0.618
k=0;                             % number of iterations

q1=a+(1-tau)*(b-a);              % computing x values
q2=a+tau*(b-a); 

% Calcular f(x)
f_x1 = custo(pars,x,q1,x_custo,mes,P);
f_x2 = custo(pars,x,q2,x_custo,mes,P);

while ((abs(b-a)>epsilon) && (k<iter))
    k=k+1;
    if(f_x1<f_x2)
        b=q2;
        q2=q1;
        q1=a+(1-tau)*(b-a);
        f_x1 = custo(pars,x,q1,x_custo,mes,P);        
        f_x2 = custo(pars,x,q2,x_custo,mes,P);   
    
    else
        a=q1;
        q1=q2;
        q2=a+tau*(b-a);
        f_x1 = custo(pars,x,q1,x_custo,mes,P);        
        f_x2 = custo(pars,x,q2,x_custo,mes,P);  
    end 
end

% chooses minimum point
if(f_x1<f_x2)
    x_min=q1;
    cf = custo(pars,x,x_min,x_custo,mes,P);
else
    x_min=q2;
    cf = custo(pars,x,x_min,x_custo,mes,P);
end


