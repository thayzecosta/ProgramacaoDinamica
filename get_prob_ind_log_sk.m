% Calcula as probabilidades independentes de ocorrencia de vaz?es
% 
% primeira coluna: vaz?es
% segunda coluna: probabilidades
% distribui??o normal

function P=get_prob_ind_log_sk(VAZ,mes,n_faixas)


P = zeros(n_faixas,2);

%tau = min(VAZ(:,mes));

min1 = min(VAZ(:,mes));
max1 = max(VAZ(:,mes));
% %[media desvio] = lognfit(VAZ(:,mes));

media = mean(VAZ(:,mes));
desvio = std(VAZ(:,mes));
sk = skewness(VAZ(:,mes)); %coeficiente de assimetria

a = 1 + sk^2/2;
b = sk^2 + sk^4/4;
phi = (a + sqrt(b))^(1/3) + (a - sqrt(b))^(1/3) -1;

media_y = 0.5*log(desvio^2/(phi^2 - phi));
desvio_y = sqrt(log(phi));
tau = media - sqrt(desvio^2/(phi -1));

    %% PLOT DAS FUNCOES DENSIDADE DE PROBABILIDADE INDEPENDENTE
    figure(2)
    subplot(3,4,mes)
    m1=media_y; s1=desvio_y;
    
     %log_min = log(min1 - tau);
     %log_max = log(max1 - tau);
     log_min = m1 - 2*s1;
     log_max = m1 + 2*s1;

    novo_min = exp(log_min) + tau;
    novo_max = exp(log_max) + tau;
    
    x = (novo_min:1:novo_max);
    y = lognpdf(x-tau,media_y,desvio_y);
    pl=plot(x,y,'k');
    pl.LineWidth = 1.5;
    title(['Mes: ', int2str(mes)]);
    hold on;
    yL = get(gca,'YLim');
    line([media media],yL,'Color','r');

    max1=novo_max; min1=novo_min;


delta = (max1-min1)/(n_faixas);
matriz_lim(1,1) = min1;


for i=1:n_faixas
    matriz_lim(i,2) = matriz_lim(i,1) + delta; 
    if i<n_faixas 
        matriz_lim(i+1,1) = matriz_lim(i,2);
    end
end    

P(:,1) = (matriz_lim(:,1) + matriz_lim(:,2))/2;

%p = normcdf([-Inf log(matriz_lim(1,2)-tau)], media_y, desvio_y);
p = logncdf([-Inf matriz_lim(1,2)-tau], media_y, desvio_y);
P(1,2) = p(2)-p(1);

%p = normcdf([log(matriz_lim(n_faixas-1,2)-tau) Inf], media_y, desvio_y);
p = logncdf([matriz_lim(n_faixas-1,2)-tau Inf], media_y, desvio_y);
P(n_faixas,2) = p(2)-p(1);

for j=2:(n_faixas-1)
    %p = normcdf([log(matriz_lim(j,1)-tau)  log(matriz_lim(j,2)-tau)], media, desvio);
    p = logncdf([matriz_lim(j,1)-tau matriz_lim(j,2)-tau], media_y, desvio_y);
    P(j,2) = p(2)-p(1);
end

