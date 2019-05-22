% Testa se duas matrizes coluna são quase iguais
% M,N: matrizes colunas a serem testadas
% tol

function T=teste(M,N,tol)
K = abs(M-N);
tam = length(K);
T=false;
if all(K(1:tam)<tol*ones(tam,1))
    T = true;
end
