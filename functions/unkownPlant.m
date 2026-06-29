function k = unkownPlant()
addpath("data")
addpath("functions/")
x_white = load("whites.mat").X;
y_white = load("whites.mat").Y;
x_corr = load("correl.mat").X;
y_corr = load("correl.mat").Y;

getOrder(x_white, y_white) % Gráfico para questão a)

% Ordem da planta é de 30, assim o vetor tem 31 casas
w_white = adaptiveFilter(x_white, y_white, 31)
w_corr = adaptiveFilter(x_corr, y_corr, 31)
end