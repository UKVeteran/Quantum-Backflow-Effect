% Test ob Norm von Eigenvec konvergiert
close all;
clear all;

%load meth_2_eigenvec3.mat % genaustes Ergebnis -> Var y
load eigenvec.mat
k = linspace(0,50,10001); % passende x-achse

max = length(y); % Länge

% sin(x^2)/x ist L^2
y2 = sin(k.^2)./k; % mit k-Werten
y2(1) = 0; % NaN sonst
% normieren
y2 = y2 / norm(y2);

norms = zeros(max, 1);
norms2 = zeros(max, 1);

for i = 1:1:max
    norms(i) = norm(y(1:i))^2;
    norms2(i) = norm(y2(1:i))^2;
end

figure(1)
hold on

plot(k, norms,'k-')
plot(k, norms2,'k--')

% zum Vergleich passenden log und sqrt zeichnen
%plot(x, sqrt(x)/80,'--')
%plot(x, log(x)/8,'-.')

%legend('Eigenvec','sin(k^2)/k','Sqrt','Log');
