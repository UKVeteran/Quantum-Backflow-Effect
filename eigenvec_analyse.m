close all;
clear all;

% Eigenvec als y einfuettern
% Achse als Vec k

% mit folgenden Werten gerechnet:
N = 10000; % groesse des Rechenintervalls [-N, N]
q = 50; % nicht zu gross sonst wird Oszillation von U zu stark, kommt auch auf N an
% iter = 1000
% Backflow-Wert: -0.02973327208522

load eigenvec.mat % Eigenvec einlesen -> y
k = linspace(-q,q,2*N+1)';

% y auffuellen
y = [zeros(N,1);y];

subplot(211)
plot(k,y,'k')
axis([-10 10 -0.05 0.14])

% Ortsdarstellung
subplot(212)
[x,Y] = symfft(k,y,-1);
plot(x,real(Y),'k')
hold on
plot(x,imag(Y),'k--')
hold off
axis([-15 15 -0.025 0.025])

% Polarplot
figure(2)
theta = angle(Y);
rho = abs(Y);
polar(theta, rho, 'k');

% Strom
% hbar/(2m) t in Einheiten [m]
% 2m/hbar j in Einh [1/m]

% fuer verschiede Zeiten
tmax = 5;
tN = 5000;
tspan = linspace(-tmax,tmax,tN);
m = 1;

% neu einlesen (nur pos. gebraucht)
load eigenvec.mat % Eigenvec einlesen -> y
k = linspace(0,q,N+1)';

for t=tspan,
    U = exp(-i*k.^2*t);
    invU = conj(U);
    
    % 4 Integrale, als Summen, Schrittweite k0/N
    I1 = sum(invU.*conj(y));
    I2 = sum(U.*k.*y);
    I3 = sum(U.*y);
    I4 = sum(invU.*k.*conj(y));
    j(m) = q/N*(I1*I2+I3*I4)/2/pi * 2; % Faktor 1/2 weg
    
    % aufintegrierter Strom
    jint(m) = sum(j);
    
    m = m+1;
end

figure(3)
plot(tspan,j,'k')
hold on
plot(tspan,zeros(length(tspan),1),'k--')
axis([-5 5 -0.28 1.02])

% Schraffieren
% Schraffur anbringen
% ersten Nulldurchgang nach Mitte suchen
for m = round(tN/2):tN
    if j(m-1)*j(m)<0,
        break
    end
end

% endpunkte
m1 = 2*round(tN/2)-m+1;
m2 = m;

% Randpunkte korrigieren (fuer Schraffierung)
j(m1)=0;j(m2)=0;

% Farbe
obj = patch(tspan(m1:m2),j(m1:m2),'w');
hatch(obj,45,'k','-',3,0.5)

% integrierter Strom
plot(tspan, jint*tmax/tN + (1-jint(length(jint))*tmax/tN)/2, 'r--')

% Normierungs-Check
sum(j)*tmax/tN % stimmt so halbwegs? < 1 nicht zu klein

% Backflow
sum(j(m1:m2))*tmax/tN