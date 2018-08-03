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
y = y'; % passt sonst nicht k zusammen
k = linspace(0,q,N+1);

plot(k,y,'b')
axis([0 10 -0.05 0.14])
hold on;

% abstand zw. z und y minimieren ueber
% variation von vorfaktor f fuer k
% und add fuer k
n = 1;
m = 1;

fspan = 0.985:0.0005:0.988;
gspan = 0.30:0.001:0.33;
d = zeros(length(fspan), length(gspan));

for g = gspan
    n = 1

    for f = fspan
        % approx mit erfi
        z = -imag(erfi((1+i)/sqrt(2)*(f*k + g)))+1;

        % normieren
        z = z/norm(z);
    
        % abstand
        d(n,m) = norm(z-y);
    
        % index 1
        n = n+1
    end
    
    % index 2
    m = m+1
end

[mindrow, minn] = min(d);
minn = minn(1);
[mind, minm] = min(mindrow);

f = fspan(minn) % gehoert zu min
g = gspan(minm)
d(minn,minm)

% damit neu berechnen
z = -imag(erfi((1+i)/sqrt(2)*(f*k + g)))+1;
z = z/norm(z);
plot(k,z,'r')

% auch graph fuer d zeichnen
figure(2)
contourf(gspan,fspan,d)
colorbar