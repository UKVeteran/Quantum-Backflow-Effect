% Calculates something
h=10;
q0=60.0;
N0=4000;
q=q0*sqrt(h);
np=N0*h;
p=linspace(0.0,q,np);
m=1.0;
c=1.0;
hbar=1.0;
T=1.0;
t=linspace(0,T,1000);
nt=length(t);
eps=sqrt(4*hbar/(m*c^2*T))
eps2=eps^2.0;
r=p./(m*c*eps);
eta=sqrt(1+eps^2*r.^2);
for ip=1:np
    p1(ip)=sqrt((eta(ip)+1.0)/eta(ip));
    p2(ip)=sqrt((eta(ip)-1.0)/eta(ip));
end
fileID = fopen('bf.txt','r');
for ip=1:np
[r(ip)]=fscanf(fileID,'%f, %f');
end
for ip=1:np
[V(ip)]=fscanf(fileID,'%f, %f');
end
fclose(fileID);
for ip=1:np
    f(ip)=exp(2*1i*eta(ip)/eps^2)*V(ip)/sqrt(m*c);
end
pre=np/(q*pi*eps2);
for it=1:nt
    for ip=1:np
        z1(ip)=p1(ip)*conj(f(ip))*exp(4.0*1i*eta(ip)*t(it)/(eps2*T));
        z2(ip)=p2(ip)*f(ip)*exp(-4.0*1i*eta(ip)*t(it)/(eps2*T));
    end
    j1(it)=trapz(r,z1);
    j2(it)=trapz(r,z2);
    j(it)=2*pre*real(j1(it)*j2(it));
end
zz=trapz(t,j)
plot(t,j);

    
        
    

    
