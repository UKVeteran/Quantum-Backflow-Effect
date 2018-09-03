% Calculates the backflow current using equation (38.245) of my notes. The
% input is h, q0, N0 which must be the same as was used in rpen to
% calculate the eigenfunction. The eigenfunction is read in through the
% file bf_xxx.txt. We also need to define here m, c, hbar and T which 
% define epsilon the relativistic backflow parameter. The current involves 
% four integrals. If we use the function f in the calculation of z1->z4 
% we get the maximum possible value of the backflow for this value of 
% epsilon. If we use fa instead we get the backflow for thatb approximation
% to the eigenfunction. This code uses an Airy function for epsilon=1.0;
h=4;
q0=20.0;
N0=4000;
q=q0*sqrt(h);
np=N0*h;
p=linspace(0.0,q,np);
m=1.0;
c=1.0;
hbar=1.0;
T=2.0;
t=linspace(0.0,T,1000);
nt=length(t);
eps=sqrt(4*hbar/(m*c^2*T))
eps2=eps^2.0;
r=p./(m*c*eps);
eta=sqrt(1+eps^2*r.^2);
for ip=1:np
    p1(ip)=sqrt((eta(ip)+1.0)/eta(ip));
    p2(ip)=sqrt((eta(ip)-1.0)/eta(ip));
end
fileID = fopen('bf_airy2.txt','r');
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
za=-1.130*(1.03*r+1.08).^(0.796);
fa1=2.14*airy(0,za)./(1.25*r+0.76).^(2.0/3.0);

V2=conj(V).*V;
fa2=conj(fa1).*fa1;
zzg=trapz(r,V2);
zzh=trapz(r,fa2);
V=V/sqrt(zzg);
fa1=fa1/sqrt(zzh);
for ip=1:np
    fa(ip)=exp(2*1i*eta(ip)/eps^2)*fa1(ip)/sqrt(m*c);
    f(ip)=exp(2*1i*eta(ip)/eps^2)*V(ip)/sqrt(m*c);
end
for it=1:nt
    for ip=1:np
        z1(ip)=p1(ip)*conj(fa(ip))*exp(4.0*1i*eta(ip)*t(it)/(eps2*T));
        z2(ip)=p2(ip)*fa(ip)*exp(-4.0*1i*eta(ip)*t(it)/(eps2*T));
        z3(ip)=p1(ip)*fa(ip)*exp(-4.0*1i*eta(ip)*t(it)/(eps2*T));
        z4(ip)=p2(ip)*conj(fa(ip))*exp(4.0*1i*eta(ip)*t(it)/(eps2*T));
    end
    j1(it)=trapz(r,z1);
    j2(it)=trapz(r,z2);
    j3(it)=trapz(r,z3);
    j4(it)=trapz(r,z4);
    j(it)=(j1(it)*j2(it)+j3(it)*j4(it))/(pi*eps2);
end
zz=trapz(t,j)
plot(r,V,r,fa1);

    
        
    

    
