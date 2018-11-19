% Code by Ehsan Mirzakhalili mirzakh@umich.edu
% https://doi.org/10.1371/journal.pone.0201302
function y=Run(t,Par)

K1=0.13;
K2=1.05;
K3=0.943;
K4=0.145;
k2_1=0.21;
k4_1=0.029;
k42=k4_1/k2_1;

Ca=0.1;


Par(25)=Par(3)*Ca^2/(Ca^2+Par(4)^2);

CaE=220;

Par(24)=Par(1)*Ca^2/(Ca^2+Par(2)^2)/(CaE-Ca);

IP3=0;
phi1=k2_1*(k42*K2*K1+K4*IP3)*Ca/(K4*K2*(K1+IP3));
phi2=k2_1*(IP3+k42*K3)/(K3+IP3);
Y=phi1/(phi1+phi2);

x0=zeros(9,1);

x0(1)=Ca;
x0(2)=CaE;
x0(7,1)=1;
x0(10,1)=-64.4;
x0(9,1)=Y;

opts = odeset('Abstol',1e-3,'Reltol',1e-3);
[~,x] = ode23(@(T,x) Model(T,x,Par),t,x0,opts);

Ca=x(:,1);

Rmin=4.963093110035229;
Rmax=32.282101799695049;
R0=Rmax*Ca(1)^1.7/(Ca(1)^1.7+2.5^1.7)/100*(Rmax-Rmin)+Rmin;
R=Rmax*Ca.^1.7./(Ca.^1.7+2.5^1.7)/100*(Rmax-Rmin)+Rmin;
y=(R-R0)/R0*100;

end
