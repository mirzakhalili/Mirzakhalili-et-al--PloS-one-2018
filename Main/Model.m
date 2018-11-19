% Code by Ehsan Mirzakhalili mirzakh@umich.edu
% https://doi.org/10.1371/journal.pone.0201302
function f = Model(t,x,Par)
Ca=x(1);
CaE=x(2);
p0=x(3);
p1=x(4);
p2=x(5);
O=x(6);
I=x(7);
IP3=x(8);
Y=x(9);
V=x(10);

K1=0.13;
K2=1.05;
K3=0.943;
K4=0.145;
K5=0.082;
k2_1=0.21;
k4_1=0.029;
k42=k4_1/k2_1;

RT=8.3145*293.1;
F=96485.3;

if t<10
    Stim=0;
elseif t<40
    Stim=1;
else
    Stim=0;
end
JPump=Par(3)*Ca^2/(Ca^2+Par(4)^2);
JSerca=Par(1)*Ca^2/(Ca^2+Par(2)^2);
JChannels=Par(6)*O*I*(Ca-2000);

phi1=k2_1*(k42*K2*K1+K4*IP3)*Ca/(K4*K2*(K1+IP3));
phi2=k2_1*(IP3+k42*K3)/(K3+IP3);
OIP3=(IP3*Ca*(1-Y)/(IP3+K1)/(Ca+K5)).^3;
%%

aV=1*26*(Ca-4/26)./(1-exp(-26*(Ca-4/26)));
bV=-10*26*(Ca-2/26)./(1-exp(0.5*26*(Ca-2/26)));
aCa=1.6/(1+exp(-0.072*(V-5)));
bCa=0.02*(V-1.31)/(exp((V-1.31)/5.36)-1);
mCa=aCa/(aCa+bCa);
JV=mCa^2*Par(8)*V*(-Ca*exp(-2*F*V/RT*1e-3)+2000)/(1-exp(-2*F*V/RT*1e-3));

%%
f(1,1)=(Par(7)*OIP3+Par(24))*(CaE-Ca)-JSerca-JPump-JChannels+JV+Par(25);
f(2,1)=Par(26)*(JSerca-(Par(7)*OIP3+Par(24))*(CaE-Ca));
f(3,1)=Par(9)*Stim/(Stim+Par(5))-Par(10)*p0;
f(4,1)=Par(11)*p0-Par(12)*p1;
f(5,1)=Par(13)*p0/(p0+Par(23))-Par(14)*p2;
f(6,1)=Par(15)*(1-O)*(p1)-Par(16)*O;
f(7,1)=Par(17)*(1-I)-Par(18)*I*p2;
f(8,1)=Par(19).*p1/(1+(Par(22)*p2)^4)*Ca.^2/(Ca.^2+Par(21).^2)-Par(20)*IP3;
f(9,1)=phi1*(1-Y)-phi2*Y;
f(10,1)=aV*(20-V)-bV*(V-Par(27));