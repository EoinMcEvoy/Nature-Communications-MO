clear all
close all

% Parameters 
r0     = 7.1e-6;   % reference cell radius 
Ag     = 63.33e-12;  % area of membrane adhered between two cells 
Ec     = 6e3;  %effective stiffness of cortical layer %original 3e3
siga   = 100; %active cortical stress
h0     = 0.6e-6;  % thickness of cortical layer 
RT     = 8.3145*310;   % qas constant x temperature
F      = 96485.3329;  % Faradays constant 
Cm     = 1e-2;  % unit membrane capicitance 
z      = -1.5;   % valence impermeable solutes
X      = 2e-13;  % impermeable solutes  
Nae    = 145; % external sodium conc
Ke     = 5;  % external potassium conc
Cle    = 110; % external cloride conc
Lpm    = 7e-12;   % membrane water permeability
gNa    = 0.1;    %  sodium channel conductivity
gK0    = 0.1;    %  minimum potassium channel conductivity 
gCl    = 2;      %  chloride channel conductivity
p      = 0.005;  %  Na/K pump current
w      = 1;      % GJ conductivity parameter
Lpg    = 1e-11;  % gap junction water permeability factor 
sigc  = 75;  % threshold stress of MS channel 
sigs  = 600; % Pa % saturating stress of MS channel 
betaK = 9*gK0/(sigs-sigc); % MS channel ion permeability factor 

% applied loading
ds1 = 1e3;  % Pa  % growth-associated stress
tl =5200; dt1=600;% time stress starts to be applied and time that it takes to reach peak, respectively

%% Transient solution

% Set initial conditions
tspan = [0 10000];
options = odeset('RelTol',2.5e-12);

% estimate initial conditions from physiological ranges
Xo0(1) = 10;
Xo0(2) = 10;
Xo0(3) = 140;
Xo0(4) = 140;
Xo0(5) = 10;
Xo0(6) = 10;
Xo0(7) = (4/3)*pi*(1.1*r0)^3;
Xo0(8) = (4/3)*pi*(1.1*r0)^3;


% ODE solver to compute volume and ions
[tt, xx]=ode23s(@dyn_ion,tspan,Xo0,options,F, Cm, z, X, Nae, Ke, Cle, RT, p, Lpm, gNa, gK0, gCl...
    ,r0,siga,Ec,h0,sigc,sigs,betaK,ds1,tl,dt1, w, Lpg, Ag);

V1  = xx(:,7);
V2  = xx(:,8);
Nai1 = xx(:,1);
Nai2 = xx(:,2);
Ki1  = xx(:,3);
Ki2  = xx(:,4);
Cli1 = xx(:,5);
Cli2 = xx(:,6);

disp('ODE analysis end')

%% Post operations

% Function for time dependent loading profile on cell 1
for i = 1:length(tt)
    
    if (tt(i)>tl)&&(tt(i)<=tl+dt1)
        sig_g1(i,1) = ds1*(tt(i)-tl)/dt1; 
    elseif (tt(i)>tl+dt1)
        sig_g1(i,1) = ds1;
    else
        sig_g1(i,1) = 0;
    end
 
end
sig_g2 = 0;

% Update cell radii
r1   = real((3.*V1./(4*3.14)).^(1/3));
sig1 = real((Ec/(r0)).*(r1-r0) + siga);
r2   = real((3.*V2./(4*3.14)).^(1/3));
sig2 = real((Ec/(r0)).*(r2-r0) + siga);

dP1  = 2.*sig1.*h0./r1 + sig_g1;
dP2  = 2.*sig2.*h0./r2 + sig_g2; 

% MS channel
for i = 1:length(tt)
    if sig1(i) < sigc   %1
        gK1(i,1)  = gK0;
    elseif sig1(i) < sigs 
        gK1(i,1)  = gK0 + betaK*(sig1(i)-sigc);
    elseif sig1(i) >= sigs
        gK1(i,1)  = gK0 + betaK*(sigs-sigc); 
    end

    if sig2(i) < sigc   %2
        gK2(i,1)  = gK0;
    elseif sig2(i) < sigs 
        gK2(i,1)  = gK0 + betaK*(sig2(i)-sigc);
    elseif sig2(i) >= sigs
        gK2(i,1)  = gK0 + betaK*(sigs-sigc); 
    end
end

A1  = 4*pi.*r1.^2;
A2  = 4*pi.*r2.^2;

phi1 = real((F./(Cm.*A1)).*(V1.*(Nai1+Ki1-Cli1) + z*X));
phi2 = real((F./(Cm.*A2)).*(V2.*(Nai2+Ki2-Cli2) + z*X));
dphi = -(phi1-phi2);

xNa1 = real((RT/(F)).*log(Nae./Nai1));
xK1  = real((RT/(F)).*log(Ke./Ki1));
xCl1 = real((RT/(F)).*log(Cle./Cli1));

xNa2 = real((RT/(F)).*log(Nae./Nai2));
xK2  = real((RT/(F)).*log(Ke./Ki2));
xCl2 = real((RT/(F)).*log(Cle./Cli2));
    
% cell 1
mNa1 = -A1.*(gNa./(F*V1)).*(phi1 - xNa1) - A1.*3*p./(F.*V1);
cNa = (Nai1+Nai2)./2;
xgNa1 = -(Ag*w./(F.*V1)).*(dphi - (RT/F).*(log(Nai2./Nai1))) - Ag*Lpg.*cNa.*(dP1-dP2)./V1;

mK1 = -A1.*(gK1./(F.*V1)).*(phi1 - xK1) + A1.*2*p./(F.*V1);
cK = (Ki1+Ki2)/2;
xgK1 = -(Ag*w./(F.*V1)).*(dphi - (RT/F).*(log(Ki2./Ki1))) - Ag*Lpg.*cK.*(dP1-dP2)./V1;

mCl1 = A1.*(gCl./(F.*V1)).*(phi1 + xCl1);
cCl = (Cli1+Cli2)./2;
xgCl1 = (Ag.*w./(F.*V1)).*(dphi + (RT./F).*(log(Cli2./Cli1))) - Ag.*Lpg.*cCl.*(dP1-dP2)./V1;

dPi1 = RT.*(Nai1+Ki1+Cli1 + X./V1 - (Nae+Ke+Cle));
mV1 = -A1.*Lpm.*(dP1-dPi1);
gV1 = -Ag.*Lpg.*(dP1-dP2);

% cell 2
mNa2 = -A2.*(gNa./(F.*V2)).*(phi2 - xNa2) - A2.*3.*p./(F.*V2);
xgNa2 = -(Ag.*w./(F.*V2)).*(-dphi - (RT./F).*(log(Nai1./Nai2))) - Ag.*Lpg.*cNa.*(dP2-dP1)./V2;

mK2 = -A2.*(gK2./(F.*V2)).*(phi2 - xK2) + A2.*2.*p./(F.*V2);
xgK2 = -(Ag.*w./(F.*V2)).*(-dphi - (RT./F).*(log(Ki1./Ki2))) - Ag.*Lpg.*cK.*(dP2-dP1)./V2;

mCl2 = A2.*(gCl./(F.*V2)).*(phi2 + xCl2);
xgCl2 = (Ag.*w./(F.*V2)).*(-dphi + (RT./F).*(log(Cli1./Cli2))) - Ag.*Lpg.*cCl.*(dP2-dP1)./V2;

dPi2 = RT.*(Nai2+Ki2+Cli2 + X./V2 - (Nae+Ke+Cle));
mV2 = -A2.*Lpm.*(dP2-dPi2);
gV2 = -Ag.*Lpg.*(dP2-dP1);

tt = tt-3600;  % start figure at steady state prior to loading
%% Figures
for d=1:1
% Plot cell volume
hfig1=figure(1); set(hfig1,'color','w'); hold all
plot(tt(:,1)./60,V1.*1e18,'-','linewidth',6);
plot(tt(:,1)./60,V2.*1e18,'-','linewidth',6);
haxsY=gca;
set(haxsY,'tickdir','out','linewidth',3,'fontsize',25,...
    'FontName','Arial','FontWeight','bold')
set(gca,'color','none')
pbaspect([2 1 1])
xlim([0 60])
box on

% Plot Na
hfig1=figure(2); set(hfig1,'color','w'); hold all
plot(tt(:,1)./60,Nai1,'-','linewidth',8);
plot(tt(:,1)./60,Nai2,'-','linewidth',6);
haxsY=gca;
set(haxsY,'tickdir','out','linewidth',3,'fontsize',35,...
    'FontName','Arial','FontWeight','bold')
set(gca,'color','none')
pbaspect([3 1 1])
xlim([0 60])
box on

% Plot K
hfig1=figure(3); set(hfig1,'color','w'); hold all
plot(tt(:,1)./60,Ki1,'-','linewidth',8);
plot(tt(:,1)./60,Ki2,'-','linewidth',6);
haxsY=gca;
set(haxsY,'tickdir','out','linewidth',3,'fontsize',35,...
    'FontName','Arial','FontWeight','bold')
set(gca,'color','none')
pbaspect([3 1 1])
xlim([0 60])
box on

% Plot Cl
hfig1=figure(4); set(hfig1,'color','w'); hold all
plot(tt(:,1)./60,Cli1,'-','linewidth',8);
plot(tt(:,1)./60,Cli2,'-','linewidth',6);
haxsY=gca;
set(haxsY,'tickdir','out','linewidth',3,'fontsize',35,...
    'FontName','Arial','FontWeight','bold')
set(gca,'color','none')
pbaspect([3 1 1])
xlim([0 60])
box on

% Plot membrane potential
hfig1=figure(5); set(hfig1,'color','w'); hold all
plot(tt(:,1)./60,phi1,'-','linewidth',6);
plot(tt(:,1)./60,phi2,'-','linewidth',6);
haxsY=gca;
set(haxsY,'tickdir','out','linewidth',3,'fontsize',25,...
    'FontName','Arial','FontWeight','bold')
set(gca,'color','none')
pbaspect([2 1 1])
xlim([0 60])
box on

% Plot pressure
hfig1=figure(6); set(hfig1,'color','w'); hold all
plot(tt(:,1)./60,dP1./1e3,'-','linewidth',6);
plot(tt(:,1)./60,dP2./1e3,'-','linewidth',6);
haxsY=gca;
set(haxsY,'tickdir','out','linewidth',3,'fontsize',25,...
    'FontName','Arial','FontWeight','bold')
set(gca,'color','none')
pbaspect([2 1 1])
xlim([0 60])
box on

% Plot osmotic
hfig1=figure(7); set(hfig1,'color','w'); hold all
plot(tt(:,1)./60,dPi1./1e3,'-','linewidth',6);
plot(tt(:,1)./60,dPi2./1e3,'-','linewidth',6);
haxsY=gca;
set(haxsY,'tickdir','out','linewidth',3,'fontsize',25,...
    'FontName','Arial','FontWeight','bold')
set(gca,'color','none')
pbaspect([2 1 1])
xlim([0 60])
box on
end


function dxdt = dyn_ion(time,x, F, Cm, z, X, Nae, Ke, Cle, RT, p, Lpm, gNa, gK0, gCl,...
    r0,siga,Ec,h0,sigc,sigs,betaK,ds1,tl,dt1,w, Lpg, Ag)

    dxdt = zeros(8,1);
    
    Nai1 = x(1);
    Nai2 = x(2);
    Ki1  = x(3);
    Ki2  = x(4); 
    Cli1 = x(5);
    Cli2 = x(6);
    V1   = x(7);
    V2   = x(8); 

    % Function for time dependent loading profile on cell 1
    if (time>tl)&&(time<=tl+dt1)
        sig_g1 = ds1*(time-tl)/dt1; 
    elseif (time>tl+dt1)
        sig_g1 = ds1;
    else
        sig_g1 = 0;
    end
    sig_g2 = 0;
    
    % Update cell radii
    r1   = real((3*V1/(4*pi))^(1/3));
    sig1 = real((Ec/(r0))*(r1-r0) + siga);
    r2   = real((3*V2/(4*pi))^(1/3));
    sig2 = real((Ec/(r0))*(r2-r0) + siga);
    
    dP1  = 2*sig1*h0/r1 + sig_g1;   
    dP2  = 2*sig2*h0/r2 + sig_g2;   
    
    % MS channel
    if sig1 < sigc   %1
        gK1  = gK0;
    elseif sig1 < sigs 
        gK1  = gK0 + betaK*(sig1-sigc);
    elseif sig1 >= sigs
        gK1  = gK0 + betaK*(sigs-sigc); 
    end
    
    if sig2 < sigc   %2
        gK2  = gK0;
    elseif sig2 < sigs 
        gK2  = gK0 + betaK*(sig2-sigc);
    elseif sig2 >= sigs
        gK2  = gK0 + betaK*(sigs-sigc); 
    end

    A1  = 4*pi*r1^2;
    A2  = 4*pi*r2^2;
    
    phi1 = real((F/(Cm*A1))*(V1*(Nai1+Ki1-Cli1) + z*X));
    phi2 = real((F/(Cm*A2))*(V2*(Nai2+Ki2-Cli2) + z*X));
    dphi = (phi1-phi2);

    xNa1 = real((RT/(F))*log((Nae)/Nai1));
    xK1  = real((RT/(F))*log(Ke/Ki1));
    xCl1 = real((RT/(F))*log(Cle/Cli1));
    
    xNa2 = real((RT/(F))*log(Nae/Nai2));
    xK2  = real((RT/(F))*log(Ke/Ki2));
    xCl2 = real((RT/(F))*log(Cle/Cli2));
%     
    % cell 1
    mNa1 = -A1*(gNa/(F*V1))*(phi1 - xNa1) - A1*3*p/(F*V1);
    cNa = (Nai1+Nai2)/2;
    xgNa1 = -(Ag*w/(F*V1))*(dphi - (RT/F)*(log(Nai2/Nai1))) - Ag*Lpg*cNa*(dP1-dP2)/V1;
    dxdt(1) = mNa1 + xgNa1;
    
    mK1 = -A1*(gK1/(F*V1))*(phi1 - xK1) + A1*2*p/(F*V1);
    cK = (Ki1+Ki2)/2;
    xgK1 = -(Ag*w/(F*V1))*(dphi - (RT/F)*(log(Ki2/Ki1))) - Ag*Lpg*cK*(dP1-dP2)/V1;
    dxdt(3) = mK1 + xgK1;
    
    mCl1 = A1*(gCl/(F*V1))*(phi1 + xCl1);
    cCl = (Cli1+Cli2)/2;
    xgCl1 = (Ag*w/(F*V1))*(dphi + (RT/F)*(log(Cli2/Cli1))) - Ag*Lpg*cCl*(dP1-dP2)/V1;
    dxdt(5) = mCl1 + xgCl1;

    dPi1 = RT*(Nai1+Ki1+Cli1 + X/V1 - (Nae+Ke+Cle));
    mV1 = -A1*Lpm*(dP1-dPi1);
    gV1 = -Ag*Lpg*(dP1-dP2);
    dxdt(7) = mV1 + gV1;
    
    % cell 2
    mNa2 = -A2*(gNa/(F*V2))*(phi2 - xNa2) - A2*3*p/(F*V2);
    xgNa2 = -(Ag*w/(F*V2))*(-dphi - (RT/F)*(log(Nai1/Nai2))) - Ag*Lpg*cNa*(dP2-dP1)/V2;
    dxdt(2) = mNa2 + xgNa2;
    
    mK2 = -A2*(gK2/(F*V2))*(phi2 - xK2) + A2*2*p/(F*V2);
    xgK2 = -(Ag*w/(F*V2))*(-dphi - (RT/F)*(log(Ki1/Ki2))) - Ag*Lpg*cK*(dP2-dP1)/V2;
    dxdt(4) = mK2 + xgK2;
    
    mCl2 = A2*(gCl/(F*V2))*(phi2 + xCl2);
    xgCl2 = (Ag*w/(F*V2))*(-dphi + (RT/F)*(log(Cli1/Cli2))) - Ag*Lpg*cCl*(dP2-dP1)/V2;
    dxdt(6) = mCl2 + xgCl2;

    dPi2 = RT*(Nai2+Ki2+Cli2 + X/V2 - (Nae+Ke+Cle));
    mV2 = -A2*Lpm*(dP2-dPi2);
    gV2 = -Ag*Lpg*(dP2-dP1);
    dxdt(8) = mV2 + gV2;

end