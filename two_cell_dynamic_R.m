clear all
close all

% Parameters listed as per SI Table S1
wg       = 5e-7;  
Lpg      = 1e-11; 
r0       = 7.1e-6; 
Ag       = 63.33e-12; 
K        = 6e3;
sigma_a  = 100; 
h        = 0.6e-6;  
beta     = 2.0e-11; 
sigma_c  = 75; 
sigma_s  = 600; 
wl       = 1.5e-9; 
gam      = 2.25*(10^-17); 
dPic     = 4e10;  
Piext    = 0.67e6; 
Lpm      = .7e-11; 
RT       = 8.3145*310;  

% applied loading
ds1 = 150;   % growth-associated stress
tl  = 60; 
dt1 = 600;

%% Determine stress free initial conditions
sigma_g1 = 0;
sigma_g2 = 0;
xout = fsolve(@(x) solve_r(x, gam,K,h,r0,sigma_a,sigma_g1,sigma_g2,dPic,sigma_c,sigma_s,beta,wl, Ag,wg),[r0,r0]);

r1 = xout(1); % unloaded steady state radius 1
r2 = xout(2); % unloaded steady state radius 2

s1 = (K/(r0))*(r1-r0) + sigma_a; % sigma 1
s2 = (K/(r0))*(r2-r0) + sigma_a; % sigma 2

p1 = 2*s1*h/r1 + sigma_g1; % Delta P1
p2 = 2*s2*h/r2 + sigma_g2; % Delta P2

V10 = (4*pi/3)*r1^3; % cell volume
V20 = (4*pi/3)*r2^3;

n10 = (p1+Piext)*V10/RT; % no. internal ions; note: (Pin = DeltaP1 + Pext)
n20 = (p2+Piext)*V20/RT;

%% Transient solution

% Set initial conditions
tspan = 0:10:4000;
xabs1 = 1e-22;
xabs2 = 1e-16;
options = odeset('RelTol',2.5e-14,'AbsTol',[xabs2 xabs2 xabs1 xabs1]);
Xo0(1) = V10;
Xo0(2) = V20;
Xo0(3) = n10;
Xo0(4) = n20;

% ODE solver to compute volume and ions
[tt, xx]=ode23s(@dyn_ion,tspan,Xo0,options,Lpg,Lpm,wg,wl,gam,beta,K,r0,sigma_a,h,sigma_c,sigma_s,RT,Ag,dPic,Piext,ds1,tl,dt1);

V1 = xx(:,1);
V2 = xx(:,2);
n1 = xx(:,3);
n2 = xx(:,4);

%% Calculate additional values for plotting

% Define vectors
sig_g1 = zeros(length(tt),1);
Piin1 = zeros(length(tt),1);
Piin2 = zeros(length(tt),1);
dPi1 = zeros(length(tt),1);
dPi2 = zeros(length(tt),1);
sig1 = zeros(length(tt),1);
sig2 = zeros(length(tt),1);
dP1 = zeros(length(tt),1);
dP2 = zeros(length(tt),1);
wms1 = zeros(length(tt),1);
wms2 = zeros(length(tt),1);
Jvg1 = zeros(length(tt),1);
Jvg2 = zeros(length(tt),1);
Jvm1 = zeros(length(tt),1);
Jvm2 = zeros(length(tt),1);
ng1 = zeros(length(tt),1);
ng2 = zeros(length(tt),1);
nms1 = zeros(length(tt),1);
nms2 = zeros(length(tt),1);
nl1 = zeros(length(tt),1);
nl2 = zeros(length(tt),1);
np1 = zeros(length(tt),1);
np2 = zeros(length(tt),1);

for i=1:length(tt)
      
    % Function for time dependent loading profile on cell 1
    if (tt(i)>tl)&&(tt(i)<=tl+dt1)
        sig_g1(i) = ds1*(tt(i)-tl)/dt1;
    elseif (tt(i)>tl+dt1)
        sig_g1(i) = ds1;
    else
        sig_g1(i) = 0;
    end
    
    sig_g2 = 0;
    
    % Update osmotic pressure
    Piin1(i) = n1(i)*RT/V1(i);
    dPi1(i)  = Piin1(i) - Piext;
    
    Piin2(i) = n2(i)*RT/V2(i);
    dPi2(i)  = Piin2(i) - Piext;
    
    % Update cell radii
    r1(i) = (3*V1(i)/(4*pi))^(1/3);
    r2(i) = (3*V2(i)/(4*pi))^(1/3);
    
   % Update membrane stress and hydrostatic pressure (cell 1)
    sig1(i) = (K/(r0))*(r1(i)-r0) + sigma_a;
    dP1(i)  = 2*sig1(i)*h/r1(i) + sig_g1(i);   
    if sig1(i) < sigma_c
        wms1(i) = 0.;
    elseif sig1(i) < sigma_s
        wms1(i) = beta*(sig1(i)-sigma_c);
    elseif sig1(i) >= sigma_s
        wms1(i) = beta*(sigma_s-sigma_c);
    end
    
   % Update membrane stress and hydrostatic pressure (cell 2)
    sig2(i) = (K/(r0))*(r2(i)-r0) + sigma_a;
    dP2(i)  = 2*sig2(i)*h/r2(i) + sig_g2;    
    if sig2(i) < sigma_c
        wms2(i) = 0;
    elseif sig2(i) < sigma_s
        wms2(i) = beta*(sig2(i)-sigma_c);
    elseif sig2(i) >= sigma_s
        wms2(i) = beta*(sigma_s-sigma_c);
    end

   % Water flux across gap junctions
    Jvg1(i) = -Lpg*(dP1(i) -(dP2(i)));
    Jvg2(i) = -Lpg*(dP2(i) -(dP1(i)));
    
   % Water flux across cell membrane
    Jvm1(i) = -Lpm*(dP1(i)-dPi1(i));
    Jvm2(i) = -Lpm*(dP2(i)-dPi2(i));
    
   % Ion flux across gap junctions
    ng1(i)  = -wg*(dPi1(i)-dPi2(i)) -((Piext)/RT)*Lpg*(dP1(i) -dP2(i));
    ng2(i)  = -wg*(dPi2(i)-dPi1(i)) -((Piext)/RT)*Lpg*(dP2(i) -dP1(i));
    
   % Ion flux across MS channels
    nms1(i)  = -wms1(i)*dPi1(i);
    nms2(i)  = -wms2(i)*dPi2(i);
    
   % Ion flux across leak channels
    nl1(i)  = -wl*dPi1(i);
    nl2(i)  = -wl*dPi2(i);
    
   % Ion flux due to active transport
    np1(i)  = gam*(dPic - dPi1(i));
    np2(i)  = gam*(dPic - dPi2(i));
    
end

% Plot cell volume
hfig1=figure(1); set(hfig1,'color','w'); hold all
plot(tt(:,1)./60,xx(:,1).*1e18,'-','linewidth',6,'color',[0.4660, 0.6740, 0.1880]);
plot(tt(:,1)./60,xx(:,2).*1e18,'-','linewidth',6,'color',[1 .5 0]);
haxsY=gca;
set(haxsY,'tickdir','out','linewidth',3,'fontsize',25,...
    'FontName','Arial','FontWeight','bold')
set(gca,'color','none')
pbaspect([1.5 1 1])
box on

% GJ ion flow rate
hfig1=figure(2); set(hfig1,'color','w'); hold all
plot(tt(:,1)./60,ng1.*Ag.*1e15,'-','linewidth',6,'color',[0.4660, 0.6740, 0.1880]);
plot(tt(:,1)./60,ng2.*Ag.*1e15,'-','linewidth',6,'color',[1 .5 0]);
haxsY=gca;
set(haxsY,'tickdir','out','linewidth',4,'fontsize',32.5,...
    'FontName','Arial','FontWeight','bold')
set(gca,'color','none')
pbaspect([2 1 1])
box on

% membrane ion flow rate
hfig1=figure(3); set(hfig1,'color','w'); hold all
plot(tt(:,1)./60,(nms1+nl1+np1).*(4*pi*r1'.^2).*1e15,'-','linewidth',6,'color',[0.4660, 0.6740, 0.1880]);
plot(tt(:,1)./60,(nms2+nl2+np2).*(4*pi*r2'.^2).*1e15,'-','linewidth',6,'color',[1 .5 0]);
haxsY=gca;
set(haxsY,'tickdir','out','linewidth',4,'fontsize',32.5,...
    'FontName','Arial','FontWeight','bold')
set(gca,'color','none')
set(gca, 'XTick', [0 30 60])
ylim([-.5 .5])
% xlim([0 120])
pbaspect([2 1 1])
box on

% membrane water flow rate
hfig1=figure(4); set(hfig1,'color','w'); hold all
plot(tt(:,1)./60,Jvm1.*(4*pi*r1'.^2).*1e18,'-','linewidth',6,'color',[0.4660, 0.6740, 0.1880]);
plot(tt(:,1)./60,Jvm2.*(4*pi*r2'.^2).*1e18,'-','linewidth',6,'color',[1 .5 0]);
haxsY=gca;
set(haxsY,'tickdir','out','linewidth',4,'fontsize',32.5,...
    'FontName','Arial','FontWeight','bold')
set(gca,'color','none')
set(gca,'linewidth',4)
pbaspect([2 1 1])
box on

% GJs water flow rate
hfig1=figure(5); set(hfig1,'color','w'); hold all
plot(tt(:,1)./60,Jvg1.*Ag.*1e18,'-','linewidth',6,'color',[0.4660, 0.6740, 0.1880]);
plot(tt(:,1)./60,Jvg2.*Ag.*1e18,'-','linewidth',6,'color',[1 .5 0]);
haxsY=gca;
set(haxsY,'tickdir','out','linewidth',4,'fontsize',32.5,...
    'FontName','Arial','FontWeight','bold')
set(gca,'color','none')
pbaspect([2 1 1])
box on

% MS channel ion flux
hfig1=figure(9); set(hfig1,'color','w'); hold all
plot(tt(:,1)./60,nms1,'-','linewidth',6,'color',[0.4660, 0.6740, 0.1880]);
plot(tt(:,1)./60,nms2,'-','linewidth',6,'color',[1 .5 0]);
haxsY=gca;
set(haxsY,'tickdir','out','linewidth',3,'fontsize',25,...
    'FontName','Arial','FontWeight','bold')
set(gca,'color','none')
set(gca,'linewidth',3)
pbaspect([1.1 1 1])
box on

% leak channel ion flux
hfig1=figure(10); set(hfig1,'color','w'); hold all
plot(tt(:,1)./60,nl1,'-','linewidth',6,'color',[0.4660, 0.6740, 0.1880]);
plot(tt(:,1)./60,nl2,'-','linewidth',6,'color',[1 .5 0]);
haxsY=gca;
set(haxsY,'tickdir','out','linewidth',3,'fontsize',25,...
    'FontName','Arial','FontWeight','bold')
set(gca,'color','none')
pbaspect([1.1 1 1])
box on

% pump ion flux
hfig1=figure(11); set(hfig1,'color','w'); hold all
plot(tt(:,1)./60,np1,'-','linewidth',6,'color',[0.4660, 0.6740, 0.1880]);
plot(tt(:,1)./60,np2,'-','linewidth',6,'color',[1 .5 0]);
haxsY=gca;
set(haxsY,'tickdir','out','linewidth',3,'fontsize',25,...
    'FontName','Arial','FontWeight','bold')
set(gca,'color','none')
set(gca, 'XTick', [0 30 60])
ylim([0 1.5e-6])
% xlim([0 120])
pbaspect([1.1 1 1])
box on

% hydrostatic pressure
hfig1=figure(6); set(hfig1,'color','w'); hold all
plot(tt(:,1)./60,dP1,'-','linewidth',6,'color',[0.4660, 0.6740, 0.1880]);
plot(tt(:,1)./60,dP2,'-','linewidth',6,'color',[1 .5 0]);
haxsY=gca;
set(haxsY,'tickdir','out','linewidth',3,'fontsize',25,...
    'FontName','Arial','FontWeight','bold')
% set(xlbl,'fontsize',25,'FontName','Arial','FontWeight','bold');
% set(ylbl,'fontsize',25,'FontName','Arial','FontWeight','bold'); %title('r17 v')
% set(gca,'TickLabelInterpreter','latex')
set(gca,'color','none')
pbaspect([1.1 1 1])
box on

% osmotic pressure
hfig1=figure(7); set(hfig1,'color','w'); hold all
plot(tt(:,1)./60,dPi1,'-','linewidth',6,'color',[0.4660, 0.6740, 0.1880]);
plot(tt(:,1)./60,dPi2,'-','linewidth',6,'color',[1 .5 0]);
haxsY=gca;
set(haxsY,'tickdir','out','linewidth',3,'fontsize',25,...
    'FontName','Arial','FontWeight','bold')
pbaspect([1.1 1 1])
box on

% applied cell 1 stress
hfig1=figure(8); set(hfig1,'color','w'); hold all
plot(tt(:,1)./60,sig_g1,'k-','linewidth',6);
haxsY=gca;
set(haxsY,'tickdir','out','linewidth',3,'fontsize',25,...
    'FontName','Arial','FontWeight','bold')
set(gca,'color','none')
set(gca,'linewidth',3)
pbaspect([1.1 1 1])
box on

%% Functions
function F = solve_r(x, gamma,K,h,r0,sigma_a,sigma_g1,sigma_g2,Pi_c,sigma_c,sigma_s,beta,wl, A_g,wg)

    s1 = K*(x(1) - r0)/(r0)  + sigma_a; 
    if s1 < sigma_c
        wms1 = wl;
    elseif s1 < sigma_s
        wms1 = beta*(s1-sigma_c) + wl;
    elseif s1 >= sigma_s
        wms1 = beta*(sigma_s-sigma_c) + wl;
    end
    
    s2 = K*(x(2) - r0)/(r0)  + sigma_a; 
    if s2 < sigma_c
        wms2 = wl;
    elseif s2 < sigma_s
        wms2 = beta*(s2-sigma_c) + wl;
    elseif s2 >= sigma_s
        wms2 = beta*(sigma_s-sigma_c) + wl;
    end
    
    % Note for wms1 and wms2, MS and leak channels both included as wms*Pi + wl*Pi = (wms+wl)*Pi
    
    r1 = x(1);
    r2 = x(2);
%     
    % Eqns F obtained by solving Eqns 1-6 at steady state (direct
    % implementation would also work)
    F(1) = sigma_a - (K*(r0 - r1))/r0 + (r1*(sigma_g1 - (4*Pi_c*pi*gamma^2*r1^2*r2^2 + 4*Pi_c*wms2*pi*gamma*r1^2*r2^2 + A_g*Pi_c*wg*gamma*r1^2 + A_g*Pi_c*wg*gamma*r2^2)/(4*pi*gamma^2*r1^2*r2^2 + A_g*wms1*wg*r1^2 + A_g*wms2*wg*r2^2 + A_g*gamma*wg*r1^2 + A_g*gamma*wg*r2^2 + 4*pi*wms1*wms2*r1^2*r2^2 + 4*pi*wms1*gamma*r1^2*r2^2 + 4*pi*wms2*gamma*r1^2*r2^2)))/(2*h);
    F(2) = sigma_a - (K*(r0 - r2))/r0 + (r2*(sigma_g2 - (Pi_c*gamma*(A_g*wg*r1^2 + A_g*wg*r2^2 + 4*pi*wms1*r1^2*r2^2 + 4*pi*gamma*r1^2*r2^2))/(4*pi*gamma^2*r1^2*r2^2 + A_g*wms1*wg*r1^2 + A_g*wms2*wg*r2^2 + A_g*gamma*wg*r1^2 + A_g*gamma*wg*r2^2 + 4*pi*wms1*wms2*r1^2*r2^2 + 4*pi*wms1*gamma*r1^2*r2^2 + 4*pi*wms2*gamma*r1^2*r2^2)))/(2*h);
    
end

function dxdt = dyn_ion(time,x,Lpg,Lpm,wg,wl,gam,beta,K,r0,siga,h0,sigma_c,sigma_s,RT,Ag,dPic,Piext,ds1,tl,dt1)

    dxdt = [0 0 0 0]';

    V1 = x(1);
    V2 = x(2);
    n1 = x(3);
    n2 = x(4);
    
    sig_g2 = 0; % No applied stress on cell 2
       
    % Function for time dependent loading profile on cell 1
    if (time>tl)&&(time<=tl+dt1)
        sig_g1 = ds1*(time-tl)/dt1;
    elseif (time>tl+dt1)
        sig_g1 = ds1;
    else
        sig_g1 = 0;
    end
    
    % Update osmotic pressures
    Piin1 = n1*RT/V1;
    dPi1  = Piin1 - Piext;
    
    Piin2 = n2*RT/V2;
    dPi2  = Piin2 - Piext;
    
    % Update cell radii
    r1 = (3*V1/(4*pi))^(1/3);
    r2 = (3*V2/(4*pi))^(1/3);
    
   % Update membrane stress and hydrostatic pressure (cell 1)
    sig1 = (K/(r0))*(r1-r0) + siga;
    dP1  = 2*sig1*h0/r1 + sig_g1;   
    if sig1 < sigma_c
        wms1 = 0;
    elseif sig1 < sigma_s
        wms1 = beta*(sig1-sigma_c);
    elseif sig1 >= sigma_s
        wms1 = beta*(sigma_s-sigma_c);
    end
    
   % Update membrane stress and hydrostatic pressure (cell 2)
    sig2 = (K/(r0))*(r2-r0) + siga;
    dP2  = 2*sig2*h0/r2 + sig_g2;    
    if sig2 < sigma_c
        wms2 = 0.0;
    elseif sig2 < sigma_s
        wms2 = beta*(sig2-sigma_c);
    elseif sig2 >= sigma_s
        wms2 = beta*(sigma_s-sigma_c);
    end

   % Water flux across gap junctions 
    Jvg1 = -Lpg*(dP1 -(dP2));
    Jvg2 = -Lpg*(dP2 -(dP1));
    
   % Water flux across cell membrane
    Jvm1 = -Lpm*(dP1-dPi1);
    Jvm2 = -Lpm*(dP2-dPi2);
    
   % Ion flux across gap junctions
    ng1  = -wg*(dPi1-dPi2) -((Piext)/RT)*Lpg*(dP1 -dP2);
    ng2  = -wg*(dPi2-dPi1) -((Piext)/RT)*Lpg*(dP2 -dP1);
    
   % Ion flux across MS channels
    nms1  = -wms1*dPi1;
    nms2  = -wms2*dPi2;
    
   % Ion flux across leak channels 
    nl1  = -wl*dPi1;
    nl2  = -wl*dPi2;
    
   % Ion flux due to active transport
    np1  = gam*(dPic - dPi1);
    np2  = gam*(dPic - dPi2);
    
    % cell volume change
    dxdt(1) = Ag*Jvg1 + (4*pi*r1^2)*Jvm1;
    dxdt(2) = Ag*Jvg2 + (4*pi*r2^2)*Jvm2;
    
    % cell number of ions change
    dxdt(3) = Ag*ng1 + (4*pi*r1^2)*(nms1 + nl1  + np1);
    dxdt(4) = Ag*ng2 + (4*pi*r2^2)*(nms2 + nl2  + np2); 

end
