set(0,'defaultaxesfontsize',20)
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultLineLineWidth',2);
set(0,'Defaultaxeslinewidth',2)

set(0,'DefaultFigureWindowStyle','docked')

%defining constants
c_c = 299792458;
c_eps_0 = 8.8542149e-12;
c_eps_0_cm = c_eps_0/100;
c_mu_0 = 1/c_eps_0/c_c^2;
c_q = 1.60217653e-19;
c_hb = 1.05457266913e-34;
c_h = c_hb*2*pi;

RL = 0i;
RR = 0i;

%beta_i = 0;
beta_r = 0;
% gain_z = gain.*(N-Ntr)./v_g;
% beta_i = (gain_z - alpha)./2;
% beta = beta_r + 1i*beta_i;
beta_spe = .3e-5;
gamma = 1.0;
SPE = 0;


kappa0 = 0;
kappaStart = 0;
kappaStop = 0;
% kappaStart1 = 0.3;
% kappaStop1 = 0.5;
% kappaStart2 = 0.6;
% kappaStop2 = 0.7;
% kappaStart3 = 0.8;
% kappaStop3 = 1;


g_fwhm = 3.53e+12/10;
LGamma = g_fwhm*2*pi;
Lw0 = 0.0;
LGain = 0.05;


%Input parameters 
InputParasL.E0 = 1e5;
InputParasL.we = 0;
InputParasL.t0 = 50e-12;
InputParasL.wg = 50e-13;
InputParasL.phi = 0;
InputParasL.rep = 500e-12;
InputParasR = 0;
%InputParasR.E0 = 1e5;
%InputParasR.we = 10e10;
%InputParasR.t0 = 2e-12;
%InputParasR.wg = 10e-10;
%InputParasR.phi = 0;

%defining the refractive index and wavelength
n_g = 3.5;
vg = c_c/n_g*1e2;
Lambda = 1550e-9;
f0 = c_c/Lambda;

Ntr = 1e18;
gain = vg*2.5e-16;
eVol = 1.5e-10*c_q;
Ion = 0.25e-9;
Ioff = 3e-9;
I_off = 0.024;
I_on = 0.1;
taun = 1e-9;
Zg = sqrt(c_mu_0/c_eps_0)/n_g;
EtoP = 1/(Zg*f0*vg*1e-2*c_hb);
alpha = 0;


%simulation setup
plotN = 100;

%boundry functions
L = 1000e-6*1e2;
XL = [0,L];
YL = [0,InputParasL.E0];


Nz = 50;
dz = L/(Nz-1);
dt = dz/vg;
fsync = dt*vg/dz;

Nt = floor(400*Nz);
tmax = Nt*dt;
t_L = dt*Nz;

z = linspace(0,L,Nz).';
time = nan(1,Nt);
InputL = nan(1,Nt);
InputR = nan(1,Nt);
OutputL = nan(1,Nt);
OutputR = nan(1,Nt);
Nave = nan(1,Nt);

%Initializing electric field vectors
Ef =zeros(size(z));
Er = zeros(size(z));

Pf = zeros(size(z));
Pr = zeros(size(z));

Efp = Ef;
Erp = Er;
Pfp = Pf;
Prp = Pr;

%defining the electric field vectors
Ef1 = @SourceFct;
ErN = @SourceFct;

%Input functions
t = 0;
time(1) = t;
InputL(1) = Ef1(t,InputParasL);
InputR(1) = ErN(t,InputParasR);

%output functions

OutputR(1) = Er(1);
OutputL(1) = Ef(Nz);


%plotting the vectors
Ef(1) = InputL(1);
Er(Nz) = InputR(1);
% 
% figure('name','Fields')
% subplot(3,2,1)
% plot(z*10000,real(Ef),'r');
% hold off
% xlabel('z(\mum)')
% ylabel('E_f')
% subplot(3,2,3)
% plot(z*10000,real(Er),'b');
% xlabel('z(\mum)')
% ylabel('E_r')
% hold off
% subplot(3,2,5)
% plot(time*1e12,real(InputL),'r'); hold on
% plot(time*1e12,real(OutputR),'r--');
% plot(time*1e12,real(InputR),'b'); hold on
% plot(time*1e12,real(OutputL),'b--');
% xlabel('time(ps)')
% ylabel('E')
% 
% hold off
% 
N = ones(size(z))*Ntr;

Nave(1) = mean (N);
for i = 2:Nt
    t = dt*(i-1);
    time(i) = t;

    InputL(i) = Ef1(t,InputParasL);%setting input equal to source fucntion
    InputR(i) = ErN(t,InputParasR);%setting input equal to source fucntion

    
    Ef(1) = InputL(i) + RL*Er(1);%setting The electric field vectors equal to the Inputs with the reflectio added
    Er(Nz) = InputR(i) + RR*Ef(Nz);%setting The electric field vectors equal to the Inputs with the reflectio added
    
    %beta = ones(size(z)).*(beta_r+1i*beta_i);
    exp_det = exp(-1i*dz*beta);
    
    kappa = kappa0*ones(size(z));
    kappa(z<L*kappaStart) = 0;
    kappa(z>L*kappaStop) = 0;
    % kappa(z<L*kappaStart1) = 80;
    % kappa(z>L*kappaStop1) = 80;
    % kappa(z<L*kappaStart2) = 50;
    % kappa(z>L*kappaStop2) = 50;
    % kappa(z<L*kappaStart3) = 80;
    % kappa(z>L*kappaStop3) = 80;


    Ef_temp = Ef(1: Nz-1);

    Pf(1) = 0;
    Pf(Nz) = 0;
    Pr(1) = 0;
    Pr(Nz) = 0;
    Cw0 = -LGamma + 1i*Lw0;

    Ef(2:Nz) = fsync*exp_det(1:Nz-1).*Ef(1:Nz-1) + 1i*dz*kappa(2:Nz).*Er(2:Nz);
    Er(1:Nz-1) = fsync*exp_det(2:Nz).*Er(2:Nz) + 1i*dz*kappa(1:Nz-1).*Ef_temp(1:Nz-1);


    Tf = LGamma*Ef(1:Nz-2) + Cw0*Pfp(2:Nz-1) + LGamma*Efp(1:Nz-2);
    Pf(2:Nz-1) = (Pfp(2:Nz-1) + 0.5*dt*Tf)./(1-0.5*dt*Cw0);
    Tr = LGamma*Er(3:Nz) + Cw0*Prp(2:Nz-1) + LGamma*Erp(3:Nz);
    Pr(2:Nz-1) = (Prp(2:Nz-1) + 0.5*dt*Tr)./(1-0.5*dt*Cw0);

    Ef(2:Nz-1) = Ef(2:Nz-1) - LGain*(Ef(2:Nz-1)-Pf(2:Nz-1));
    Er(2:Nz-1) = Er(2:Nz-1) - LGain*(Er(2:Nz-1)-Pr(2:Nz-1));


    gain_z = gain.*(N-Ntr)./vg;
    beta_i = (gain_z - alpha)./2;
    beta = beta_r + 1i*beta_i;
    
    A = sqrt(gamma*beta_spe*c_hb*f0*L*1e-2/taun)/(2*Nz);
    if SPE>0
        Tf = (randn(Nz,1)+1i*randn(Nz,1))*A;
        Tr = (randn(Nz,1)+1i*randn(Nz,1))*A;
    else
        Tf = (ones(Nz,1))*A;
        Tr = (ones(Nz,1))*A;
    end

    EsF = Tf*abs(SPE).*sqrt(N.*1e6);
    EsR = Tr*abs(SPE).*sqrt(N.*1e6);

    Ef = Ef +EsF;
    Er = Er +EsR;

    Efp = Ef;
    Erp = Er;
    Pfp = Pf;
    Prp = Pr;
    
 
    % if GenGifs
    %     system(['rm' gifFile]);
    % end
    

    S = (abs(Ef).^2 + abs(Er).^2).*EtoP*1e-6;

    if t < Ion || t > Ioff
        I_injv = I_off;
    else 
        I_injv = I_on;
    end
    Stim = gain.*(N - Ntr).*S;
    N = (N + dt*(I_injv/eVol - Stim))./(1+dt/taun);
    Nave(i) = mean (N);

 

    OutputR(i) = Ef(Nz)*(1-RR);
    OutputL(i) = Er(1)*(1-RL);

    %plotting Ef
    if mod(i,plotN) == 0
        subplot(3,4,1)
        plot(z*10000,real(Ef),'r'); hold on
        plot(z*10000,imag(Ef),'r--'); hold off
        xlim(XL*1e4)
        % ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_f')
        legend('\Re','\Im')
        hold off
        subplot(3,4,2)
        plot(z*10000, real(N),'r'); hold on
         xlabel('z(\mum)')
        ylabel('N')
        hold off

        subplot(3,4,[5,6])
        plot(time*1e12,Nave,'b'); hold on
         xlabel('time(ps)')
        ylabel('Nave')
        hold off

        subplot(3,4,[9,10])
        plot(time*1e12,real(InputL),'r'); hold on
        plot(time*1e12,real(OutputR),'b');
        plot(time*1e12,real(InputR),'g');
        plot(time*1e12,real(OutputL),'m');
        %plot(time*1e12,real(Nave),'b'); hold on
        % xlim([0,Nt*dt*1e12])
        % ylim(YL)
        xlabel('time(ps)')
        ylabel('O')
        legend('Left Input', 'Right Output', 'Right Input', 'Left Output' ...
            ,'Location','east')
        hold off
        % subplot(3,2,1)
        % plot(z*10000,real(Ef),'r'); hold on
        % plot(z*10000,imag(Ef),'r--'); hold off
        % % xlim(XL*1e4)
        % % ylim(YL)
        % xlabel('z(\mum)')
        % ylabel('E_f')
        % legend('\Re','\Im')
        % hold off
        % subplot(3,2,2)
        % % plot(z*10000,real(Er),'b'); hold on
        % % plot(z*10000,imag(Er),'b--'); hold off
        % plot(z*10000,real(N),'r'); hold on
        % % xlim(XL*1e4)
        % % ylim(YL)
        % xlabel('z(\mum)')
        % ylabel('E_r')
        % legend('\Re','\Im')
        % 
        % hold off
        % %plotting Er
        % subplot(3,2,[3,4])
        % % plot(time*1e12,real(InputL),'r'); hold on
        % % plot(time*1e12,real(OutputR),'b');
        % % plot(time*1e12,real(InputR),'g');
        % % plot(time*1e12,real(OutputL),'m');
        % plot(time*1e12,real(Nave),'b'); hold on
        % % xlim([0,Nt*dt*1e12])
        % % ylim(YL)
        % xlabel('time(ps)')
        % ylabel('0')
        % legend('Left Input', 'Right Output', 'Right Input', 'Left Output' ...
        %     ,'Location','east')
        % hold off
        pause(0.01)
    end
end
fftOutput = fftshift(fft(OutputR));
fftInput = fftshift(fft(InputL));
%fftOutputL = fftshift(fft(OutputL));
omega = fftshift(wspace(time));
 subplot(3,4,[3,4])
        plot(time*1e12,real(InputL),'r'); hold on
        plot(time*1e12,real(OutputR),'g');
        plot(time*1e12,real(InputR),'b');
        plot(time*1e12,real(OutputL),'m');

        xlabel('time(s)')
        ylabel(' output')
        legend('Right','Left')
        hold off
  subplot(3,4,[7,8])
        plot( omega, 20*log(fftOutput),'b'); hold on
         plot(omega, 20*log(fftInput),'r');
         xlim([-1,1]*1e13)
        xlabel('Omega(Hz)')
        ylabel('20log|E|')
        legend('Output','Input')
        hold off
  subplot(3,4,[11,12])
        plot( omega, unwrap(angle(fftOutput)),'r'); hold on
        plot(omega ,unwrap(angle(fftInput)),'b');
        xlabel('Omega(Hz)')
        ylabel('phase(E)')
        hold off
