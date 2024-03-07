%Milestone 6
set(0, 'defaultaxesfontsize',20)
set(0,'DefaultFigureWindowStyle','normal')
set(0,'DefaultLineLineWidth',2);
set(0,'Defaultaxeslinewidth',2)

set(0, 'DefaultFigurePosition', get(0, 'Screensize'))
set(0,'DefaultFigureWindowStyle','normal')

c_c = 299792458;            % m/s TWM speed of light
c_eps_0 = 8.8542149e-12;    % F/m vacuum permittivity
c_eps_0_cm = c_eps_0/100;   % F/cm
c_mu_0 = 1/c_eps_0/c_c^2;
c_q = 1.60217653e-19;
c_hb = 1.05457266913e-34;                % Dirac constant
c_h = c_hb*2*pi;

%Milestone 1
RL = 0.0i;  
RR = 0.0i;  %Reflective Efficiency 0.9i

%Milestone 2
beta_i = 0;
beta_r = 0;

InputParasL.E0=100e5;     %Amplitude? 100e5
InputParasL.we = 0;   %Frequency for modulation (1e13)
InputParasL.t0 = 2e-12;
InputParasL.wg = 5e-13; %5e-13 
InputParasL.phi = 0;
InputParasL.rep = 5e-10;
InputParasR = 0;

n_g = 3.5; 
vg = c_c/n_g*1e2;       % TWM cm/s group velocity
Lambda = 1550e-9;

plotN = 100;

L = 1000e-6*1e2;    %cm
XL = [0,L];
YL =[1*-InputParasL.E0,1*InputParasL.E0]; %vertical scale
%YL = [-1,1];

Nz =500;            
dz =10*L/(Nz-1);
dt = dz/vg;
fsync = dt*vg/dz;

Nt =floor(40*Nz);        %designates length of simulation
tmax = Nt*dt;
t_L = dt*Nz;               % time to travel length

z = linspace(0,L,Nz).';    % Nz points, nz-1 segments
time = nan(1,Nt);
InputL = nan(1,Nt);        % create arrays of "NaN" 1 x Nt
InputR = nan(1,Nt);
OutputL = nan(1,Nt);
OutputR = nan(1,Nt);

%Milestone 3
kappa0 = 0; %100
kappaStart = 1/3;
kappaStop = 2/3;
kappa = kappa0*ones(size(z)); %fill a matrix of size z with the value of kappa0
kappa(z<L*kappaStart) = 0;
kappa(z>L*kappaStop) = 0;

Ef = zeros(size(z));       % create array of 0 
Er = zeros(size(z));

%Milestone 4
Pf = zeros(size(z));
Pr = zeros(size(z));
Efp = Ef;
Erp = Er;
Pfp = Pf;
Prp = Pr;
g_fwhm = 5*3.53e+012/10;
LGamma = g_fwhm*2*pi;
Lw0 = 0.0;
LGain = 0;

%Milestone 6
Ntr = 1e18;
n_g = 3.5;
v_g = c_c/n_g*1e2; %TWM cm/s group velocity
Lambda = 1550e-9; %Cm
f0 = c_c/Lambda;
Nave = nan(1,Nt);
N = ones(size(z))*Ntr;
Nave(1)=mean(N);
% if GenGifs
%     system(['rm' gifFile]);
% end
gain = v_g*2.5e-16;
eVol = 1.5e-10*c_q;
Ion = 0.25e-9;
Ioff = 3e-9;
I_off = 0.024;
I_on = 0.1;
taun = 1e-9;
Zg = sqrt(c_mu_0/c_eps_0)/n_g;
EtoP = 1/(Zg*f0*v_g*1e-2*c_hb);
alpha = 0;

%Milestone 7
beta_r = 0;
gain_z = gain.*(N-Ntr)./v_g;
beta_i = (gain_z - alpha)./2;
beta = beta_r + 1i*beta_i;
beta_spe = .3e-5;
gamma = 1.0;
SPE = 7;

Ef1 = @SourceFct; %Handle creation
ErN = @SourceFct;

t = 0;                      %  start of time 
time(1) = t;

InputL(1) = Ef1(t,InputParasL);    %Calling SourceFct
InputR(1) = ErN(t,InputParasR);

OutputR(1) = Ef(Nz);        %Replacing location is OutputR with location in Ef
OutputL(1) = Er(1);

Ef(1) = InputL(1);          %Replacing location in Ef with location InputL
Er(Nz) = InputR(1);  

%Milestone 2
%beta = ones(size(z))*(beta_r+1i*beta_i); %Initializing Beta
exp_det = exp(-1i*dz*beta);

%Milestone 6-7 Figure
figure('name', 'Spontainious')
subplot(3,4,1)
plot(1);
xlabel('z(\mum)')
ylabel('Ef')
hold off
subplot(3,4,2)
plot(z*10000,real(N));
hold off
xlabel('z(\mum)')
ylabel('y')
subplot(3,4,[3,4])
plot(time*1e12,real(OutputR),'r');hold on
plot(time*1e12,real(InputL),'g');
xlabel('time(ps)')
ylabel('0utput')
hold off
subplot(3,4,[5,6])
plot(1);
xlabel('time(ps)')
ylabel('Nave')
hold off
subplot(3,4,[7,8])
plot(1);
xlabel('GHz')
ylabel('20 log|E|')
hold off
subplot(3,4,[9,10])
plot(time*1e12,real(InputL),'r'); hold on
plot(time*1e12,real(OutputR),'r--'); 
plot(time*1e12,real(InputR),'b'); hold on
plot(time*1e12,real(OutputL),'b--');
xlabel('time(ps)')
ylabel('E')
hold off
subplot(3,4,[11,12])
plot(1);
xlabel('GHz')
ylabel('phase (E)')
hold off

%Create all initial graphs
% figure('name', 'Fields')
% subplot(3,1,1)
% plot(z*10000,real(Ef),'r');
% hold off
% xlabel('z(\mum)')
% ylabel('E_f')
% subplot(3,1,2)
% plot(z*10000,real(Er),'b');
% xlabel('z(\mum)')
% ylabel('E_r')
% hold off
% subplot(3,1,3)
% plot(time*1e12,real(InputL),'r'); hold on
% plot(time*1e12,real(OutputR),'r--'); 
% plot(time*1e12,real(InputR),'b'); hold on
% plot(time*1e12,real(OutputL),'b--');
% xlabel('time(ps)')
% ylabel('E')
% 
% hold off

%Keep updating graphs and recalculating propegation 
for i = 2:Nt
    t = dt*(i-1);
    time(i) = t;
    Pf(1) = 0;
    Pf(Nz) = 0;
    Pr(1) = 0;
    Pr(Nz) = 0;
    Cw0 = -LGamma + 1i*Lw0;
    if i > 12000
        InputParasL.rep = 1;
    end
    %Milestone 6
    S = (abs(Ef).^2 +abs(Er).^2).*EtoP*1e-6;
    if t < Ion || t > Ioff
        I_injv = I_off;
    else
        I_injv = I_on;
    end
    Stim = gain.*(N-Ntr).*S;
    N = (N + dt*(I_injv/eVol - Stim))./(1+ dt/taun);
    Nave(i) = mean(N);
    

    %Milestone7
    A = sqrt(gamma*beta_spe*c_hb*f0*L*1e-2/taun)/(2*Nz);
    if SPE > 0
        Tf = (randn(Nz,1)+1i*randn(Nz,1))*A;
        Tr = (randn(Nz,1)+1i*randn(Nz,1))*A;
    else
        Tf = (ones(Nz,1))*A;
        Tr = (ones(Nz,1))*A;
    end
    EsF = Tf*abs(SPE).*sqrt(N.*1e6);
    EsR = Tr*abs(SPE).*sqrt(N.*1e6);



    InputL(i) = Ef1(t,InputParasL);
    InputR(i) = ErN(t,0);

    Ef(1) = InputL(i) + RL*Er(1); %Milestone 1 RL
    Er(Nz) = InputR(i) + RR*Ef(Nz); %Milestone 1 RR
    
    Ef_temp = Ef(1:Nz-1);
    %Milestone 3 inlcude kappa and couple
    Ef(2:Nz) = fsync*Ef(1:Nz-1) + 1i*dz*kappa(2:Nz).*Er(2:Nz);
    %Milestone 3 include kappa and couple 
    Er(1:Nz-1) = fsync*Er(2:Nz) + 1i*dz*kappa(1:Nz-1).*Ef(1:Nz-1);

    Tf = LGamma*Ef(1:Nz-2) + Cw0*Pfp(2:Nz-1) + LGamma*Efp(1:Nz-2);
    Pf(2:Nz-1) = (Pfp(2:Nz-1)+0.5*dt*Tf)./(1-0.5*dt*Cw0);
    Tr = LGamma*Er(3:Nz) + Cw0*Prp(2:Nz-1) + LGamma*Erp(3:Nz);
    Pr(2:Nz-1) = (Prp(2:Nz-1) + 0.5*dt*Tr)./(1-0.5*dt*Cw0);

    Ef(2:Nz-1) = Ef(2:Nz-1) - LGain*(Ef(2:Nz-1)-Pf(2:Nz-1));
    Er(2:Nz-1) = Er(2:Nz-1) - LGain*(Er(2:Nz-1)-Pr(2:Nz-1));
    
    %Milestone 7
    Ef = Ef + EsF;
    Er = Er + EsR;

    %Milestone 1
    OutputR(i) = Ef(Nz)*(1-RR);   
    OutputL(i) = Er(1)*(1-RL);    %Reflecting

    if mod(i,plotN) == 0
        % subplot(3,1,1)
        % plot(z*10000,real(Ef),'r'); hold on
        % plot(z*10000,imag(Ef),'r--'); hold off
        % xlim(XL*1e4)
        % ylim(YL)
        % xlabel('z(\mum)')
        % ylabel('E_f')
        % legend('\Re','\Im')
        % hold off
        % subplot(3,1,2)
        % plot(z*10000,real(Er),'b'); hold on
        % plot(z*10000,imag(Er),'b--'); hold off
        % xlim(XL*1e4)
        % ylim(YL)
        % xlabel('z(\mum)')
        % ylabel('E_r')
        % legend('\Re','\Im')
        % 
        % hold off
        % subplot(3,1,3);
        % plot(time*1e12,real(InputL),'r'); hold on
        % plot(time*1e12,real(OutputR),'g'); 
        % plot(time*1e12,real(InputR),'b');
        % plot(time*1e12,real(OutputL),'m');
        % xlim([0,Nt*dt*1e12])
        % ylim(YL)
        % xlabel('time(ps)')
        % ylabel('0')
        % legend('Left Input','Right Output', 'Right Input', 'Left Output', 'Location', 'east')
        % hold off
        
        %Milestone 6 plots
        subplot(3,4,1)
        plot(real(Ef));hold on
        plot(imag(Ef));
        xlabel('z(\mum)')
        ylabel('Ef')
        legend('R','i','east')
        hold off
        subplot(3,4,2)
        plot(z*10000,real(N));
        xlim(XL*1e4)
        ylim([1e18,5e18])
        xlabel('z(\mum)')
        ylabel('N')
        hold off
        subplot(3,4,[5,6])
        plot(time*1e12,real(Nave));
        xlim([0,Nt*dt*1e12])
        ylim()
        xlabel('time(ps)')
        ylabel('Nav')
        hold off
        subplot(3,4,[9,10]);
        plot(time*1e12,real(InputL),'r'); hold on
        plot(time*1e12,real(OutputR),'g'); 
        plot(time*1e12,real(InputR),'b');
        plot(time*1e12,real(OutputL),'m');
        xlim([0,Nt*dt*1e12])
        ylim(YL)
        xlabel('time(ps)')
        ylabel('0')
        legend('Left Input','Right Output', 'Right Input', 'Left Output', 'Location', 'east')
        hold off
        pause(0.01)
    end
    Efp = Ef;
    Erp = Er;
    Pfp = Pf;
    Prp = Pr;
end

%Milestone 7 plots
fftEf = fftshift(fft(Ef));
fftEr = fftshift(fft(Er));
omega = fftshift(wspace(time));

subplot(3,4,[3,4])
plot(time*1e12,real(OutputR),'r');hold on
plot(time*1e12,real(InputL),'g');
xlim([0,Nt*dt*1e12])
ylim(YL)
xlabel('time(ps)')
ylabel('0utput')
legend('Output','Input','east')
hold off
subplot(3,4,[7,8])
%plot(omega, unwrap(angle(fftEf))); hold on
%plot(omega, unwrap(angle(fftEr)));
% xlim()
% ylim()
xlabel('GHz')
ylabel('20 log|E|')
hold off
subplot(3,4,[11,12])
plot(1);
% xlim()
% ylim()
xlabel('GHz')
ylabel('phase (E)')
hold off

% %Milestone 2 Taking FFT of the output and input 
% fftOutput = fftshift(fft(OutputR));
% fftInput = fftshift(fft(InputL));
% %Getting the vector of frequencies that are based of time
% omega = fftshift(wspace(time));
% figure('name', 'FFT')
% subplot(3,1,1)
% plot(time*1e12,real(InputL),'r'); hold on
% plot(time*1e12,real(OutputR),'g'); 
% plot(time*1e12,real(InputR),'b'); hold on
% plot(time*1e12,real(OutputL),'m');
% legend('Left Input','Right Output', 'Right Input', 'Left Output', 'Location', 'east')
% xlabel('time(ps)')
% ylabel('E')
% hold off
% subplot(3,1,2)
% plot(omega, 20*log(abs(fftOutput))); hold on 
% plot(omega, 20*log(abs(fftInput)));
% legend('Output', 'Input','east');
% xlim([-0.5E14,0.5E14])
% xlabel('THz*10')
% ylabel('20log(|E|)')
% hold off
% subplot(3,1,3)
% plot(omega, unwrap(angle(fftOutput))); hold on%Unwrap goes from 0-360-0 to 0-360-720
% plot(omega, unwrap(angle(fftInput)));
% xlabel('THz')
% ylabel('phase (E)')
% legend('Output', 'Input','east');
% hold off
