%Milestone 4
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
beta_i = 8;
beta_r = 0;

InputParasL.E0=100e5;     %Amplitude?
InputParasL.we = 0;   %Frequency for modulation (1e13)
InputParasL.t0 = 2e-12;
InputParasL.wg = 5e-13; %5e-13 
InputParasL.phi = 0;
InputParasR = 0;

n_g = 3.5; 
vg = c_c/n_g*1e2;       % TWM cm/s group velocity
Lambda = 1550e-9;

plotN = 10;

L = 1000e-6*1e2;    %cm
XL = [0,L];
YL =[1*-InputParasL.E0,1*InputParasL.E0]; %vertical scale

Nz =500;            
dz =L/(Nz-1);
dt = dz/vg;
fsync = dt*vg/dz;

Nt =floor(2*Nz);        %designates length of simulation
tmax = Nt*dt;
t_L = dt*Nz;               % time to travel length

z = linspace(0,L,Nz).';    % Nz points, nz-1 segments
time = nan(1,Nt);
InputL = nan(1,Nt);        % create arrays of "NaN" 1 x Nt
InputR = nan(1,Nt);
OutputL = nan(1,Nt);
OutputR = nan(1,Nt);

%Milestone 3
kappa0 = 100; %100
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
LGain = 0.1;

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
beta = ones(size(z))*(beta_r+1i*beta_i); %Initializing Beta
exp_det = exp(-1i*dz*beta);

%Create all initial graphs
figure('name', 'Fields')
subplot(3,1,1)
plot(z*10000,real(Ef),'r');
hold off
xlabel('z(\mum)')
ylabel('E_f')
subplot(3,1,2)
plot(z*10000,real(Er),'b');
xlabel('z(\mum)')
ylabel('E_r')
hold off
subplot(3,1,3)
plot(time*1e12,real(InputL),'r'); hold on
plot(time*1e12,real(OutputR),'r--'); 
plot(time*1e12,real(InputR),'b'); hold on
plot(time*1e12,real(OutputL),'b--');
xlabel('time(ps)')
ylabel('E')

hold off

%Keep updating graphs and recalculating propegation 
for i = 2:Nt
    t = dt*(i-1);
    time(i) = t;
    Pf(1) = 0;
    Pf(Nz) = 0;
    Pr(1) = 0;
    Pr(Nz) = 0;
    Cw0 = -LGamma + 1i*Lw0;
    
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

    %Milestone 1
    OutputR(i) = Ef(Nz)*(1-RR);   
    OutputL(i) = Er(1)*(1-RL);    %Reflecting

    if mod(i,plotN) == 0
        subplot(3,1,1)
        plot(z*10000,real(Ef),'r'); hold on
        plot(z*10000,imag(Ef),'r--'); hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_f')
        legend('\Re','\Im')
        hold off
        subplot(3,1,2)
        plot(z*10000,real(Er),'b'); hold on
        plot(z*10000,imag(Er),'b--'); hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_r')
        legend('\Re','\Im')

        hold off
        subplot(3,1,3);
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

%Milestone 2 Taking FFT of the output and input 
fftOutput = fftshift(fft(OutputR));
fftInput = fftshift(fft(InputL));
%Getting the vector of frequencies that are based of time
omega = fftshift(wspace(time));
figure('name', 'FFT')
subplot(3,1,1)
plot(time*1e12,real(InputL),'r'); hold on
plot(time*1e12,real(OutputR),'g'); 
plot(time*1e12,real(InputR),'b'); hold on
plot(time*1e12,real(OutputL),'m');
legend('Left Input','Right Output', 'Right Input', 'Left Output', 'Location', 'east')
xlabel('time(ps)')
ylabel('E')
hold off
subplot(3,1,2)
plot(omega, 20*log(abs(fftOutput))); hold on 
plot(omega, 20*log(abs(fftInput)));
legend('Output', 'Input','east');
xlim([-0.5E14,0.5E14])
xlabel('THz*10')
ylabel('20log(|E|)')
hold off
subplot(3,1,3)
plot(omega, unwrap(angle(fftOutput))); hold on%Unwrap goes from 0-360-0 to 0-360-720
plot(omega, unwrap(angle(fftInput)));
xlabel('THz')
ylabel('phase (E)')
legend('Output', 'Input','east');
hold off
