%Milestone 2
set(0, 'defaultaxesfontsize',20)
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultLineLineWidth',2);
set(0,'Defaultaxeslinewidth',2)

set(0,'DefaultFigureWindowStyle','docked')

c_c = 299792458;            % m/s TWM speed of light
c_eps_0 = 8.8542149e-12;    % F/m vacuum permittivity
c_eps_0_cm = c_eps_0/100;   % F/cm
c_mu_0 = 1/c_eps_0/c_c^2;
c_q = 1.60217653e-19;
c_hb = 1.05457266913e-34;                % Dirac constant
c_h = c_hb*2*pi;

RL = 0.9i;  %Milestone 1
RR = 0.9i;  %Reflective Efficiency

InputParasL.E0=1e5;     %Amplitude?
InputParasL.we = 0;
InputParasL.t0 = 2e-12;
InputParasL.wg = 5e-13;
InputParasL.phi = 45;
InputParasR = 0;

n_g = 3.5; 
vg = c_c/n_g*1e2;       % TWM cm/s group velocity
Lambda = 1550e-9;

plotN = 10;

L = 1000e-6*1e2;    %cm
XL = [0,L];
YL =[-InputParasL.E0,InputParasL.E0];

Nz =500;            
dz =L/(Nz-1);
dt = dz/vg;
fsync = dt*vg/dz;

Nt =floor(6*Nz);        %designates length of simulation
tmax = Nt*dt;
t_L = dt*Nz;               % time to travel length

z = linspace(0,L,Nz).';    % Nz points, nz-1 segments
time = nan(1,Nt);
InputL = nan(1,Nt);        % create arrays of "NaN" 1 x Nt
InputR = nan(1,Nt);
OutputL = nan(1,Nt);
OutputR = nan(1,Nt);

Ef = zeros(size(z));       % craete array of 0 
Er = zeros(size(z));

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

    InputL(i) = Ef1(t,InputParasL);
    InputR(i) = ErN(t,0);

    Ef(1) = InputL(i) + RL*Er(1); %Milestone 1
    Er(Nz) = InputR(i) + RR*Ef(Nz);

    Ef(2:Nz) = fsync*Ef(1:Nz-1);
    Er(1:Nz-1) = fsync*Er(2:Nz);

    OutputR(i) = Ef(Nz)*(1-RR);   %Milestone 1
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
end