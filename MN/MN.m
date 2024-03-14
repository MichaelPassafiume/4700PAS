G1 = 1/1;
G2 = 1/2;
G3 = 1/10;
G4 = 1/0.1;
G0 = 1/1000;
G5 = 1/1000.1;
a = 100;

%capacitor
%c = 0.25;
c = normrnd(0.25,0.05);


l = 0.2;
V1 = 0;
V2 = 0;
V3 = 0;
V4 = 0;
V5 = 0;
Iin = 0;
IL = 0;
IV = 0;
w = 1;

%DC
%Vin = -10;
%AC
Vin = 1;


%V = [V1 V2 V3 V4 V5 Iin IL IV];

G = [1 0 0 0 0 0 0 0;
    -G1 (G1+G2) 0 0 0 0 1 0;
     G1 -G1 0 0 0 1 0 0;
     0 0 -a*G3 1 0 0 0 0;
     0 -1 1 0 0 0 0 0;
     0 0 G3 0 0 0 1 0;
     0 0 0 G4 -G4 0 0 1;
     0 0 0 -G4 (G4+G0) 0 0 0;
     0 0 (-a*G3*G5) 0 0 0 0 1];
C = [0 0 0 0 0 0 0 0;
    -c c 0 0 0 0 0 0;
    c -c 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 l 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0];
F = [Vin;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0];

Vfive = [];
Vone = [];
Vthree = [];
imatrix = [];
cmatrix = [];

% %DC
% for i = 1:20
%     F(1) = F(1) + 1;
%     Vop = G\F;
%     Vfive(i) = Vop(5);
%     imatrix(i) = F(1);
% end
% figure 
% plot(imatrix,Vfive);

% %AC
% for i = 1:100
%     w = w + 1;
%     Vop = (G + 1i*w*C)\F;
%     Vone(i) = Vop(1);
%     Vthree(i) = Vop(3);
%     Vfive(i) = Vop(5);
%     imatrix(i) = w;
% end
% figure
% subplot(2,1,1)
% plot(imatrix,real(Vthree)); hold on
% plot(imatrix,real(Vfive));
% legend ('three','five','east')
% hold off
% gain = real(Vfive)./real(Vone);
% subplot (2,1,2)
% histogram(gain,14);

%AC 2
for i = 1:100
    c = normrnd(0.25,0.05);
    C = [0 0 0 0 0 0 0 0;
    -c c 0 0 0 0 0 0;
    c -c 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 l 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0];
    cmatrix(i) = c;
    w = 3.14;
    Vop = (G + 1i*w*C)\F;
    Vone(i) = Vop(1);
    Vfive(i) = Vop(5);
    imatrix(i) = i;
end
figure
subplot(2,2,[1,2])
plot(imatrix,imag(Vfive)); hold on
plot(imatrix,real(Vfive));
legend ('imag','real','east')
hold off
gain = real(Vfive)./real(Vone);
subplot(2,2,3)
histogram(cmatrix,14);
subplot (2,2,4)
histogram(gain,14);

