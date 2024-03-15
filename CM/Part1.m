I_s = 0.01e-12;
I_b = 0.1e-12;
V_b = 1.3;
G_p = 0.1;
j = (0.7+1.95)/200;
V = -1.95:j:0.7;

I = I_s*(exp(1.2*V/0.025)-1) + G_p*V - I_b*(exp(-1.2*(V+V_b)/0.025)-1);

I_rand = I_s*(exp(1.2*V/0.025)-1) + G_p*V - I_b*(exp(-1.2*(V+V_b)/0.025)-1);
I_rand = I_rand + 0.2*I.*randn(1,201);

figure(1)
title('polyfit()')
subplot(2,4,[1,2])
plot(V,I,'b');hold on
fit1 = polyfit(V,I,4);
fit2 = polyfit(V,I,8);
plot(V,polyval(fit1,V),'r--');
plot(V,polyval(fit2,V)),'m--';
ylabel("I");
xlabel("V");
legend('Without Noise','Fit 4th','Fit 8th');
hold off
subplot(2,4,[3,4])
semilogy(V,abs(I)); hold on
semilogy(V,abs(polyval(fit1,V)),'r--');
semilogy(V,abs(polyval(fit2,V)),'m--');
legend('Without Noise','Fit 4th','Fit 8th');
hold off
subplot(2,4,[5,6])
plot(V,I_rand,'b'); hold on
fit5 = polyfit(V,I,4);
fit6 = polyfit(V,I,8);
plot(V,polyval(fit5,V),'r--'); 
plot(V,polyval(fit6,V),'m--');
ylabel("I");
xlabel("V");
legend('With Noise','Fit 4th','Fit 8th');
hold off
subplot(2,4,[7,8]) 
semilogy(V,abs(I_rand)); hold on
semilogy(V,abs(polyval(fit5,V)),'r--');
semilogy(V,abs(polyval(fit6,V)),'m--');
legend('With Noise','Fit 4th','Fit 8th');
hold off

fo1 = fittype('A.*(exp(1.2*x/25e-3)-1) + 0.1.*x - C*(exp(1.2*(-(x+1.3))/25e-3)-1)');
ff1 = fit(V',I',fo1);
fo2 = fittype('A.*(exp(1.2*x/25e-3)-1) + B.*x - C*(exp(1.2*(-(x+1.3))/25e-3)-1)');
ff2 = fit(V',I',fo2);
fo3 = fittype('A.*(exp(1.2*x/25e-3)-1) + B.*x - C*(exp(1.2*(-(x+D))/25e-3)-1)');
ff3 = fit(V',I',fo3);
figure(2)
title('polyfit()')
subplot(1,3,1)
plot(V,ff1(V),'m');hold on
plot(V,I,'b--');
ylim([-4 4]);
ylabel("I");
xlabel("V");
legend('Fit1','I','east');
hold off
subplot(1,3,2)
plot(V,ff2(V),'g'); hold on
plot(V,I,'b--');
ylim([-4 4]);
ylabel("I");
xlabel("V");
legend('Fit2','I','east');
hold off
subplot(1,3,3)
plot(V,ff3(V),'r'); hold on
plot(V,I,'b');
ylim([-4 4]);
ylabel("I");
xlabel("V");
legend('Fit3','I','east');
hold off

%Nural Network
inputs = V.';
targets = I.';
hiddenLayerSize = 10;
net = fitnet(hiddenLayerSize);
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;
[net,tr] = train(net,inputs,targets);
outputs = net(inputs);
errors = gsubtract(outputs,targets);
performance = perform(net,targets,outputs);
view(net)
Inn = outputs;
figure(3)
title('network')
plot(V,I,'b'); hold on
plot(V,Inn,'r--');
ylabel("I");
xlabel("V");
legend('I','Inn','east');
hold off

