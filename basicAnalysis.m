%% Calculate Integrals
saveGData;
%% Load Integrals
gData = load('G_r512_k256_l4_s2.mat');
%% Load Data
imagePath = 'sample_highcounts.dat'; % 'sample_lowcounts.dat'
data = load(imagePath);
%% Do Inversion
x0 = 512;
y0 = 512;
fold = resizeFolded(foldQuadrant(data,x0,y0),512);
out = pbasex(fold,gData,1);
%% Plot Inversion Data
alpha = 3.8e-5;
E = alpha*out.r.^2;
figure;
plot(E,out.Ir/max(out.Ir))
hold on
plot(E,out.betas)
xlabel('Energy (eV)')
ylabel('Intensity (normalized), Beta-value')
legend('Intensity','beta-2','beta-4')
title('Sample CPBASEX Inversion Data')
figure;
subplot(1,3,1)
imagesc(data);
title('Raw Image')
subplot(1,3,2)
imagesc(out.recon)
title('Non-Inverted Fit')
subplot(1,3,3)
imagesc(out.inv)
title('Inverted Fit')