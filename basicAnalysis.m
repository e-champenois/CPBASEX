%% Calculate Integrals
saveGData;
%% Load Integrals
<<<<<<< HEAD
gData = loadG('/Users/elio/Google Drive/Python/CPBASEX/G_r512_k256_l4.0_s1.4.h5');
%% Load Data
high_counts = 'sample_highcounts.dat';
low_counts = 'sample_lowcounts.dat';
data = load(high_counts); % load one image into a 2D array
data(:,:,2) = load(low_counts); % stack second image to make a 3D array
=======
gData = load('G_r512_k256_l4_s2.mat');
%% Load Data
imagePath = 'sample_highcounts.dat'; % 'sample_lowcounts.dat'
data = load(imagePath);
>>>>>>> 629229caabfed460d2a5fd6b8e5d81782bcc3051
%% Do Inversion
x0 = 512;
y0 = 512;
fold = resizeFolded(foldQuadrant(data,x0,y0),512);
<<<<<<< HEAD
%%
out = pbasex(fold,gData,1);
%% Plot Inversion Data
figure;
subplot(2,1,1)
plot(E,out.IE*diag(1./mean(out.IE))) % normalize by image and plot
xlabel('Energy (eV)')
ylabel('Intensity (normalized)')
legend('High Counts','Low Counts')
title('Sample CPBASEX Inversion Data')
subplot(2,1,2)
plot(E,out.betas(:,:,1)) % plot betas versus energy for the first image
xlabel('Energy (eV)')
ylabel('Beta-value')
legend('beta-2','beta-4')
figure;
subplot(1,3,1)
imagesc(data(:,:,1)); % first image raw data
title('Raw Image')
subplot(1,3,2)
imagesc(out.recon(:,:,1)) % first image fit data
title('Non-Inverted Fit')
subplot(1,3,3)
imagesc(out.inv(:,:,1)) % first image inverted fit data
=======
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
>>>>>>> 629229caabfed460d2a5fd6b8e5d81782bcc3051
title('Inverted Fit')