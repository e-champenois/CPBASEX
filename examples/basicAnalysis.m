%% Load Integrals
gData = loadG('G_r512_k128_l4.h5',1);
%% Load Data
high_counts = 'high_counts.bin';
low_counts = 'low_counts.bin';
h_fid = fopen(high_counts,'r');
l_fid = fopen(low_counts,'r');
data = fread(h_fid,[1023,1023],'double'); % load one image into a 2D array
data(:,:,2) = fread(l_fid,[1023,1023],'double'); % stack second image to make a 3D array
fclose(h_fid);
fclose(l_fid);
%% Do Inversion
x0 = 512;
y0 = 512;
fold = resizeFolded(foldQuadrant(data,x0,y0),512);
%%
out = pbasex(fold,gData,1);
%% Plot Inversion Data
figure;
subplot(2,1,1)
plot(out.E,out.IE*diag(1./mean(out.IE))) % normalize by image and plot
xlabel('Energy (eV)')
ylabel('Intensity (normalized)')
legend('High Counts','Low Counts')
title('Sample CPBASEX Inversion Data')
subplot(2,1,2)
plot(out.E,out.betas(:,:,1)) % plot betas versus energy for the first image
xlabel('Energy (eV)')
ylabel('Beta-value')
legend('beta-2','beta-4')
%%
figure;
subplot(1,3,1)
imagesc(data(:,:,1)); % first image raw data
title('Raw Image')
subplot(1,3,2)
imagesc(out.recon(:,:,1)) % first image fit data
title('Non-Inverted Fit')
subplot(1,3,3)
imagesc(out.inv(:,:,1)) % first image inverted fit data
title('Inverted Fit')
caxis([0,1e3])