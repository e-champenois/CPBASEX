%% Load Integrals
gData = loadG('../gData/G_r512_k64_l4.h5',1);
%% Load Data
h_fid = fopen('high_counts.bin','r');
l_fid = fopen('low_counts.bin','r');
f_fid = fopen('fancy.bin','r');
data = fread(h_fid,[1023,1023],'double'); % load one image into a 2D array
data(:,:,2) = fread(l_fid,[1023,1023],'double'); % stack second image to make a 3D array
data(:,:,3) = fread(f_fid,[1023,1023],'double');
fclose(h_fid);
fclose(l_fid);
fclose(f_fid);
%% Do Inversion
x0 = 512;
y0 = 512;
fold = resizeFolded(foldQuadrant(data,x0,y0,[1,1,1,1]),512);
%%
out = pbasex(fold,gData,1);
%% Plot Inversion Data
figure;
subplot(2,1,1)
plot(out.E,out.IE(:,1:2)*diag(1./mean(out.IE(:,1:2)))) % normalize by image and plot
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
imagesc(data(:,:,3)); % image raw data
title('Raw Image')
subplot(1,3,2)
imagesc(out.recon(:,:,3)) % image fit data
title('Non-Inverted Fit')
subplot(1,3,3)
imagesc(out.inv(:,:,3)) % image inverted fit data
title('Inverted Fit')
caxis([0,1e2])