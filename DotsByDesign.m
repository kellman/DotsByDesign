clear all; close all; clc;
warning off;
%% loading image data
mkdir('Results');

filename = 'Rooster';
img = importdata([filename '.jpg']);

% img = phantom();
% img = importdata('cameraman.jpg');
img = mean(img,3);
img = 1 - img./max(img(:));
[Nx,Ny] = size(img);

% size of averaging kernel
nb = 36;

% lineary maps distinct averaged blocks to dots with radius between rmin
% and rmax.
rmax = 10;
rmin = 2;
%% Compute the radii

R = blkproc(img, [nb nb], 'mean2');
R = conv2(R,1/4*ones(2,2),'same');
R = R./max(R(:));
R = (rmax - rmin)*R + rmin;

figure,
subplot 131, imagesc(img);
axis image;
title('Original Image');

subplot 132, imagesc(R);
axis image; colormap gray;
title('Averaged & Smoothed');

%% Drawing dots

subplot 133;
for ii = 1:size(R,2)
    for jj = 1:size(R,1)
        plot(ii, size(R,1) - jj, '.r', 'MarkerSize',R(jj,ii));
        hold on;
    end
end
axis off; axis([1 size(R,2) 1 size(R,1)]); axis image; 
title('Dot - ify');

%% Final
% on white background

figure;
for ii = 1:size(R,2)
    for jj = 1:size(R,1)
        plot(ii, size(R,1) - jj, '.r', 'MarkerSize',R(jj,ii));
        hold on;
    end
end
axis off; axis([1 size(R,2) 1 size(R,1)]); axis image;
set(gcf,'color','w');
savefilename = ['Results/dbd_' filename '.png'];
saveas(gcf,savefilename);

%% extension (1)
% random colors on white background

figure;
for ii = 1:size(R,2)
    for jj = 1:size(R,1)
        p = plot(ii,size(R,1) - jj, '.', 'MarkerSize',R(jj,ii));
        set(p,'Color',rand(1,3));
        hold on;
    end
end
axis off; axis([1 size(R,2) 1 size(R,1)]); axis image;
set(gcf,'color','w');
savefilename = ['Results/dbdr_' filename '.png'];
saveas(gcf,savefilename);


%% extension (2)
% red and blue color offset gridded dots

figure;
hold on;
for ii = 1:size(R,2)
    for jj = 1:size(R,1)
        plot(ii, size(R,1) - jj, '.r', 'MarkerSize',R(jj,ii));
        plot(ii + 1/4,size(R,1) - jj + 1/4, '.b', 'MarkerSize',R(jj,ii));
    end
end
axis off; axis([1 size(R,2) 1 size(R,1)]); axis image;
set(gcf,'color','w');
savefilename = ['Results/dbd2c_' filename '.png'];
saveas(gcf,savefilename);

%% extension (3)
% edges are emphasized
kernel = [[1 -1];[-1 1]];
gradimg = conv2(img,kernel,'same');

R = blkproc(gradimg, [nb nb], 'mean2');
R = conv2(R,1/4*ones(2,2),'same');
R = R - min(R(:));
R = R./max(R(:));
R = (rmax - rmin)*R + rmin;

figure;
for ii = 1:size(R,2)
    for jj = 1:size(R,1)
        plot(ii, size(R,1) - jj, '.b', 'MarkerSize',R(jj,ii));
        hold on;
    end
end
axis off; axis([1 size(R,2) 1 size(R,1)]); axis image;
set(gcf,'color','w');
savefilename = ['Results/dbde_' filename '.png'];
saveas(gcf,savefilename);