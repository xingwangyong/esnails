clear

underSampling = 3; % for fully sampled (underSampling=1) data, espirit in bart would fail (segmentation fault), while espirit matlab toolbox always works.
crop_sz = 32;
llrMaxIter = 100;
mat_fname = 'kdata_ktraj.mat';



%%
wind1 = hamming(crop_sz);
wind2 = wind1(:)*wind1(:).';


load(mat_fname)


tic
%%%%%%%%%%%%%%%%%%%%%%% 0 gridding %%%%%%%%%%%%%%%%%%%%%%%          
ksp_crop_cart = zeros(crop_sz,crop_sz ,Nshots/underSampling,numCha);
ksp_crop_cart_wind = ksp_crop_cart;
ksp_full_cart = zeros(N,N,Nshots/underSampling,numCha);
for iter_shots = 1:Nshots/underSampling
    dcf = voronoidens( ksp_traj_bart_tmp(1,:,iter_shots)+1i*ksp_traj_bart_tmp(2,:,iter_shots) );
    dcf = dcf.';
    kgridded = zeros(N,N,numCha);
    for iter_chan = 1:numCha
        kgridded(:,:,iter_chan) = gridmat(ksp_traj_0p5_tmp(1,:,iter_shots)+1i*ksp_traj_0p5_tmp(2,:,iter_shots),ksp_bart_tmp(:,:,iter_shots,iter_chan),dcf,N);
    end
    ksp_crop_cart_wind(:,:,iter_shots,:) = wind2.*bart(sprintf('resize -c 0 %d 1 %d',crop_sz,crop_sz),kgridded);
    ksp_crop_cart(:,:,iter_shots,:) = 1.*bart(sprintf('resize -c 0 %d 1 %d',crop_sz,crop_sz),kgridded);
    ksp_full_cart(:,:,iter_shots,:) = kgridded;
end


%%%%%%%%%%%%%%%%%%%%%%% 1 snails %%%%%%%%%%%%%%%%%%%%%%%
ksp_zeropad_cart_wind = bart(sprintf('resize -c 0 %d 1 %d',N,N),ksp_crop_cart_wind);
img_zpad_wind = bart('fft -i -u 3 ', ksp_zeropad_cart_wind);
tmpimg = reshape(img_zpad_wind, N,N,1,Nshots/underSampling,numCha);
tmpimg = permute(tmpimg, [1 2 5 4 3]);
sens_from_lowres = tmpimg;

params = [];
params.N = N;
params.Nshots = Nshots/underSampling;
params.sens = img_zpad_wind; % no squeeze(), use bart format here!
params.traj = ksp_traj_bart_tmp;

encoding_func = @apply_sense_spi_shots_combined;
A_for = @(in)encoding_func(in,'notransp',params);
A_adj = @(in)encoding_func(in,'transp',params);
AHA = @(in) A_adj(A_for(in));
ksp_adj = A_adj(ksp_bart_tmp);
res2 = symmlq(AHA,ksp_adj(:));

img_sn = mean(reshape(res2,N,N,1/1),3);
im_old_snails = img_sn;


%%%%%%%%%%%%%%%%%%%%%%% 2 esnails %%%%%%%%%%%%%%%%%%%%%%%
% estimate combined sensitivity map, esnail
ksp_zeropad_cart = bart(sprintf('resize -c 0 %d 1 %d',N,N),ksp_crop_cart);
img_zpad = bart('fft -i -u 3 ', ksp_zeropad_cart);
ksp_zeropad_cart2 = reshape(ksp_zeropad_cart,N,N,1,Nshots/underSampling*numCha);
ecal_spiral = bart('ecalib -d5  -m 1 ',ksp_zeropad_cart2);
sens = ecal_spiral;
sens_from_ecalib = sens;

params = [];
params.N = N;
params.Nshots = Nshots/underSampling;
params.sens = reshape(sens,N,N,Nshots/underSampling,numCha); % no squeeze(), use bart format here!
params.traj = ksp_traj_bart_tmp;

encoding_func = @apply_sense_spi_shots_combined;
A_for = @(in)encoding_func(in,'notransp',params);
A_adj = @(in)encoding_func(in,'transp',params);
AHA = @(in) A_adj(A_for(in));
ksp_adj = A_adj(ksp_bart_tmp);
res2 = symmlq(AHA,ksp_adj(:));

img_new = mean(reshape(res2,N,N,1/1),3);
im_esnails = img_new;


%%%%%%%%%%%%%%%%%%%%%%% 3 llr %%%%%%%%%%%%%%%%%%%%%%%
ksp_shotllr = permute(ksp_bart_tmp, [1 2 6 4 5 3]);
ksp_traj_shotllr =   permute(ksp_traj_bart_tmp, [1 2 6 4 5 3]);
cmd = sprintf('pics -S -R L:7:7:%f -s 1e-2  -i %d -t', 0.001, llrMaxIter);
tmp_img = bart(cmd,ksp_traj_shotllr,ksp_shotllr,sens_flip21);
img = mean(abs(squeeze(tmp_img)),3);
im_llr = img;

toc

%% 
im_esnails = rot90(im_esnails,3);
im_old_snails = rot90(im_old_snails,3);
im_llr = rot90(im_llr,3);
sens_from_ecalib = rot90(sens_from_ecalib,3);

figure;imshow(abs(im_old_snails),[]);
figure;imshow(abs(im_esnails),[]);
figure;imshow(abs(im_llr),[]);
