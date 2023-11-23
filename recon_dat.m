clear

underSampling = 3; % for fully sampled (underSampling=1) data, espirit in bart would fail (segmentation fault), while espirit matlab toolbox always works.
crop_sz = 32;
llrMaxIter = 100;

sens_fname = '/autofs/cluster/berkin/xingwang/Syncbb/epti_proj/data/sens/20231103_bay4_sens.mat';
data_file_path = '/autofs/space/marduk_001/users/berkin/2023_11_03_bay4_spider_invivo/meas_MID00282_FID03826_b1k_ms6_dds_snails_20slc.dat';
seq_file_path = '/autofs/cluster/berkin/xingwang/Syncbb/seq-for-experiment/20231103_bay4_invivo/b1k_ms6_dds_snails_20slc.seq';
traj_file_path = '/autofs/cluster/berkin/xingwang/Syncbb/seq-for-experiment/20231103_bay4_invivo/b1k_ms6_dds_snails_20slc.mat';


%% sanity check
twix_obj_spiral = mapVBVD(data_file_path);
seq = mr.Sequence(mr.opts('B0',2.89));
seq.read(seq_file_path,'detectRFuse');
twix_sig = twix_obj_spiral{end}.hdr.Dicom.tSequenceVariant;
seq_sig = seq.signatureValue;
if ~isequal(twix_sig,seq_sig) && ~contains(twix_sig,seq_sig) % ideally these two should be equal, but somehow sometimes the former has some extra trailing chars
    error('Signature mismatch')
end


%% load
Nshots = seq.getDefinition('NumShots');
Nslices = seq.getDefinition('NumSlices');
N = seq.getDefinition('MatrixSize');
num_bval = seq.getDefinition('num_bval');
num_bdir = seq.getDefinition('num_bdir');
numSamplesPerShot = seq.getDefinition('numSamplesPerShot');
bTable = seq.getDefinition('bTable');
bTable = reshape(bTable,[],3);
Nechoes = seq.getDefinition('Nechoes');
num_volumes = seq.getDefinition('num_volumes');

tmpStuct = load(traj_file_path);
ktraj_adc_splited = tmpStuct.ktraj_adc_splited;
ktraj_adc_splited(1,:) = rescale(ktraj_adc_splited(1,:),-0.5*N,0.5*N);
ktraj_adc_splited(2,:) = rescale(ktraj_adc_splited(2,:),-0.5*N,0.5*N);
ktraj_adc_splited(3,:) = 0;
ktraj_bart_all_slices = repmat(ktraj_adc_splited,1,1,1,Nslices,1,num_volumes);

ktraj_bart_all_slices_scale0p5 = ktraj_adc_splited;
ktraj_bart_all_slices_scale0p5 = rescale(ktraj_bart_all_slices_scale0p5,-0.5,0.5);
ktraj_bart_all_slices_scale0p5 = repmat(ktraj_bart_all_slices_scale0p5,1,1,1,Nslices,1,num_volumes);    


%% load pre-calculated sensitivity map
tmpStruct = load(sens_fname);
sens_sequential = tmpStruct.sens;
sens_sequential_flip21 = flip(flip(sens_sequential,2),1);


%% calculate slice indices to undo slice interleaving, after undo, the 1st to last slice is from feet to head
% for spiral
slice_indices_odd = 1:2:Nslices;
slice_indices_even = 2:2:Nslices;
slice_indices = [slice_indices_odd,slice_indices_even];
[~,ind_undo_interleave] = sort(slice_indices);
ind_undo_interleave = flip(ind_undo_interleave); % TODO, I don't know why the spiral acquisition is from head to feet
ind_do_interleave = slice_indices;


%% 
if ~exist('twix_obj_spiral','var')
    twix_obj_spiral = mapVBVD(data_file_path);
end
if iscell(twix_obj_spiral)
    data = twix_obj_spiral{end}.image;
else
    data = twix_obj_spiral.image;
end
data.flagRemoveOS = false;
data.flagDoAverage = true;
data.flagIgnoreSeg = true;

kspace = data.unsorted();
kspace = permute(kspace,[1,3,2]);
kspace = reshape(kspace, [size(kspace, 1), size(kspace, 2), 1, size(kspace, 3)]);
kspace_svd = squeeze(kspace);

numCha = size(kspace_svd,3);
kspace_splitted = reshape(kspace_svd,numSamplesPerShot,Nechoes,Nslices,Nshots,num_volumes,numCha);
kspace_sequential = kspace_splitted(:,:,ind_undo_interleave,:,:,:,:);
im_esnails = zeros(N,N,Nslices,num_volumes,Nechoes);
im_llr = im_esnails; % llr, locally low rank
im_old_snails = im_esnails;
sens_from_ecalib = zeros(N,N,Nslices,numCha*Nshots/underSampling,num_volumes);
sens_from_lowres = zeros(N,N,Nslices,numCha,Nshots/underSampling,num_volumes);
wind1 = hamming(crop_sz);
wind2 = wind1(:)*wind1(:).';


%%
tic
for iter_echo = 1:Nechoes
    parfor iter_slice = 1:Nslices
        fprintf('%d/%d\n',iter_slice,Nslices);
        for iter_bvol = 1:num_volumes
            ksp_traj_bart = ktraj_bart_all_slices(:,:,iter_echo,iter_slice,:,iter_bvol);
            ksp_traj_bart = squeeze(ksp_traj_bart);
            ksp_traj_0p5 = ktraj_bart_all_slices_scale0p5(:,:,iter_echo,iter_slice,:,iter_bvol);
            ksp_traj_0p5 = squeeze(ksp_traj_0p5);
            ksp_bart = squeeze( kspace_sequential(:,iter_echo,iter_slice,:,iter_bvol,:) );
            ksp_bart = reshape(ksp_bart, [1,size(ksp_bart)]);

            %%%%%%%%%%%%%%%%%%%%%%% 0 gridding %%%%%%%%%%%%%%%%%%%%%%%
            ksp_bart_tmp = ksp_bart(:,:,1:underSampling:Nshots,:);
            ksp_traj_bart_tmp = ksp_traj_bart(:,:,1:underSampling:Nshots);
            ksp_traj_0p5_tmp = ksp_traj_0p5(:,:,1:underSampling:Nshots);
            ksp_crop_cart = zeros(crop_sz,crop_sz ,Nshots/underSampling,numCha);
            ksp_crop_cart_wind = ksp_crop_cart;
            img_no_phase_shots = zeros(N,N,Nshots/underSampling,numCha);
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
            sens_from_lowres(:,:,iter_slice,:,:,iter_bvol) = tmpimg;

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
            im_old_snails(:,:,iter_slice,iter_bvol,iter_echo) = img_sn;


            %%%%%%%%%%%%%%%%%%%%%%% 2 esnails %%%%%%%%%%%%%%%%%%%%%%%
            % estimate combined sensitivity map, esnail
            ksp_zeropad_cart = bart(sprintf('resize -c 0 %d 1 %d',N,N),ksp_crop_cart);
            img_zpad = bart('fft -i -u 3 ', ksp_zeropad_cart);
            ksp_zeropad_cart2 = reshape(ksp_zeropad_cart,N,N,1,Nshots/underSampling*numCha);
            ecal_spiral = bart('ecalib -d5  -m 1 ',ksp_zeropad_cart2);
            sens = ecal_spiral;
            sens_from_ecalib(:,:,iter_slice,:,iter_bvol) = sens;

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
            im_esnails(:,:,iter_slice,iter_bvol,iter_echo) = img_new;


            %%%%%%%%%%%%%%%%%%%%%%% 3 llr %%%%%%%%%%%%%%%%%%%%%%%
            ksp_shotllr = permute(ksp_bart_tmp, [1 2 6 4 5 3]);
            ksp_traj_shotllr =   permute(ksp_traj_bart_tmp, [1 2 6 4 5 3]);
            cmd = sprintf('pics -S -R L:7:7:%f -s 1e-2  -i %d -t', 0.001, llrMaxIter);
            tmp_img = bart(cmd,ksp_traj_shotllr,ksp_shotllr,sens_sequential_flip21(:,:,iter_slice,:));
            img = mean(abs(squeeze(tmp_img)),3);
            im_llr(:,:,iter_slice,iter_bvol,iter_echo) = img;
        end
    end
end
toc

%% 
im_esnails = rot90(im_esnails,3);
im_old_snails = rot90(im_old_snails,3);
im_llr = rot90(im_llr,3);
sens_from_ecalib = rot90(sens_from_ecalib,3);

figure;imshow( abs(im_llr(:,:,1,1,1)),[] );