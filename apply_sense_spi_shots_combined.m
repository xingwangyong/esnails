function res = apply_sense_spi_shots_combined(in, tflag, params )
% apply sense for each shot, i.e. each shot gives an image


N = params.N;
N = N(1);
Nshots = params.Nshots;
sens = params.sens;
traj = params.traj;

sz = size(sens);
num_chan = sz(end);
readout_len = size(traj,2);


if strcmp(tflag,'transp')
    
    % Transposed SENSE operator:
    % IFFT coil k-space, multiply by conjugate of coil sensitivities, then
    % sum across channels

    
    dim = sprintf('-d%d:%d:%d',N,N,1);
    kspace = reshape(in, [1,readout_len, Nshots, num_chan]);
    img = zeros(N,N, Nshots,num_chan);
    for iter_shot = 1:Nshots
        img(:,:,iter_shot,:) = bart(['nufft ',dim,' -a'], traj(:,:,iter_shot),kspace(:,:,iter_shot,:));
    end

    Res = sum( sum(img .* conj(sens), 4) , 3);
    
    res = Res(:);
    % res = Res;
    
elseif strcmp(tflag,'notransp')
    
    % Forward SENSE operator:
    % multiply by coil sensitivities, take undersampled FFT

    img = reshape(in, N,N,1) .* sens;
    kspace = zeros([1,readout_len,Nshots, num_chan]);
    for iter_shot = 1:Nshots
        kspace(:,:,iter_shot,:) = bart('nufft', traj(:,:,iter_shot), img(:,:,iter_shot,:));
    end
    
    res = kspace(:);
    % res = kspace;
else
    error('Syntax error, un-recognized 2nd param: %s', tflag);
end
end