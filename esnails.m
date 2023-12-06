% before running this sequence on scanner, please check the acoustic
% resonance frequency carefully, see 
% pulseq\matlab\demoUnsorted\gradSpectrum.m for detail

clear 
%% input
% system
sys_type                  = 'prisma'; % prisma, skyra, trio
slew_safety_magrin        = 0.7;
grad_safety_magrin        = 0.9; 
spiral_slew_safety_margin = 0.5; % to reduce PNS. decrease this would not lengthen readout too much
spiral_grad_safety_margin = 0.8;
diff_slew_safety_margin   = 0.4; % to reduce PNS. decrease this would not lengthen TE too much
diff_grad_safety_margin   = 0.97;

% geometry
fov                    = 224e-3; 
Nx                     = 224;
low_res_img_resolution = 6e-3; % determines many portion of k-space center is fully sampled
sliceThickness         = 3e-3;
Nslices                = 1;

%
Narms             = 6; % number of spiral interleaves
Nshots            = Narms;
TR                = 4000e-3;
TE                = 38e-3;
bvalues           = [1]*1e3; % b0 is always acquired
adcDwell          = 2e-6; % the final dwell may slightly deviate from this
fsSpoilerDuration = 5000e-6; % not use the shortest possible duration, to avoid PNS
bTable            = [0 0 0;0 0 1];

% misc
rot_angle = 2*pi/Nshots;
flipAngle = 90;
tRFex     = 3e-3;
tRFref    = 3e-3;
fspS      = 0.5;
MAX_ADCS  = 100; % should not be changed, this is from interpreter cpp file
result_traj_fname = 'esnails_traj.mat';

% not used
Ndirections = 1;
Nechoes     = 1;
Ndummies    = 0;


%%
if strcmp(sys_type,'prisma')
    physical_slew_max = 200;
    physical_grad_max = 80;
    B0=2.89; % 1.5 2.89 3.0
elseif strcmp(sys_type,'skyra')
    physical_slew_max = 180;
    physical_grad_max = 43;
    B0=2.89; % 1.5 2.89 3.0
elseif strcmp(sys_type,'trio')
    physical_slew_max = 170;
    physical_grad_max = 38;
    B0=2.89; % 1.5 2.89 3.0
else
    error('Undefined')
end
slew_max = physical_slew_max*slew_safety_magrin; % T/m/s
grad_max = physical_grad_max*grad_safety_magrin; % mT/m
sys = mr.opts('MaxGrad',grad_max,'GradUnit','mT/m',...
    'MaxSlew',slew_max,'SlewUnit','T/m/s',...
    'rfRingdownTime', 30e-6, 'rfDeadtime', 100e-6, 'adcDeadTime', 10e-6,'B0',B0);
sys_diff = mr.opts('MaxGrad',physical_grad_max*diff_grad_safety_margin,'GradUnit','mT/m',...
    'MaxSlew',physical_slew_max*diff_slew_safety_margin,'SlewUnit','T/m/s',...
    'rfRingdownTime', 30e-6, 'rfDeadtime', 100e-6, 'adcDeadTime', 10e-6,'B0',B0);

tRefwd    = tRFref+sys.rfRingdownTime+sys.rfDeadTime;
seq=mr.Sequence(sys);      % Create a new sequence object


%%
%%%%%%%%%%%%%%%%%%%% VDS design %%%%%%%%%%%%%%%%%%%%
res_mm                = 1e3*fov/Nx; % millimetter
low_res_img_kmax      = 1/low_res_img_resolution/2;
full_img_kax          = 1/(res_mm*1e-3)/2;
k_ratio_full_sampling = low_res_img_kmax./full_img_kax;
k_ratio2              = max([k_ratio_full_sampling + 0.01, 0.25]);
sampling_locs         = [0,k_ratio_full_sampling,k_ratio2,1];
sampling_dens         = [Narms*fov,Narms*fov,fov,fov];
[k_riv,g_riv,s_riv,time_riv,Ck_riv] = vdSpiralDesign(Narms,0, res_mm,100*sampling_dens,sampling_locs,physical_grad_max*spiral_grad_safety_margin/10,spiral_slew_safety_margin*physical_slew_max/10,1e3*sys.gradRasterTime,[],'pchip');
fprintf('K-space sampling locations: %s\n', num2str(sampling_locs,'%.2f '));
fprintf('K-space sampling densities: %s\n', num2str(sampling_dens/fov,'%.2f '));
seq.setDefinition('k_ratio_full_sampling',k_ratio_full_sampling);

g_mT_m = g_riv(:,1:2)*10; % G/cm to mT/m
grad_std = mr.convert(g_mT_m,'mT/m');
s_T_m_s = s_riv(:,1:2)*10;
slew_std = mr.convert(s_T_m_s,'T/m/s');

L = length(s_riv);
time_array = linspace(0,time_riv,size(g_riv,1));
figure
subplot(2,2,1)
plot(k_riv(:,1), k_riv(:,2)); 
title('k-space');% axis([-6 6 -6 6]);
subplot(2,2,2)
plot(time_array,g_riv(:,1)); 
%axis([0,L,-4.5,4.5]); 
title('gradient waveforms')
hold on
plot(time_array,g_riv(:,2), 'r');
xlabel 'ms'
legend('gx', 'gy', 'Location', 'NorthEast');
subplot(2,2,3)
plot(time_array,(g_riv(:,1).^2 + g_riv(:,2).^2).^0.5, 'r');
xlabel 'ms'
%axis([0 L 0 6]);
title('gradient magnitude')
subplot(2,2,4), plot(time_array(1:end-1),(s_riv(:,1).^2 + s_riv(:,2).^2).^0.5, 'r'); 
xlabel 'ms'
title('slew-rate magnitude');
%%%%%%%%%%%%%%%%%%%% VDS design %%%%%%%%%%%%%%%%%%%%

% gradients and slew rates
% ka = k_std;
% ka = [real(ka); imag(ka)];
ga = grad_std.';
sa = slew_std.';

% [gos, sos]=mr.traj2grad(kopt_smooth);
gos = ga;
sos = sa;

g_for_plot = g_mT_m.';
s_for_plot = s_T_m_s.';
figure;plot([g_for_plot;abs(g_for_plot(1,:)+1i*g_for_plot(2,:))]');title('gradient with abs(grad) constraint');ylabel 'mT/m';
figure;plot([s_for_plot;abs(s_for_plot(1,:)+1i*s_for_plot(2,:))]');title('slew rate with abs(slew) constraint');ylabel 'mT/m/ms';

spiral_grad_shape=gos;


%% fat sat and 180 crusher
warning('OFF', 'mr:restoreShape'); % restore shape is not compatible with spirals and will throw a warning from each plot() or calcKspace() call

% Create fat-sat pulse 
% (in Siemens interpreter from January 2019 duration is limited to 8.192 ms, and although product EPI uses 10.24 ms, 8 ms seems to be sufficient)
sat_ppm=-3.45;
sat_freq=sat_ppm*1e-6*B0*sys.gamma;
rf_fs = mr.makeGaussPulse(110*pi/180,'system',sys,'Duration',8e-3,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq);
gz_fs = mr.makeTrapezoid('z',sys,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm
if (gz_fs.riseTime+gz_fs.fallTime+gz_fs.flatTime) < fsSpoilerDuration
    gz_fs = mr.makeTrapezoid('z',sys,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4,'duration',fsSpoilerDuration);
end

% Create 90  degree slice selection pulse and gradient
[rf90, gs_ex, gzReph] = mr.makeSincPulse(flipAngle*pi/180,'system',sys,'Duration',tRFex,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4);

% Create 180 degree slice refocusing pulse and gradients
[rf180, gz180] = mr.makeSincPulse(pi,'system',sys,'Duration',tRFref,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,'use','refocusing');
GSref = mr.makeTrapezoid('z',sys,'amplitude',gs_ex.amplitude,'FlatTime',tRefwd);

% GS7 is left crusher of 180, GS5 is right crusher for 180. 
% GS4 is slice selection gradient for 180 (only flat)
% These gradients are copied from writeTSE.m
AGSex=gs_ex.area/2;
GSspr = mr.makeTrapezoid('z',sys,'area',AGSex*(1+fspS));

GS4times=[0 GSref.flatTime];
GS4amp=[GSref.amplitude GSref.amplitude];
GS4 = mr.makeExtendedTrapezoid('z','times',GS4times,'amplitudes',GS4amp);

if GSspr.flatTime==0
    GS5times=[0 GSspr.riseTime  GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
    GS5amp=[GSref.amplitude GSspr.amplitude  0];
    GS5 = mr.makeExtendedTrapezoid('z',sys,'times',GS5times,'amplitudes',GS5amp);

    GS7times=[0 GSspr.riseTime  GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
    GS7amp=[0 GSspr.amplitude  GSref.amplitude];
    GS7 = mr.makeExtendedTrapezoid('z',sys,'times',GS7times,'amplitudes',GS7amp);
else
    GS5times=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
    GS5amp=[GSref.amplitude GSspr.amplitude GSspr.amplitude 0];
    GS5 = mr.makeExtendedTrapezoid('z',sys,'times',GS5times,'amplitudes',GS5amp);

    GS7times=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
    GS7amp=[0 GSspr.amplitude GSspr.amplitude GSref.amplitude];
    GS7 = mr.makeExtendedTrapezoid('z',sys,'times',GS7times,'amplitudes',GS7amp);
end


%% change slice index, for interleaving slices: 1,3,5; 2,4,6
slice_indices_odd = 1:2:Nslices;
slice_indices_even = 2:2:Nslices;
slice_indices = [slice_indices_odd,slice_indices_even];


%% calculate ADC
adcTime = sys.gradRasterTime*size(spiral_grad_shape,2);
adcSamplesDesired = round( adcTime/adcDwell );
adcSamplesPerSegment = ceil( adcSamplesDesired/MAX_ADCS );
adcSamplesPerSegment = ceil(adcSamplesPerSegment/100)*100; % must be integer multiple of 100ns
if adcSamplesDesired/adcSamplesPerSegment > MAX_ADCS
    %warning('Current #ADC in a block is %.1f, exceeding maximal #ADC events allowed, which is %d', adcSamplesDesired/adcSamplesPerSegment,MAX_ADCS);
    error('Current #ADC in a block is %.1f, exceeding maximal #ADC events allowed, which is %d', adcSamplesDesired/adcSamplesPerSegment,MAX_ADCS);
end
adcSegments=fix(adcSamplesDesired/adcSamplesPerSegment);
adcSamples=adcSegments*adcSamplesPerSegment;
adcDwell=round(adcTime/adcSamples/100e-9)*100e-9; % on Siemens adcDwell needs to be aligned to 100ns (if my memory serves me right)
if adcDwell < 1e-6
    error('ADC dwell < 1us, it is highly likely you would fail on Siemens scanner!')
end
adcSegmentDuration=adcSamplesPerSegment*adcDwell; % with the 100 samples above and the 100ns alignment we automatically fullfill the segment alignment requirement
% if mod(adcSegmentDuration, sys.gradRasterTime)>eps 
if abs(   adcSegmentDuration/sys.gradRasterTime - round(adcSegmentDuration/sys.gradRasterTime)    ) > 1e-10
    error('ADC segmentation model results in incorrect segment duration');
end
% update segment count
adcSegments=floor(adcTime/adcSegmentDuration);
adcSamples=adcSegments*adcSamplesPerSegment;
adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',sys.adcDeadTime);
fprintf('<strong>ADC duration is %.1f ms\n</strong>', adc.duration*1e3);


%%
% extend spiral_grad_shape by repeating the last sample
% this is needed to accomodate for the ADC tuning delay
spiral_grad_shape = [spiral_grad_shape spiral_grad_shape(:,end)];

% readout grad 
gx_spi_out = mr.makeArbitraryGrad('x',spiral_grad_shape(1,:),sys,'Delay',sys.adcDeadTime);
gy_spi_out = mr.makeArbitraryGrad('y',spiral_grad_shape(2,:),sys,'Delay',sys.adcDeadTime);
gx_spi_out.first=0;
gy_spi_out.first=0;

% spoil
gz_spoil=mr.makeTrapezoid('z',sys,'Area',1/fov*Nx*4);
gx_spoil=mr.makeExtendedTrapezoid('x','times',[0 mr.calcDuration(gz_spoil)],'amplitudes',[spiral_grad_shape(1,end),0]); %todo: make a really good spoiler
gy_spoil=mr.makeExtendedTrapezoid('y','times',[0 mr.calcDuration(gz_spoil)],'amplitudes',[spiral_grad_shape(2,end),0]); %todo: make a really good spoiler


%% timing
rfCenterInclDelay    = rf90.delay + mr.calcRfCenter(rf90);
rf180centerInclDelay = rf180.delay + mr.calcRfCenter(rf180);
delayTE1     = TE/2   -   (mr.calcDuration(gs_ex)-rfCenterInclDelay)   -  mr.calcDuration(gzReph) - rf180centerInclDelay;
delayTE1     = ceil(delayTE1/sys.gradRasterTime)*sys.gradRasterTime;
delayTE2     = TE/2   -  (mr.calcDuration(GSref)-rf180centerInclDelay);
delayTE2     = ceil(delayTE2/sys.gradRasterTime)*sys.gradRasterTime;
delayTE1_b0  = TE/2   -   (mr.calcDuration(gs_ex)-rfCenterInclDelay)   -  mr.calcDuration(gzReph) - mr.calcDuration(GS7)  -  rf180centerInclDelay;
delayTE1_b0  = ceil(delayTE1_b0/sys.gradRasterTime)*sys.gradRasterTime;
delayTE2_b0  = TE/2   -  mr.calcRfCenter(rf180) - sys.rfRingdownTime  - mr.calcDuration(GS5);
delayTE2_b0  = ceil(delayTE2_b0/sys.gradRasterTime)*sys.gradRasterTime;

time_1_slice = mr.calcDuration(gz_fs) + rfCenterInclDelay + 0.5*TE + 0.5*TE + mr.calcDuration(gx_spi_out,gy_spi_out,adc)  +  mr.calcDuration(gx_spoil,gy_spoil,gz_spoil);
delayTR      = ceil(  (TR- Nslices*time_1_slice)  /  seq.gradRasterTime  )*seq.gradRasterTime;

fprintf('<strong>Minimal TR is %.2f, current TR is %.2f\n</strong>', time_1_slice*Nslices, TR);
delay_after_1slice = TR/Nslices - time_1_slice;
assert(delay_after_1slice>=0, 'At least, you should increase TR to %.2f seconds',time_1_slice*Nslices);
delay_after_1slice = ceil(delay_after_1slice/sys.blockDurationRaster)*sys.blockDurationRaster;
fprintf('<strong>Delay after each shot is %.2f\n</strong>',delay_after_1slice)
assert(delayTE1>=0);
assert(delayTE2>=0);
assert(delayTE1_b0>=0);
assert(delayTE2_b0>=0);
assert(all(delayTR>=0));


%%
small_delta = min([delayTE1,delayTE2])-ceil(sys_diff.maxGrad/sys_diff.maxSlew/sys_diff.gradRasterTime)*sys_diff.gradRasterTime;
big_delta   = delayTE1 + mr.calcDuration(GSref);
maxAvailBvalue = calc_bval_trap(physical_grad_max*diff_grad_safety_margin,small_delta,big_delta,ceil(sys_diff.maxGrad/sys_diff.maxSlew/sys_diff.gradRasterTime)*sys_diff.gradRasterTime);
assert(max(bvalues)<=maxAvailBvalue,'Required b-value too large, max achievable is %.1f, required is %.1f',maxAvailBvalue,max(bvalues));

num_bdir = size(bTable,1);
num_bval = numel(bvalues); % b0 is a special dir, instead of a b-value
num_volumes = num_bdir+(num_bval-1)*(num_bdir-Ndummies-1);
for iter_bval = 1:num_bval

for iter_bdir = 1:num_bdir
    if all(bTable(iter_bdir,:)==0)
        isb0 = true;
    else
        isb0 = false;
    end    

    if isb0 && iter_bval>1
        % from 2nd bvalue, skip if b0 , we don't want to acquire too many
        % b0. This continue keyword is ugly, TODO, fix it
        continue;
    end

    g = sqrt(bvalues(iter_bval)*1e6/bFactCalc_trap(1,small_delta,big_delta,ceil(sys_diff.maxGrad/sys_diff.maxSlew/sys_diff.gradRasterTime)*sys_diff.gradRasterTime));
    g_x=g.*bTable(iter_bdir,1);
    g_y=g.*bTable(iter_bdir,2);
    g_z=g.*bTable(iter_bdir,3);
    g_xr=ceil(abs(g_x)/sys_diff.maxSlew/sys_diff.gradRasterTime)*sys_diff.gradRasterTime;
    g_yr=ceil(abs(g_y)/sys_diff.maxSlew/sys_diff.gradRasterTime)*sys_diff.gradRasterTime;
    g_zr=ceil(abs(g_z)/sys_diff.maxSlew/sys_diff.gradRasterTime)*sys_diff.gradRasterTime;
    g_diff_rise = max([g_xr,g_yr,g_zr]);

    if ~isb0
        gDiff_x=mr.makeTrapezoid('x','amplitude',g_x,'riseTime',g_diff_rise,'flatTime',small_delta-g_diff_rise,'system',sys_diff);
        gDiff_y=mr.makeTrapezoid('y','amplitude',g_y,'riseTime',g_diff_rise,'flatTime',small_delta-g_diff_rise,'system',sys_diff);
        gDiff_z=mr.makeTrapezoid('z','amplitude',g_z,'riseTime',g_diff_rise,'flatTime',small_delta-g_diff_rise,'system',sys_diff);
        assert(mr.calcDuration(gDiff_x,gDiff_y,gDiff_z)<=delayTE1);
        assert(mr.calcDuration(gDiff_x,gDiff_y,gDiff_z)<=delayTE2);
    end

    % Define sequence blocks
    for shot_ind = 1:Nshots
        phi = (shot_ind-1) * rot_angle;
        for slc_ind=slice_indices
            seq.addBlock(rf_fs,gz_fs); % fat-sat

            rf90.freqOffset=gs_ex.amplitude*sliceThickness*(slc_ind-1-(Nslices-1)/2);
            rf90.phaseOffset=pi/2-2*pi*rf90.freqOffset*mr.calcRfCenter(rf90); % compensate for the slice-offset induced phase
            rf180.freqOffset=gz180.amplitude*sliceThickness*(slc_ind-1-(Nslices-1)/2);
            rf180.phaseOffset=-2*pi*rf180.freqOffset*mr.calcRfCenter(rf180); % compensate for the slice-offset induced phase

            % excitation
            seq.addBlock(rf90,gs_ex);
            % rephasing
            seq.addBlock(gzReph);

            if isb0
                seq.addBlock(mr.makeDelay(delayTE1_b0));
                seq.addBlock(GS7);
                seq.addBlock(rf180,GS4);
                seq.addBlock(GS5);
                seq.addBlock(mr.makeDelay(delayTE2_b0));
            else % no crusher needed for bvalue~=0
                seq.addBlock(mr.makeDelay(delayTE1),gDiff_x,gDiff_y,gDiff_z);
                seq.addBlock(rf180,GSref);
                seq.addBlock(mr.makeDelay(delayTE2),gDiff_x,gDiff_y,gDiff_z);
            end

            % acq
            seq.addBlock(mr.rotate_yxw('z',phi,sys,gx_spi_out,gy_spi_out,adc));
            seq.addBlock(mr.rotate_yxw('z',phi,sys,gx_spoil,gy_spoil,gz_spoil));

            if delay_after_1slice > 0
                seq.addBlock(mr.makeDelay(delay_after_1slice));
            end
        end
    end
end
end


%%
% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end
scanTime = sum(seq.blockDurations);
fprintf('<strong>Sequence duration is %02d:%02d\n</strong>', floor(scanTime/60),ceil(rem(scanTime, 60)));

%
seq.setDefinition('FOV', [fov fov sliceThickness*Nslices]);
seq.setDefinition('Name', 'spider');
seq.setDefinition('MaxAdcSegmentLength', adcSamplesPerSegment); % this is important for making the sequence run automatically on siemens scanners without further parameter tweaking
seq.setDefinition('NumShots',Nshots);
seq.setDefinition('NumArms',Narms);
seq.setDefinition('NumSlices',Nslices);
seq.setDefinition('MatrixSize',Nx);
seq.setDefinition('bValues',bvalues);
% seq.setDefinition('TR',TR); % may error if Siemens interpreter assumes a fixed TR
% seq.setDefinition('TE',TE);
seq.setDefinition('RepetitionTime',TR); 
seq.setDefinition('EchoTime',TE); 
seq.setDefinition('Resolution',fov/Nx);
seq.setDefinition('bTable',bTable);
seq.setDefinition('numSamplesPerShot',adc.numSamples);
seq.setDefinition('num_bdir',num_bdir);
seq.setDefinition('num_bval',numel(bvalues));
seq.setDefinition('slice_indices',slice_indices);
seq.setDefinition('num_volumes',num_volumes);
seq.setDefinition('Ndummies',Ndummies);
seq.setDefinition('Ndirections',Ndirections);
seq.setDefinition('Nechoes',Nechoes);
seq.setDefinition('sys_type',sys_type);
seq.setDefinition('grad_safety_magrin',grad_safety_magrin);
seq.setDefinition('slew_safety_magrin',slew_safety_magrin);
seq.setDefinition('spiral_slew_safety_margin',spiral_slew_safety_margin);
seq.setDefinition('spiral_grad_safety_margin',spiral_grad_safety_margin);
seq.setDefinition('diff_grad_safety_margin',diff_grad_safety_margin);
seq.setDefinition('diff_slew_safety_margin',diff_slew_safety_margin);


seq.write('esnails.seq');   % Output sequence for scanner

% the sequence is ready, so let's see what we got 
seq.plot('timeRange',[0,2*TR]);             % Plot sequence waveforms


%% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

ktraj_adc_splited = reshape(ktraj_adc, 3, adc.numSamples, Nechoes,Nslices,Nshots,num_volumes);
save(result_traj_fname,'ktraj_adc_splited') 
t_adc_splited = reshape(t_adc, 1, adc.numSamples, Nechoes,Nslices,Nshots,num_volumes);
figure;
plot(ktraj(1,:),ktraj(2,:),'b');
hold on;
ind_echo = Nechoes;
ind_slice = Nslices;
ind_volume = num_volumes;
for ind_shot = 1%:Nshots   
    plot(ktraj_adc_splited(1,:,ind_echo,ind_slice,ind_shot,ind_volume),ktraj_adc_splited(2,:,ind_echo,ind_slice,ind_shot,ind_volume),'r.')
end
legend('All trajectory','ADC samples of 1st shot')

xlim(1.2*Nx/2*1/fov*[-1 1])
ylim(1.2*Nx/2*1/fov*[-1 1])


%% calculate PNS and check whether we are still within limits
try
    % change the .asc file path accordingly!
    if strcmp(sys_type,'prisma')
        [pns,tpns]=seq.calcPNS('idea/asc/MP_GPA_K2309_2250V_951A_AS82.asc');
    elseif strcmp(sys_type,'skyra')
        [pns,tpns]=seq.calcPNS('idea/asc/MP_GPA_K2309_2250V_793A_GC99.asc');
    else
        error('Undefined system')
    end
catch
    warning('Please setup the path of the .asc file properly and do PNS simulation. Large PNS will prevent sequence from running')
end


%% forbbiden frequency calculation
try
    gradSpectrum;    
catch
    warning('This is important! Please add gradSpectrum to path and check forbidden frequency carefully to avoid hardware damage: pulseq\matlab\demoUnsorted\gradSpectrum.m')
end



%% helper functions
function b=bFactCalc(g, delta, DELTA)
% see DAVY SINNAEVE Concepts in Magnetic Resonance Part A, Vol. 40A(2) 39–65 (2012) DOI 10.1002/cmr.a
% b = gamma^2  g^2 delta^2 sigma^2 (DELTA + 2 (kappa - lambda) delta)
% in pulseq we don't need gamma as our gradinets are Hz/m
% however, we do need 2pi as diffusion equations are all based on phase
% for rect gradients: sigma=1 lambda=1/2 kappa=1/3 
% for trapezoid gradients: TODO
sigma=1;
%lambda=1/2;
%kappa=1/3;
kappa_minus_lambda=1/3-1/2;
b= (2*pi * g * delta * sigma)^2 * (DELTA + 2*kappa_minus_lambda*delta);
end

function b = bFactCalc_trap(g, delta, DELTA, rup)
% handbook, page 287, example 9.2
b = (2*pi*g)^2 * (   delta^2*(DELTA-delta/3)  +  rup^3/30-delta*rup^2/6  );
end

function bvalue = calc_bval_trap(grad_max,small_delta,big_delta,rup)
% grad unit is mT/m, time unit is seconds
% handbook, page 287, example 9.2
gamma = 2*pi*42.57e3;           % rad/mT
bvalue = (gamma.*grad_max/1000).^2*( small_delta^2*(big_delta-small_delta/3) + rup^3/30-small_delta*rup^2/6  ); 
end