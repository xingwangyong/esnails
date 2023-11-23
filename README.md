# eSNAILS

Enhancing Self-Navigated Interleaved Spiral with ESPIRiT (eSNAILS). Use ESPIRiT to improve the composite sensitivity map estimation in multi-shot spiral dMRI.

## Installation
Download this repo first, then download 3rd party toolboxes accordingly.

### For sequence generation
1. Download and add the following toolboxes to Matlab path
    - [pulseq](https://github.com/pulseq/pulseq), this repo uses [v1.4.1](https://github.com/pulseq/pulseq/releases/tag/v1.4.1)
    - [Time Optimal Gradient Design](http://people.eecs.berkeley.edu/~mlustig/Software.html), this repo uses [v0.2](http://people.eecs.berkeley.edu/~mlustig/software/tOptGrad_V0.2.tar.gz)

2. copy rotate_yxw.m to the folder \pulseq\matlab\+mr\
In pulseq v.1.4.1, the mr.rotate() does not take `system` as an input, thus have no idea about system properties and uses default gradient and slew. This would fail if strong gradient and slew are used.
### For reconstruction
Download and add the following toolbox and files to Matlab path
- [bart](https://github.com/mrirecon/bart)
- [gridmat.m](https://github.com/mribri999/MRSignalsSeqs/blob/master/Matlab/gridmat.m)
- [voronoidens.m](https://github.com/zenghf/NMRTool/blob/aa1199937a7252887f5ae64612e8d131e10c93e0/csMRI/utils/voronoidens.m)
- (optional) [mapVBVD](https://github.com/pehses/mapVBVD), needed when you want to reconstrcut from Siemens rawdata




## Usage
### Sequence file generation
run the following command, which would generate a `esnails.seq` that can be used on scanner.
```
>>esnails
```
### Recon
Two similar recon scripts are provided.
- recon_mat.m: data are saved in .mat file
- recon_dat.m: read data from Siemens rawdata (.dat) file, will need pulseq on matlab path.