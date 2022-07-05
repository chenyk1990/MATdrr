**MATdrr**
======

## Description

**MATdrr** is a Matlab package for the damped rank reduction (DRR) method and its several variants. The DRR method has a variety of applications in both exploration and earthquake seismology, including but not limited to seismic denoising, seismic reconstruction, seismic diffraction separation, constrained LSRTM, constrained FWI, etc.

## Reference
    Huang, W., Wang, R., Chen, Y., Li, H., & Gan, S. (2016). Damped multichannel singular spectrum analysis for 3D random noise attenuation. Geophysics, 81(4), V261-V270.

    Chen, Y., Huang, W., Zhang, D., & Chen, W. (2016). An open-source Matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction. Computers & Geosciences, 95, 59-66.
    
    Chen, Y., Zhang, D., Jin, Z., Chen, X., Zu, S., Huang, W., & Gan, S. (2016). Simultaneous denoising and reconstruction of 5-D seismic data via damped rank-reduction method. Geophysical Journal International, 206(3), 1695-1717.
    
BibTeX:

	@article{huang2016dmssa,
	  title={Damped Multichannel Singular Spectrum Analysis for 3{D} Random Noise Attenuation},
	  author={Weilin Huang and Runqiu Wang and  Yangkang Chen and Huijian Li and Shuwei Gan},
	  journal={Geophysics},
	  volume={81},
	  number={4},
	  issue={4},
	  pages={V261-V270},
	  year={2016},
	  publisher={Society of Exploration Geophysicists}
	}

	@article{chen2016drr5d,
	  title={Simultaneous denoising and reconstruction of 5{D} seismic data via damped rank-reduction method},
	  author={Yangkang Chen and Dong Zhang and Zhaoyu Jin and Xiaohong Chen and Shaohuan Zu and Weilin Huang and Shuwei Gan},
	  journal={Geophysical Journal International},
	  volume={206},
	  number={3},
	  issue={3},
	  pages={1695-1717},
	  year={2016}
	}
	
	@article{chen2016drr3d,
	  title={An open-source Matlab code package for improved rank-reduction 3{D} seismic data denoising and reconstruction},
	  author={Yangkang Chen and Dong Zhang and Weilin Huang and Wei Chen},
	  journal={Computers \& Geosciences},
	  volume={95},
	  pages={59-66},
	  year={2016}
	}

-----------
## Copyright
    Initial version: Yangkang Chen (chenyk2016@gmail.com), 2021-2022
	Later version: MATdrr developing team, 2022-present
-----------

## License
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)   

-----------

## Install
Using the latest version

    git clone https://github.com/chenyk1990/MATdrr
    cd MATdrr
    addpath(genpath('./')); #in Matlab command line
    
-----------
## Examples
    The "demo" directory contains all runable scripts to demonstrate different applications of MATdrr. 
    
-----------
## Gallery
The gallery figures of the MATdrr package can be found at
    https://github.com/chenyk1990/gallery/tree/main/matdrr
Each figure in the gallery directory corresponds to a DEMO script in the "demo" directory with the exactly the same file name.

-----------
## Dependence Packages
* Matlab 2015 and later versions

-----------
## Naming criteria
    drr3d_xxx.m corresponds to all functions.
    drr_xxx.m corresponds to all subroutines.
    localxxx.m corresponds to key functions in local windows.
    test_xxx.m corresponds to all DEMO scripts that are runnable.
    
-----------
## Development
    The development team welcomes voluntary contributions from any open-source enthusiast. 
    If you want to make contribution to this project, feel free to contact the development team. 

-----------
## Contact
    Regarding any questions, bugs, developments, collaborations, please contact  
    Yangkang Chen
    chenyk2016@gmail.com

-----------


