# OAM-Ctrl

OAM-Ctrl is the open source code of the paper "Xin Liu, Yifan Peng, Shijie Tu, Jun Guan, Cuifang Kuang, Xu Liu, and Xiang Hao, [Generation of arbitrary longitudinal polarization vortices by pupil function manipulation](https://doi.org/10.1002/adpr.202000087). Advanced Photonics Research 2(1), 2000087 (2021)."

- **Please cite this paper if you use this code for the related calculation in your work.**
- The current code is for proof of concept only. Correspondence and requests for a faster version (1000-time acceleration) should be addressed to Xiang Hao (haox@zju.edu.cn).

## Getting Started

- The code is named with the index of the figures in the paper.  
- Run the code, and you can obtain the corresponding result in the paper.

### Prerequisites

MATLAB 2020a.

### Example

1. Open the script named 'Figure3a_d_singleEzTunable.m'.
2. Set the radial polarization ratio value in line 14.
```
A = 0;
```
- The pupil function is calculated and shown as follows:  
![image](https://github.com/Hao-Laboratory/OAM-Ctrl/blob/master/OAM-Ctrl/data/Pupil%20Function.png)

## Author

* **Xin Liu** - *Initial work* - [LiuX2018](https://github.com/LiuX2018).

## Note!!!

- This version is only for the basic **demo** purposes.  It does not include the performance optimization modules, e.g., Parallel Computing, and Graphical User Interface.
- Please contact us if you want the **complete** version (fully functional) of this code.
- Email: liuxin.optics@gmail.com; haox@zju.edu.cn.

## License

This project is licensed under the GNU GPLv3 License - see the [LICENSE.md](LICENSE.md) file for details.
