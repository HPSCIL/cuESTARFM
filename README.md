cuESTARFM
========
Version 1.0

Overview
========
MODIS and Landsat surface reflectance products have complementary characteristics in terms of spatial and temporal resolutions. To fully exploit these datasets, the Spatial and Temporal Adaptive Reflectance Fusion Model (STARFM) was developed by Gao et al. (2006). The STARFM approach blends the high-frequency temporal information from MODIS and the high-resolution spatial information from Landsat to generate synthetic surface reflectance products at 30m spatial resolution and daily temporal resolution. STARFM uses one or more pairs of Landsat-MODIS images collected on the same dates to predict surface reflectance at Landsat resolution on other MODIS observation dates. In order to better predict the reflectance of sub-pixel consisting of heterogeneous landscapes, an enhanced STARFM (ESTARFM) was developed by Zhu et al. (2010), which is based on the spectral unmixing theory and uses a “conversion coefficient” to help the prediction. However, the computational performance of ESTARFM has been a bottleneck for mass production.

To overcome the computational barrier and support mass production of large-size images, we designed and implemented a GPU-enabled ESTARFM program based on the Compute Unified Device Architecture (CUDA), called cuESTARFM. By taking advantages of the large amount of concurrent computing threads of a GPU, cuESTARFM can greatly reduce the computing time and improve the computational performance. Experiments showed that cuESTARFM achieved a speedup of 75 using a Nvidia Tesla K40 GPU, compared with a sequential ESTARFM program running on an Intel Xeon E3-1226 CPU.

Key features of cuESTARFM:
========
+ Supports a wide range of CUDA-enabled GPUs (https://developer.nvidia.com/cuda-gpus)  
  - Automatic setting of the numbers of threads and thread blocks according to the GPU’s available computing resources (e.g., memory, streaming multiprocessors, and warp)  
  - Adaptive cyclic task assignment to achieve better load balance
  - Optimized use of registers to improve the computational performance
  - Adaptive data decomposition when the size of images exceeds the GPU’s memory  
  -	All above are completely transparent to users
+	Intakes any number of pairs of Landsat-MODIS images as the input
+ Outputs any number of prediction images
+ Supports a wide range of image formats (see http://gdal.org/formats_list.html)
+ Supports both Windows and Linux/Unix operating systems

References
========
+ Gao, F.; Masek, J.; Schwaller, M. and Hall, F. On the Blending of the Landsat and MODIS Surface Reflectance: Predict Daily Landsat Surface Reflectance, IEEE Transactions on Geoscience and Remote Sensing. 2006, 44(8):2207-2218.   
+ Zhu, X.; Chen, J.; Gao, F.; Chen, X. and Masek, J. An enhanced spatial and temporal adaptive reflectance fusion model for complex heterogeneous regions, Remote Sensing of Environment, 2010, 114(11): 2610–2623.

To Cite cuESTARFM in Publications
========
+ A paper describing cuESTARFM will be submitted to a scientific journal for publication soon
+	For now, please cite the following reference:  
Gao, H., Zhu, X., Guan, Q., Yang, X., Yao, Y., Zeng, W., Peng, X., 2021. cuFSDAF: An Enhanced Flexible Spatiotemporal Data Fusion Algorithm Parallelized Using Graphics Processing Units. IEEE Transactions on Geoscience and Remote Sensing. https://doi.org/10.1109/TGRS.2021.3080384
Compilation
========
+ Requirements:
  -	A computer with a CUDA-enabled GPU (https://developer.nvidia.com/cuda-gpus)
  -	A C/C++ compiler (e.g., Microsoft Visual Studio for Windows, and gcc/g++ for Linux/Unix) installed and tested
  -	Nvidia CUDA Toolkit (https://developer.nvidia.com/cuda-downloads) installed and tested
  -	Geospatial Data Abstraction Library (GDAL, http://gdal.org) installed and tested
+ For the Windows operating system (using MS Visual Studio as an example)
  1. Open all the source codes in Visual Studio
  2. Click menu Project -> Properties -> VC++ Directories -> Include Directories, and add the “include” directory of GDAL (e.g., C:\GDAL\include\)
  3. Click menu Project -> Properties -> VC++ Directories -> Lib Directories, and add the “lib” directory of GDAL (e.g., C:\GDAL\lib\)
  4. Click menu Build -> Build Solution  
  Once successfully compiled, an executable file, cuESTARFM.exe, is created.
+ For the Linux/Unix operating system (using the CUDA compiler --- nvcc)  
In a Linux/Unix terminal, type in: 
  - $ cd /the-directory-of-source-codes/
  - $ nvcc -o cuESTARFM kernel.cu cuLayer.cpp cuESTARFM_util.cpp trans.cpp -lgdal  
  Once successfully compiled, an executable file, cuESTARFM, is created.
  
Usage 
========
+ Before running the program, make sure that all Landsat and MODIS images have been pre-processed and co-registered. They must have:
  - the same spatial resolution (i.e., Landsat resolution --- 30m)
  - the same image size (i.e., numbers of rows and columns)
  - the same map projection
+ A text file must be manually created to specify the input and output images, and other parameters for the ESTARFM model.  
Example (# for comments):

>ESTARFM_PARAMETER_START
>
>#The number of input pairs of Landsat-MODIS images (>=2)  
>  NUM_IN_PAIRS = 2
>
>#The input MODIS images   
>#File names are separated by space  
>  IN_PAIR_MODIS_FNAME = D:\cuda\shikong\testdata\M_2002_01_04.tif  D:\cuda\shikong\testdata\M_2002_02_21.tif
>
>#The input Landsat images  
>#File names are separated by space  
>  IN_PAIR_LANDSAT_FNAME = D:\cuda\shikong\testdata\L_2002_01_04.tif D:\cuda\shikong\testdata\L_2002_02_21.tif
>
>#The MODIS images for the prediction dates  
>#Multiple images can be given  
>#File names are separated by space  
>  IN_PDAY_MODIS_FNAME = D:\cuda\shikong\testdata\M_2002_02_12.tif
>
>#The output synthetic prediction images  
>#Multiple images can be given  
>#File names are separated by space  
>  OUT_PDAY_LANDSAT_FNAME = D:\cuda\shikong\testdata\estr.tif
>
>#The_width of searching_window  
>  The_width_of_searching_window = 51
>
>#Assumed_number of classifications  
>  Assumed_number_of_classifications = 6
>
>#Landsat sensor error   
>  sensor_uncertain = 0.0028
>
>#Output image format (optional)  
>#Will be used when the extension of the output files  
>#is not given  
>  G_Type = GTIff
>
>ESTARFM_PARAMETER_END

+ The program runs as a command line. You may use the Command (i.e., cmd) in Windows, or a terminal in Linux/Unix. 
   - For the Windows version:    
   $ cuESTARFM.exe parameters.txt 
   - For the Linux/Unix version:   
   $ ./cuESTARFM parameters.txt 

+ Note: The computational performance of cuESTARFM largely depends on the GPU. The more powerful is the GPU, the better performance. 
