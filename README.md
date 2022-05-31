# k-t-TMNN

This is the code of the paper 'Dynamic Cardiac MRI Reconstruction Using Combined Tensor Nuclear Norm and Casorati Matrix Nuclear Norm Regularizations'.
```
Dynamic Cardiac MRI Reconstruction Using Combined Tensor Nuclear Norm and Casorati Matrix Nuclear Norm Regularizations.
Zhang, Yinghao, and Yue Hu. 
2022 IEEE 19th International Symposium on Biomedical Imaging (ISBI). IEEE, 2022.
```


* For the k-t TMNN algorithm, run the `Main_dmri_TMnn.m` to test it. The file uses the fast algorithm in frequency domain. The file also contains the fast algorithms for only tensor nuclear norm (TNN) and only Casorati matrix nuclear norm (MNN).  

<br/>

* For the comparation of the classical algorithm in spatiotemporal domain and the fast algorithm in frequency domain, run the `Main_dmri_TMnn_vs_fast.m`   

<br/>

* The data used in this code can be found in [yhao-z/dMRI-Data: dMRI data that I always use (github.com)](https://github.com/yhao-z/dMRI-Data)

* The files in the path `algs` :

`mnnAlg.m` and `mnnAlg_fast.m` denotes the classical algorithm using Casorati matrix nuclear norm based on ADMM and the fast algorithm of that.

so is the `tnnAlg.m` and `tnnAlg_fast.m` , where tnn denotes tensor nuclear norm.

`tmnnAlg.m` and `tmnnAlg_fast.m` is the same relationship.

`prox_nuclear.m` and `prox_tnn.m` implement the singular value thresholding on Casorati and the tensor singular value thresholding.

`calc_tnn.m` calculate the tensor nuclear norm of a tensor, e.g., <a href="https://www.codecogs.com/eqnedit.php?latex=\mathcal{A}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathcal{A}" title="\mathcal{A}" /></a>   

<br/>

* the `mnnAlg.m` and the `tnnAlg.m` is not used in the `Main` codes, but it has been tested.

