# Model independent tests of Cosmology
Repository of the codes and data files I have been developing to perform model-independent tests of Cosmology using real H(z) measurements, as well as Euclid and SKA simulations.

The python code 'model_ind_tests.py' carries out Gaussian Processes reconstructions with GaPP (https://github.com/astrobengaly/GaPP) in order to obtain model-independent cosmological measurements such as: 
  - Hubble Constant H0 
  - deceleration parameter q0
  - null tests of the LCDM model, namely the Om and Lm tests as in [https://arxiv.org/abs/0807.3548], [https://arxiv.org/abs/0807.4304] (see also [https://arxiv.org/abs/0712.3457])

Data files in this repository consist of: 
  - Real H(z) measurements from galaxy ages, aka Cosmic Chromoneters, and from the radial BAO mode obtained by the SDSS redshift survey
  - Simulated radial BAO measurements from Euclid galaxy survey, and SKA band1 and band2 21cm Intensity Mapping. Fiducial Cosmology assumes Planck 2018 (TT, TE, EE+lowE+lensing) best-fit

The code receives the number of data points and name of the survey as inputs
  - For real data: nz1 refers to the number of CC data points, nz2 to the BAO ones. Put 0 in nz2 in case you just want the former 
  - For simulations: nz1 refers to the number of Euclid or SKA band1 data points, nz2 to the SKA band2. Put 0 in nz2 if you do not want to include band 2 simulations
  - survey refers to 'euclid', 'ska', or 'cc'

Plots of the reconstructions assuming SKA band1 + band2 and observational CC + SDSS radial BAO are also provided in this repository. 

I hope this code is helpful for students and anyone interested in performing this analysis in the future. 
If you use it, please cite the following papers

  - Null tests of the concordance model in the era of Euclid and the SKA; Phys. Dark. Univ. accepted [https://arxiv.org/abs/2007.04879]
  - Evidence for cosmic acceleration with next-generation surveys: A model-independent approach; Mon.Not.Roy.Astron.Soc. 499 (2020) 1, L6-L10 
[https://arxiv.org/abs/1912.05528]
  - The Hubble constant tension with next-generation galaxy surveys; JCAP05(2020)053 [https://arxiv.org/abs/1908.04619] 
  
Please contact carlosbengaly@on.br or carlosap87@gmail.com for further enquiries. Suggestions are always welcome. 
