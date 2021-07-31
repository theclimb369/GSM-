# GSM 

**Description** 
  Grid search modelling method developed to model 1-D seismic structure of the lowermost mantle 
**Prerequisite** 
  WKBJ code (Chapman and Orcutt, 1985) & Python installed
  Preferable package - GMT, SAC

**Data Processing**
1. Data processing Supported software - SAC
2. Data grouping 
3. Data alignment & stacking - Adaptive alignment stacking (Rawlinson and Kennett, 2004)
		Code: Adaptive_alignment.py; Stack_process_data.py;
 
**GSM**
1. Build synthetic waveform library
	  Code: Model_generator_5p.py; 
2. Calculate misfits between synthetics and data
		Code: Misfitcalculation.py; 
3. Conduct Likelihood Ratio Test - 1D & 2D
		Code: Likelihood_ratio_test_1D_2D.py; 
