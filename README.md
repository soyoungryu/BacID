# BacID 
Author: So Young Ryu <br />
Reference: Ryu S. Likelihood-based bacterial identification approach for bimicrobial mass spectrometry data. (Under Review) <br />

## Codes Folder
bacID_functions.R: R functions needed for BacID Prob and BacID Score. <br />
BacID_Prob.R: An example R file to run BacID Prob using an example dataset. <br />
BacID_Score.R: An example R file to run BacID Score using an example dataset. <br />

## RawData
Data.RData: This contains both Standard Bacterial Mixture Dataset and Co-Cultured Bacterial Mixture Dataset <br />
There are three objects in this file: Mix.test, Mix.train, and ref. <br /> 

### Mix.test 
Mix.test has the following information: peak.list (m/z and intensity values before normalization at the first column), true bacterial IDs, and dataset.name (whether a mass spectrum is from Standard Bacterial Mixture or Co-Cultured Bacterial Mixture.). Mix.test contains a total of 127 mass spectra, thus there are a total of 254 (=127*2) species in these mass spectra. <br />

### Mix.train
Mix.train contains the same information as Mix.test. Mix.train contains a total of 54 mass spectra, thus there are a total of 108 (=54*2) species in these mass spectra. <br />

### ref
ref is a reference database including  reference mass spectra of Bacillus subtilis (Bs), Enterobacter cloacae (El), Escherichia coli (Ec), Klebsiella oxytoca (Ko), Klebsiella pneumoniae (Kp), Pseudomonas aeruginosa (Pa), Pseudomonas fluorescens (Pf), and Staphylococcus aureus (Sa) in addition to 1,000 decoy mass spectra. <br />

## Results
### BacID_Prob_Results.RData
It contains BacID Prob results of Mix.test data which contains both Standard and Co-cultured Bacterial Mixture Datasets. There are analysis results of 127 mass spectra (107 mass spectra from Standard dataset and 20 mass spectra from Co-Cultured dataset.) <br />

### BacID_Score_Results.RData
It contains BacID Score results of Mix.test data which contains both Standard and Co-cultured Bacterial Mixture Datasets. There are analysis results of 127 mass spectra (107 mass spectra from Standard dataset and 20 mass spectra from Co-Cultured dataset.)



