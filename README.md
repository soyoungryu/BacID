# BacID 
**Author: So Young Ryu <br />**
**Reference: Ryu S. (In press) Likelihood-based bacterial identification approach for bimicrobial mass spectrometry data. *Annals of Applied Statistics*. <br />**

## File Description
### Codes Folder
bacID_functions.R: R functions needed for BacID. <br />
BacID.R: An example R file to run BacID using an example dataset. <br />

### RawData Folder
Data.RData: This contains both Standard Bacterial Mixture Dataset and Co-Cultured Bacterial Mixture Dataset. There are two objects in this file: Mix and ref. <br /> 
Data Reference: Yang, Y., Lin, Y. and Qiaq, L. (2018). Direct MALDI-TOF MS identification of bacterial mixtures. *Analytical Chemistry* 90 10400–10408. <br /> 

**Mix** <br /> 
Mix has the following information: peak.list (m/z and intensity values before normalization at the first column), true bacterial IDs, and dataset.name (whether a mass spectrum is from Standard Bacterial Mixture or Co-Cultured Bacterial Mixture.). Mix contains a total of 127 mass spectra, thus there are a total of 254 (=127*2) species in these mass spectra. <br />

**ref** <br />
ref is a reference database including  reference mass spectra of *Bacillus subtilis* (Bs), *Enterobacter cloacae* (El), *Escherichia coli* (Ec), *Klebsiella oxytoca* (Ko), *Klebsiella pneumoniae* (Kp), *Pseudomonas aeruginosa* (Pa), *Pseudomonas fluorescens* (Pf), and *Staphylococcus aureus* (Sa) in addition to 1,000 decoy mass spectra. <br />

### Results Folder
**BacID_Results.RData**
It contains BacID results of Mix data which contains both Standard and Co-cultured Bacterial Mixture Datasets. There are analysis results of 127 mass spectra (107 mass spectra from Standard dataset and 20 mass spectra from Co-Cultured dataset.) <br />

## How to Run Codes
Prerequisite: R or R studio, and *doParallel* R package. <br />
1. Please download codes from Codes folder and data from RawData folder in one folder of user's computer. <br />
2. Open BacID.R. <br />
3. Specify direcotry at the first line of BacID.R. <br />
For example, <br />
```
directory = 'C:/Documents/'
```
4. Modify *no.codes* in BacID.R. <br />
```
no.cores=500
```
5. Run the R codes. <br />

**This will produce an output file named as BacID_Results.RData.**


