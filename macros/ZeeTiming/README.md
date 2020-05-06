# ZeeTiming README
### Summer 2017
General structure of code:
1. Include statements
2. A lot of helper function definitions
3. Main function called `Zee_DG()`
   - Specify ntuple
   - Specify which processes to run (`bool do...`)
   - Define variables and set branch address
   - Loop through ntuple once to find min/max values
   - Initialize then declare histograms
   - For each process:
     - Verify process is turned on
     - Loop through ntuple again, filling histograms with the specific cuts/values
     - Call helper functions to generate plots

Things to do before running code:
- Ensure input ntuple file is in same directory as Zee_DG.C
- Update ntuple filename in code (`string filename = "All2016.root"`).
- Select the analyses you wish to perform, e.g. set `doRun = 1` instead of 0.
- With Zee_DG.C in ./ we also need to create directories ./plots/c, ./plots/pdf, ./plots/png

Adjusting cuts/parameters:
- Search for `float crystalEnergyMin`. In that area, I define ranges for crystal energy, transparency, etc.
- Many bins are determined by the line `int nSteps = 50`. This divides all the ranges into 50 bins.
- For the 2D ieta vs iphi vs z_variable plots, `int BinSize` sets the number of ieta and iphi within 1 bin. Note that there are 170 ieta and 360 iphi, so BinSize must be 1, 2, 5, or 10.
- In some helper functions, there is `bool singleGauss`, which is 0 if you want to perform a double Gaussian fit, and 1 to do a single Gaussian.
- In `crystalTOF2()` and `crystalTOF3()`, the cuts are hardcoded because there was not much of a need to change them. Specifically, `energyCut_Rel` and `deltaRcut`.

## Ask Daniel Gawerc if clarification is needed
