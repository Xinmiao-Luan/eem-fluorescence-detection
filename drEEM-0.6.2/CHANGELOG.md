### v0.6.2
- Modifications:
  - `outliertest`: Fixed R2020b compatibility and spelling mistake.
  - `handlescatter`: Fixed behaviour for selective scatter excision & interpolation.
  - `checkdataset`: Fixed displaying message when no models are in the structure.
  - `splitanalysis`: Fixed a bug in the non-parallel function execution (thanks to Morimaru Kida)

### v0.6.1
- Modifications:
  - `splitanalysis`: Bugfix in input argument parsing. Now works as documented.
  - `handlescatter`: Bugfix. Smoothing now works as expected (previously interpolation was done for all types and ignored the user input in opt.interpolate.
  - `modelout` & `modeloutmac`: Bugfix in drEEM version reporting. Support for `handlescatter` added (no error message anymore)
  - `assembledataset`: Added .i as a standard field to the function output.
- New:
  - `modelexport`: Based on `modelout`, but works faster and is platform independent.
     This function may require Matlab R2019a or newer to work properly.

### v0.6.0
- Modifications
  - `normeem`: Automated detection of blanks / low intensity samples. User is asked if sample-exclusion should occur.
  - `normeem`: New option allows automated choice between leaving and removing low-intensity samples.
  - `nwayparafac`: Minor changes to improve forward compatibility & speed.
  - `fastnnls`: Minor changes to improve speed.
  - `nmodel`: Minor changes to improve speed.
  - `outliertest`: Parallelization support & improved graphical output.
  - `eemreview`: Changed graphical layout. Now includes drop-down menu for sample selection and checkboxes for options.
  - `fingerprint`: New option to show scatter as lines (superimposed)

- Deletions
  - `prandinitanal`: Function replaced with new `randinitanal`. Multitreading capabilities retained.
  - `psplitanalysis`: Function replaced with new version of `splitanalysis`. Multitreading capabilities retained.
  - `multiwaitbar`: Replaced with improved function (`blockbar`)

- New
  - **Documentation was completely overhauled.** Each function now has a html help file. These can be called
    by the command `doc functionname` or `help functionname` (& click on link).
    The old help in each m-file still exists, but **is no longer updated**.
  - `openfluormatches`: Plot exported results from database matches (openfluor.org)
  - `randinitanal`. Completely rewritten. Automatic detection parallel computing capabilities. Changes:
    - Allow user to cancel unfinished models (parallel mode).
    - Updated help section (minor fixes).
    - Fixed error in calculation of component size (field "compsize").
    - User can specify whether PLS_toolbox or N-way toolbox should be used for PARAFAC.
    - Switched to input parsing of input arguments, but old input is also recognized.
    - Plots existing models if dataset with models is provided and no other inputs are given
  - `splitanalysis`. Completely rewritten. Automatic detection parallel computing capabilities. Changes:
    - User can specify whether PLS_toolbox or N-way toolbox should be used for PARAFAC.
    - Switched to input parsing of input arguments, but old input is also recognized.
	- User can track all models with detailed progress bar and cancel components if desired.
  - `checkdataset`. Troubleshooting function to validate a dataset structure and identify issues.
    - Function runs quietly in other functions and reports problems if they are found.
  - `handlescatter`: `smootheem` with new iterative interpolation method (inpainting).
  - `dreemfig`: A interal function that makes all of drEEM's plots look much nicer.
  - Support for datacursor in various plotting functions. When using the cursor, more information is
    now shown, such as filename of the sample, wavelength, index of a sample or wavelength etc.

### v0.5.1
- Modifications:
  - `prandinitanal`: 
    - Fixed calculation of component contribution
    - Fixed overwriting issue that resulted in cancellation of user input to convergence criterion (`CC`)
    - Other minor fixes & improvements in function help and function code.
  - `coreandvar`: Minor updates to account for the possibility of > 10 components.
  - `psplitanalysis`: Updated function help section and added console output.
  - `eemreview`: Minor fixes of spelling mistakes.
  - multiple minor improvements.

### v0.5.0
- New
  - `eemreview`: Visualize EEMs and select emission and excitation slices for detailed inspection
  - `scanview`: Visualize EEMs one excitation measurement at a time
  - `coreandvar`: Plot core consistency and % expl. variance for different number of components
  - `spectralvariance`: Visualize the variance within a given CDOM and FDOM dataset by means of the spectral standard deviation
  - `slopefit`: Fit and diagnose CDOM spectral slopes (Exponential + log-linear CDOM slopes)
  - `pickpeaks`: Pick FDOM peaks and calculate fluorescence indicies
  - `prandinitanal`: Utilize all available physical CPU's to speed up PARAFAC fitting
  - `psplitanalysis`: Utilize all available physical CPU's to speed up PARAFAC fitting in dataset splits
  - `diffeem`: Calculate & visualize the difference between a reference EEM and other samples.
  - `errorsandleverages`: Inspect leverage of samples and wavelengths agains errors to identify problematic variables / samples.
  - `dreeminstall`: Forget addpath(...), just let us make the setup of drEEM for you.
  

- Modifications: 
  - `spectralloadings`: Changed default color to new Matlab default and added cursor callbacks with detail info.
  - `ini.m` (nway): replaced random number seed function to ensure forward compatibility 
  - `nwayparafac` (nway): The function will not create temp.mat files anymore
  - `specsse`: Visualize the sum-of-squared errors in all three measurement modes (sample x emission x excitation)

- Modified tutorials
  - `drEEM_parallelcomputing`: 				 Short demonstration of how parallel computing works in MATLAB
  - `drEEM_shorttutorial`:     				 drEEM tutorial, slightly shorter version
  - `drEEM_shorttutorial_parallelcomputing`: drEEM tutorial, slightly shorter version, with parallel computing
  - `drEEM_fulltutorial_parallelcomputing`:	 drEEM tutorial with parallel computing

- Comments:
  - PLS_toolbox as PARAFAC engine now supported. Add PLS_toolbox beneath drEEM to the path, and the functions will chose to model your dataset with PLS_toolbox instead of nway. Note: PLS_toolbox PARAFAC takes a somewhat longer to initialize, this is not due to drEEM.

### v0.4.0 (Kathleen Murphy)
- Several minor bug fixes since previous releases. variable / function name conflict in `nprocess` / `normeem` caused issues with new Matlab version.

### v0.2.0 (Kathleen Murphy)

- What is new

  (1) Now compatible with MAC/Linux operating systems:
    ```
	readineems
	readinscans
	fdomcorrect
	modeloutmac
    ```

  (2) Now compatible with AquaLog fluorometer files:
    ```
	readineems
	readinscans
    ```

	For example, create a cube of AquaLog EEMs (X) and an EEM dataset (DS):
	```
    [X,Emmat,Exmat,fl,outdata]=readineems(3,'*Waterfall Plot Sample.dat','A2..BV126',[0 1],0,2);
	DS = assembledataset(X,Exmat(1,:),Emmat(:,1),'AU','filelist',fl,[]);
    ```

  (3) Other changes/enhancements:	
	 (a) `fdomcorrect`
 	   - outputs to .dat by default (Excel is optional)
	   - exports sample names if you use DS not X as the input variable (see above example)
	   - faster and more accurate calculation of Raman Areas
	   
	 (b) `randinitanal`
	   - selection of the best LS model excludes models that did not converge
	  
	 (c) `outliertest`
	   - fixed bug that produces wrong nSample value in the case that nSample< max(nEx,nEm)
	   
	 (d) `splitds`
	   - fixed bug that causes failure to recognise split types in newer MATLAB versions

  (4) New EEM correction demo mfile

	drEEM_demo_020
	   
  (4) New functionality:	
	`openfluor` - export model loadings ready for direct upload to OpenFluor (www.openfluor.org)
	`readlogfile` - import a log file (*.csv) into a data structure (Windows/Mac/Linux)
	`alignds` - align samples on the basis of the above log file (Windows/Mac/Linux)
	`ABAife` - correct for inner filter effects by absorbance method in 1 cm cell
	`CDAife` - correct for inner filter effects by controlled dilution approach (Luciani et al. 2009)
	`IFEtest_demo` - comparison of inner filter effect correction methods (Kothawala et al. 2013)
	
  (5) Obselete functions
    `matchsamples` - Use `readlogfile` and `alignds` instead of matchsamples
