# Usage
Written and ran using: MATLAB R2019b

Each folder contains `testAnalyze.m` that parses the output FCS files using [`fca_readfcs.m` ](https://www.mathworks.com/matlabcentral/fileexchange/9608-fca_readfcs) and plots the populations either as a scatter plot or histogram. [`tightfig.m`](https://www.mathworks.com/matlabcentral/fileexchange/34055-tightfig-hfig?s_tid=srchtitle) removes excess white space from figures. FCS files are included in each folder to reproduce the plots.

`%% Load data` section assigns the selected populations into two variables `[name,nameHdr]`, wherein `name` contains the fluorescence intensity values for each experiment and `nameHdr`is a [structured array](https://www.mathworks.com/help/matlab/ref/struct.html) which contains the `par` variable to indicate the row number of the fluorescence channel of interest.

**Minimum working example:**
* Go to `figure2_flying` folder → Run `testAnalyze.m`
* Double-clicking the `plane1Hdr` opens up the content of the structured array. Row 7 contains the `FJComp-APC-A` which corresponds to the column of intensity values in `plane1`.
