Make any updates based on changes to the code or research in general 
23/03/2026 - Kornel: Created the R document and added the following libraries (tidyverse, tidyfinance, scales, frenchdata (main dataset), dplyr, moments, sandwich, lmtest). Added code to keep track of the timeframe, download data from the Kenneth French Library, created the excess returns column for all except the FF5 and rf, and then the main returns series is stored in managed_portfolios

24/03/2026 - Yadhu (approx 11:00): Created the Git repo connecting our Overleaf, R code, and the place to note updates.
For those reading, you need to install the following within VSCode to make this work
For R: R extension, R debugger; within R (do it in R console for ease): install packages languageserver and httpgd
For Latex: Install MiKTeX
You have to run your files locally, so I don't know how you have saved it, but just Chat to figure out Git
At the end of it all, the network-analysis-volatility-management folder in VSCode should have all the files visible, and you should be able to edit and run them locally
Please push your commits regularly, and for development notes, just edit them in GitHub, cause it is easier than having to isolate the file and push commits (cause it might fail sometimes due to some other files being untracked)
But yeah, worst case, just use Chat to get the steps to configure the things on your system

24/03/2026 - Yadhu (12:30): Added library lubridate and created dataframe trading_days to keep track of the number of trading days per month per year
Yadhu (13:40): Removed trading days, cause not needed to calculate RV. Added rv_monthly calculation

24/03/2026 - Kornel (13:00): Added the beggining of Benchmark section with defining estimation sample. SR and Optimal MVE wieghts functions added.
Kornel (17:30): Chenged the section of Benchmarks for 'Functions' with function definitions. Added benchmarks section, where both the MVE strategy as well as equally weighted buy-and-hold are implemented. Still need to verify the code and make sure of few Moreira and Muir approaches, but largerly it should be it.

24/03/2026 - Yadhu (19:00): Added a function that estimates GARCH(1,1) for a single return series. Added the rugarch library (and two others which we do not need any longer, and will remove in next commit).

24/03/2026 - Zhang (19:24): Just downloaded everything and cloned. Failed to download miktek tho so maybe wont be able to code overleaf here but i can do it on overleaf itself. 

24/03/2026 - Yadhu (20:39): Changed the function for the return series to monthly restimation, and added the function to estimate GARCH across all factors. Ran a smaller case for the first 1100 observations (so 596 forecasts) to check if speed is reasonable and it took 5 minutes. For the actual runtime, we should be looking upwards of an hour, but the code does run. Too tired to write the full explanation of the function right now. Will do it in my next update (or not)

24/03/2026 - Yadhu (21:27): Added a function to construct matrices of standardised residual, diagonal matrices and ran check for the matrices and the mean and standard deviation measures. All checks have passed so far (do note for a pretty limited set, but should be working as normal later, will just take longer)

24/03/2026 - Yadhu (23:41): Added most of the DCC functionality manually. We cannot use a DCC package and then apply non-linear shrinkage, because we need to apply the shrinkage to uncondtional covariance/correlation matrix, but using the DCC package does not allow for that, as it treats the covariance structure as part of the likelihood estimation, so we cant really apply it when we need. The functions I have added do all the DCC steps with a placeholder for the non-linear shrinkage (which I will work on tomorrow). 
Also, the code is really crowded right now, I have added as many comments as I think we need right now, but tomorrow, after I have confirmed everything with the functions and checked if they actually generate real results, I will do a writeup on the total functionality.
DISCLAIMER: I did use Chat quite a bit for this step, but most of it is because we are just hard-coding the steps outlined in the Engel, Ledoit, and Wolf (2019) paper. I have checked to ensure that the formulas suggested are the actual formulas used, and the estimation windows are all set up properly. At this point, I can explain what the code does, but I still need to conduct some checks before I can confidently move onto the next step. 
This is definitely the most intense part of the code, and probably will take the longest finish but the rest of the methodology can honestly be finished in one evening.
So hopefully, we can get this part done tomorrow, and I can confidently say that I understand every single part of it. I can guarantee and write up for this part after the non-linear shrinkage part and finalising all functions

24/03/2025 - Yadhu (23:57): Nvm, added the nonlinear shrinkage part. Everything runs for now, but it could all be fake, and I might be hallucinating

25/03/2026 - Kornel (13:10): Created two additional functions (defined in functions section) for summary statistics table. Description available at the beggining of the corresponding functions.