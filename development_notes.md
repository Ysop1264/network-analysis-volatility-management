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
