RCA of Cliff &amp; Caruso 1998 <img src="man/figures/200px-Rti-logo.png" align="right" />
========================================================

RCA was developed by RTI International to implement the Reliable Component 
Analysis (RCA) of Cliff &amp; Caruso (1998).

## Installation

Go to [https://github.com/ICTatRTI/RCA/tags](https://github.com/ICTatRTI/RCA/tags)
and click on the top most `tar.gz` link to download the newest version of the `RCA` package.


In R, use

```
install.packages("path/to/RCA-0.1.0.tar.gz", repos = NULL, type = "source")
```

substituting the path where you saved the `tar.gz` file. 

## Simulated Data Example

```
# Replicate the analyses of Cliff and Caruso (1998)

library(RCA)

# simulate data using the Cliff and Caruso (1998) correlation matrix
data(corr)
data <- corrSim(corr)
reliab <- c(.90, .89, .96, .86, .84, .82, .85, .76, .89, .68, .86)

# view the RelComp documentation
?RelComp

# run RelComp on the simulated data
output <- RelComp(data, names(data), reliab)

# view the reliabilities, compare to the first row in Table 8 of 
# Cliff and Caruso (1998)
round(output$reliabilities[1:5], 3)
```

## Real Data Example

-- pending --
