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
# simulate example data 
data(corr)
data <- corrSim(corr)
reliab <- c(.90, .89, .96, .86, .84, .82, .85, .76, .89, .68, .86)

# view the RelComp documentation
?RelComp

# run RelComp on the simulated data
RelComp(data, names(data), reliab)
```

## Real Data Example

-- pending --