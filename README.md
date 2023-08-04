# Picewise Linear
The code available in this repository works with R and RStudio. It is designed to produce two types of analysis. First, it considers the cumulative number of events of a particular type (for example, the number of protests in some city). Then, based on the assumption of a time-varying rate, detect the fluctuations of the expected number of daily events.

One of that types of events are terrorism attacks. Here, data from ACLED is used as an example to analyse the daily rate of Boko Haram events.
More data is available at: https://acleddata.com/
Raleigh, Clionadh, et al. "Introducing ACLED: An armed conflict location and event dataset." Journal of peace research 47.5 (2010): 651-660.

## Code structure
The code produces two functions: piecewise_linear and piecewise_linear_plotter

### piecewise_linear 
The function takes three inputs:
dates - Date of each event in dd/mm/yyyy format
counts - The count for each date to generate cumulative sum
max_pieces - The maximum number of breakpoints to consider in analysis, be wary of using too many break points

the function gives a list containing
- Aggregated dataset
- breakpoints found (date)
- gradients between break points
- piecewise linear models generated

### piecewise_linear_plotter
The outputs are:
A time series with the time-varying rate of events
An estimate of the number of shocks
Identification of the significant shocks, where the rate changes more than 20% in less than one month.

The repository contains a set of events in which the picewise linear functions can be tested and the number of shocks observed.


