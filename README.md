# Picewise Linear
The code available in this repository works with R and RStudio. It is designed to produce two types of analysis. First, it considers the cumulative number of events of a particular type (for example, the number of protests in some city). Then, based on the assumption of a time-varying rate, detect the fluctuations of the expected number of daily events.

One of that types of events is terrorist attacks. Here, data from ACLED is used as an example to analyse the daily rate of Boko Haram events.
More data is available at: https://acleddata.com/
Raleigh, Clionadh, et al. "Introducing ACLED: An armed conflict location and event dataset." Journal of peace research 47.5 (2010): 651-660.

The code output is a time series with the time-varying rate of events, an estimate of the number of shocks, and identifying the significant shocks, where the rate changes by more than 20% in less than one month.

## Code structure
The code produces two functions: piecewise_linear and piecewise_linear_plotter
The repository contains a set of events in which the functions can be tested and the number of shocks observed.

### piecewise_linear 
The function takes three inputs:
dates - Date of each event in dd/mm/yyyy format
counts - The count for each date to generate a cumulative sum
max_pieces - The maximum number of breakpoints to consider in analysis, be wary of using too many breakpoints

the function gives a list containing
- Aggregated dataset
- breakpoints found (date)
- gradients between breakpoints
- piecewise linear models generated

### piecewise_linear_plotter
The function takes inputs:
pwlf, the output from piecewise_linear()
breakpoints_to_plot=2, which number of breakpoints to plot
legend.position=c(.85,.25) - where to position legend, 'none' for no legend
breakpoint_lines=TRUE

The outputs are:
A figure with the cumulative curve of events and the curve where the breakpoints are identified.

