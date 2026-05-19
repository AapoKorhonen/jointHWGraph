# jointHWGraph, joint Hierachical Wishart prior for graphical models




# Hyperparameters
## Nu
We recommend a default value of p*10, p is the number of variables
## Delta
We recommend a default value of 1, but increasing it to p*10 can also be used.
## Target matrix B
On default, an indentity matrix I is used. 

# Edge selection

## Optimal threshold based on approximated F1 value
limit_FN controls, if the estimated FN is limimted to be at least 0. If limit_FN = F, then estimated FN can be negative. When limit_FN is used, the expected_number_of_connections is more informative and behaves more like upper limit for connectios to be allowed for the estimated network.  
## fdrtool, minimal false non discovery rate

## Target FDR

# Missing data values

