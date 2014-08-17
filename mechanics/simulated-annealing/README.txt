How the scripts be run:
The main script that should be run is mainprog.m
In the MATLAB Command Window, change directory to
where the MATLAB scripts are located by:

>> cd <some_directory>/simulated-annealing

in the command window.

Run mainprog by typing the following in the command window:

>>  mainprog(<choice>)

where <choice> is:
1: For analysis on the bump function

2: For analysis on the 10 bar structure using SA

For choice = 1:
A figure will be generated with average optimal f vs. cooling
parameter for different epsilon. The goal is to find the best
cooling parameter and epsilon for this particular SA.

For choice = 2:
The average minimum weight of the structure will be displayed.
A figure will appear showing the convergence trend using one-pass,
and quadratic penalty for different values for P. Statistical values
are also shown.

NOTE: The run time of this program is very slow.