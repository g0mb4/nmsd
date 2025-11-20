set datafile separator ";"

set multiplot layout 2,1

plot 'fixed.csv' using 1:2 with lines title 'AN', \
     'fixed.csv' using 1:3 with lines title 'EE1', \
     'fixed.csv' using 1:5 with lines title 'VV2'

plot 'fixed.csv' using 1:4 with lines title '|E_E_E_1|', \
     'fixed.csv' using 1:6 with lines title '|E_V_V_2|'
