set datafile separator ";"

set multiplot layout 2,1

set xlabel "t, s"
set ylabel "x, m"
set grid
plot 'fixed.csv' using 1:2 with lines linewidth 2 title 'AN', \
     'fixed.csv' using 1:3 with lines linewidth 2 title 'EE1', \
     'fixed.csv' using 1:5 with lines linewidth 2 title 'VV2', \
     'fixed.csv' using 1:7 with lines linewidth 2 title 'RK4', \
     'dp54.csv' using 1:3 with lines linewidth 2 title 'DP54', \

set xlabel "t, s"
set ylabel "|E|, m"
set grid
plot 'fixed.csv' using 1:4 with lines linewidth 2 title '|E_E_E_1|', \
     'fixed.csv' using 1:6 with lines linewidth 2 title '|E_V_V_2|', \
     'fixed.csv' using 1:8 with lines linewidth 2 title '|E_R_K_4|', \
     'dp54.csv' using 1:4 with lines linewidth 2 title '|E_D_P_5_4|'
