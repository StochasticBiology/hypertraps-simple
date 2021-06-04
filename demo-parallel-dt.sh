# this script runs and analyses three parallel chains so that convergence can (very roughly) be visually assessed. consider Gelman-Rubin or similar for formal work!

date
gcc -o3 hypertraps-dt.c -lm -o hypertraps-dt.ce
./hypertraps-dt.ce fake-cross-samples-1.txt 1 1 5 0
./hypertraps-dt.ce fake-cross-samples-1.txt 2 1 5 0
./hypertraps-dt.ce fake-cross-samples-1.txt 3 1 5 0
date
gcc -o3 posteriors.c -lm -o posteriors.ce
./posteriors.ce 1 fake-cross-samples-1.txt-posterior-0-1-1-5.txt
./posteriors.ce 1 fake-cross-samples-1.txt-posterior-0-2-1-5.txt
./posteriors.ce 1 fake-cross-samples-1.txt-posterior-0-3-1-5.txt
date
gnuplot -e "set term svg size 400,400; unset key; set xrange [-0.9:4.9]; set yrange [-0.9:4.9]; set xtics 1; set ytics 1; set output \"demo-parallel-output.svg\"; set xlabel \"Feature\"; set ylabel \"Order\"; myps = 0.5; plot \"fake-cross-samples-1.txt-posterior-0-1-1-5.txt.process\" u 3:1:(sqrt(\$5)*myps):(0):(120) w circles, \"fake-cross-samples-1.txt-posterior-0-2-1-5.txt.process\" u 3:1:(sqrt(\$5)*myps):(120):(240) w circles, \"fake-cross-samples-1.txt-posterior-0-3-1-5.txt.process\" u 3:1:(sqrt(\$5)*myps):(240):(360) w circles; quit;"
