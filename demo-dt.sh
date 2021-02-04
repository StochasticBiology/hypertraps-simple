date
gcc -o3 hypertraps-dt.c -lm -o hypertraps-dt.ce
./hypertraps-dt.ce fake-cross-samples-1.txt 1 1 5 0
date
gcc -o3 posteriors.c -lm -o posteriors.ce
./posteriors.ce 1 fake-cross-samples-1.txt-posterior-0-1-1-5.txt
date
gnuplot -e "set term svg size 400,400; unset key; set output \"demo-output.svg\"; set xlabel \"Feature\"; set ylabel \"Order\"; plot \"fake-cross-samples-1.txt-posterior-0-1-1-5.txt.process\" u 3:1:(\$5*10) ps variable pt 7; quit;"
