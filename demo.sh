date
gcc -o3 hypertraps.c -lm -o hypertraps.ce
./hypertraps.ce fake-cross-samples-1.txt 1 1 5 0
date
gcc -o3 posteriors.c -lm -o posteriors.ce
./posteriors.ce 1 fake-cross-samples-1.txt-posterior-0-1-1-5.txt
date
Rscript plot-process.R fake-cross-samples-1.txt-posterior-0-1-1-5.txt.process
