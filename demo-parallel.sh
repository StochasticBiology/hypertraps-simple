# this script runs and analyses three parallel chains so that convergence can (very roughly) be visually assessed. consider Gelman-Rubin or similar for formal work!

date
gcc -o3 hypertraps.c -lm -o hypertraps.ce
./hypertraps.ce fake-cross-samples-1.txt 1 1 5 0
./hypertraps.ce fake-cross-samples-1.txt 2 1 5 0
./hypertraps.ce fake-cross-samples-1.txt 3 1 5 0
date
gcc -o3 posteriors.c -lm -o posteriors.ce
./posteriors.ce 1 fake-cross-samples-1.txt-posterior-0-1-1-5.txt
./posteriors.ce 1 fake-cross-samples-1.txt-posterior-0-2-1-5.txt
./posteriors.ce 1 fake-cross-samples-1.txt-posterior-0-3-1-5.txt
date
