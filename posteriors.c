// posteriors.c
// analyses samples from the posteriors produced by hypertraps-ct code and summarises ordering and time distributions
// note time distributions are only meaningful when the inference was performed using the continuous time options!

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND drand48()

// number of trajectories to simulate for each parameterisation
#define NTRAJ 100

// number of times to call GetRoutes. each call runs NSAMP trajectories and records one sampled route, so the balance of the two controls the number of explicit recorded routes vs the number of sampled trajectories.
#define NSAMP 10

// maximum continuous-time value above which results are truncated
#define MAXCT 1000

// just used in assigning ordinal labels to different features
#define FLEN 15

// simulate trajectories on a given hypercube parameterisation, and store a bunch of summary data about those trajectories
// mean[i] stores the mean acquisition ordering for feature i
// ctrec[MAXCT*i + ref] stores a histogram of acquisitions of feature i at continuous time reference ref
// times[t] stores the continuous time at which feature t is acquired in the first simulated run
// betas[t] stores the exit propensity after feature t is acquired in the first simulated run
// route[t] is the feature changed at step t
void GetRoutes(int *matrix, int len, int ntarg, double *ntrans, int *rec, double *mean, double *ctrec, double *times, double *betas, int *route)
{
  int run, t;
  double time1;
  int *state;
  double totrate;
  double *rate;
  double *cumsum;
  double r;
  int i, j;
  int startt;
  int *checker;
  int numhits;
  double tgap;
  double continuoustime;

  state = (int*)malloc(sizeof(int)*len);
  rate = (double*)malloc(sizeof(double)*len);
  cumsum = (double*)malloc(sizeof(double)*len);
  checker = (int*)malloc(sizeof(int)*ntarg);
  
  for(i = 0; i < ntarg; i++)
    checker[i] = 0;

  for(i = 0; i < len; i++)
    mean[i] = 0;
  
  /* loop through NTRAJ simulated trajectories */
  for(run = 0; run < NTRAJ; run++)
    {
      startt = 0; time1 = 0;

      // start at initial state
      for(i = 0; i < len; i++)
	state[i] = 0;

      // track the (continuous) time elapsed
      // (but continuous time is not interpretable unless the posteriors have been produced in the continuous time paradigm)
      continuoustime = 0;

      // loop through feature acquisitions
      for(t = 0; t < len; t++)
	{
	  totrate = 0;
	  // compute the rate for feature i given the current set of features
	  for(i = 0; i < len; i++)
	    {
	      // ntrans stores the transition matrix: ntrans[0]-ntrans[LEN] are the bare rates; ntrans[LEN+j*LEN+i] is the modifier for i from j
	      if(state[i] == 0)
		{
	          rate[i] = ntrans[i];
	          for(j = 0; j < len; j++)
		    rate[i] += state[j]*ntrans[len+j*len+i];
		  rate[i] = exp(rate[i]);
		}
	      else // we've already lost this gene
		rate[i] = 0;

	      // roulette wheel calculations as normal
	      cumsum[i] = (i == 0 ? 0 : rate[i-1]+cumsum[i-1]);
	      totrate += rate[i];
	    }

	  // choose a step
	  for(i = 0; i < len; i++)
	    cumsum[i] /= totrate;
	  r = RND;
	  continuoustime += (1./totrate)*log(1./r);

#ifdef VERBOSE
	  for(i = 0; i < len; i++)
	    printf("%.2f ", cumsum[i]);
	  printf("\n");
#endif

	  r = RND;
	  for(i = 0; i < len-1; i++)
	    {
	      if(cumsum[i] < r && cumsum[i+1] > r) { break; }
	    }

#ifdef VERBOSE
	  printf("Rolled %f, chose %i\n", r, i);
#endif

	  // we've chosen feature i, at ordering t, and a timescale continuoustime
	  state[i] = 1;
	  mean[i] += t;

	  // rec[t*len + i] increments if we acquire feature i at ordering t
	  // ctrec[MAXCT*i + ref] increments if we acquire feature i at ct-reference ref
	  // pay attention here! we scale continuous times by x100 to produce a reference that allows sensible storage in an integer-referenced histogram, bounded by 0 and MAXCT (element MAXCT-1 stores the number of cases that exceed this)

  	  rec[t*len+i]++;
	  if(continuoustime*100. < MAXCT)
	    ctrec[MAXCT*i+((int)(100.*continuoustime))]++;
	  else
	    ctrec[MAXCT*i + MAXCT-1]++;

	  // sample the statistics of the first simulated run. 
	  if(run == 0)
	    {
	      times[t] = continuoustime;
	      betas[t] = totrate;
	      route[t] = i;
	    }

#ifdef VERBOSE
	  for(i = 0; i < len; i++)
	    printf("%i", state[i]);
	  printf(" (%i)\n", t);
#endif

	}
    }

  for(i = 0; i < len; i++)
    mean[i] /= NTRAJ;

  free(state);
  free(rate);
  free(cumsum);
  free(checker);
}

// construct labels for different features
// for different specific studies this can be adapted to help output
void Label(char *names, int len)
{
  int i;
  for(i = 0; i < len; i++)
    {
      sprintf(&names[i*FLEN], "feature_%i", i);
    }
}

int main(int argc, char *argv[])
{
  FILE *fp;
  int *matrix;
  int len, ntarg;
  double *trans, *ntrans;
  int t;
  int i, j;
  char ch;
  double lik, nlik;
  int *rec, *order;
  double *drec, *sortdrec, *mean;
  int maxt, allruns;
  int seed = 0;
  char str[200];
  char shotstr[200];
  double tmp;
  int change;
  char names[2000];
  int expt;
  int count;
  double *meanstore, *fmeanstore;
  int specref;
  double *ctrec, ctnorm;
  double *times, *betas;
  int *route;
  FILE *fp1, *fp2, *fp3;
  char fstr[200];
  int tlen;
  int verbose;
  int fref;

  // process command-line arguments
  if(argc < 3)
    {
      printf("Usage:\n posteriors.ce [verbose flag] [posterior sample file(s)]\n");
      return 0;
    }
  verbose = atoi(argv[1]);
  printf("Verbose flag is %i\n", verbose);

  // open posterior file to assess format
  fp = fopen(argv[2], "r");
  if(fp == NULL)
    {
      printf("Couldn't open file %i: %s\n", 1, argv[2]);
      return 0;
    }
  tlen = 0;
  do{
    ch = fgetc(fp);
    if(ch == ' ') tlen++;
  }while(ch != '\n' && ch != EOF);
  fclose(fp);

  if(ch == EOF)
    {
      printf("Couldn't find appropriate samples in file %s\n", argv[2]);
      return 0;
    }

  // figure out if posterior file is presented in L*(L+1) format; get L if so
  len = 0;
  for(i = 1; i < 200; i++)
    {
      if(tlen == i*(i+1))
	{
	  len = i;
	  break;
	}
    }
  if(len == 0)
    {
      printf("Couldn't determine number of features from %s\n", argv[2]);
      return 0;
    }

  printf("Based on %s , there are %i features\n", argv[2], len);


  // initialise and allocate a lot of different arrays to compute and store statistics
  srand48(121+seed);
  allruns  =0;
  ntarg = 0;
  Label(names, len);
  
  matrix = (int*)malloc(sizeof(int)*10000);
  ctrec = (double*)malloc(sizeof(double)*MAXCT*len);
  times = (double*)malloc(sizeof(double)*len);
  betas = (double*)malloc(sizeof(double)*len);
  route = (int*)malloc(sizeof(int)*len);

  trans = (double*)malloc(sizeof(double)*len*(len+1)); /* transition matrix */
  ntrans = (double*)malloc(sizeof(double)*len*(len+1));
  rec = (int*)malloc(sizeof(int)*len*len); /* stores step ordering, modified by getlikelihood */
  mean = (double*)malloc(sizeof(double)*len);
  meanstore = (double*)malloc(sizeof(double)*len);
  fmeanstore = (double*)malloc(sizeof(double)*len);
  order = (int*)malloc(sizeof(int)*len);
  drec = (double*)malloc(sizeof(double)*len*len);
  sortdrec = (double*)malloc(sizeof(double)*len*len);

  // initialise
  for(i = 0; i < MAXCT*len; i++)
    ctrec[i] = 0;
  ctnorm = 0;

  for(i = 0; i < len*len; i++)
    rec[i] = 0;

  for(i = 0; i < len; i++)
    fmeanstore[i] = 0;

  // set up file outputs
  if(verbose)
    {
      sprintf(fstr, "%s-routes.txt", argv[2]);
      fp1 = fopen(fstr, "w");
      sprintf(fstr, "%s-betas.txt", argv[2]);
      fp2 = fopen(fstr, "w");
      sprintf(fstr, "%s-times.txt", argv[2]);
      fp3 = fopen(fstr, "w");
    }

  // loop through specified posterior files
  printf("Trying %i files\n", argc-2);
  for(fref = 2; fref < argc; fref++)
    {
      // try to open this file
      fp = fopen(argv[fref], "r");
      if(fp == NULL)
	{
	  printf("Couldn't open file %i: %s\n", fref, argv[fref]);
	  return 0;
	}
      count = 0;
      printf("Working on file %i: %s\n", fref-1, argv[fref]);

      // loop through posterior samples in this file
      while(!feof(fp))
	{
	  // read in single posterior sample
	  for(i = 0; i < len*(len+1); i++)
	    fscanf(fp, "%lf", &ntrans[i]);

	  // this if statement controls which samples get processed
	  // if we want to include burn-in or subsampling, can put it here
	  if(!feof(fp))// && count > 100 && count % 5 == 0)
	    {
	      // loop through iterations
	      for(j = 0; j< NSAMP; j++)
		{
		  for(i = 0; i < len; i++)
		    meanstore[i] = 0;
		  // simulate behaviour on this posterior and add statistics to counts and histograms
		  GetRoutes(matrix, len, ntarg, ntrans, rec, meanstore, ctrec, times, betas, route);
		  for(i = 0; i < len; i++)
		    fmeanstore[i] += meanstore[i];
		  ctnorm += NTRAJ;
		  allruns++;

		  if(verbose)
		    {
		      for(i = 0; i < len; i++)
 		        fprintf(fp1, "%i ", route[i]);
		      for(i = 0; i < len; i++)
			fprintf(fp2, "%.15f ", betas[i]);
		      for(i = 0; i < len; i++)
			fprintf(fp3, "%.3f ", times[i]);
		      fprintf(fp1, "\n");
		      fprintf(fp2, "\n");
		      fprintf(fp3, "\n");
		    }
		}
	    }
	  count++;
	}
      fclose(fp);
      if(verbose)
	{
	  fclose(fp1);
	  fclose(fp2);
	  fclose(fp3);
	} 
    }

  // output various summaries
  for(i = 0; i < len; i++)
    printf("%i %f\n", i, fmeanstore[i]/allruns);

  // compute mean orderings
  // rec[t*len+i] is prob of obtaining i at time t

  for(i = 0; i < len*len; i++)
    drec[i] = (double)rec[i]/(allruns*NTRAJ);

  for(i = 0; i < len; i++)
    {
      mean[i] = 0;
      order[i] = i;
      for(t = 0; t < len; t++)
	mean[i] += t*drec[t*len+i];
    }

  // simple bubble sort orders features by mean acquisition order
  do{
    change = 0;
    for(i = 0; i < len-1; i++)
      {
	if(mean[i] > mean[i+1])
	  {
	    tmp = mean[i]; mean[i] = mean[i+1]; mean[i+1] = tmp;
	    tmp = order[i]; order[i] = order[i+1]; order[i+1] = tmp;
	    change = 1;
	  }
      }
  }while(change == 1);
  seed--;

  // output the set of summary statistics
  // rec[t*len+i] is prob of obtaining i at time t

  // this produces the heatmap of acquisition probability by feature and order
  // outputs both the original feature ordering and the above mean-sorted references
  sprintf(str, "%s.process", argv[2]);
  fp = fopen(str, "w");
  for(t = 0; t < len; t++)
    {
      for(i = 0; i < len; i++)
	fprintf(fp, "%i %i %i %s %.15f\n", t, i, order[i], &names[FLEN*order[i]], drec[t*len+order[i]]);
      fprintf(fp, "\n");
    }

  // these appended comments give gnuplot commands for axis labels if required: both in original and mean-sorted orderings
  fprintf(fp, "# set xtics (");
  for(i = 0; i < len; i++)
    fprintf(fp, "\"%s\" %i%c", &names[FLEN*order[i]], i, (i == len-1 ? ')' : ','));
  fprintf(fp, "\n");
  fprintf(fp, "# default-order set xtics (");
  for(i = 0; i < len; i++)
    fprintf(fp, "\"%s\" %i%c", &names[FLEN*i], i, (i == len-1 ? ')' : ','));
  fprintf(fp, "\n");

  fprintf(fp, "# (");
  for(i = 0; i < len; i++)
    fprintf(fp, "%i, ", order[i]);
  fprintf(fp, ")\n");

  fprintf(fp, "# ");
  for(i = 0; i < len; i++)
    fprintf(fp, "%s %.4f, ", &names[FLEN*order[i]], mean[i]);
  fprintf(fp, ")\n");

  fclose(fp);

  // this stores the time histograms associated with acquisition times for each feature
  // remember here that we've scaled by x100 to store in an integer-referenced array (see GetRoutes())
  sprintf(str, "%s.ctrec.process", argv[2]);
  fp = fopen(str, "w");
  for(i = 0; i < len; i++)
    {
      tmp = 0;
      for(j = 0; j < MAXCT; j++)
	{
	  fprintf(fp, "%i %i %.6f\n", i, j, ctrec[MAXCT*i+j]/ctnorm);
	  tmp += ctrec[MAXCT*i+j]*j;
	}
      printf("%i %.4f\n", i, tmp/ctnorm);
      fprintf(fp, "\n");
    }

  return 0;
}
