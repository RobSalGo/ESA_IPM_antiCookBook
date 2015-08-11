OPTIONS FORMCHAR="|----|+|---+=|-/\<>*"; * fixes output format;

proc iml; title 'subroutine construction';
* starts IML programming language;

*------CREATE SUBROUTINES--------------------------------------;

* INSERT YOUR OWN, CORRECT PARAMETER VALUES, CUT-AND-PASTED FROM THE STATISTICAL ANALYSES;

* SURVIVAL SUBROUTINE (named surv);
 
* probability of surviving (psurv) as a function of size;
* size1 is input, psurv is output;
* target & adj are inputs for the elasticity analysis;
* a binomial distribution and logit link function are used;
* recall that the logit function a = logit(p) = log(p/(1-p));
* recall that the back-transformation is p = exp(a)/((exp(a)+1);  

start surv(size1,psurv,target,adj);           * begins subroutine;
  logsize1 = log(size1);
  * predict logit(survival probability);
  logitpsurv = -2.5 + 1.7*logsize1;  * replace with your own values;
  * backtransform predicted value to a probability;
  explogit = exp(logitpsurv);
  psurv = explogit / (1 + explogit);    
  if (target='surv') then psurv = psurv*adj;  * for elasticity;
finish;                                       * ends subroutine;

* GROWTH SUBROUTINE (named growth);

* this subroutine calculates the probability (pgrow) of a plant of       size1 that survives becoming a plant of size2;
* the sum of probabilities for a given value of size1 must add to 1;
* i.e., pgrow = transition probability from size1 to size2;
* tmax is the maximum size (and the dimension of the IPM matrix)
* size1 and size2 are input, pgrow is output;
* target & adj are inputs for the elasticity analysis;
* distribution is truncated negative binomial, with parameters mu and k (terminology follows Bolker);
* values of 0 are not possible, so all the probabilities have to be adjusted accordingly - this is the truncation;
* both mu and k are functions of initial size;
* I use SAS built-in negative binomial cdf (cumulative probability), named probnegb;
* probnegb is parameterized in n & p, not mu and k, so I have to re-parameterize from mu and k to n and p;
* note that in effect the function uses log(size 1) and log(size2);

start growth(size1,size2,pgrow,target,adj,tmax);
  logsize1 = log(size1);
  mu = exp(0.8 + 0.8*logsize1);  * replace with your own values; 
  if (target='grow') then mu = mu * adj; 
  k = 2.5 + 0.1*size1;            * replace with your own values;     
  var = mu + (mu**2)/k;
  p = mu/var; n = mu*p/(1-p);
  if size2 < tmax then do;
  	prob = probnegb(p,n,size2)- probnegb(p,n,size2-1);
                       end;
  * any predicted sizes > tmax are reset to tmax, with corresponding adjustment in the probability of entering size class tmax;
  if size2 = tmax then do; 
	prob = 1 - probnegb(p,n,tmax-1); 
                       end;
  * this is the adjustment for absent 0 values;
  * to do it, divide by (1 - P(zero));
  prob0 = probnegb(p,n,0);       * probablity of getting 0;
  probtrunc = prob/(1-prob0);    * adjusts all other probabliities;
  pgrow = probtrunc; 
finish;

* PROBABILITY OF REPRODUCING SUBROUTINE (named repro);

* probability of reproducing is a function of initial size;
* size1 is input, prep is output;
* target & adj are inputs for the elasticity analysis;

start repro(size1,prep,target,adj);
  logsize1 = log(size1); 
  logitprep = -7 + 2*logsize1; * replace with your own values; 
  explogit = exp(logitprep);
  prep = explogit / (1 + explogit);
  if (target='prep') then prep = prep*adj;
finish;

* SEEDSET SUBROUTINE (named seedset);

* seedset/reproductive plant is a function of initial size;
* size1 is input, fec is output;
* target & adj are inputs for the elasticity analysis;

start seedset(size1,fec,target,adj);
  logsize1 = log(size1);
  mu = 0.5 + 0.6*logsize1;  * replace with your own values; 
  fec = exp(mu);
  if (target='seed') then fec = fec*adj;
finish;

* SEED-TO-NEW-RECRUIT SUBROUTINE (named seedtorec);

* probability of survival from seed to recruit;
* a fixed value is used, output as psr;
* target & adj are inputs for the elasticity analysis;

start seedtorec(psr,target,adj);
  psr = 0.3;                * replace with your own values; 
  if (target='stor') then psr = psr*adj;
finish;

* SIZE DISTRIBUTION OF RECRUITS SUBROUTINE (named recsize);

* recsize is input (a given size of plants newly recruited to the population);
* output is precsize = probability of a new recruit being the give size; 
* target & adj are inputs for the elasticity analysis;
* a truncated negative binomial distribution is used;
* terminology again follows Bolker, as does re-parameterization from mu and k to n and p;
* instead of the cdf function used in the growth subroutine, I use the pdf (probability density function) for the negative binomial;
* truncation to adjust for no 0 values as in the growth subroutine;

start recsize(recsize,precsize,target,adj);
  k = 0.25;  mu = 0.44;            * replace with your own values;  
* mu could be a function of treatment, etc. in another model;
  if (target='srec') then mu = mu*adj;
  * convert to n,p parameterization, following Bolker;
  n=k;  p = k / (k+mu);
  rawp = pdf('NEGBIN',recsize,p,n);
  prob0 = pdf('NEGBIN',0,p,n);
  probtrunc = rawp / (1 - prob0);
  precsize = probtrunc;
finish;

* SAVE THE SUBROUTINES IN TEMPORARY STORAGE;

* the temporary storage catalog is called WORK.IMLSTOR;

store module = (surv growth repro seedset recsize seedtorec);

* summary of subroutines

subroutine name		requires					returns
surv					size1,target,adj			psurv
growth				size1,size2,target,adj		pgrow
repro				size1,target,adj			prep
seedset				size1,target,adj			fec
recsize				trecsize,target,adj			precsize
seedtorec				target,adj				psr;

* leave IML, returning to regular SAS;
quit;

proc iml; title 'IPM construction, elastiticies';
* starts IML programming language;

* load the subroutines;

load module = (surv growth repro seedset recsize seedtorec growth);

* summary of subroutines

subroutine name		requires					returns
surv					size1,target,adj			psurv
growth				size1,size2,target,adj		pgrow
repro				size1,target,adj			prep
seedset				size1,target,adj			fec
recsize				trecsize,target,adj			precsize
seedtorec				target,adj				psr;


*----------SET THE OVERALL PARAMETERS--------------------------;
tmax = 200;              * maximum size, also matrix dimensions;
print 'matrix dimensions', tmax;

* INSERT YOUR OWN, CORRECT BASELINE LAMBDA, CUT-AND-PASTED FROM THE STATISTICAL ANALYSES;
baselambda =  1.05723;   * insert correct value of baseline lambda      here, obtained from lambda when adj=1.0 in preliminary output;

* parameters for the elasticity analysis;

* changes in variable values are stored in adjvector;
* 1.05 is a 5% increase, 0.95 is a 5% decrease, etc;
adjvector = {0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2};  	
nadj = nrow(adjvector);       * # of elements of adjvector;


* # of different variables whose elasticities will be calculated;
nparms = 6;
finalnrows = nparms * nadj;

* alternative: baseline model only; * set nparms = 1;
target = 'null';
	
* elasticities are calculated for only one variable at a time;
* the variable being used is determined by targnumber; 

do targnumber = 1 to nparms;

  if targnumber = 1 then target = 'surv';  * P(survival);
  if targnumber = 2 then target = 'grow';  * mean predicted size;
  if targnumber = 3 then target = 'prep';  * P(reproducing);
  if targnumber = 4 then target = 'seed';  * seedset;
  if targnumber = 5 then target = 'stor';  * P(seed to recruit);
  if targnumber = 6 then target = 'srec';  * mean size-new recruits;

* if you want only the baseline model, set nparms = 1, above, and make all 6 "if targnumber = ...." lines comments by putting * in front of each one;  

* create an output matrix to save the lambdas of this target;

elasmat = j(nadj,4,0);

do elascounter = 1 to nadj;         * select elasticity adjustment;

  adj = adjvector[elascounter];

* GROWTH MATRIX, G ;

G = j(tmax,tmax,0);
do i = 1 to tmax;     * final size = matrix row;
  do j = 1 to tmax;   * initial size = matrix column;
    size1 = j; size2 = i; 
    call surv(size1,psurv,target,adj);  
       if (psurv > 1.0) then psurv = 1;
	call growth(size1,size2,pgrow,target,adj,tmax);
	   if (size2 < 0) then size2 = 0.1;
	gij = psurv*pgrow;  
	G[i,j] = round(gij, 0.0000001);
end; end;
* print G;

* REPRODUCTION MATRIX, R;

R = j(tmax,tmax,0);
do i = 1 to tmax;     * final size = matrix row;
  do j = 1 to tmax;   * initial size = matrix column;
    size1 = j; recsize = i; 
	call repro(size1,prep,target,adj); 
	  if (prep > 1.0) then prep = 1;
	call seedset(size1,fec,target,adj); 
	  if (fec < 0) then fec = 0.1;
     call seedtorec(psr,target,adj);     
	call recsize(recsize,precsize,target,adj);
	rij = prep*fec*psr*precsize;
	R[i,j] = round(rij, 0.0000001);      
end; end;
* print R;

* COMBINE G & R TO MAKE L;

L = G + R;                    * element-by-element matrix addition;
* print L;

* CALCULATE LAMBDA & CHANGE IN LAMBDA;

call eigen(evals,evecs,L);
* print evals; 
lambda = evals[1,1];    * print elascounter,lambda;
* one row of evals is one eigenvalue;
* eigenvalues listed in order of size, real, then imaginary parts;

* principalReigenvector = evecs[,1]
* it's not used here, but column 1 of evecs contains the principal 
right eigenvector, which, once divided by its sum, gives the stable age distribution;

lchange = (lambda-baselambda)/baselambda;

elasmat[elascounter,1] = elascounter;
elasmat[elascounter,2] = adj;
elasmat[elascounter,3] = lambda;
elasmat[elascounter,4] = lchange;

end;                                  * end adjustment loop;

* print elasmat;

* save the eigenvalue, eigenvalue change, etc.;

* matrix of numerical values;
if (targnumber = 1) then numout = elasmat;
if (targnumber > 1) then numout = numout//elasmat;

* matrix of character values, holding the target names;
dummymat = j(nadj,1,target);    
if (targnumber = 1) then charmat = dummymat;
if (targnumber > 1) then charmat = charmat//dummymat;

end;                                 * end target loop;

* CREATE REGULAR SAS DATA SETS;

* numerical output;
cols1 = {elascount, adj, lambda, lchange}; * new variable names;  
create datnumout from numout [colname=cols1];  
append from numout; 

cols2 = {target};
create datcharout from charmat [colname=cols2];
append from charmat;

quit;                                * leave iml;

* proc print data=datnumout; title 'numerical output';
* proc print data=datcharout; title 'character output';

data datelasticities; 
	merge datnumout datcharout;  * combine the two data sets;
proc print data=datelasticities; 
	title 'elasticities';
run;


