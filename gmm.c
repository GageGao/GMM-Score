#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include "gmm.h"


int InitGMM(char *fn,GMM *pgmm)
{
	FILE *fp;
	char *p;
	int i,j;

	fp = fopen(fn,"rb");
	if ( NULL == fp )
		return ERROR_OPENFILE;

	fread(&(pgmm->m),1,sizeof(int),fp);
	fread(&(pgmm->l),1,sizeof(int),fp);

	p = (char *)malloc((2*pgmm->m+2*pgmm->m*pgmm->l)*sizeof(float));
	if ( NULL == p )
		return ERROR_MALLOC;
	fread(p,1,(2*pgmm->m+2*pgmm->m*pgmm->l)*sizeof(float),fp);
	pgmm->wt = (float *)p;
	pgmm->mean = (float *)(p+pgmm->m*sizeof(float));
	pgmm->var = (float *)(p+pgmm->m*sizeof(float)+pgmm->m*pgmm->l*(sizeof(float)));
	pgmm->gconst = (float *)(p+pgmm->m*sizeof(float)+2*pgmm->m*pgmm->l*(sizeof(float)));

	for ( i = 0; i < pgmm->m; i++ )
	{
		pgmm->wt[i] = log(pgmm->wt[i]);
	}

	for ( i = 0; i < pgmm->m; i++ )
	{
		for ( j = 0; j < pgmm->l; j++ )
		{
			pgmm->var[i*pgmm->l+j] = 1.0/pgmm->var[i*pgmm->l+j];
		}
	}

	fclose(fp);
	return SUCCESS_END;
}

int UInitGMM(GMM *pgmm)
{
	free(pgmm->wt);
}

float log_wgd(const int m,const float *dat,GMM *pgmm)
{
	int l, ll;
	float sum, *diff = NULL, tmp, lwgd;

	sum = pgmm->gconst[m];

	for (l = 0; l < pgmm->l; l++) {
		tmp = dat[l] - pgmm->mean[m*pgmm->l+l];
		sum += (tmp * tmp) * pgmm->var[m*pgmm->l+l];
	}

	lwgd = pgmm->wt[m] - 0.5 * sum;
	return (lwgd);
}

float log_add(float logx, float logy)
{
	float swap, diff, minLogExp, z;

	if (logx < logy) {
		swap = logx;
		logx = logy;
		logy = swap;
	}

	diff = logy - logx;
	minLogExp = -log(-LZERO);

	if (diff < minLogExp)
		return ((logx < LSMALL) ? LZERO : logx);
	else {
		z = exp(diff);
		return (logx + log(1.0 + z));
	}
}

int ScoreGMM(float *dat,int l,float *s,GMM *pgmm)
{
	int i;
	float logwgd, logb;

	if ( l != pgmm->l )
		return ERROR_INVALID;
	
	for (i = 0, logb = LZERO; i < pgmm->m; i++) {
		logwgd = log_wgd(i,dat,pgmm);
		logb = log_add(logb, logwgd);
	}
	*s = logb;

	return SUCCESS_END;
}