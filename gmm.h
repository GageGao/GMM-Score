#ifndef GMM_H
#define GMM_H

#define ERROR_OPENFILE	(0)
#define ERROR_MALLOC	(1)
#define ERROR_INVALID	(2)
#define SUCCESS_END		(3)

#define LZERO		(-1.0e+10)
#define LSMALL		(-0.5e+10)

typedef struct{
	int m;
	int l;
	float *wt;
	float *mean;
	float *var;
	float *gconst;
}GMM;

int InitGMM(char *fn,GMM *pgmm);
int UInitGMM(GMM *pgmm);

int ScoreGMM(float *dat,int l,float *s,GMM *pgmm);

#endif	//GMM_H