#include <stdio.h>
#include "gmm.h"

void main()
{
	FILE *fp;
	float mfc[39];
	float s;
	GMM gmm;
	
	InitGMM("gmm.bin",&gmm);
	fp = fopen("test","rb");
	while(!feof(fp))
	{
		fread(mfc,1,39*sizeof(float),fp);
		ScoreGMM(mfc,39,&s,&gmm);
		printf("%f\n",s);
	}
	fclose(fp);
	UInitGMM(&gmm);

	return;
}