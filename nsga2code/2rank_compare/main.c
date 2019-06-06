/*compare fitness of last generation of two operations 3-4、3-5...3-9*/
/*4-5、4-6...4-9...8-9 with non dominated sorting and make new rank */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define population 50
#define minvar 3
#define maxvar 9
#define nfunc 3
int twopop = population * 2;
typedef struct individual
{
    double fitness[3];
    int varn;
    int flag;
    int rank;
	int rankno[100];
} indiv;
int indcmp(double p1f[3],double p2f[3]);

int main()
{
    int now = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, q = 0, tpa = 0 , tpb = 0,
        rnk,
        val,
        nondom,
        maxrank1,
        rankarr[twopop];


    float *ptr1,*ptr2;
    double s1=0.0,s2=0.0,per1=0.0,per2=0.0;

    FILE *f1 = fopen("C:/Users/Josh/Desktop/nsga2/nsga2code/opt3_o.txt","r");
    FILE *f2 = fopen("C:/Users/Josh/Desktop/nsga2/nsga2code/opt4_o.txt","r");
    FILE *f3 = fopen("C:/Users/Josh/Desktop/nsga2/nsga2code/opt5_o.txt","r");
    FILE *f4 = fopen("C:/Users/Josh/Desktop/nsga2/nsga2code/opt6_o.txt","r");
    FILE *f5 = fopen("C:/Users/Josh/Desktop/nsga2/nsga2code/opt7_o.txt","r");
    FILE *f6 = fopen("C:/Users/Josh/Desktop/nsga2/nsga2code/opt8_o.txt","r");
    FILE *f7 = fopen("C:/Users/Josh/Desktop/nsga2/nsga2code/opt9_o.txt","r");
	FILE *f8 = fopen("compare.txt","w+");

    indiv opt[maxvar-minvar+1][population];

    while((!feof(f1))&& now < population)
    {
        fscanf(f1,"%lf %lf %lf %d ",&opt[0][now].fitness[0],&opt[0][now].fitness[1],&opt[0][now].fitness[2],&opt[0][now].varn);
        fscanf(f2,"%lf %lf %lf %d ",&opt[1][now].fitness[0],&opt[1][now].fitness[1],&opt[1][now].fitness[2],&opt[1][now].varn);
        fscanf(f3,"%lf %lf %lf %d ",&opt[2][now].fitness[0],&opt[2][now].fitness[1],&opt[2][now].fitness[2],&opt[2][now].varn);
        fscanf(f4,"%lf %lf %lf %d ",&opt[3][now].fitness[0],&opt[3][now].fitness[1],&opt[3][now].fitness[2],&opt[3][now].varn);
        fscanf(f5,"%lf %lf %lf %d ",&opt[4][now].fitness[0],&opt[4][now].fitness[1],&opt[4][now].fitness[2],&opt[4][now].varn);
        fscanf(f6,"%lf %lf %lf %d ",&opt[5][now].fitness[0],&opt[5][now].fitness[1],&opt[5][now].fitness[2],&opt[5][now].varn);
        fscanf(f7,"%lf %lf %lf %d ",&opt[6][now].fitness[0],&opt[6][now].fitness[1],&opt[6][now].fitness[2],&opt[6][now].varn);
        opt[0][now].flag = -1;//var ==3
        opt[0][now].rank = -1;
        opt[1][now].flag = -1;//var ==4
        opt[1][now].rank = -1;
        opt[2][now].flag = -1;//var ==5
        opt[2][now].rank = -1;
        opt[3][now].flag = -1;//var ==6
        opt[3][now].rank = -1;
        opt[4][now].flag = -1;//var ==7
        opt[4][now].rank = -1;
        opt[5][now].flag = -1;//var ==8
        opt[5][now].rank = -1;
        opt[6][now].flag = -1;//var ==9
        opt[6][now].rank = -1;
        now+=1;
    }


    //start comparing
    /*Initializing the ranks to zero*/
    rnk = 0 ;

    nondom = 0 ;
    maxrank1 = 0;
    indiv pop[twopop];

    for(m=minvar; m<maxvar; m++) //m:3~8
    {
        for(n=m+1; n<=maxvar; n++) //n:4.5.6.7.8~9
        {
			rnk = 0;
			tpa = 0;
			tpb = 0;
			s1 = 0.0;
			s2 = 0.0;
			per1 = 0.0;
			per2 = 0.0;
            for(now=0; now<twopop; now++)
            {
                //combine two var indiv
                if(now>=population){
					pop[now].fitness[0] = opt[n-minvar][now-population].fitness[0];
					pop[now].fitness[1] = opt[n-minvar][now-population].fitness[1];
					pop[now].fitness[2] = opt[n-minvar][now-population].fitness[2];
					pop[now].varn = opt[n-minvar][now-population].varn;
					pop[now].flag = 2;
					pop[now].rank = -1;
                }
                else{
					pop[now].fitness[0] = opt[m-minvar][now].fitness[0];
					pop[now].fitness[1] = opt[m-minvar][now].fitness[1];
					pop[now].fitness[2] = opt[m-minvar][now].fitness[2];
					pop[now].varn = opt[m-minvar][now].varn;
					pop[now].flag = 2;
					pop[now].rank = -1;
                }


            }
            q = 0;
            for(k=0; k<twopop; k++,q=0)
            {
                for(j=0; j<twopop; j++)
                {
                    if(pop[j].flag!=1)break;
                    /*Break if all the individuals are assigned a rank*/
                }
                if(j==twopop)break;

                rnk = rnk + 1 ;

                for(j=0; j<twopop; j++)
                {
                    if(pop[j].flag==0)pop[j].flag = 2;
                    /*Set the flag of dominated(looser) individuals to 2*/
                }
                for(i=0; i<twopop; i++)
                {
                    //printf("i=%d\n",i);
                    /*Select an individual which rank to be assigned*/
                    if(pop[i].flag!=1 && pop[i].flag!=0)
                    {

                        ptr1 = &(pop[i].fitness[0]);

                        for(j=0; j<twopop; j++)
                        {
                            //printf("j=%d\n",j);
                            /*Select the other individual which has not got a rank*/
                            if(i != j)
                            {
                                if(pop[j].flag != 1)
                                {
                                    ptr2 = &(pop[j].fitness[0]);
									//printf("ptr1 = %lf , ptr2 = %lf \n",*(ptr1),*(ptr2));
                                    /*Compare the two individuals for fitness*/
                                    //val = indcmp(ptr1,ptr2);
                                    val = indcmp(pop[i].fitness,pop[j].fitness);

									/*VAL = 2 for dominated individual which rank to be given*/
                                    /*VAL = 1 for dominating individual which rank to be given*/

                                    /*VAL = 3 for non comparable individuals*/

                                    if( val == 2)
                                    {
                                        pop[i].flag = 0;/* individual 1 is dominated */
										//printf("i = %d , j = %d , val = %d\n",i,j,val);
                                        break;
                                    }

                                    if(val == 1)
                                    {
                                        pop[j].flag = 0;/* individual 2 is dominated */
										//printf("i = %d , j = %d , val = %d\n",i,j,val);
                                    }

                                    if(val == 3)
                                    {
                                        nondom++;/* individual 1 & 2 are non dominated */
										//printf("i = %d , j = %d , val = %d\n",i,j,val);
                                        if(pop[j].flag != 0)
                                            pop[j].flag = 3;
                                    }


                                }


                            }
                        }
                        if(j == twopop)
                        {
                            pop[i].rank = rnk;
                            pop[i].flag = 1;
                            rankarr[q] = i;
                            q++;
                        }
                    }
                }
				//pop[i]
            }
			for(i=0;i<twopop;i++){
				if(pop[i].rank==1 && pop[i].varn==m){
					tpa+=1;
				}
				else if(pop[i].rank==1 && pop[i].varn==n){
					tpb+=1;
				}
			}
			s1 = (double)tpa/((double)tpa+(double)tpb);
			s2 = (double)tpb/((double)tpa+(double)tpb);
			per1 = s1 * 100.0;
			per2 = s2 * 100.0;
			printf("level 1 in %d variables : %d, level 1 in %d variables : %d \n",m,tpa,n,tpb);
			fprintf(f8,"var = %d , Percentage: %.2f%% , var = %d , Percentage: %.2f%% \n",m,per1,n,per2);
            //printf("var = %d , percent: %lf , var = %d , percent: %lf \n",m,s1,n,s2);




        }
    }
	/*
    for(i=0; i<twopop; i++)
    {
        printf("%d rank : %d source %d fitness: %lf %lf %lf \n",i,pop[i].rank,pop[i].varn,pop[i].fitness[0],pop[i].fitness[1],pop[i].fitness[2]);
    }
	*/
	/*
	for(i=0;i<twopop;i++){
		printf("rankarr %d\n",rankarr[i]);
	}
    */


    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);
    fclose(f5);
    fclose(f6);
    fclose(f7);
	fclose(f8);
	system("pause");
    return 0;
}
int indcmp(double p1f[3],double p2f[3])
{
    float fit1[nfunc],fit2[nfunc];
    int i,value,m,n;
    for(i = 0; i < nfunc ; i++)
    {
        fit1[i] = p1f[i];
        fit2[i] = p2f[i];
    }
    m = 0;
    n = 0;
    while(m < nfunc && fit1[m] <= fit2[m])
    {
        if((fit2[m] -  fit1[m]) < 1e-7) n++;
        m++;
    }
    if(m == nfunc)
    {
        if(n == nfunc) value = 3;
        else value = 1;             /*value = 1 for dominationg*/
    }
    else
    {
        m = 0;
        n = 0;
        while(m < nfunc && fit1[m] >= fit2[m])
        {
            if((fit1[m] - fit2[m]) < 1e-7) n++;
            m++;
        }
        if(m == nfunc)
        {
            if(n != nfunc)
                value = 2;                       /*value =  2 for dominated */
            else value =3;
        }
        else value = 3;                   /*value = 3 for incomparable*/
    }

    return value;
}
