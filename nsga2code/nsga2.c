/* This is a Multi-Objective GA program.
**********************************************************************
*  This program is the implementation of the NSGA-2 proposed by      *
*                                                                    *
*  Prof. Kalyanmoy Deb and his students .                            *
*                                                                    *
*  copyright Kalyanmoy Deb
**********************************************************************

18.08.2003: The keepaliven.h file is modified to have normalized
            crowding distance calculation. The previous version of
            the code did not have this feature. This way, maintaining
            a good distribution of solutions in problems having quite
            a different range of objective functions were difficult.
            Hopefully, with this modification, such difficulties will
            not appear. --  K. Deb
18.08.2003: Also the dfit.h file is deleted. It was not needed any way.

The user have to give the input manualy or through a data file.

The user needs to enter objective functions in func-con.h
The code can also take care of the constraints. Enter the constraints
in the space provided in the func-con.h file.
Constraints must be of the following type:
g(x) >= 0.0
Also normalize all constraints (see the example problem in func-con.h)

If your program asks you to increase the values of some parameters in the
program come to main program and accordingly changed the values which are
defined against #define ...

The program generates few output files. These are described as
1.output.out
*           This file has the detailed record for all the variables,
*           the fitness values, constraint values, overall constraint
            violation (penalty)  and their ranks for all the members
*           of old population in the left hand side of the |**|
*           and of new population in the right hand side.

2.all_fitness.out
*         This file prints the record of all the fitness values for
*         different individual of new popultion created at all
*         generations.

3.g_rank_record.out
*        This file maintains the record of individuals in global pop-
*        -ulation at different ranks for all the generations.

4.ranks.out
*         This file prints the number of individual at different ranks
*          in old and new population and finds rank ratios

5.final_fitness.out
*                 This file has the fitness value of all feasible and
                  non-dominated individuals at the final generation

6.final_var.out
*                 This file has the all the variables of the feasible
                  and non-dominated individuals at the final generation.
                  The i-th solutions here corresponds to the i-th solution
                  in the final_fitness.out file.

7.plot.out        This file contains gnuplot-based file for plotting
                  the non-dominated feasible solutions obtained by the code.
*************************************************************************
*         This is recommended to delete or rename all the *.out files
*         obtained from the previous runs as some files are opened in
*         append mode so they give false resemblence of data if the
*         user is not careful

Compilation procedure:  gcc nsga2.c -lm
Run ./a.out with or without an input file

Input data files: Three files are included, but at one time one is needed
depending on the type of variables used:
inp-r (template file input-real)  : All variables are real-coded
inp-b (template file input-binary): All variables are binary-coded
inp-rb(template file input-rl+bin): Some variables are real and some are binary
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define square(x) ((x)*(x))
#define maxpop   500  /*Max population */
#define maxchrom 200  /*Max chromosome length*/
#define maxvar    20  /*Max no. of variables*/
#define maxfun    10  /*Max no. of functions */
#define maxcons   20  /*Max no. of Constraints*/
//0312
#define pi 3.14159265358979
typedef unsigned int bool;
#define false 0
#define true 1
double d1,d2,del,gam,veg;
double r = 0.01035;
double todate = 14.0;//距到期日(11/21) 交易日11/7
double totday = 252.0;
double t = 14.0/252.0;//11/2: 19.0， 11/3: 18.0,11/6: 15,11/7:14
double q = 0.0;
double S0 = 9824.95;//11/2:9844.74，11/5:9906.59, 11/6:9889.81,11/7:9824.95
//0401 因增加short運算 由[16] 改為 [32]，xpg:(買權8賣權8)x2，mpg:(買權結算價8賣權結算價8)x2


//double xpg[32] = {9600.0,9700.0,9800.0,9900.0,10000.0,10100.0,10200.0,10300.0,9000.0,9100.0,9200.0,9300.0,9400.0,9500.0,9600.0,9700.0,9600.0,9700.0,9800.0,9900.0,10000.0,10100.0,10200.0,10300.0,9000.0,9100.0,9200.0,9300.0,9400.0,9500.0,9600.0,9700.0};
double xpg[32] = {9700.0,9800.0,9900.0,10000.0,10100.0,10200.0,10300.0,10400.0,9200.0,9300.0,9400.0,9500.0,9600.0,9700.0,9800.0,9900.0,9700.0,9800.0,9900.0,10000.0,10100.0,10200.0,10300.0,10400.0,9200.0,9300.0,9400.0,9500.0,9600.0,9700.0,9800.0,9900.0};

//double mpg[32] = {302.0,232.0,173.0,123.0,80.0,51.0,31.0,18.0,31.5,37.5,45.0,56.0,70.0,91.0,114.0,145.0,302.0,232.0,173.0,123.0,80.0,51.0,31.0,18.0,31.5,37.5,45.0,56.0,70.0,91.0,114.0,145.0};//11/1
//double mpg[32] = {359.0,283.0,215.0,157.0,105.0,67.0,41.0,23.5,24.0,28.5,35.5,42.5,53.0,67.0,86.0,111.0,359.0,283.0,215.0,157.0,105.0,67.0,41.0,23.5,24.0,28.5,35.5,42.5,53.0,67.0,86.0,111.0};//11/2
//double mpg[32] = {289.0,220.0,158.0,108.0,69.0,41.5,24.0,13.0,30.0,36.0,44.0,54.0,68.0,89.0,113.0,144.0,289.0,220.0,158.0,108.0,69.0,41.5,24.0,13.0,30.0,36.0,44.0,54.0,68.0,89.0,113.0,144.0};//11/6
double mpg[32] = {245.0,173.0,115.0,71.0,40.0,21.0,11.0,5.4,26.5,33.5,42.0,55.0,72.0,94.0,126.0,167.0,245.0,173.0,115.0,71.0,40.0,21.0,11.0,5.4,26.5,33.5,42.0,55.0,72.0,94.0,126.0,167.0};//11/7

#define numInOnePortfolio 3 //一個交易組合所持有之不同選擇權檔數

//0312

int gener,       /*No of generations*/
    nvar,nchrom,          /*No of variables*/
    ncons,         /*No of Constraints*/
    vlen[maxvar],  /*Array to store no of bits for each variable*/
    nmut,          /* No of Mutations */
    ncross,        /*No of crossovers*/
    ans;
float seed,      /*Random Seed*/
      pcross,        /*Cross-over Probability*/
      pmut_b, pmut_r,          /*Mutation Probability*/
      lim_b[maxvar][2], lim_r[maxvar][2];/*Limits of variable in array*/
float di,        /*Distribution Index for the Cross-over*/
      dim,           /*Distribution Index for the Mutation*/
      delta_fit,     /* variables required forfitness for fitness sharing */
      min_fit,
      front_ratio;
int optype,      /*Cross-over type*/
    nfunc,         /*No of functions*/
    sharespace;    /*Sharing space (either parameter or fitness)*/

double coef[maxvar]; /*Variable used for decoding*/

static int popsize,  /*Population Size*/
       chrom;             /*Chromosome size*/

typedef struct       /*individual properties*/
{
    int genes[maxchrom], /*bianry chromosome*/
        rank,              /*Rank of the individual*/
        flag;              /*Flag for ranking*/
    float xreal[maxvar], /*list of real variables*/
          xbin[maxvar];      /*list of decoded value of the chromosome */
    float fitness[maxfun],/*Fitness values */
          constr[maxcons],     /*Constraints values*/
          cub_len,             /*crowding distance of the individual*/
          error;              /* overall constraint violation for the individual*/
} individual;       /*Structure defining individual*/


typedef struct
{
    int maxrank;            /*Maximum rank present in the population*/
    float rankrat[maxpop];  /*Rank Ratio*/
    int rankno[maxpop];     /*Individual at different ranks*/
    individual ind[maxpop], /*Different Individuals*/
               *ind_ptr;
} population ;            /*Popuation Structure*/

#include "random.h"       /*Random Number Generator*/

#include "input.h"        /*File Takes Input from user*/

#include "realinit.h"     /*Random Initialization of the populaiton*/

#include "init.h"         /*Random Initialization of the population*/

#include "decode.h"       /*File decoding the binary dtrings*/

#include "ranking.h"      /*File Creating the Pareto Fronts*/

#include "rancon.h"       /*File Creating the Pareto Fronts when
			    Constraints are specified*/

#include "func-con.h"     /*File Having the Function*/

#include "select.h"       /*File for Tournament Selection*/

#include "crossover.h"    /*Binary Cross-over*/

#include "uniformxr.h"    /*Uniform Cross-over*/

#include "realcross2.h"   /*Real Cross-over*/

#include "mut.h"          /*Binary Mutation*/

#include "realmut1.h"     /*Real Mutation*/

#include "keepaliven.h"   /*File For Elitism and Sharing Scheme*/

#include "report.h"       /*Printing the report*/

population oldpop,
           newpop,
           matepop,
           *old_pop_ptr,
           *new_pop_ptr,
           *mate_pop_ptr;
/*Defining the population Structures*/
//double mpg[32] = {289.0,220.0,158.0,108.0,69.0,41.5,24.0,13.0,30.0,36.0,44.0,54.0,68.0,89.0,113.0,144.0,289.0,220.0,158.0,108.0,69.0,41.5,24.0,13.0,30.0,36.0,44.0,54.0,68.0,89.0,113.0,144.0};//11/6
main()
{

    /*Some Local variables to this Problem (Counters And some other pointers*/
    printf("main start\n");
    int i,j,l,f,maxrank1;
    float *ptr,tot;
    FILE
    *rep_ptr,
    *gen_ptr,
    *rep2_ptr,
    *end_ptr,
    *g_var,
    *lastit,
	//0531 new 
	*fts_ptr,
	*opt3,*opt4,*opt5,*opt6,*opt7,*opt8,*opt9;
    /*File Pointers*/

    rep_ptr = fopen("output.out","w");
    gen_ptr =fopen("all_fitness.out","w");
    rep2_ptr = fopen("ranks.out","w");
    end_ptr = fopen("final_fitness.out","w");
    g_var = fopen("final_var.out","w");
    lastit = fopen("plot.out","w");
	//fts_ptr = fopen("fitness_variable.txt","w");
	
	//0601 new
	opt3 = fopen("opt3_o.txt","a");
	opt4 = fopen("opt4_o.txt","a");
	opt5 = fopen("opt5_o.txt","a");
	opt6 = fopen("opt6_o.txt","a");
	opt7 = fopen("opt7_o.txt","a");
	opt8 = fopen("opt8_o.txt","a");
	opt9 = fopen("opt9_o.txt","a");
    /*Opening the files*/

    old_pop_ptr = &(oldpop);

    nmut = 0;
    ncross = 0;

    /*Get the input from the file input.h*/
    input(rep_ptr);

    fprintf(rep_ptr,"Results in a file\n");
    fprintf(end_ptr,"# Last generation population (Feasible and non-dominated)\n");
    fprintf(end_ptr,"# Fitness_vector (first %d)  Constraint_violation (next %d)  Overall_penalty\n",nfunc,ncons);
    fprintf(g_var,"#Feasible Variable_vectors for non-dominated solutions at last generation\n");
    fprintf(g_var,"# Real (first %d)  Binary (next %d)\n",nvar,nchrom);
    fprintf(lastit,"# Feasible and Non-dominated Objective Vector\n");

    /*Initialize the random no generator*/
    warmup_random(seed);
    //printf("warmup_random(seed) \n");
    /*Binary Initializaton*/
    if (nchrom > 0)
        init(old_pop_ptr);
    if (nvar > 0)
        realinit(old_pop_ptr);

    old_pop_ptr = &(oldpop);

    // decode binary strings
    decode(old_pop_ptr);

    old_pop_ptr = &(oldpop);
    new_pop_ptr = &(newpop);

    for(j = 0; j < popsize; j++)
    {
        /*Initializing the Rank array having different individuals
        at a particular  rank to zero*/
        old_pop_ptr->rankno[j] = 0;
        new_pop_ptr->rankno[j] = 0;
    }

    old_pop_ptr = &(oldpop);
    printf("decode binary strings before func \n");
    func(old_pop_ptr);
    printf("decode binary strings after func \n");
    /*Function Calculaiton*/

    fprintf(rep_ptr,"----------------------------------------------------\n");
    fprintf(rep_ptr,"Statistics at Generation 0 ->\n");
    fprintf(rep_ptr,"--------------------------------------------------\n");

    /********************************************************************/
    /*----------------------GENERATION STARTS HERE----------------------*/
    for (i = 0; i < gener; i++)
    {
        printf("Generation = %d\n",i+1);
        old_pop_ptr = &(oldpop);
        mate_pop_ptr = &(matepop);
        fprintf(rep_ptr,"Population at generation no. -->%d\n",i+1);
        fprintf(gen_ptr,"#Generation No. -->%d\n",i+1);
        fprintf(gen_ptr,"#Variable_vector  Fitness_vector Constraint_violation Overall_penalty\n");

        /*--------SELECT----------------*/
        nselect(old_pop_ptr,mate_pop_ptr );

        new_pop_ptr = &(newpop);
        mate_pop_ptr = &(matepop);

        /*CROSSOVER----------------------------*/
        if (nchrom > 0)
        {

            if(optype == 1)
            {
                crossover(new_pop_ptr,mate_pop_ptr );
                /*Binary Cross-over*/
            }

            if(optype == 2)
            {
                unicross(new_pop_ptr,mate_pop_ptr );
                /*Binary Uniform Cross-over*/
            }
        }
        if (nvar > 0)
            realcross(new_pop_ptr,mate_pop_ptr );
        /*Real Cross-over*/


        /*------MUTATION-------------------*/
        new_pop_ptr = &(newpop);

        if (nchrom > 0)
            mutate(new_pop_ptr );
        /*Binary Mutation */

        if (nvar > 0)
            real_mutate(new_pop_ptr );
        /*Real Mutation*/

        new_pop_ptr = &(newpop);

        /*-------DECODING----------*/
        if(nchrom > 0)
            decode(new_pop_ptr );
        /*Decoding for binary strings*/

        /*----------FUNCTION EVALUATION-----------*/
        new_pop_ptr = &(newpop);
        func(new_pop_ptr );

        /*-------------------SELECTION KEEPING FRONTS ALIVE--------------*/
        old_pop_ptr = &(oldpop);
        new_pop_ptr = &(newpop);
        mate_pop_ptr = &(matepop);

        /*Elitism And Sharing Implemented*/
        keepalive(old_pop_ptr,new_pop_ptr,mate_pop_ptr,i+1);

        mate_pop_ptr = &(matepop);
        if(nchrom > 0)
            decode(mate_pop_ptr );

        mate_pop_ptr = &(matepop);
        /*------------------REPORT PRINTING--------------------------------*/
        report(i,old_pop_ptr,mate_pop_ptr,rep_ptr,gen_ptr, lastit );

        /*==================================================================*/

        /*----------------Rank Ratio Calculation------------------------*/
        new_pop_ptr = &(matepop);
        old_pop_ptr = &(oldpop);

        /*Finding the greater maxrank among the two populations*/

        if(old_pop_ptr->maxrank > new_pop_ptr->maxrank)
            maxrank1 = old_pop_ptr->maxrank;
        else
            maxrank1 = new_pop_ptr->maxrank;

        fprintf(rep2_ptr,"--------RANK AT GENERATION %d--------------\n",i+1);
        fprintf(rep2_ptr,"Rank old ranks   new ranks     rankratio\n");

        for(j = 0; j < maxrank1 ; j++)
        {
            /*Sum of the no of individuals at any rank in old population
              and the new populaion*/

            tot = (old_pop_ptr->rankno[j])+ (new_pop_ptr->rankno[j]);

            /*Finding the rank ratio for new population at this rank*/

            new_pop_ptr->rankrat[j] = (new_pop_ptr->rankno[j])/tot;

            /*Printing this rank ratio to a file called ranks.dat*/

            fprintf(rep2_ptr," %d\t  %d\t\t %d\t %f\n",j+1,old_pop_ptr->rankno[j],new_pop_ptr->rankno[j],new_pop_ptr->rankrat[j]);

        }

        fprintf(rep2_ptr,"-----------------Rank Ratio-------------------\n");
        /*==================================================================*/

        /*=======Copying the new population to old population======*/

        old_pop_ptr = &(oldpop);
        new_pop_ptr = &(matepop);

        for(j = 0; j < popsize; j++)
        {
            old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[j]);
            new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[j]);
            if(nchrom > 0)
            {
                /*For Binary GA copying of the chromosome*/

                for(l = 0; l < chrom; l++)
                    old_pop_ptr->ind_ptr->genes[l]=new_pop_ptr->ind_ptr->genes[l];

                for(l = 0; l < nchrom; l++)
                    old_pop_ptr->ind_ptr->xbin[l] = new_pop_ptr->ind_ptr->xbin[l];
            }
            if(nvar > 0)
            {
                /*For Real Coded GA copying of the chromosomes*/
                for(l = 0; l < nvar; l++)
                    old_pop_ptr->ind_ptr->xreal[l] = new_pop_ptr->ind_ptr->xreal[l];
            }

            /*Copying the fitness vector */
            for(l = 0 ; l < nfunc ; l++)
                old_pop_ptr->ind_ptr->fitness[l] = new_pop_ptr->ind_ptr->fitness[l];

            /*Copying the dummy fitness*/
            old_pop_ptr->ind_ptr->cub_len = new_pop_ptr->ind_ptr->cub_len;

            /*Copying the rank of the individuals*/
            old_pop_ptr->ind_ptr->rank = new_pop_ptr->ind_ptr->rank;

            /*Copying the error and constraints of the individual*/

            old_pop_ptr->ind_ptr->error = new_pop_ptr->ind_ptr->error;
            for(l = 0; l < ncons; l++)
            {
                old_pop_ptr->ind_ptr->constr[l] = new_pop_ptr->ind_ptr->constr[l];
            }

            /*Copying the flag of the individuals*/
            old_pop_ptr->ind_ptr->flag = new_pop_ptr->ind_ptr->flag;
        }   // end of j

        maxrank1 = new_pop_ptr->maxrank ;

        /*Copying the array having the record of the individual
        at different ranks */
        for(l = 0; l < popsize; l++)
        {
            old_pop_ptr->rankno[l] = new_pop_ptr->rankno[l];
        }

        /*Copying the maxrank */
        old_pop_ptr->maxrank = new_pop_ptr->maxrank;

        /*Printing the fitness record for last generation in a file last*/
        if(i == gener-1)
        {
            // for the last generation
            old_pop_ptr = &(matepop);
            for(f = 0; f < popsize ; f++) // for printing
            {
                old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[f]);

                if ((old_pop_ptr->ind_ptr->error <= 0.0) && (old_pop_ptr->ind_ptr->rank == 1))  // for all feasible solutions and non-dominated solutions
                {
                    for(l = 0; l < nfunc; l++)
					{
						fprintf(end_ptr,"%f\t",old_pop_ptr->ind_ptr->fitness[l]);
						//0601
						if(nchrom==3){
							fprintf(opt3,"%f\t",old_pop_ptr->ind_ptr->fitness[l]);
						}
						else if(nchrom==4){
							fprintf(opt4,"%f\t",old_pop_ptr->ind_ptr->fitness[l]);
						}
						else if(nchrom==5){
							fprintf(opt5,"%f\t",old_pop_ptr->ind_ptr->fitness[l]);
						}
						else if(nchrom==6){
							fprintf(opt6,"%f\t",old_pop_ptr->ind_ptr->fitness[l]);
						}
						else if(nchrom==7){
							fprintf(opt7,"%f\t",old_pop_ptr->ind_ptr->fitness[l]);
						}
						else if(nchrom==8){
							fprintf(opt8,"%f\t",old_pop_ptr->ind_ptr->fitness[l]);
						}
						else if(nchrom==9){
							fprintf(opt9,"%f\t",old_pop_ptr->ind_ptr->fitness[l]);
						}
						
					}
                        
                    for(l = 0; l < ncons; l++)
                    {
                        fprintf(end_ptr,"%f\t",old_pop_ptr->ind_ptr->constr[l]);
                    }
                    if (ncons > 0)
                        fprintf(end_ptr,"%f\t",old_pop_ptr->ind_ptr->error);
                    fprintf(end_ptr,"\n");

                    if (nvar > 0)
                    {
                        for(l = 0; l < nvar ; l++)
                        {
                            fprintf(g_var,"%f\t",old_pop_ptr->ind_ptr->xreal[l]);
                        }
                        fprintf(g_var,"  ");
                    }

                    if(nchrom > 0)
                    {
                        for(l = 0; l < nchrom; l++)
                        {
                            fprintf(g_var,"%f\t",old_pop_ptr->ind_ptr->xbin[l]);
                        }
                    }
                    fprintf(g_var,"\n");
					//06061
					if(nchrom==3){
							fprintf(opt3,"%d\n",nchrom);
						}
						else if(nchrom==4){
							fprintf(opt4,"%d\n",nchrom);
						}
						else if(nchrom==5){
							fprintf(opt5,"%d\n",nchrom);
						}
						else if(nchrom==6){
							fprintf(opt6,"%d\n",nchrom);
						}
						else if(nchrom==7){
							fprintf(opt7,"%d\n",nchrom);
						}
						else if(nchrom==8){
							fprintf(opt8,"%d\n",nchrom);
						}
						else if(nchrom==9){
							fprintf(opt9,"%d\n",nchrom);
						}
					
                }  // feasibility check
            } // end of f (printing)

        } // for the last generation
    }  // end of i

    /*                   Generation Loop Ends                                */
    /************************************************************************/

    fprintf(rep_ptr,"NO. OF CROSSOVER = %d\n",ncross);
    fprintf(rep_ptr,"NO. OF MUTATION = %d\n",nmut);
    fprintf(rep_ptr,"------------------------------------------------------------\n");
    fprintf(rep_ptr,"---------------------------------Thanks---------------------\n");
    fprintf(rep_ptr,"-------------------------------------------------------------\n");
    printf("NOW YOU CAN LOOK IN THE FILE OUTPUT2.DAT\n");

    /*Closing the files*/
    fclose(rep_ptr);
    fclose(gen_ptr);
    fclose(rep2_ptr);
    fclose(end_ptr);
    fclose(g_var);
    fclose(lastit);
	
	fclose(opt3);
	fclose(opt4);
	fclose(opt5);
	fclose(opt6);
	fclose(opt7);
	fclose(opt8);
	fclose(opt8);
}





