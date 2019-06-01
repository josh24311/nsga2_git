/*This is the program used to evaluate the value of the function & errors
************************************************************************/
double value(double x,int w);
static double CND(double d);
void func(population *pop_ptr);

void func(population *pop_ptr)
{
    /*File ptr to the file to store the value of the g for last iteration
      g is the parameter required for a particular problem
      Every problem is not required*/

    float *realx_ptr, /*Pointer to the array of x values*/
          *binx_ptr,      /* Pointer to the binary variables */
          *fitn_ptr,      /*Pointer to the array of fitness function*/
          x[2*maxvar],     /* problem variables */
          f[maxfun],     /*array of fitness values*/
          *err_ptr,      /*Pointer to the error */
          cstr[maxcons];

    int i,j,k;
    float error, cc;

    pop_ptr->ind_ptr= &(pop_ptr->ind[0]);
    //printf("In func before for %d\n",popsize);
    /*Initializing the max rank to zero*/
    pop_ptr->maxrank = 0;
    for(i = 0; i < popsize; i++)
    {
        //printf("In func pop for i = %d\n",i);
        pop_ptr->ind_ptr = &(pop_ptr->ind[i]);
        realx_ptr = &(pop_ptr->ind_ptr->xreal[0]);
        binx_ptr = &(pop_ptr->ind_ptr->xbin[0]);

        for(j = 0; j < nvar; j++)
        {
            // Real-coded variables
            x[j] = *realx_ptr++;
        }

        for(j = 0; j < nchrom; j++)
        {
            // Binary-codced variables
            x[nvar+j] = *binx_ptr++;
        }

        fitn_ptr = &(pop_ptr->ind_ptr->fitness[0]);
        err_ptr = &(pop_ptr->ind_ptr->error);



        /*   DO NOT CHANGE ANYTHING ABOVE   */
        /*----------------------CODE YOUR OBJECTIVE FUNCTIONS HERE------------*/
        /*All functions must be of minimization type, negate maximization
              functions */
        /*============Start Coding Your Function From This Point=============*/


        // First fitness function  //f[0]:delta, f[1];gamma, f[2]:vega


        if(nchrom==3)
        {
            //3opt
            f[0] = fabs(value((double)x[0],0)+value((double)x[1],0)+value((double)x[2],0));
            f[1] = fabs(value((double)x[0],1)+value((double)x[1],1)+value((double)x[2],1));
            f[2] = fabs(value((double)x[0],2)+value((double)x[1],2)+value((double)x[2],2));
        }
        else if(nchrom==4)
        {
            //4opt
            f[0] = fabs(value((double)x[0],0)+value((double)x[1],0)+value((double)x[2],0)+value((double)x[3],0));
            f[1] = fabs(value((double)x[0],1)+value((double)x[1],1)+value((double)x[2],1)+value((double)x[3],1));
            f[2] = fabs(value((double)x[0],2)+value((double)x[1],2)+value((double)x[2],2)+value((double)x[3],2));
        }
        else if(nchrom==5)
        {
            //5opt
            f[0] = fabs(value((double)x[0],0)+value((double)x[1],0)+value((double)x[2],0)+value((double)x[3],0)+value((double)x[4],0));
            f[1] = fabs(value((double)x[0],1)+value((double)x[1],1)+value((double)x[2],1)+value((double)x[3],1)+value((double)x[4],1));
            f[2] = fabs(value((double)x[0],2)+value((double)x[1],2)+value((double)x[2],2)+value((double)x[3],2)+value((double)x[4],2));
        }
        else if(nchrom==6)
        {
            //6opt
            f[0] = fabs(value((double)x[0],0)+value((double)x[1],0)+value((double)x[2],0)+value((double)x[3],0)+value((double)x[4],0)+value((double)x[5],0));
            f[1] = fabs(value((double)x[0],1)+value((double)x[1],1)+value((double)x[2],1)+value((double)x[3],1)+value((double)x[4],1)+value((double)x[5],1));
            f[2] = fabs(value((double)x[0],2)+value((double)x[1],2)+value((double)x[2],2)+value((double)x[3],2)+value((double)x[4],2)+value((double)x[5],2));
        }
        else if(nchrom==7)
        {
            //7opt
            f[0] = fabs(value((double)x[0],0)+value((double)x[1],0)+value((double)x[2],0)+value((double)x[3],0)+value((double)x[4],0)+value((double)x[5],0)+value((double)x[6],0));
            f[1] = fabs(value((double)x[0],1)+value((double)x[1],1)+value((double)x[2],1)+value((double)x[3],1)+value((double)x[4],1)+value((double)x[5],1)+value((double)x[6],1));
            f[2] = fabs(value((double)x[0],2)+value((double)x[1],2)+value((double)x[2],2)+value((double)x[3],2)+value((double)x[4],2)+value((double)x[5],2)+value((double)x[6],2));
        }
        else if(nchrom==8)
        {
            //8opt
            f[0] = fabs(value((double)x[0],0)+value((double)x[1],0)+value((double)x[2],0)+value((double)x[3],0)+value((double)x[4],0)+value((double)x[5],0)+value((double)x[6],0)+value((double)x[7],0));
            f[1] = fabs(value((double)x[0],1)+value((double)x[1],1)+value((double)x[2],1)+value((double)x[3],1)+value((double)x[4],1)+value((double)x[5],1)+value((double)x[6],1)+value((double)x[7],1));
            f[2] = fabs(value((double)x[0],2)+value((double)x[1],2)+value((double)x[2],2)+value((double)x[3],2)+value((double)x[4],2)+value((double)x[5],2)+value((double)x[6],2)+value((double)x[7],2));
        }
        else if(nchrom==9)
        {
            //9opt
            f[0] = fabs(value((double)x[0],0)+value((double)x[1],0)+value((double)x[2],0)+value((double)x[3],0)+value((double)x[4],0)+value((double)x[5],0)+value((double)x[6],0)+value((double)x[7],0)+value((double)x[8],0));
            f[1] = fabs(value((double)x[0],1)+value((double)x[1],1)+value((double)x[2],1)+value((double)x[3],1)+value((double)x[4],1)+value((double)x[5],1)+value((double)x[6],1)+value((double)x[7],1)+value((double)x[8],1));
            f[2] = fabs(value((double)x[0],2)+value((double)x[1],2)+value((double)x[2],2)+value((double)x[3],2)+value((double)x[4],2)+value((double)x[5],2)+value((double)x[6],2)+value((double)x[7],2)+value((double)x[8],2));
        }
		else{
			printf("Variable number error\n");
		}


















        /*=========End Your Coding Upto This Point===============*/

        /******************************************************************/
        /*              Put The Constraints Here                          */
        /******************************************************************/
        // g(x) >= 0 type (normalize g(x) as in the cstr[1] below)
        /*===========Start Coding Here=============*/

        //cstr[0] = x[0]*x[0]+x[1]*x[1]-1.0-0.1*cos(16.0*atan(x[0]/x[1]));
        //cstr[1] = (-square(x[0]-0.5) - square(x[1]-0.5) + 0.5)/0.5;

        //0403 若x相差16 則 此組合不符合限制條件，令其為-1
        //int varn = 4;//變數數量可調整 改為 nchrom 變數數量
        int m = 0,u = 1,smm = 0;
        //printf("Variables num:%d\n",nchrom);
        for(m=0; m<nchrom-1; m++)
        {

            //printf("m= %d\n",m);
            //0-1,0-2,1-2(變數3個)
            for(u=m+1; u<nchrom; u++)
            {
                //printf("u= %d\n",u);
                if(fabs(x[m]-x[u])==16.0)
                {
                    smm = 1;
                    break;
                }
            }
            if(smm==1) break;
        }
        //printf("after constrain  for i = %d\n",i);
        cstr[0] = smm==1? -1 : (x[0]*x[0]+x[1]*x[1]-1.0-0.1*cos(16.0*atan(x[0]/x[1])));


        /*===========Constraints Are Coded Upto Here=============*/
        /*   DO NOT CHANGE ANYTHING BELOW  */



        for(k = 0 ; k < nfunc ; k++)
        {
            *fitn_ptr++  = f[k];
        }

        for (k = 0; k < ncons; k++)
        {
            pop_ptr->ind_ptr->constr[k] = cstr[k];
        }

        error = 0.0;
        for (k = 0; k < ncons; k++)
        {
            cc = cstr[k];
            if(cc < 0.0)
                error = error - cc;
        }

        *err_ptr = error;
    }

    /*---------------------------* RANKING *------------------------------*/

    if(ncons == 0)
        ranking(pop_ptr);
    else
        rankcon(pop_ptr);

    return;
}
//0312 added
static double CND(double d)
{
    const double       A1 = 0.31938153;
    const double       A2 = -0.356563782;
    const double       A3 = 1.781477937;
    const double       A4 = -1.821255978;
    const double       A5 = 1.330274429;
    const double RSQRT2PI = 0.39894228040143267793994605993438;

    double K = 1.0 / (1.0 + 0.2316419 * fabs(d));
    double cnd = RSQRT2PI * exp(- 0.5 * d * d) * (K * (A1 + K * (A2 + K * (A3 + K * (A4 + K * A5)))));
    if (d > 0)
        cnd = 1.0 - cnd;
    return cnd;
}
double value(double x,int w)
{
    //0527 fprintf
    /*
    FILE *c_ptr;
    c_ptr = fopen("C_approach.txt","w");

    */
    double tp1 = -1, tp2 = -1, tp3 = -1, sum = -1;
    int now = -1;
    double X = 0.0;
    double mktp = 0.0;
    int xnew = (int)x;
    int isCall = 0;//Call 為 1 ， Put 為 0
    double sigma, C, P, upper, lower,lstm,ttp;
    sum = 0.0;
    d1 = 0.0;
    d2 = 0.0;
    C = 0.0;
    P = 0.0;
    sigma = 0.3;
    lstm = 0.0;//0430-last minus : fabs(C - mktp)
    ttp = 0.0;
    X = xpg[xnew];
    mktp = mpg[xnew];
    //printf("X=: %lf, mktp=:  %lf\n",X,mktp);
    upper = 1.000000;
    lower = 0.000000;
    if(x<8||(x>=16&&x<=23))// long/short call
    {


        while(fabs(C - mktp)>1e-6)
        {

            //fprintf(c_ptr,"C=%lf,mktp=%lf,minus=%lf,SIGMA=%lf,UPPER=%lf,LOWER=%lf\n",C,mktp,fabs(C-mktp),sigma,upper,lower);

            //if(fabs(C - mktp) < 1) break;//提早結束迴圈測試
            d1 = (log(S0/X)+(r-q+pow(sigma,2)/2)*t)/(sigma*sqrt(t));
            d2 = d1 - sigma * sqrt(t);
            C = S0 * exp(-q*t)*CND(d1)-X*exp(-r*t)*CND(d2);//C 更新有問題，random seed 更換?
            //fprintf(c_ptr,"new C = %lf now SIGMA = %lf, UPPER = %lf,d1=%lf, d2=%lf\n",C,sigma,upper,d1,d2);

            if(C - mktp>0.0)
            {
                //fprintf(c_ptr,"C > mktp move LEFT,SIGMA now = %lf\n",sigma);
                upper = sigma;//problem upper 值變超大
                sigma = (sigma+lower)/2;
                //fprintf(c_ptr,"new LOWER=%lf, new SIGMA=%lf, new UPPER=%lf\n");
            }
            else
            {
                //fprintf(c_ptr,"C <= mktp move RIGHT\n");
                lower = sigma;
                sigma = (sigma+upper)/2;
                //fprintf(c_ptr,"new LOWER=%lf, new SIGMA=%lf, new UPPER=%lf\n");
            }


            if(sigma==0.0)
            {

                break;
            }
            //printf("sigma AFTER = %lf lower = %lf upper = %lf \n",sigma,lower,upper);
            //printf("sigma found(call)\n");
            /*test value*/
        }
        //fclose(c_ptr);//關檔
        d1 = (log(S0/X)+(r-q+pow(sigma,2)/2)*t)/(sigma*sqrt(t));
        if(x<8)
        {
            del = CND(d1);
            gam = (exp(-(d1*d1)/2))/sqrt(2*pi)/(S0*sigma*sqrt(t));
            veg = S0 * sqrt(t) * (exp(-(d1*d1)/2)) / sqrt(2*pi);
        }
        else
        {
            //short 變號
            del = -CND(d1);
            gam = -(exp(-(d1*d1)/2))/sqrt(2*pi)/(S0*sigma*sqrt(t));
            veg = -S0 * sqrt(t) * (exp(-(d1*d1)/2)) / sqrt(2*pi);
        }

        tp3 = del;//暫存
    }
    else //put
    {
        while((fabs(P - mktp)>1e-6))
        {

            d1 = (log(S0/X)+(r-q+pow(sigma,2)/2)*t)/(sigma*sqrt(t));
            d2 = d1 - sigma * sqrt(t);
            //P = S0 * exp(-q*t)*CND(d1)-X*exp(-r*t)*CND(d2);
            P =X*exp(-r*t)*CND(-d2)-S0*exp(-q*t)*CND(-d1);
            if(P - mktp>0)
            {
                upper = sigma;
                sigma = (sigma+lower) / 2;
            }
            else
            {
                lower = sigma;
                sigma = (sigma+upper) / 2;
            }
            /*test value*/
            //printf("sigma found(put)\n");
        }
        d1 = (log(S0/X)+(r-q+pow(sigma,2)/2)*t)/(sigma*sqrt(t));
        if(x>=8&&x<=15)
        {
            del = CND(d1)-1;
            gam = (exp(-(d1*d1)/2))/sqrt(2*pi)/(S0*sigma*sqrt(t));
            veg = S0 * sqrt(t) * (exp(-(d1*d1)/2)) / sqrt(2*pi);
        }
        else //short put
        {
            del = -(CND(d1)-1);
            gam = -(exp(-(d1*d1)/2))/sqrt(2*pi)/(S0*sigma*sqrt(t));
            veg = -S0 * sqrt(t) * (exp(-(d1*d1)/2)) / sqrt(2*pi);
        }

        tp3 = del;
    }
    if(w==0)
    {
        return del;
    }
    else if(w==1)
    {
        return gam;
    }
    else
    {
        return veg;
    }
    return 0.0;

}





