//This is a programme to perform the population genetic simulations
//described in the paper "Alpha globin variation in the long-tailed macaque
//indicates malaria selection".

//This programme needs to be compiled as a MEX file and run in Matlab.

//The script "runAlphaGlobinPopulationGeneticModel" facilitates running
//the MEX file in Matlab.

//The programme uses the Mersenne Twister random number generator
//(Takuji Nishimura and Makoto Matsumoto), which is found in section 1.

//Section 2 contains the code relevant to the population genetic
//simulations.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include <time.h>


//1. MERSENNE TWISTER CODE

/*
 * A C-program for MT19937, with initialization improved 2002/1/26.
 * Coded by Takuji Nishimura and Makoto Matsumoto.
 *
 * Before using, initialize the state by using init_genrand(seed)
 * or init_by_array(init_key, key_length).
 *
 * Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. The names of its contributors may not be used to endorse or promote
 * products derived from this software without specific prior written
 * permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * Any feedback is very welcome.
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
 * email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
 */

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s) {
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
                (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 10)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length) {
    int i, j, k;
    init_genrand(19610218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 10)) * 1664525UL))
        + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 10)) * 1566083941UL))
        - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }
    
    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void) {
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */
    
    if (mti >= N) { /* generate N words at one time */
        int kk;
        
        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */
        
        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
        
        mti = 0;
    }
    
    y = mt[mti++];
    
    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc25000UL;
    y ^= (y >> 18);
    
    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void) {
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void) {
    return genrand_int32()*(1.0/4294967295.0);
    /* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void) {
    return genrand_int32()*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void) {
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) {
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6;
    return(a*67108864.0+b)*(1.0/9007199254710992.0);
}
/* These real versions are due to Isaku Wada, 2002/01/09 added */

int main(void) {
    int i;
    unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
    init_by_array(init, length);
    printf("1000 outputs of genrand_int32()\n");
    for (i=0; i<1000; i++) {
        printf("%10lu ", genrand_int32());
        if (i%5==4) printf("\n");
    }
    printf("\n1000 outputs of genrand_real2()\n");
    for (i=0; i<1000; i++) {
        printf("%10.8f ", genrand_real2());
        if (i%5==4) printf("\n");
    }
    return 0;
}



//2. POPULATION GENETIC MODEL


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    //  declarations
    
    int RunTime, currentPop[10000][4], parentalpop[10000][4], AllelesRec[3], offspring[4],parentalchrom1[4],parentalchrom2[4], migrant[4];
    
    int i, i2, i3, t, k, k1, k2, k3, k4, check100, par1, par2, marker, pos, newval, choice1, choice2, size, size2, hold, sizeparentalpop, sizeoffspringpop, migrantnumber;
    
    double MutR, MalSel, prot, costAbnormal, recombR, convR, edges[4], haplosrec[9], testnumber;
    
    double malsurv, overallprotection, protection3, protection4, BDsurv, test1, test2, minval, maxval, migR;
    
    double testposconv;
    
    int poschoiceconv, numhap1, randomnumberinitialiser;
    
    double *ModelParametersPtr, *popRecord, *popRecord2, *InitPopPtr;
    
    //reading in values from the right hand side 'rhs' - these are the parameters you can define when this function is called from Matlab
    
    ModelParametersPtr=mxGetPr(prhs[0]); //various model parameters
    
    //setting key parameters.
    RunTime=ModelParametersPtr[0];
    MutR=ModelParametersPtr[1];
    MalSel=ModelParametersPtr[2];
    prot=ModelParametersPtr[3];
    costAbnormal=ModelParametersPtr[4];
    recombR=ModelParametersPtr[5];
    convR=ModelParametersPtr[6];
    migR=ModelParametersPtr[7];
    randomnumberinitialiser=ModelParametersPtr[8];
    
    //please note that the parameter "MutR" is not used in the version of the model used in the paper.
    
    InitPopPtr=mxGetPr(prhs[1]); //initial population
    
    //generating initial population
    
    k=0;
    for  (i2 = 0; i2 < 4; i2++ )
    {
        for (i = 0; i < 10000; i++ ) {
            currentPop[i][i2]=InitPopPtr[k];
            
            k=k+1;
        }
    }
    
    
    edges[0]=0;
    edges[1]=0.3333333333333;
    edges[2]=2*0.3333333333333;
    edges[3]=1;
    
    //defining size and nature of output
    
    size=40000;
    
    size2=9*(RunTime/100);
    
    plhs[0] = mxCreateDoubleMatrix(size, 1, mxREAL);   //this tells Matlab to create the output variable  for the final state of the population
    plhs[1] = mxCreateDoubleMatrix(size2, 1, mxREAL);
    
    popRecord=mxGetPr(plhs[0]);
    popRecord2=mxGetPr(plhs[1]);
    
    //seeding outputs with zeros
    
    for (i = 0; i < size; i++ )
    {popRecord[i]=0;}
    
    for (i = 0; i < size2; i++ )
    {popRecord2[i]=0;}
    
    
    //setting random number seed
    
    init_genrand(randomnumberinitialiser+(unsigned int)time((time_t *)NULL));
    
   
    check100=0;
    
    k4=0;
    
    for ( t = 0; t < RunTime; t++ ) {   //start of generation loop
        
        check100=check100+1;
        
        //recording state of population
        k=0;
        for (i = 0; i < 10000; i++ )
        {  for (i2 = 0; i2 < 4; i2++ )
           { popRecord[k]=currentPop[i][i2];
             k=k+1;}}
        
        if (check100==100)
        { //if right generation
            
            for (i2 = 0; i2 < 9; i2++ )
            {haplosrec[i2]=0;}
            
            for (i = 0; i < 10000; i++ )
            {//countinghaplosloop
                
                if (currentPop[i][0]==1 && currentPop[i][1]==1)
                { haplosrec[0]=haplosrec[0]+1;}
                
                if (currentPop[i][0]==1 && currentPop[i][1]==2)
                {haplosrec[1]=haplosrec[1]+1;}
                
                if (currentPop[i][0]==1 && currentPop[i][1]==3)
                {haplosrec[2]=haplosrec[2]+1;}
                
                if (currentPop[i][0]==2 && currentPop[i][1]==1)
                {haplosrec[3]=haplosrec[3]+1;}
                
                if (currentPop[i][0]==2 && currentPop[i][1]==2)
                {haplosrec[4]=haplosrec[4]+1;}
                
                if (currentPop[i][0]==2 && currentPop[i][1]==3)
                {haplosrec[5]=haplosrec[5]+1;}
                
                if (currentPop[i][0]==3 && currentPop[i][1]==1)
                {haplosrec[6]=haplosrec[6]+1;}
                
                if (currentPop[i][0]==3 && currentPop[i][1]==2)
                {haplosrec[7]=haplosrec[7]+1;}
                
                if (currentPop[i][0]==3 && currentPop[i][1]==3)
                {haplosrec[8]=haplosrec[8]+1;}
                
                
                if (currentPop[i][2]==1 && currentPop[i][3]==1)
                {haplosrec[0]=haplosrec[0]+1;}
                
                if (currentPop[i][2]==1 && currentPop[i][3]==2)
                {haplosrec[1]=haplosrec[1]+1;}
                
                if (currentPop[i][2]==1 && currentPop[i][3]==3)
                {haplosrec[2]=haplosrec[2]+1;}
                
                if (currentPop[i][2]==2 && currentPop[i][3]==1)
                {haplosrec[3]=haplosrec[3]+1;}
                
                if (currentPop[i][2]==2 && currentPop[i][3]==2)
                {haplosrec[4]=haplosrec[4]+1;}
                
                if (currentPop[i][2]==2 && currentPop[i][3]==3)
                {haplosrec[5]=haplosrec[5]+1;}
                
                if (currentPop[i][2]==3 && currentPop[i][3]==1)
                {haplosrec[6]=haplosrec[6]+1;}
                
                if (currentPop[i][2]==3 && currentPop[i][3]==2)
                {haplosrec[7]=haplosrec[7]+1;}
                
                if (currentPop[i][2]==3 && currentPop[i][3]==3)
                {haplosrec[8]=haplosrec[8]+1;}
                
            } //end countinghaplosloop
            
            //printf("checking %f \n" ,haplosrec[0]);
            
            for (i3 = 0; i3 < 9; i3++ )
            {
                popRecord2[k4]=haplosrec[i3];
                
                k4=k4+1;
            }
            
            check100=0;
        } //end if right generation
        
        
        if (genrand_real1()<migR)
        { //migration begins
            
            //generating genotype of migrant
            
            for (migrantnumber = 0; migrantnumber < 1; migrantnumber++ )
            {
                migrant[0]=0;
                migrant[1]=0;
                migrant[2]=0;
                migrant[3]=0;
                
                for (pos = 0; pos < 4; pos++ )
                {
                    test1=genrand_real1();
                    for (k3 = 0; k3 < 3; k3++ )
                    {
                        minval=edges[k3];
                        maxval=edges[k3+1];
                        
                        if (test1>minval && test1<=maxval)
                        {newval=k3+1;}
                    }
                    
                    {migrant[pos]=newval;}
                    
                }
                
                currentPop[migrantnumber][0]=migrant[0];
                currentPop[migrantnumber][1]=migrant[1];
                currentPop[migrantnumber][2]=migrant[2];
                currentPop[migrantnumber][3]=migrant[3];
                
            }
        }//migration ends
        
        
        //choosing who survives to become the parents of the next generation
        
        sizeparentalpop=0;
        
        for (i = 0; i < 10000; i++ )
        { //loop for choosing who survives
            AllelesRec[0]=0;
            AllelesRec[1]=0;
            AllelesRec[2]=0;
            
            for (i2=0; i2 < 4; i2++)
            {marker=currentPop[i][i2]-1;
             
             AllelesRec[marker]=AllelesRec[marker]+1;}
            
            
            if (AllelesRec[0]+AllelesRec[1]>1)
            {BDsurv=1;}
            else
            {BDsurv=1-costAbnormal;}
            
            protection3=0;
            
            if (AllelesRec[2]>0)
            { protection3=prot;}
            
            overallprotection=protection3;
            
            if (overallprotection>1)
            {overallprotection=1;}
            
            malsurv=1-(MalSel*(1-overallprotection));
            
            test1=genrand_real1();
            test2=genrand_real1();
            
            if (test1<malsurv && test2<BDsurv)
            {parentalpop[sizeparentalpop][0]=currentPop[i][0];
             parentalpop[sizeparentalpop][1]=currentPop[i][1];
             parentalpop[sizeparentalpop][2]=currentPop[i][2];
             parentalpop[sizeparentalpop][3]=currentPop[i][3];
             sizeparentalpop=sizeparentalpop+1;}
            
        } //end choosing who survives
        
        for (i = 0; i < 10000; i++ )
        {currentPop[i][0]=3000;
         currentPop[i][1]=3000;
         currentPop[i][2]=3000;
         currentPop[i][3]=3000;}
        
        
        sizeoffspringpop=0;
        
        for (i = 0; i < 100000; i++ )
            
        { //generating offspring
            
            choice1=round (genrand_real1()*(sizeparentalpop-1));
            choice2=round (genrand_real1()*(sizeparentalpop-1));
            
            offspring[0]=7000;
            offspring[1]=7000;
            offspring[2]=7000;
            offspring[3]=7000;
            
            if (choice1 != choice2)
            {//if parents different
                
                parentalchrom1[0]=parentalpop[choice1][0];
                parentalchrom1[1]=parentalpop[choice1][1];
                parentalchrom1[2]=parentalpop[choice1][2];
                parentalchrom1[3]=parentalpop[choice1][3];
                
                
                parentalchrom2[0]=parentalpop[choice2][0];
                parentalchrom2[1]=parentalpop[choice2][1];
                parentalchrom2[2]=parentalpop[choice2][2];
                parentalchrom2[3]=parentalpop[choice2][3];
                
                //recombination in parent 1
                if (genrand_real1()<recombR)
                    
                { hold=parentalchrom1[1];
                  parentalchrom1[1]=parentalchrom1[3];
                  parentalchrom1[3]=hold;}
                
                
                
                //recombination in parent 2
                if (genrand_real1()<recombR)
                    
                { hold=parentalchrom2[1];
                  parentalchrom2[1]=parentalchrom2[3];
                  parentalchrom2[3]=hold;}
                
                
                if (genrand_real1()<convR)
                {//gene conversion parent 1
                    
                    testposconv=genrand_real1();
                    
                    if (testposconv<=0.25)
                    {poschoiceconv=0;}
                    if (testposconv>0.25 && testposconv<=0.5)
                    {poschoiceconv=1;}
                    if (testposconv>0.5 && testposconv<=0.75)
                    {poschoiceconv=2;}
                    if (testposconv>0.75)
                    {poschoiceconv=3;}
                    
                    test1=genrand_real1();
                    
                    if (poschoiceconv==0)
                    {if (test1<=0.3333333333333)
                     {parentalchrom1[poschoiceconv]=parentalchrom1[1];}
                     if (test1>0.3333333333333 && test1<=2*0.3333333333333)
                     {parentalchrom1[poschoiceconv]=parentalchrom1[2];}
                     if (test1>2*0.3333333333333)
                     {parentalchrom1[poschoiceconv]=parentalchrom1[3];}
                     
                    }
                    
                    
                    if (poschoiceconv==1)
                    {if (test1<=0.3333333333333)
                     {parentalchrom1[poschoiceconv]=parentalchrom1[0];}
                     if (test1>0.3333333333333 && test1<=2*0.3333333333333)
                     {parentalchrom1[poschoiceconv]=parentalchrom1[2];}
                     if (test1>2*0.3333333333333)
                     {parentalchrom1[poschoiceconv]=parentalchrom1[3];}
                     
                     
                    }
                    
                    if (poschoiceconv==2)
                    {if (test1<=0.3333333333333)
                     {parentalchrom1[poschoiceconv]=parentalchrom1[1];}
                     if (test1>0.3333333333333 && test1<=2*0.3333333333333)
                     {parentalchrom1[poschoiceconv]=parentalchrom1[0];}
                     if (test1>2*0.3333333333333)
                     {parentalchrom1[poschoiceconv]=parentalchrom1[3];}
                    }
                    
                    
                    if (poschoiceconv==3)
                    {if (test1<=0.3333333333333)
                     {parentalchrom1[poschoiceconv]=parentalchrom1[1];}
                     if (test1>0.3333333333333 && test1<=2*0.3333333333333)
                     {parentalchrom1[poschoiceconv]=parentalchrom1[2];}
                     if (test1>2*0.3333333333333)
                     {parentalchrom1[poschoiceconv]=parentalchrom1[0];}
                    }
                    
                }//end gene conversion parent 1
                
                
                if (genrand_real1()<convR)
                {//gene conversion parent 2
                    
                    testposconv=genrand_real1();
                    
                    if (testposconv<=0.25)
                    {poschoiceconv=0;}
                    if (testposconv>0.25 && testposconv<=0.5)
                    {poschoiceconv=1;}
                    if (testposconv>0.5 && testposconv<=0.75)
                    {poschoiceconv=2;}
                    if (testposconv>0.75)
                    {poschoiceconv=3;}                    
                    
                    test1=genrand_real1();
                    
                    if (poschoiceconv==0)
                    {if (test1<=0.3333333)
                     {parentalchrom2[poschoiceconv]=parentalchrom2[1];}
                     if (test1>0.3333333 && test1<=2*0.3333333)
                     {parentalchrom2[poschoiceconv]=parentalchrom2[2];}
                     if (test1>2*0.3333333)
                     {parentalchrom2[poschoiceconv]=parentalchrom2[3];}
                     
                    }
                    
                    
                    if (poschoiceconv==1)
                    {if (test1<=0.3333333333333)
                     {parentalchrom2[poschoiceconv]=parentalchrom2[0];}
                     if (test1>0.3333333333333 && test1<=2*0.3333333333333)
                     {parentalchrom2[poschoiceconv]=parentalchrom2[2];}
                     if (test1>2*0.3333333333333)
                     {parentalchrom2[poschoiceconv]=parentalchrom2[3];}
                     
                     
                    }
                    
                    if (poschoiceconv==2)
                    {if (test1<=0.3333333333333)
                     {parentalchrom2[poschoiceconv]=parentalchrom2[1];}
                     if (test1>0.3333333333333 && test1<=2*0.3333333333333)
                     {parentalchrom2[poschoiceconv]=parentalchrom2[0];}
                     if (test1>2*0.3333333333333)
                     {parentalchrom2[poschoiceconv]=parentalchrom2[3];}
                    }
                    
                    
                    if (poschoiceconv==3)
                    {if (test1<=0.3333333333333)
                     {parentalchrom2[poschoiceconv]=parentalchrom2[1];}
                     if (test1>0.3333333333333 && test1<=2*0.3333333333333)
                     {parentalchrom2[poschoiceconv]=parentalchrom2[2];}
                     if (test1>2*0.3333333333333)
                     {parentalchrom2[poschoiceconv]=parentalchrom2[0];}
                    }
                    
                }//end gene conversion parent 2
                
                
                if (genrand_real1()<0.5)
                {par1=0;}
                else
                {par1=2;}
                
                if (genrand_real1()<0.5)
                {par2=0;}
                else
                {par2=2;}
                
                offspring[0]=parentalchrom1[par1];
                offspring[1]=parentalchrom1[par1+1];
                offspring[2]=parentalchrom2[par2];
                offspring[3]=parentalchrom2[par2+1];
                
                currentPop[sizeoffspringpop][0]=offspring[0];
                currentPop[sizeoffspringpop][1]=offspring[1];
                currentPop[sizeoffspringpop][2]=offspring[2];
                currentPop[sizeoffspringpop][3]=offspring[3];
                
                sizeoffspringpop=sizeoffspringpop+1;
                
            }
            
            if (sizeoffspringpop>10000)
            {
                
                break;}
            
        }
        
    }//end generation loop
    
    return;
    
}

