#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define mu0 (M_PI*.4e-6)
#define DELTA   .01

#define MIN_RES 0.01
#define MAX_RES 20000.
#define RES_INI 10
#define AVERAGE 4
#define MAX_ITER 5
#define DA 0
#define STOP_VAL 1.0

typedef struct {
    double coord;
    double coord_y;
    double hor_dist;
    double ver_dist;
    double relief;
    double alt;
    double prim;
    double *w;
} geometry;


double *upper;
double *lower;

double bessj0( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of first kind and order  */
/*          0 at input x                                      */
/*------------------------------------------------------------*/
{
   double ax,z;
   double xx,y,ans,ans1,ans2;

   if ((ax=fabs(x)) < 8.0) {
      y=x*x;
      ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
         +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
      ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
         +y*(59272.64853+y*(267.8532712+y*1.0))));
      ans=ans1/ans2;
   } else {
      z=8.0/ax;
      y=z*z;
      xx=ax-0.785398164;
      ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
         +y*(-0.2073370639e-5+y*0.2093887211e-6)));
      ans2 = -0.1562499995e-1+y*(0.1430488765e-3
         +y*(-0.6911147651e-5+y*(0.7621095161e-6
         -y*0.934935152e-7)));
      ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
   }
   return ans;
}

double complex PartSum(double n0,double hh,double complex n1,double r,double complex Imp)
{
    double complex s;
    if(fabs(r)>.001)
        s = bessj0(n0*r);
    else {
        s  = 1;
    }
    double complex A = exp(-n0*hh)*(n1-n0*Imp)*n0*n0*.25/(n1+n0*Imp)/M_PI;
    s = A*s;
    return s;
}

double complex Impedance(int n, double n0, double om, double complex *rho, double *dep) {
    int i, m;
    double complex ni,nim1;
    double complex Imp;
    double dpth;

    Imp = 1;
    nim1 = csqrt(n0*n0 - I*om*mu0/rho[n-1]);
    dpth = 0;
    m = n - 1;
    while(m>0 && fabs(rho[m] - rho[m-1])<0.01) m--;

    for(i=m;i>0;i--){ 
		     dpth+=dep[i-1];
		     while(i>1 && fabs(rho[i-2] - rho[i-1])<0.01){
			   dpth+=dep[i-2];
			   i--;
		     }
                     ni = nim1;
                     nim1 = csqrt(n0*n0 - I*om*mu0/rho[i-1]);
                     Imp = ctanh(nim1*dpth+catanh(nim1/ni*Imp));
		     dpth = 0;
    }
    return Imp;
}

double complex integral(int n,double hh,double r,double complex *rho,double *dep,double f)
{
    double complex PS;         // f  -  Hz!!!
    double complex intl = 0;
    double dn0;
    double n0=0;//in0;
    double complex n1,c;
    double sigma = 1./rho[0];
    double complex Imp;
    double om = f*2*M_PI;
    c = I*(-om*sigma*mu0);        //  4pi * 1.e-7 * 2pi (Hz -> 1/sec)

    #define VAL .001
    for(n0=VAL,dn0=VAL;n0<1;n0+=dn0) {
        n1 = csqrt(n0*n0+c);
        Imp = Impedance(n,n0,om,rho,dep);
        PS  = PartSum(n0,hh,n1,r,Imp);
	if(isnan(creal(PS))){
               if(isnan(cimag(PS))){
                    PS = 0 + 0*I;
	       } else{
                    PS = 0 + cimag(PS)*I;
	       }
	} else {
               if(isnan(cimag(PS)))
		    PS = creal(PS) + 0*I;
 
        }
        intl += dn0*PS;
	}
    return intl;
}


double complex ImHz(int n, double r,double z,double f,double complex *rho, double *dep) {
    return integral(n,z,r,rho,dep,f);
}


void fdfun(geometry geo, int nlay, int bfr,double *x, double *cole, int *pol_inds, int pol_num, double *y, double *freqs) {
    int i, polc, j;
    int freq_num = bfr;
    double complex refl;
    double complex rho[nlay];
    double da = 0;
    double m, tau, c, f;
    if(nlay>1) da = x[nlay];
    memset(rho, 0, sizeof(rho));
    for(i=0;i<freq_num;i++) {
	f = freqs[i];
	if(pol_num > 0){
	    polc = 0;
	    for(j = 0; j < nlay; j++){
                 if(j == pol_inds[polc]){
	             m = 1 - cole[3*polc]/x[j];
		     tau = cole[3*polc + 1];
		     c = cole[3*polc + 2];
		     rho[j] = x[j]*(1 - m*(1 - 1/(1 + pow(2*M_PI*f*tau, c)*cexp(5*I*M_PI*c/2))));
                     polc++;
                 } else {
                     rho[j] = x[j];
	         }
            }
	} else {
             for(j = 0; j < nlay; j++) rho[j] = x[j];
	}
        refl = ImHz(nlay,geo.hor_dist,2*(geo.alt + da)+geo.ver_dist,freqs[i],rho,&(x[nlay + 1]))*I/geo.prim ;
        y[2*i] = fabs(creal(refl));
        y[2*i+1] = fabs(cimag(refl));
        if(i>bfr) break;
    }
    
        for(i=0;i<bfr;i++) {
            if(i<freq_num-1){
                y[2*i+1] = y[2*(i+1)+1] - y[2*i+1];
	    }
            else{
		if(bfr > 1)
                       y[2*i+1] = y[2*(i-1)+1];
	    }
        }
     
}


// FK 1 step
void proc(int n, double *x, double *S, double z, double *h, double sg2) {
// Square Root Matrix method for LSM(Kalman)
    int i,j,k;
    double f[n + 1],e[n + 1],d[2],bk,ck,dz;
    d[0] = sg2;
    for(i=0;i<n;i++) {
        f[i] = 0;
        for(j=0;j<n;j++)
            f[i]+=S[j*n+i]*h[j]; // S^T !!!
    }
    memset(e,0,sizeof(e));
    for(k=0;k<n+1;k++) {
        d[1] = d[0] + f[k]*f[k];
        bk = sqrt(d[0]/d[1]);
        ck = f[k]/sqrt(d[0]*d[1]);
        for(i=0;i<n;i++) {
            double tmp = S[i*n+k]*f[k];
            S[i*n+k] = bk*S[i*n+k]-ck*e[i];
            e[i] += tmp;
        }
        d[0] = d[1];
    }
    dz = z;
    for(i=0;i<n;i++) dz -= h[i]*x[i];
    dz/=d[0];
    for(i=0;i<n;i++) x[i] += e[i]*dz;
}

// combined free and fixed layers inversion
// nlay shoud be the total number of layers
int flinversion(geometry geo,
                int bfr,// now it is the # of used frqss`
                int nlay,
                double *x_ini,
		double *cole_ini,
		int *pol_inds,
		int pol_num,
                double *dpth,
                double *y_ini,
                double *y_mes,
                double *residual,
                int *up,
                double *S,
		double *freqs,
		double *upper,
		double *lower) {// both x_ and y_ are in log-axes

    int lay_num = nlay;
    int freq_num = bfr;
    double y1[2*freq_num];
    double dx[2*lay_num + 1 + 3*pol_num],x0[2*lay_num + 1],x1[2*lay_num + 1],xini[2*lay_num + 1];
    double dcole[3*pol_num], cole0[3*pol_num], cole1[3*pol_num], coleini[pol_num];
    double join_vec[lay_num + 1 + 3*pol_num];
    int charge = 1;


    double Jacob[2*freq_num][lay_num + 1 + 3*pol_num];
    double res = 0;
    int i,j;

    memset(dx,0,sizeof(dx));
    memset(Jacob, 0,sizeof(Jacob));
    for(i=0;i<nlay + 1;i++) {
        xini[i] = x_ini[i];
        if(i<nlay-1)
            xini[i+nlay + 1] = x0[i+nlay + 1] = x1[i+nlay + 1] = dpth[i];
    }
    for(i = 0; i < 3*pol_num; i++) coleini[i] = cole_ini[i];


    // first forward calculation for the model
    if(*residual<0) {
        fdfun(geo,nlay, bfr, xini, coleini, pol_inds, pol_num, y_ini, freqs);
        *residual = 0;
        for(j=0;j<2*bfr;j++) {
            double val = (y_ini[j]-y_mes[j])*(y_ini[j]-y_mes[j])*geo.w[j]*geo.w[j]/(2*bfr);
            *residual += val;//(val>*residual)?val:*residual;
        }
    }

    // a good enough first approach
    if(*residual<0.01)
        return 0;

    // Jacobian matrix calculation
    for(i=0;i<nlay + 1 + 3*pol_num;i++) {
        for(j=0;j<nlay + 1 + 3*pol_num;j++) {
            double k = 10;
	    if(j < nlay) x1[j] = x_ini[j] * ( (i!=j)? 1. : (1.-k*DELTA) );

	    if(j == nlay) x1[j] = x_ini[j] - ( (i!=j)? 0. : 0.1);

	    if(j > nlay) cole1[j-nlay-1] = cole_ini[j] * ( (i!=j)? 1. : (1. - k*DELTA) );
        }
        fdfun(geo,nlay, bfr,x1, cole1, pol_inds, pol_num, y1, freqs);
	double div = 1;

	if (i < nlay){
            div = log(x1[i]/x_ini[i]);
	}
	if (i == nlay){
            div = 0.1;
        }
        if (i > nlay){
            div = log(cole1[i-nlay-1]/cole_ini[i-nlay-1]);
	}

        for(j=0;j<2*bfr;j++) {
            Jacob[j][i] = log(fabs(y1[j]/y_ini[j]));
	    if(isnan(Jacob[j][i])) Jacob[j][i] = 0;
            Jacob[j][i] = Jacob[j][i]/div;
        }
    }
//We need to make joint variable of x + cole for proc update
    // Step (by step) a-la FK
    for(j=0;j<2*bfr;j++) {
        proc(nlay+1 + 3*pol_num,dx,S,log(fabs(y_ini[j]/y_mes[j])),Jacob[j],1./(geo.w[j]*y_mes[j]*geo.w[j]*y_mes[j]));
    }

    for(i=0;i<nlay + 1 + 3*pol_num;i++) {
	if(i < nlay)
           x0[i] = xini[i]*exp(-dx[i]);
	if(i==nlay)
           x0[i] = xini[i] + dx[i];
	if(i > nlay)
	   cole0[i - nlay-1] = coleini[i - nlay-1]*exp(-dx[i]);
	if(i <= nlay){
           if(isnan(x0[i])) x0[i] = xini[i];
           if(x0[i]>upper[i]) x0[i] = upper[i];
           if(x0[i]<lower[i]) x0[i] = lower[i];
	}
	if(i > nlay){
           if(isnan(cole0[i - nlay-1])) cole0[i - nlay-1] = coleini[i - nlay-1];
	   if(cole0[i - nlay-1] > upper[i]) cole0[i - nlay-1] = upper[i];
	   if(cole0[i - nlay-1] < lower[i]) cole0[i - nlay-1] = lower[i];
	}
    }

    fdfun(geo,nlay,bfr,x0, cole0, pol_inds, pol_num, y_ini, freqs);
    for(j=0;j<2*bfr;j++){
        double val = (y_ini[j]-y_mes[j])*(y_ini[j]-y_mes[j])*geo.w[j]*geo.w[j]/(2*bfr);
        res += val;//(val>res)?val:res;
    }

    int cntr = 0;
    while(res>*residual*1.01) {
        if(cntr++ > 3) break;
        for(i=0;i<nlay + 1 + 3*pol_num;i++) {
        	if(i < nlay)
                   x0[i] = xini[i]*exp(-dx[i]);
	        if(i==nlay)
                   x0[i] = xini[i] + dx[i];
	       if(i > nlay)
	           cole0[i - nlay-1] = coleini[i - nlay-1]*exp(-dx[i]);
	       if(i <= nlay){
                   if(isnan(x0[i])) x0[i] = xini[i];
                   if(x0[i]>upper[i]) x0[i] = upper[i];
                   if(x0[i]<lower[i]) x0[i] = lower[i];
	       }
	       if(i > nlay){
                   if(isnan(cole0[i - nlay-1])) cole0[i - nlay-1] = coleini[i - nlay-1];
	           if(cole0[i - nlay-1] > upper[i]) cole0[i - nlay-1] = upper[i];
	           if(cole0[i - nlay - 1] < lower[i]) cole0[i - nlay - 1] = lower[i];
	       }
        }


        res = 0;
        fdfun(geo,nlay, bfr,x0, cole0, pol_inds, pol_num, y_ini, freqs);
        for(j=0;j<2*bfr;j++){
            double val = (y_ini[j]-y_mes[j])*(y_ini[j]-y_mes[j])*geo.w[j]*geo.w[j]/(2*bfr);
            res += val;//(val>res)?val:res;
        }
    }

    if(res>*residual) up[0]++;
    else up[0] = 0;
    *residual = res;
    for(i=0;i<nlay + 1 + 3*pol_num;i++){
	if(i <= nlay){
           x_ini[i] = x0[i];
	} else{
           cole_ini[i - nlay - 1] = cole0[i - nlay - 1];
        }
    }
    
    return 1;
}

double primField(double hd,
                 double vd)
{

    double E[10], R[4], RR[10], Hp[4], k, MR, Ampl, M;
    R[1]=hd; R[2]=0; R[3]=vd; //towed cable = 39.11
    //R[1]=33; R[2]=0; R[3]=21; //towed cable = 39.11
    //M[1]=0; M[2]=0; M[3]=1.e9;
    M = 1;
    E[1]=1; E[2]=0; E[3]=0;
    E[4]=0; E[5]=1; E[6]=0;
    E[7]=0; E[8]=0; E[9]=1;
    RR[1]=R[1]*R[1]; RR[2]=R[1]*R[2]; RR[3]=R[1]*R[3];
    RR[4]=R[2]*R[1]; RR[5]=R[2]*R[2]; RR[6]=R[2]*R[3];
    RR[7]=R[3]*R[1]; RR[8]=R[3]*R[2]; RR[9]=R[3]*R[3];
    MR=R[1]*R[1]+R[2]*R[2]+R[3]*R[3];
    k = M / MR / (sqrt(MR)) / 4 / M_PI;
    Hp[1] = (3*RR[3] / MR-E[3]) * k;
    Hp[2] = (3*RR[6] / MR-E[6]) * k;
    Hp[3] = (3*RR[9] / MR-E[9]) * k;

    Ampl = sqrt(Hp[1] * Hp[1] + Hp[2] * Hp[2] + Hp[3] * Hp[3]);
    return Ampl;
}

int main(int argc, char **argv)
{
    

    printf("inversion started\n");
    if(argc != 3) {
        printf("inversion file_in file_out\n");
        return 0;
    }
    char *full_name[1024];
    strcpy(full_name, argv[0]); 
    char *cur_name = strrchr(full_name, '/');
    char ququ[1024] = "configs"; 
    //strcat(cur_name, ".conf");
    strcat(ququ, cur_name);
    strcat(ququ, ".conf");
   printf("%s\n", ququ);
    FILE *conf = fopen(ququ, "rt");
    printf("conf file opened");
    char buf[3000], tmp[256];

    int i, freq_num, pos, is_digit, j, read_num;
    size_t n = 0;
    double *values, *args, *freqs;

    //Read freq num
    fgets(buf, 3000, conf);
    fgets(buf, 3000, conf);
    freq_num = atoi(buf);
    values = (double *)malloc((2*freq_num + 100)*sizeof(double));
    args = (double *)malloc((12 + 3*freq_num)*sizeof(double));
    freqs = (double *)malloc(freq_num*sizeof(double));
    args[0] = freq_num;

    memset(buf, 0, sizeof(buf));

    //Read lay_num, first_thick, step
    fgets(buf, 2000, conf);
    fgets(buf, 2000, conf);
    is_digit = 0;
    j = 0;
    read_num = 0;
    for(i = 0; i < 2000; i++){
        if((buf[i]>='0' && buf[i]<='9')||buf[i]=='.'||buf[i]==','){
               tmp[j] = buf[i];
	       j++;
	       is_digit = 1;
	} else{
               if(is_digit==1){
		  values[read_num] = atof(tmp);
		  read_num++;
                  is_digit = 0;
		  for(j = 0; j < 256; j++)
	          tmp[j] = (char) 0;
	       }
	       j = 0;
        }
    }
    
    args[1] = values[0];
    args[2] = values[1];
    args[3] = values[2];

    //Read freqs
    fgets(buf, 2000, conf);
    fgets(buf, 2000, conf);

    is_digit = 0;
    j = 0;
    read_num = 0;
    for(i = 0; i < 2000; i++){
        if((buf[i]>='0' && buf[i]<='9')||buf[i]=='.'||buf[i]==','){
               tmp[j] = buf[i];
	       j++;
	       is_digit = 1;
	} else{
               if(is_digit==1){
		  values[read_num] = atof(tmp);
		  read_num++;
                  is_digit = 0;
		  for(j = 0; j < 256; j++)
	          tmp[j] = (char) 0;
	       }
	       j = 0;
        }
    }

    for(i = 0; i < freq_num; i++){
	   args[4+i] = values[i];
           freqs[i] = values[i];
    }


//Read R_i_i
    fgets(buf, 2000, conf);
    fgets(buf, 2000, conf);

 
    is_digit = 0;
    j = 0;
    read_num = 0;
    for(i = 0; i < 2000; i++){
        if((buf[i]>='0' && buf[i]<='9')||buf[i]=='.'||buf[i]==','){
               tmp[j] = buf[i];
	       j++;
	       is_digit = 1;
	} else{
               if(is_digit==1){
		  values[read_num] = atof(tmp);
		  read_num++;
                  is_digit = 0;
		  for(j = 0; j < 256; j++)
	          tmp[j] = (char) 0;
	       }
	       j = 0;
        }
    }



    for(i = 0; i < freq_num; i++){
	   args[4 + freq_num + 2*i] = values[i];
           args[4 + freq_num + 2*i + 1] = values[i];
    }

//Read ERR_INI, COR_INI
    fgets(buf, 2000, conf);
    fgets(buf, 2000, conf);


    is_digit = 0;
    j = 0;
    read_num = 0;
    for(i = 0; i < 2000; i++){
        if((buf[i]>='0' && buf[i]<='9')||buf[i]=='.'||buf[i]==','){
               tmp[j] = buf[i];
	       j++;
	       is_digit = 1;
	} else{
               if(is_digit==1){
		  values[read_num] = atof(tmp);
		  read_num++;
                  is_digit = 0;
		  for(j = 0; j < 256; j++)
	          tmp[j] = (char) 0;
	       }
	       j = 0;
        }
    }



    double ERR_INI = args[4 + 3*freq_num] = values[0];
    double COR_INI = args[5 + 3*freq_num] = values[1];
//Read position of measurements, hor_dist, ver_dist, alt
    fgets(buf, 2000, conf);
    fgets(buf, 2000, conf);


    is_digit = 0;
    j = 0;
    read_num = 0;
    for(i = 0; i < 2000; i++){
        if((buf[i]>='0' && buf[i]<='9')||buf[i]=='.'||buf[i]==','){
               tmp[j] = buf[i];
	       j++;
	       is_digit = 1;
	} else{
               if(is_digit==1){
		  values[read_num] = atof(tmp);
		  read_num++;
                  is_digit = 0;
		  for(j = 0; j < 256; j++)
	          tmp[j] = (char) 0;
	       }
	       j = 0;
        }
    }

    args[6 + 3*freq_num] = values[0];
    args[7 + 3*freq_num] = values[1];
    args[8 + 3*freq_num] = values[2];
    args[9 + 3*freq_num] = values[3];
    args[10 + 3*freq_num] = values[4];
    
//Read number of pol layers and their positions
    fgets(buf, 2000, conf);
    fgets(buf, 2000, conf);


    is_digit = 0;
    j = 0;
    read_num = 0;
    for(i = 0; i < 2000; i++){
        if((buf[i]>='0' && buf[i]<='9')||buf[i]=='.'||buf[i]==','){
               tmp[j] = buf[i];
	       j++;
	       is_digit = 1;
	} else{
               if(is_digit==1){
		  values[read_num] = atof(tmp);
		  read_num++;
                  is_digit = 0;
		  for(j = 0; j < 256; j++)
	          tmp[j] = (char) 0;
	       }
	       j = 0;
        }
    }

    int pol_num = values[0];
    int pol_num_loc;
    if(pol_num > 0){
        pol_num_loc = pol_num;
    } else{
        pol_num_loc = 1;
    }

    int pol_inds[pol_num_loc];

    if(pol_num==0) pol_inds[0] = -1;

    for(i = 0; i < pol_num; i++) pol_inds[i] = values[1+i];

    printf("conf params obtained");
    geometry geo;
    char *data;
    //char time[23];
    //for double spaces
    char time[13];
    double y_mes[2*freq_num];
    int lay_num = (int)(args[1] + 0.1);
    double rho[lay_num];
    double cole[3*pol_num];
    double rho_DA_ini[4];
    double x_ini[lay_num + 1];
    double cole_ini[3*pol_num];
    double S_ini[(lay_num + 1 + 3*pol_num)*(lay_num+1+3*pol_num)];
    double S_pre[(lay_num + 1 + 3*pol_num)*(lay_num+1+3*pol_num)];
    double S0[(lay_num+1+3*pol_num)*(lay_num+1+3*pol_num)];
    double y_ini[2*freq_num];
    double upper[lay_num + 1 + 3*pol_num];
    double lower[lay_num + 1 + 3*pol_num];
    double upper_ini[2];
    double lower_ini[2];
    int hor_dist_pos = (int)(args[8 + 3*freq_num]);
    int ver_dist_pos = (int)(args[9 + 3*freq_num]);
    int alt_pos = (int)(args[10 + 3*freq_num]);
    int first_mes_pos = (int)(args[6 + 3*freq_num]);		    
    int up = 0;
    int s7 = AVERAGE, s7c = 0;
    double mesv[2*freq_num];
    double alta = 0;
    double vda = 0;
    int ft = 1;
    int start_with_rho_ini = 1;
    int data_cntr = 0;
    double dpth[lay_num],d=args[2];
    int nlay = lay_num;
    double weight = 1000.;
    for(int i=0;i<lay_num;i++)
        upper[i] = 20000;
    upper[lay_num] = d;
    for(i = 0; i < pol_num; i++){
        upper[lay_num + 1 + 3*i] = 20000;
	lower[lay_num + 1 + 3*i] = 0.01;
        upper[lay_num + 1 + 3*i+1] = 10;
	lower[lay_num + 1 + 3*i+1] = 0.000001;
        upper[lay_num + 1 + 3*i+2] = 3;
	lower[lay_num + 1 + 3*i+2] = 0.000001;
    }
    for(int i=0;i<lay_num;i++)
        lower[i] = 0.01;
    lower[lay_num] = -d;
    upper_ini[0] = 20000;
    upper_ini[1] = 0;
    lower_ini[0] = 0.01;
    lower_ini[1] = 0;

    memset(mesv,0,sizeof(mesv));

    memset(time, 0, sizeof(time));

    // Layers thicknesses are fixed here !!!
    for(int i=0;i<lay_num-1;i++,d*=args[3]) dpth[i] = d;

    data = buf + 13; // it's an offset to skip time in the string

    FILE *fin  = fopen(argv[1],"rt");
    FILE *fout = fopen(argv[2],"wt");
    printf("data files opened");
    //double primField(double hd, double vd)
    double mom = 10000/primField(28.8,27.8);
//    double mom = 10000/primField(28.28,28.33);
    geo.prim = 1./mom; //5099041845; // Equator's primary to get 10000 A/m in FD!!! - greater than 5.e9 Am2
                       //5563098240;//!!!!!! For 2200:5 in Croatia!!!! // primField(28.21,28.54)/.92;
                       //6489173050;//!!!!!! For 2110:2 in Croatia!!!! // primField(29.3,27.8)/.82;
    geo.w = (double *)malloc(2*freq_num*sizeof(double));
    memset(geo.w,0,sizeof(geo.w));
    //Setting R_i_i
    
    for(i = 0; i < 2*freq_num; i++)
	    geo.w[i] = 1/args[4 + freq_num + i];

    for(int i=0;i<freq_num-1;i++)
        geo.w[i*2+1] = sqrt(geo.w[i*2+1]*geo.w[i*2+1]+geo.w[(i+1)*2+1]*geo.w[(i+1)*2+1]);

    int lr;
    for(lr = 0;lr<lay_num;lr++)
        rho[lr] = RES_INI;

    for(lr = 0; lr < pol_num; lr++){
        cole[3*lr] = RES_INI;
	cole[3*lr + 1] = 0.001;
        cole[3*lr + 2] = 0.5;	
    }

    memset(x_ini,0,sizeof(x_ini));
    memset(S_ini,0,sizeof(S_ini));
    memset(buf,0,sizeof(buf));
    for(int i=0;i<lay_num;i++)
        x_ini[i] = rho[i];
    x_ini[lay_num] = 0;
    for(i = 0; i < 3*pol_num; i++)
	cole_ini[i] = cole[i];

    for(i=0;fgets(buf,3000,fin);i++) {
        if(buf[0]=='/') {  // reading comments
            memset(buf,0,sizeof(buf));
            continue;
        }

        if(buf[0]=='L' || buf[0]=='B' || buf[0]=='T') { // reading Lines or Base or Trends
            fputs(buf,fout);
            printf("%s",buf);
            for(int lr = 0;lr<lay_num;lr++)
                rho[lr] = RES_INI;
            for(lr = 0; lr < pol_num; lr++){
                cole[3*lr] = RES_INI;
	        cole[3*lr + 1] = 0.001;
                cole[3*lr + 2] = 0.5;	
            }

            memset(buf,0,sizeof(buf));
            ft = 1;
            continue;
        }

        memcpy(time,buf,sizeof(time)-1); // reading time

        if(strstr(data,"*")) { // all the nonmeasurements are skipped
            memset(buf,0,sizeof(buf));
            continue;
        }
        is_digit = 0;
        j = 0;
        read_num = 0;
        for(i = 0; i < 2000; i++){
            if((data[i]>='0' && data[i]<='9')||data[i]=='.'||data[i]==','){
                 tmp[j] = data[i];
	         j++;
	         is_digit = 1;
	    } else{
                 if(is_digit==1){
		    values[read_num] = atof(tmp);
		    read_num++;
                    is_digit = 0;
		    for(j = 0; j < 256; j++)
	               tmp[j] = (char) 0;
	         }
	         j = 0;
            }
        }
        ft = 1; 
	geo.hor_dist = values[hor_dist_pos];
	geo.ver_dist = values[ver_dist_pos];
	geo.alt = values[alt_pos];

	for(i = 0; i<freq_num; i++){
		y_mes[2*i] = values[first_mes_pos + freq_num + i];
                y_mes[2*i + 1] = values[first_mes_pos + i]; 
        }


        memset(buf,0,sizeof(buf));
        data_cntr++;
        //if(data_cntr<1400) continue;

        geo.alt += DA;

    // averaging for AVERAGE samples
        for(int i=0;i<2*freq_num;i++)
            mesv[i] += y_mes[i];
        alta += geo.alt;
        vda += geo.ver_dist;
        s7c++;
        if(--s7)
            continue;
        for(int i=0;i<2*freq_num;i++)
            y_mes[i] = mesv[i]/s7c;
        geo.alt = alta/s7c;
        geo.ver_dist = vda/s7c;
        s7c = 0;
        s7 = AVERAGE;
        memset(mesv,0,sizeof(mesv));
        alta = 0;
        vda = 0;

	// excluding all negatives
        for(int i=0;i<2*freq_num;i++)
           y_mes[i] = fabs(y_mes[i]);

    // It's a cricial point: we use differences of Inphase components instead of Inphase components!!!
            for(i=0;i<freq_num;i++) {
                if(i<freq_num-1){
                    y_mes[2*i+1] = fabs(y_mes[2*(i+1)+1] - y_mes[2*i+1]);
		}
                else{
		    if(freq_num > 1)
                           y_mes[2*i+1] = y_mes[2*(i-1)+1];
		}
            }
        

    // calculating of a half-space for the lowerest frequency!!!
    // itereative inversion for fixed layers
    

        double res = -1;
        double rho_ini = RES_INI;
	rho_DA_ini[0] = RES_INI;
	rho_DA_ini[1] = 0;
	rho_DA_ini[2] = 0;
	rho_DA_ini[3] = 0;
        int itr;
	double empty_cole[1];
	int empty_pol_inds[1];
	empty_pol_inds[0] = -1;
        for (itr = 0;itr < MAX_ITER; itr++) {
            double d_ini[4];
	    memset(d_ini, 0, sizeof(d_ini));
	    d_ini[0] = ERR_INI;
            d_ini[3] = 0.1;
            int up = 0;
            flinversion(geo,1,1, rho_DA_ini, empty_cole, empty_pol_inds, 0, dpth,y_ini,y_mes,&res,&up, d_ini, freqs, upper_ini, lower_ini);
            rho_ini = rho_DA_ini[0];
	    if(sqrt(res) <STOP_VAL) break;
            if(up) break;
        }

        res = sqrt(res);

        fprintf(fout, "%s %f %d %f",time, res, itr, rho_DA_ini[0]);
        fprintf(fout, "\n");


	//inversion
        int iter = 0;
        //double res1;

        double Sr[(lay_num + 1 + 3*pol_num)*(lay_num+1+3*pol_num)];
        res = -1;
        up = 0;

        if(ft) { // for each new line restart KF
            memset(S_ini,0,sizeof(S_ini));
            // ordinary resistivity
            for(int i=nlay-1;i>=0;i--) {
                if(i==nlay-1)
                    S_ini[i+(nlay+1+3*pol_num)*i] = ERR_INI;
                else {
                    S_ini[i+1+(nlay+1+3*pol_num)*i] =
                            COR_INI*ERR_INI/S_ini[i+1+(nlay+1+3*pol_num)*(i+1)];
                    S_ini[i+(nlay+1+3*pol_num)*i] =
                     sqrt(ERR_INI*ERR_INI-S_ini[i+1+(nlay+1+3*pol_num)*i]*
                            S_ini[i+1+(nlay+1+3*pol_num)*i]);
                }
		S_ini[nlay + nlay*(nlay+1+3*pol_num)] = 0.25;
		if(start_with_rho_ini)
                	x_ini[i] = rho_ini;
            }
	    for(i = nlay + 1; i < nlay+ 1 + 3*pol_num; i++){
		switch(i%3){
                    case 0:
			    cole_ini[i] = rho_ini;
			    break;
		    case 1:
			    cole_ini[i] = 0.001;
			    break;
	            case 2:
			    cole_ini[i] = 0.5;
			    break;
		}
	        S_ini[i + i*(nlay + 3*pol_num + 1)] = 0.1;
	    }
	    x_ini[nlay] = 0;
	    //start_with_rho_ini = 0;
            ft = !ft;
            memcpy(S0,S_ini,sizeof(Sr));
        } else { // restart S for the next point. x_ini is tied to the initial resistivity
            memset(S_ini,0,sizeof(S_ini));
            for(int i=nlay-1;i>=0;i--) {
                S_ini[i+(nlay+1+3*pol_num)*i] = S0[i+(nlay+1+3*pol_num)*i];
                if(i<nlay-1)
                    S_ini[i+(nlay+1+3*pol_num)*i+1] = S0[i+(nlay+1+3*pol_num)*i+1];
                x_ini[i] = (rho_ini*(1-1./weight)+x_ini[i]/weight);
            }
	    x_ini[nlay] = 0;
	    for(i = nlay+1; i < nlay+1+3*pol_num; i++){
                switch(i%3){
                    case 0:
			    cole_ini[i] = (rho_ini*(1-1./weight)+cole_ini[i]/weight);;
			    break;
		    case 1:
			    cole_ini[i] = (0.001*(1-1./weight)+cole_ini[i]/weight);;
			    break;
	            case 2:
			    cole_ini[i] = (0.5*(1-1./weight)+x_ini[i]/weight);;
			    break;
		}

                S_ini[i+(nlay+1+3*pol_num)*i] = S0[i+(nlay+1+3*pol_num)*i];
	    }
        }
        memcpy(Sr,S_ini,sizeof(Sr));

	// itereative inversion for fixed layers
        for (iter = 0;iter < MAX_ITER; iter++) {
            memcpy(S_ini,Sr,sizeof(Sr));
            flinversion(geo,freq_num,nlay,x_ini,cole_ini, pol_inds, pol_num, dpth,y_ini,y_mes,&res,&up,S_ini, freqs, upper, lower);
	/*
        for(int i=0;i<nlay;i++) {
            double v = 0, v1 = 0;
            for (int j=i;j<nlay;j++) {
                v+= S_ini[i*nlay+j]*S_ini[i*nlay+j];
                v1+= Sr[i*nlay+j]*Sr[i*nlay+j];
            }
            printf("%f ",1 - sqrt(v/v1));
        }
	printf("\n");

	for(int i = 0; i < nlay; i++) printf("%f ", x_ini[i]);
	printf("\n");
        */
            if(sqrt(res) <STOP_VAL) break;
            if(up) break;
        }
        
        res = sqrt(res);
        weight = (res>1)?res:1.;

        printf("%s %f %d ",time, res, iter);
        printf("\n");

        fprintf(fout, "%s %f %d ",time, res, iter);
        for(int i=0;i<nlay + 1;i++) {
                fprintf(fout,"%.3f ",x_ini[i]);
                printf("%f ",x_ini[i]);
        }
        double depp = 0;
        for(int i=0;i<nlay-1;i++) {
            fprintf(fout,"%.3f ",depp+dpth[i]*.5);
            depp += dpth[i];
        }
        printf("%f ",dpth[0]);

        fprintf(fout,"%.3f %.3f ",depp+dpth[nlay-2]*.5,rho_ini);
        printf("\n");

        for(int i=0;i<nlay+1;i++) {
            double v = 0, v1 = 0;
            for (int j=i;j<nlay+1;j++) {
                v+= S_ini[i*(nlay+1+3*pol_num)+j]*S_ini[i*(nlay+1+3*pol_num)+j];
                v1+= Sr[i*(nlay+1+3*pol_num)+j]*Sr[i*(nlay+1+3*pol_num)+j];
            }
            fprintf(fout,"%.7f ",1 - sqrt(v/v1));
            printf("%f ",1 - sqrt(v/v1));
        }

        fprintf(fout,"\n");
        printf("\n");
        continue;
        // here ends the routine
    }
    printf("\n end");
    return 0;
}


