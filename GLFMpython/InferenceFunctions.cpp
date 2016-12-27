#include "GeneralFunctions.h"

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"

// Functions
int AcceleratedGibbs (int maxK,int bias, int N, int D, int K, char *C,  int *R, double alpha, double s2B, double s2Y, gsl_matrix **Y, gsl_matrix *Z, int *nest, gsl_matrix *P, gsl_matrix *Pnon, gsl_matrix **lambda, gsl_matrix **lambdanon){
   int TK=2;
   gsl_matrix_view Zn;
   gsl_matrix_view Ydn;
   gsl_matrix_view Pnon_view;
   gsl_matrix_view Lnon_view;
   //gsl_matrix_view P_view;
   //gsl_matrix_view L_view;
   gsl_matrix *muy;
   gsl_matrix *s2y_p= gsl_matrix_alloc(1,1);
   gsl_matrix *aux;
   gsl_matrix *Snon;
   double s2y_num;
   
   gsl_matrix_memcpy (Pnon, P);
   for (int d =0; d<D; d++){
        gsl_matrix_memcpy (lambdanon[d], lambda[d]);
    }
   
   for (int n=0; n<N; n++){
       double p[TK];       
       // Pnon, LambdaNon
       Zn = gsl_matrix_submatrix (Z, 0, n, K, 1);
       Pnon_view = gsl_matrix_submatrix (Pnon, 0, 0, K, K);
       Snon = gsl_matrix_calloc(K,K);
       gsl_matrix_memcpy (Snon, &Pnon_view.matrix);
       inverse(Snon, K);
       matrix_multiply(&Zn.matrix,&Zn.matrix,&Pnon_view.matrix,-1,1,CblasNoTrans,CblasTrans);
       for (int d =0; d<D; d++){
               Lnon_view= gsl_matrix_submatrix (lambdanon[d], 0, 0, K, R[d]);
               Ydn= gsl_matrix_submatrix (Y[d], 0, n, R[d], 1);
               matrix_multiply(&Zn.matrix,&Ydn.matrix,&Lnon_view.matrix,-1,1,CblasNoTrans,CblasTrans);
       }
       // Sampling znk for k=1...K
       for (int k=bias; k<K; k++){
           if (gsl_matrix_get(&Zn.matrix,k,0)==1){nest[k]--;}
          
           if (nest[k]>0){
               aux = gsl_matrix_alloc(1,K);
               // z_nk=0
               gsl_matrix_set (&Zn.matrix, k, 0, 0);
               matrix_multiply(&Zn.matrix,Snon,aux,1,0,CblasTrans,CblasNoTrans);
               gsl_matrix_set(s2y_p,0,0,s2Y);
               matrix_multiply(aux,&Zn.matrix,s2y_p,1,1,CblasNoTrans,CblasNoTrans);
               s2y_num=gsl_matrix_get(s2y_p,0,0);
               double lik0=0;
               for (int d =0; d<D; d++){
                   Ydn= gsl_matrix_submatrix (Y[d], 0, n, R[d], 1);
                   Lnon_view= gsl_matrix_submatrix (lambdanon[d], 0, 0, K, R[d]);
                   muy= gsl_matrix_alloc(1,R[d]);
                   matrix_multiply(aux,&Lnon_view.matrix,muy,1,0,CblasNoTrans,CblasNoTrans);   
                   if (C[d]=='c'){
                       for (int r=0; r<R[d]-1; r++){
                           lik0-=0.5/s2y_num*pow((gsl_matrix_get(&Ydn.matrix,r,0)- gsl_matrix_get(muy,0,r)),2)+0.5*gsl_sf_log(2*M_PI*s2y_num);
                           }
                       }else{
                           lik0-=0.5/s2y_num*pow((gsl_matrix_get(&Ydn.matrix,0,0)- gsl_matrix_get(muy,0,0)),2)+0.5*gsl_sf_log(2*M_PI*s2y_num);
                       }
                   gsl_matrix_free(muy);
               }

               // z_nk=1
               gsl_matrix_set (&Zn.matrix, k, 0, 1);
               matrix_multiply(&Zn.matrix,Snon,aux,1,0,CblasTrans,CblasNoTrans);
               gsl_matrix_set(s2y_p,0,0,s2Y);
               matrix_multiply(aux,&Zn.matrix,s2y_p,1,1,CblasNoTrans,CblasNoTrans);
               s2y_num=gsl_matrix_get(s2y_p,0,0);
               double lik1=0;
               for (int d =0; d<D; d++){               
                   Ydn= gsl_matrix_submatrix (Y[d], 0, n, R[d], 1);
                   Lnon_view= gsl_matrix_submatrix (lambdanon[d], 0, 0, K, R[d]);
                   muy= gsl_matrix_alloc(1,R[d]);
                   matrix_multiply(aux,&Lnon_view.matrix,muy,1,0,CblasNoTrans,CblasNoTrans); 
                   if (C[d]=='c'){
                       for (int r=0; r<R[d]-1; r++){
                           lik1-=0.5/s2y_num*pow((gsl_matrix_get(&Ydn.matrix,r,0)- gsl_matrix_get(muy,0,r)),2)+0.5*gsl_sf_log(2*M_PI*s2y_num);
                           }
                   }else{
                       lik1-=0.5/s2y_num*pow((gsl_matrix_get(&Ydn.matrix,0,0)- gsl_matrix_get(muy,0,0)),2)+0.5*gsl_sf_log(2*M_PI*s2y_num);
                   }
                   gsl_matrix_free(muy);
               }

               //printf("lik0=%f , lik1=%f \n", lik0, lik1);
               double p0= gsl_sf_log(N-nest[k])+lik0;
               double p1= gsl_sf_log(nest[k])+lik1;
               double p1_n, p0_n;
               //printf("p1=%f p0=%f \n", p1,p0);
               if (p0>p1){               
                   p1_n=expFun(p1-p0);
                   p0_n=1;
               }else{ 
                   p0_n=expFun(p0-p1);
                   p1_n=1;
               }
               p1_n=p1_n/(p1_n+p0_n);
               //printf("p1_n=%f ", p1_n);
                
               //sampling znk
               if (drand48()>p1_n){
                   gsl_matrix_set (&Zn.matrix, k, 0, 0);
                   p[0]=lik0;
               }else{
                   nest[k]+=1;
                   p[0]=lik1;}
               gsl_matrix_free(aux);
           }else{
               gsl_matrix_set (&Zn.matrix, k, 0, 0);
           } 
           
        
       }
       gsl_matrix_free(Snon);

       // remove empty features
       int flagDel=0;
       int Kdel=0;
       for (int k=0; k<K;k++){
           if ((nest[k]==0) & (K>1)){
                Kdel++;
                flagDel=1;
                for (int kk=k; kk<K-1; kk++){
                    gsl_vector_view Zrow= gsl_matrix_row (Z,kk+1);
                    gsl_matrix_set_row (Z, kk, &Zrow.vector);
                    nest[kk]=nest[kk+1];
                }
            }
               
       }
       for (int k=K-Kdel; k<K; k++){
           gsl_vector_view Zrow= gsl_matrix_row (Z,k);
           gsl_vector_set_zero (&Zrow.vector);
           nest[k]=0;
           }       
       K-=Kdel;
       
       
       if (flagDel){ 
           gsl_matrix_set_identity (P);
           matrix_multiply(Z,Z,P,1,s2B,CblasNoTrans,CblasTrans);
           gsl_matrix_memcpy (Pnon, P);
           Pnon_view = gsl_matrix_submatrix (Pnon, 0, 0, K, K);
           Zn = gsl_matrix_submatrix (Z, 0, n, K, 1);
           matrix_multiply(&Zn.matrix,&Zn.matrix,&Pnon_view.matrix,-1,1,CblasNoTrans,CblasTrans);
           for (int d =0; d<D; d++){
                matrix_multiply(Z,Y[d],lambda[d],1,0,CblasNoTrans,CblasTrans);
                gsl_matrix_memcpy (lambdanon[d], lambda[d]);
                Lnon_view= gsl_matrix_submatrix (lambdanon[d], 0, 0, K, R[d]);
                Ydn= gsl_matrix_submatrix (Y[d], 0, n, R[d], 1);
                matrix_multiply(&Zn.matrix,&Ydn.matrix,&Lnon_view.matrix,-1,1,CblasNoTrans,CblasTrans);
           }
       
       }
       
       // Adding new features
       gsl_matrix_view Znew= gsl_matrix_submatrix (Z, K, 0, TK, N);
       gsl_matrix_set_zero (&Znew.matrix);
       
       if (K+TK<maxK){
       double pmax=p[0];
       double kk=0;
       for (int j=1; j<TK; j++){  
           gsl_vector_view Pnon_colum= gsl_matrix_column (Pnon, K+j-1);
           gsl_vector_set_zero (&Pnon_colum.vector);
           Pnon_colum= gsl_matrix_row (Pnon, K+j-1);
           gsl_vector_set_basis (&Pnon_colum.vector, K+j-1);
           
           aux = gsl_matrix_alloc(1,K+j);
           Pnon_view = gsl_matrix_submatrix (Pnon, 0, 0, K+j, K+j);
           Snon = gsl_matrix_calloc(K+j,K+j);
           gsl_matrix_memcpy (Snon, &Pnon_view.matrix);
           inverse(Snon, K+j);
           Zn = gsl_matrix_submatrix (Z, 0, n, K+j, 1);
           gsl_matrix_set (&Zn.matrix, K+j-1, 0, 1);
           matrix_multiply(&Zn.matrix,Snon,aux,1,0,CblasTrans,CblasNoTrans);
           gsl_matrix_set(s2y_p,0,0,s2Y);
           matrix_multiply(aux,&Zn.matrix,s2y_p,1,1,CblasNoTrans,CblasNoTrans);
           s2y_num=gsl_matrix_get(s2y_p,0,0);
           double lik=0;
           for (int d =0; d<D; d++){               
               Ydn= gsl_matrix_submatrix (Y[d], 0, n, R[d], 1);
               Lnon_view= gsl_matrix_submatrix (lambdanon[d], 0, 0, K+j, R[d]);
               muy= gsl_matrix_alloc(1,R[d]);
               matrix_multiply(aux,&Lnon_view.matrix,muy,1,0,CblasNoTrans,CblasNoTrans); 
               if (C[d]=='c'){
                   for (int r=0; r<R[d]-1; r++){
                       lik-=0.5/s2y_num*pow((gsl_matrix_get(&Ydn.matrix,r,0)- gsl_matrix_get(muy,0,r)),2)+0.5*gsl_sf_log(2*M_PI*s2y_num);
                       }
               }else{
                   lik-=0.5/s2y_num*pow((gsl_matrix_get(&Ydn.matrix,0,0)- gsl_matrix_get(muy,0,0)),2)+0.5*gsl_sf_log(2*M_PI*s2y_num);
               }
               gsl_matrix_free(muy);
           }
           p[j]=lik+j*gsl_sf_log(alpha/N)-gsl_sf_log(factorial(j));
           gsl_matrix_free(aux);
           gsl_matrix_free(Snon);
     
           if(pmax<p[j]){pmax=p[j];}
           kk=j;
       }
   double den=0;
   for (int k=0;k<=kk; k++){
       p[k]-=pmax;
       p[k]=expFun(p[k]);
       den+=p[k];
       
   }
   for (int k=0;k<=kk; k++){
       p[k]=p[k]/den;
   }
   int Knew=mnrnd(p, kk);
   if (Knew>0){
       for (int k=K; k<K+Knew; k++){nest[k]=1;}
       }
   K+=Knew;
   }
   //Adding Zn
   Zn = gsl_matrix_submatrix (Z, 0, n, K, 1);
   Pnon_view = gsl_matrix_submatrix (Pnon, 0, 0, K, K);
   matrix_multiply(&Zn.matrix,&Zn.matrix,&Pnon_view.matrix,1,1,CblasNoTrans,CblasTrans);
   gsl_matrix_memcpy (P, Pnon);
   for (int d =0; d<D; d++){
        Ydn= gsl_matrix_submatrix (Y[d], 0, n, R[d], 1);
        Lnon_view= gsl_matrix_submatrix (lambdanon[d], 0, 0, K, R[d]);
        matrix_multiply(&Zn.matrix,&Ydn.matrix,&Lnon_view.matrix,1,1,CblasNoTrans,CblasTrans);
        gsl_matrix_memcpy (lambda[d], lambdanon[d]);
   }
   
   }
   gsl_matrix_free(s2y_p);    
   
   return K;
}



//Sample Y
void SampleY (double missing, int N, int d, int K, char Cd,  int Rd, double wd, double s2Y, double s2theta, gsl_matrix *X, gsl_matrix *Z, gsl_matrix *Yd,  gsl_matrix *Bd, gsl_vector *thetad, const gsl_rng *seed){
    double s2u=1;
    gsl_matrix_view Zn;
    gsl_matrix_view Bd_view;
    gsl_matrix *muy;
    double xnd;
    switch(Cd){
        case 'g':
            muy= gsl_matrix_alloc(1,1);
            for (int n=0; n<N; n++){
                xnd=gsl_matrix_get (X, d, n);
                if (xnd==missing || gsl_isnan(xnd)){
                    Zn = gsl_matrix_submatrix (Z, 0, n, K, 1);
                    Bd_view = gsl_matrix_submatrix (Bd, 0, 0, K, 1);
                    
                    matrix_multiply(&Zn.matrix,&Bd_view.matrix,muy,1,0,CblasTrans,CblasNoTrans); 
                    gsl_matrix_set (Yd, 0, n, gsl_matrix_get(muy,0,0)+  gsl_ran_gaussian (seed, s2Y));

                }
                
            }
            gsl_matrix_free(muy);    
            
            break;
            
        case 'p': 
            muy= gsl_matrix_alloc(1,1);
            for (int n=0; n<N; n++){
                xnd=gsl_matrix_get (X, d, n);
                Zn = gsl_matrix_submatrix (Z, 0, n, K, 1);
                Bd_view = gsl_matrix_submatrix (Bd, 0, 0, K, 1);
                matrix_multiply(&Zn.matrix,&Bd_view.matrix,muy,1,0,CblasTrans,CblasNoTrans);
                if (xnd==missing || gsl_isnan(xnd)){
                     gsl_matrix_set (Yd, 0, n, gsl_matrix_get(muy,0,0)+  gsl_ran_gaussian (seed, s2Y));
                }else{
                    gsl_matrix_set (Yd, 0, n, (f_1(xnd, wd)/s2u + gsl_matrix_get(muy,0,0)/s2Y)/(1/s2Y+1/s2u) +  gsl_ran_gaussian (seed, 1/(1/s2Y+1/s2u)));
                }
            }
            gsl_matrix_free(muy); 
            break;
            
        case 'n': 
            muy= gsl_matrix_alloc(1,1);
            for (int n=0; n<N; n++){
                xnd=gsl_matrix_get (X, d, n);
                Zn = gsl_matrix_submatrix (Z, 0, n, K, 1);
                Bd_view = gsl_matrix_submatrix (Bd, 0, 0, K, 1);
                matrix_multiply(&Zn.matrix,&Bd_view.matrix,muy,1,0,CblasTrans,CblasNoTrans);
                if (xnd==missing || gsl_isnan(xnd)){
                     gsl_matrix_set (Yd, 0, n, gsl_matrix_get(muy,0,0)+  gsl_ran_gaussian (seed, s2Y));
                }else{
                    gsl_matrix_set (Yd, 0, n, truncnormrnd(gsl_matrix_get(muy,0,0), s2Y, f_1(xnd,wd),f_1(xnd+1, wd)));
                }
            }
            gsl_matrix_free(muy); 
            break;
            
        case 'c': 
            muy= gsl_matrix_alloc(1,Rd);
            for (int n=0; n<N; n++){
                xnd=gsl_matrix_get (X, d, n);   
                Zn = gsl_matrix_submatrix (Z, 0, n, K, 1);
                Bd_view = gsl_matrix_submatrix (Bd, 0, 0, K, Rd);
                matrix_multiply(&Zn.matrix,&Bd_view.matrix,muy,1,0,CblasTrans,CblasNoTrans);
                if (xnd==missing || gsl_isnan(xnd)){
                    for(int r=0; r<Rd; r++){
                        gsl_matrix_set (Yd, r, n, gsl_matrix_get(muy,0,r)+gsl_ran_gaussian (seed, s2Y));
                    }
                }else{
                    double maxY=0;
                    double ytrue=gsl_matrix_get (Yd, xnd-1, n);
                    for(int r=0; r<Rd; r++){
                        double ydr= gsl_matrix_get (Yd, r, n);
                        if ((ydr!=ytrue) & (ydr>maxY)){maxY=ydr;}
                    }
                    gsl_matrix_set (Yd, xnd-1, n, truncnormrnd(gsl_matrix_get(muy,0,xnd-1), s2Y, maxY, GSL_POSINF));
                    for(int r=0; r<Rd; r++){
                        if (r!=xnd-1){
                            gsl_matrix_set (Yd, r, n, truncnormrnd(gsl_matrix_get(muy,0,r), s2Y, GSL_NEGINF, gsl_matrix_get (Yd, xnd-1, n)));
                        }
                    }
                }
            }
            gsl_matrix_free(muy);
            break;
          case 'o': 
                // Sample Y
                gsl_vector *Ymax= gsl_vector_calloc (Rd); 
                gsl_vector *Ymin= gsl_vector_alloc (Rd); 
                gsl_vector_set_all (Ymin, GSL_POSINF);
                muy= gsl_matrix_alloc(1,1);
                for (int n=0; n<N; n++){
                    xnd=gsl_matrix_get (X, d, n);
                    Zn = gsl_matrix_submatrix (Z, 0, n, K, 1);
                    Bd_view = gsl_matrix_submatrix (Bd, 0, 0, K, 1);
                    
                    matrix_multiply(&Zn.matrix,&Bd_view.matrix,muy,1,0,CblasTrans,CblasNoTrans);
                    if (xnd==missing || gsl_isnan(xnd)){                        
                         gsl_matrix_set (Yd, 0, n, gsl_matrix_get(muy,0,0)+gsl_ran_gaussian (seed, s2Y));
                    }else if (xnd==1){
                         gsl_matrix_set(Yd, 0, n, truncnormrnd(gsl_matrix_get(muy,0,0), s2Y, GSL_NEGINF, gsl_vector_get (thetad, xnd-1)));
                         if (gsl_matrix_get(Yd, 0, n)>gsl_vector_get(Ymax,xnd-1)){gsl_vector_set(Ymax,xnd-1,gsl_matrix_get(Yd, 0, n));}
                         if (gsl_matrix_get(Yd, 0, n)<gsl_vector_get(Ymin,xnd-1)){gsl_vector_set(Ymin,xnd-1,gsl_matrix_get(Yd, 0, n));}
                    }else{
                         gsl_matrix_set(Yd, 0, n, truncnormrnd(gsl_matrix_get(muy,0,0), s2Y, gsl_vector_get (thetad, xnd-2), gsl_vector_get (thetad, xnd-1)));
                         if (gsl_matrix_get(Yd, 0, n)>gsl_vector_get(Ymax,xnd-1)){gsl_vector_set(Ymax,xnd-1,gsl_matrix_get(Yd, 0, n));}
                         if (gsl_matrix_get(Yd, 0, n)<gsl_vector_get(Ymin,xnd-1)){gsl_vector_set(Ymin,xnd-1,gsl_matrix_get(Yd, 0, n));}
                    }
                    
                }
                break; 
                gsl_matrix_free(muy);
                
                //Sample Theta
                for(int r=1; r<Rd-1; r++){
                    double xlo;
                    double xhi;
                    if( gsl_vector_get (thetad, r)>gsl_vector_get(Ymax,r)){xlo=gsl_vector_get (thetad, r);}
                    else{xlo=gsl_vector_get(Ymax,r);}
                    if( gsl_vector_get (thetad, r+1)<gsl_vector_get(Ymin,r+1)){xhi=gsl_vector_get (thetad, r+1);}
                    else{xhi=gsl_vector_get(Ymin,r+1);}
                    gsl_vector_set (thetad, r, truncnormrnd(0, s2theta, xlo, xhi));
                }
            
            
        
   
    }
    
}

int IBPsampler_func (double missing, gsl_matrix *X, char *C, gsl_matrix *Z, gsl_matrix **B, gsl_vector **theta, int *R, double *w, int maxR, int bias, int N, int D, int K, double alpha, double s2B, double s2Y,int maxK,int Nsim){
//Starting C function
       //.....INIZIALIZATION........//
    double s2theta=1;
    // random numbers
    srand48(time(NULL));
    gsl_rng *seed = gsl_rng_alloc(gsl_rng_taus);
    time_t clck=time(NULL);
    gsl_rng_set(seed, clck);
    
    // auxiliary variables
    int Kest=K;
    gsl_matrix *P= gsl_matrix_alloc(maxK,maxK);
    gsl_matrix_set_identity (P);
    matrix_multiply(Z,Z,P,1,s2B,CblasNoTrans,CblasTrans);
    gsl_matrix *Pnon= gsl_matrix_alloc(maxK,maxK);
    
    // initialize counts
    int nest[maxK];
    for (int k=0; k<Kest; k++){
        nest[k]=(gsl_matrix_get (P, k, k)-1);
        }
    
    gsl_matrix **Y=(gsl_matrix **) calloc(D,sizeof(gsl_matrix*));
    gsl_matrix **lambda=(gsl_matrix **) calloc(D,sizeof(gsl_matrix*));
    gsl_matrix **lambdanon=(gsl_matrix **) calloc(D,sizeof(gsl_matrix*));
    for (int d =0; d<D; d++){
        //Initialize Y !!!!!!
         switch(C[d]){
            double xnd;
            case 'g':
                Y[d] =gsl_matrix_alloc(1,N); 
                for (int n=0; n<N; n++){
                    xnd=gsl_matrix_get (X, d, n);
                    if (xnd==missing || gsl_isnan(xnd)){
                        gsl_matrix_set (Y[d], 0, n, gsl_ran_gaussian (seed, s2Y));
                   }else{
                        gsl_matrix_set(Y[d], 0, n, xnd);
                    }

                }

                break;

            case 'p': 
                Y[d] =gsl_matrix_alloc(1,N); 
                for (int n=0; n<N; n++){
                    xnd=gsl_matrix_get (X, d, n);
                    
                    if (xnd==missing || gsl_isnan(xnd)){
                         gsl_matrix_set (Y[d], 0, n, gsl_ran_gaussian (seed, s2Y));
                    }else{
                         gsl_matrix_set(Y[d], 0, n, f_1(xnd,w[d])+gsl_ran_gaussian (seed, s2Y));
                    }
                }
                break;

            case 'n':
                Y[d] =gsl_matrix_alloc(1,N); 
                lambda[d] =gsl_matrix_calloc(maxK,1);
                matrix_multiply(Z,Y[d],lambda[d],1,0,CblasNoTrans,CblasTrans);
                for (int n=0; n<N; n++){
                    xnd=gsl_matrix_get (X, d, n);
                    
                    if (xnd==missing || gsl_isnan(xnd)){
                         gsl_matrix_set (Y[d], 0, n, gsl_ran_gaussian (seed, s2Y));
                    }else{
                         gsl_matrix_set(Y[d], 0, n, f_1(xnd+gsl_ran_beta (seed, 5,1),w[d]));
                    }
                }
                break;
                
            case 'c': 
                Y[d] =gsl_matrix_alloc(R[d],N);
                for (int n=0; n<N; n++){
                    xnd=gsl_matrix_get (X, d, n);                   
                    if (xnd==missing || gsl_isnan(xnd)){
                        for(int r=0; r<R[d]; r++){
                            gsl_matrix_set (Y[d], r, n, gsl_ran_gaussian(seed, s2Y));
                        }
                    }else{
                        gsl_matrix_set (Y[d], xnd-1, n, truncnormrnd(0, s2Y, 0, GSL_POSINF));
                        for(int r=0; r<R[d]; r++){
                            if (r!=xnd-1){
                                gsl_matrix_set (Y[d], r, n, truncnormrnd(0, s2Y, GSL_NEGINF, gsl_matrix_get (Y[d], xnd-1, n)));
                            }
                        }
                    }
                }
                break;
                
             case 'o': 
                Y[d] =gsl_matrix_alloc(R[d],N); 
                gsl_vector_set (theta[d], 0, gsl_ran_gaussian(seed, s2theta));
                for(int r=1; r<R[d]-1; r++){
                    gsl_vector_set (theta[d], r, truncnormrnd(0, s2theta, gsl_vector_get (theta[d], r-1), GSL_POSINF));
                }
                gsl_vector_set (theta[d], R[d]-1, GSL_POSINF);
                for (int n=0; n<N; n++){
                    xnd=gsl_matrix_get (X, d, n);
                    
                    if (xnd==missing || gsl_isnan(xnd)){
                         gsl_matrix_set (Y[d], 0, n, gsl_ran_gaussian (seed, s2Y));
                    }else if (xnd==1){
                         gsl_matrix_set(Y[d], 0, n, truncnormrnd(0, s2Y, GSL_NEGINF, gsl_vector_get (theta[d], xnd-1)));
                    }else{
                         gsl_matrix_set(Y[d], 0, n, truncnormrnd(0, s2Y, gsl_vector_get (theta[d], xnd-2), gsl_vector_get (theta[d], xnd-1)));
                    }
                }
                 break;

        }
            lambda[d] =gsl_matrix_calloc(maxK,R[d]);
            matrix_multiply(Z,Y[d],lambda[d],1,0,CblasNoTrans,CblasTrans);
            lambdanon[d]= gsl_matrix_calloc(maxK,R[d]);        
    }
    
    //....Body functions....//      
    for (int it=0; it<Nsim; it++){
        Kest=AcceleratedGibbs (maxK,bias,N, D, Kest, C, R, alpha, s2B, s2Y, Y, Z, nest, P, Pnon, lambda, lambdanon);
        
        gsl_matrix_view P_view = gsl_matrix_submatrix (P, 0, 0, Kest, Kest);
        gsl_matrix *S= gsl_matrix_calloc(Kest,Kest);
        gsl_matrix_memcpy (S, &P_view.matrix);
        inverse(S, Kest);
        gsl_matrix *MuB = gsl_matrix_alloc(Kest,1);
        for (int d =0; d<D; d++){
        //Sample Bs
             if (C[d]=='c'){
                 gsl_vector_view Bd_view;
                 for (int r=0; r<R[d]-1; r++){
                    gsl_matrix_view L_view= gsl_matrix_submatrix (lambda[d], 0, r, Kest, 1);                       
                    matrix_multiply(S,&L_view.matrix,MuB,1,0,CblasNoTrans,CblasNoTrans);                       
                    Bd_view =  gsl_matrix_subcolumn (B[d], r, 0, Kest);
                    gsl_vector_view MuB_view =  gsl_matrix_column (MuB, 0);
                    mvnrnd(&Bd_view.vector, S, &MuB_view.vector, Kest, seed);
                    
                 }  
                 Bd_view =  gsl_matrix_subcolumn (B[d], R[d]-1, 0, Kest);
                 gsl_vector_set_zero (&Bd_view.vector);
             }else{
                gsl_matrix_view Lnon_view= gsl_matrix_submatrix (lambda[d], 0, 0, Kest, 1);
                matrix_multiply(S,&Lnon_view.matrix,MuB,1,0,CblasNoTrans,CblasNoTrans);
                
                gsl_vector_view Bd_view =  gsl_matrix_subcolumn (B[d], 0, 0, Kest);
                gsl_vector_view MuB_view =  gsl_matrix_subcolumn (MuB, 0, 0, Kest);
                mvnrnd(&Bd_view.vector, S, &MuB_view.vector, Kest, seed);   
             }
             
         //Sample Y  
         SampleY (missing, N, d, Kest, C[d],  R[d], w[d], s2Y, s2theta, X, Z, Y[d],  B[d], theta[d], seed);
         
         //Update lambda
         matrix_multiply(Z,Y[d],lambda[d],1,0,CblasNoTrans,CblasTrans);     
            
         }
        gsl_matrix_free(S);
        }

    for (int d=0; d<D; d++){    
        gsl_matrix_free(Y[d]);
        gsl_matrix_free(lambda[d]);
        gsl_matrix_free(lambdanon[d]);
        //gsl_vector_free(theta[d]);
    }
    gsl_matrix_free(P);
    gsl_matrix_free(Pnon);
    free(lambda);
    free(lambdanon);
    free(Y);
    
    return Kest;
}


