#include <iostream>
#include <fstream>
#include <ctime>            
#include <boost/random/mersenne_twister.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <sys/time.h>

#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <queue>

#include <math.h>
#include <matrix.h>
#include <mex.h>

boost::random::mt19937 gen;

using namespace std;

int collisions (int &overlap, int &dimx_xp, deque<int> &xist_gene, deque<int> &tsix_gene, 
        int &sum_xp, int &sum_tp, boost::random::uniform_int_distribution<> &dist)
{
    int q;
    for (q=0; q<overlap; q++)
            {
            if ((xist_gene[q+(dimx_xp-overlap)] + tsix_gene[q])>1) {
                if (dist(gen)==0) {
                    xist_gene[q+(dimx_xp-overlap)]=0;
                    sum_xp--;
                } else {
                    tsix_gene[q]=0;
                    sum_tp--;
                }
            }
        }
}

int elongation (int &overlap, int &dimx_xp, deque<int> &xist_gene, deque<int> &tsix_gene, int &sum_xp, int &sum_tp, 
        bool &xpol_ini, bool &tpol_ini, int &xist_rna, int &tsix_rna, 
        vector<int> &all_xist, int &temp_sum_xr, int &temp_sum_tr, int &temp_sum_xp, int &temp_sum_tp,
        boost::random::uniform_int_distribution<> &dist, int &step, double &p8)
{
    xist_gene.push_front(xpol_ini);
    sum_xp = sum_xp - xist_gene.back() + xist_gene.front();
    xist_rna = xist_rna + (int)xist_gene.back();
    xist_gene.pop_back();
    if (p8>0) {
        collisions (overlap, dimx_xp, xist_gene, tsix_gene, sum_xp, sum_tp, dist);
    }   
    tsix_gene.push_back(tpol_ini);
    sum_tp = sum_tp - tsix_gene.front() + tsix_gene.back();
    tsix_rna = tsix_rna + (int)tsix_gene.front();
    tsix_gene.pop_front();
    if (p8>0) {
        collisions (overlap, dimx_xp, xist_gene, tsix_gene, sum_xp, sum_tp, dist);
    }              
    all_xist[step] = xist_rna;
                      
    temp_sum_xr = temp_sum_xr + xist_rna;
    temp_sum_tr = temp_sum_tr + tsix_rna;
    temp_sum_xp = temp_sum_xp + sum_xp;
    temp_sum_tp = temp_sum_tp + sum_tp;
}

void mexFunction(int nhls, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    //initiate random number generator
    struct timeval t1;
    gettimeofday(&t1, NULL);
    gen.seed(static_cast<unsigned int>(t1.tv_usec * t1.tv_sec));
    
	//declare variables // pointers that point at input and output variables
	mxArray *xist_pol, *tsix_pol, *xist_pol2, *tsix_pol2, *par, *const_par; // *xist_smooth;
    mxArray *time_out, *xist_pol_out, *tsix_pol_out, *xist_rna_out, *tsix_rna_out, 
            *xist_pol_out2, *tsix_pol_out2, *xist_rna_out2, *tsix_rna_out2, *test; //, *xist_prom_out;
	const mwSize *dims;
	double *xp, *tp, *xr, *tr, *xpr, *xp2, *tp2, *xr2, *tr2, *xpr2, *p, *cp, 
            t, sil_threshold, t_max;
    double *xpo, *tpo, *xro, *tro, *xpro, *xpo2, *tpo2, *xro2, *tro2, *xpro2, 
            *time, v = 1.0/1440.0, out_step, *to, k_adv_sil;
	int dimx_xp, dimx_tp, xist_rna, tsix_rna, xist_prom, xist_rna2, tsix_rna2, xist_prom2, overlap_X_T;
	int n_steps, n_out_steps, n_small_steps, sec_chr, step, sel_rx, tsix_silenced;
	
	//associate inputs: first chromosome
    par = mxDuplicateArray(prhs[0]);
    const_par = mxDuplicateArray(prhs[1]);
	xist_pol = mxDuplicateArray(prhs[2]);
    tsix_pol = mxDuplicateArray(prhs[3]);  
    xist_rna = mxGetScalar(prhs[4]);
    tsix_rna = mxGetScalar(prhs[5]);
    xist_prom = mxGetScalar(prhs[6]);
              
	//associate pointers
    xp = mxGetPr(xist_pol);
    tp = mxGetPr(tsix_pol);
    p = mxGetPr(par);
    cp = mxGetPr(const_par);
    t = cp[0];
    t_max =  cp[1];
    sil_threshold = cp[2];
    out_step = cp[3];
    sec_chr = cp[4];
    k_adv_sil = cp[5];
    tsix_silenced = cp[6];
    overlap_X_T = cp[7]; 
    
    //associate inputs: second chromosome
   if (sec_chr==1) {
        xist_pol2 = mxDuplicateArray(prhs[7]);
        tsix_pol2 = mxDuplicateArray(prhs[8]);  
        xist_rna2 = mxGetScalar(prhs[9]);
        tsix_rna2 = mxGetScalar(prhs[10]);
        xist_prom2 = mxGetScalar(prhs[11]);
        
        xp2 = mxGetPr(xist_pol2);
        tp2 = mxGetPr(tsix_pol2);
   }
    
    n_out_steps = 1+ceil(t_max/out_step);
    n_small_steps = round(out_step/v);
    n_steps = (n_out_steps-1)*n_small_steps;
	
	//figure out dimensions and associate output
	dims = mxGetDimensions(prhs[2]);
	dimx_xp = (int)dims[0];
    dims = mxGetDimensions(prhs[3]);
	dimx_tp = (int)dims[0];
    
    //output
    time_out = plhs[0] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
    xist_pol_out = plhs[1] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
    tsix_pol_out = plhs[2] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
    xist_rna_out = plhs[3] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
    tsix_rna_out = plhs[4] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
    xpo = mxGetPr(xist_pol_out);
    tpo= mxGetPr(tsix_pol_out);
    xro = mxGetPr(xist_rna_out);
    tro = mxGetPr(tsix_rna_out);
    time =   mxGetPr(time_out);
        
    xist_pol_out2 = plhs[5] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
    tsix_pol_out2 = plhs[6] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
    xist_rna_out2 = plhs[7] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
    tsix_rna_out2 = plhs[8] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
    test = plhs[9] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
        
    xpo2 = mxGetPr(xist_pol_out2);
    tpo2= mxGetPr(tsix_pol_out2);
    xro2 = mxGetPr(xist_rna_out2);
    tro2 = mxGetPr(tsix_rna_out2);
    to = mxGetPr(test);
    
    //write pol mxarray in deque and sum up polymerases
    deque<int> xist_gene, tsix_gene, xist_gene2, tsix_gene2;
    int q, sum_xp=0, sum_tp=0, sum_xp2 = 0, sum_tp2 = 0;
    
    for (q=0; q<dimx_xp; q++) 
        {
            xist_gene.push_back (xp[q]);
            sum_xp = sum_xp + xp[q];
             if (sec_chr==1) {
                xist_gene2.push_back (xp2[q]);
                sum_xp2 = sum_xp2 + xp2[q];
             }
        }
    for (q=0; q<dimx_tp; q++) 
        {
            tsix_gene.push_back (tp[q]);
            sum_tp = sum_tp + tp[q];
            if (sec_chr==1) {
                tsix_gene2.push_back (tp2[q]);
                 sum_tp2 = sum_tp2 + tp2[q];
            }
        }  
        
    bool xpol_ini=0, tpol_ini=0, xpol_ini2=0, tpol_ini2=0;
    boost::random::uniform_int_distribution<> dist(0, 1);
    boost::random::uniform_real_distribution<> dist_real(0.0, 1.0);
    double t_next, r, sum_rx;
    vector<double> rx(16), cum_rx(16);
    vector<int> all_xist(n_steps), all_xist2(n_steps);
    
    int  sil_tsix = 0, sil_tsix2 = 0, sil_state_tsix = 0, sil_state_tsix2 = 0; 
    int m = 0, b = 0, temp_sum_xr = 0, temp_sum_tr = 0, temp_sum_xp = 0, temp_sum_tp = 0,
            temp_sum_xr2 = 0, temp_sum_tr2 = 0, temp_sum_xp2 = 0, temp_sum_tp2 = 0;
    all_xist[0] = xist_rna;
    xro[b] = xist_rna;
    tro[b] = tsix_rna;
    xpo[b] = sum_xp;
    tpo[b] = sum_tp;
    if (sec_chr==1) {
        all_xist2[0] = xist_rna2;
        xro2[b] = xist_rna2;
        tro2[b] = tsix_rna2;
        xpo2[b] = sum_xp2;
        tpo2[b] = sum_tp2;
    }
    
    //Set initial condition for Xi: Tsix should be silenced from beginning
    if (tsix_silenced==1){sil_state_tsix=p[12]; sil_tsix=1;}
    //while (act_time<t_max)
    for (step=0; step<n_steps; step++)
        { 
        
	//If Tsix silencing has reached final state, Tsix is silenced //p[12]=parameter that determines # of intermediate states until silencing of tsix occurs
     	if (sil_state_tsix>=p[12]&xist_rna>=sil_threshold) {sil_tsix=1;} 
	//If silencing has not reached the final state but Xist rna level is below threshold system falls back into initial tsix silencing state
	else {
	if (xist_rna<sil_threshold&sil_state_tsix<p[12]) {sil_state_tsix=0; sil_tsix=0;}
	}  
	if (p[6]==0&tsix_gene.front()==1) 
        	{
            	xist_prom = 2;
		}
        if (sec_chr==1) {
     	if (sil_state_tsix2>=p[12]&xist_rna2>=sil_threshold) {sil_tsix2=1;} 
	else {
	if (xist_rna2<sil_threshold&sil_state_tsix2<p[12]) {sil_state_tsix2=0; sil_tsix2=0;}
	} 
	if (p[6]==0&tsix_gene2.front()==1) 
		{
                	xist_prom2 = 2;
		}
        }
        
//Gillespie reactions
        t_next = 0;
        
        xpol_ini=0, tpol_ini=0, xpol_ini2=0, tpol_ini2=0;  
        while (t_next<v) {
			//without XA influence => set k1 to the effective value it would have with one XA dosage
            rx[0] = p[0]*1*(xist_prom==1)*(tsix_gene.front()==0); //Xist initiation
            rx[1] = p[1]*(sil_tsix<1); //Tsix initiation
            rx[2] = (p[2]*(xist_prom==0)+p[7]*p[2]*(xist_prom==2))*(tsix_gene.front()==0); // turning on the Xist promoter: p[7]=percentage of p[2] with which transition OFF state 2 -> ON state 1 happens (p[2]=rate with which transition OFF state 0 -> ON state 1 happens)
            rx[3] = p[3]*xist_rna; //Xist degradation
            rx[4] = p[4]*tsix_rna; //Tsix degradation
            /////////NEW!!!!!
            rx[5] = p[9]*(xist_prom==1); //turning off the Xist promoter - basal rate
            
            //second chr
            //without XA influence => set k1 to the effective value it would have with one XA dosage
            rx[6] = (sec_chr==1)*p[0]*1*(xist_prom2==1)*(tsix_gene2.front()==0); //Xist initiation
            rx[7] = (sec_chr==1)*p[1]*(sil_tsix2<1);   //Tsix initiation
            rx[8] = (sec_chr==1)*(p[2]*(xist_prom2==0)+p[7]*p[2]*(xist_prom2==2))*(tsix_gene2.front()==0);// turning on the Xist promoter
            rx[9] = (sec_chr==1)*p[3]*xist_rna2; //Xist degradation
            rx[10] = (sec_chr==1)*p[4]*tsix_rna2; //Tsix degradation
            rx[11] = (sec_chr==1)*p[9]*(xist_prom2==1); //turning off the Xist promoter

	    //New!!! Transition between silencing states as Gillespie reaction
	    ////chromosome 1: now k_adv_sil= rate for transition between silencing states
	    rx[12] = k_adv_sil*(xist_rna>=sil_threshold)*(sil_state_tsix<p[12]); //Tsix silencing advance one state
	    //if system is in last silencing state of tsix but xist rna fall below the threshold the silencing state can go back to the initial state with rate k[14]
	    rx[13] = p[14]*(xist_rna<sil_threshold)*(sil_state_tsix>=p[12]);
	    ////second chromosome
	    rx[14] = (sec_chr==1)*k_adv_sil*(xist_rna2>=sil_threshold)*(sil_state_tsix2<p[12]);
	    rx[15] = (sec_chr==1)*p[14]*(xist_rna2<sil_threshold)*(sil_state_tsix2>=p[12]);
	    
            
            
            sum_rx = 0, sel_rx = 100;
            for (q=0; q<rx.size(); q++) {
                sum_rx = sum_rx + rx[q];
                cum_rx[q] = sum_rx;
            }
         
            if (sum_rx>0) {
                t_next = t_next + (-log(dist_real(gen)))/(sum_rx);
            } else {t_next=1;}
            
            if (t_next<v){
                // decide which reaction occurs
                r = dist_real(gen);
                for (q=0; q<rx.size(); q++) {
                    if (r<(cum_rx[q]/sum_rx)){
                        sel_rx = q;
                        break;
                    }
                }
                //execute reaction
                
                switch (sel_rx){
                    case 0:
                        xpol_ini = 1; break;
                    case 1:
                        tpol_ini = 1; break;
                    case 2:
                        xist_prom = 1; break;
                    case 3:
                        xist_rna--; break;
                    case 4:
                        tsix_rna--; break;
                    case 5:
                        xist_prom = 0; break;
                                       
                    //second chr                        
                      
                   case 6:
                        xpol_ini2 = 1; break;
                    case 7:
                        tpol_ini2 = 1; break;
                    case 8:
                        xist_prom2 = 1; break;
                    case 9:
                        xist_rna2--; break;
                    case 10:
                        tsix_rna2--; break;
                    case 11:
                        xist_prom2 = 0; break;
		    case 12:
			if (sil_state_tsix<p[12]) {sil_state_tsix++;}
			break;
		    case 13:
			sil_state_tsix=0; sil_tsix=0; break;
		    case 14:
			if (sil_state_tsix2<p[12]) {sil_state_tsix2++;}
			break;
		    case 15:
			sil_state_tsix2=0; sil_tsix2=0; break;
                    default:
                        break;
                }
            }
        }

        //elongation and collisions       
        elongation (overlap_X_T, dimx_xp, xist_gene, tsix_gene, sum_xp, sum_tp, 
            xpol_ini, tpol_ini, xist_rna, tsix_rna, 
            all_xist, temp_sum_xr, temp_sum_tr, temp_sum_xp, temp_sum_tp, dist, step, p[8]);
        
               
        //second chromosome
        if (sec_chr==1) {
            elongation (overlap_X_T, dimx_xp, xist_gene2, tsix_gene2, sum_xp2, sum_tp2, 
            xpol_ini2, tpol_ini2, xist_rna2, tsix_rna2, 
            all_xist2, temp_sum_xr2, temp_sum_tr2, temp_sum_xp2, temp_sum_tp2, dist, step, p[8]);
            
        } 
        
        //write out smoothed variable
        m++;          
        if (m==n_small_steps) {
            b++;
            xro[b] = temp_sum_xr/n_small_steps;
            tro[b] = temp_sum_tr/n_small_steps;
            xpo[b] = temp_sum_xp/n_small_steps;
            tpo[b] = temp_sum_tp/n_small_steps;

            time[b] = step*v;
            temp_sum_xr = 0;
            temp_sum_tr = 0;
            temp_sum_xp = 0;
            temp_sum_tp = 0;

            if (sec_chr==1) {
                xro2[b] = temp_sum_xr2/n_small_steps;
                tro2[b] = temp_sum_tr2/n_small_steps;
                xpo2[b] = temp_sum_xp2/n_small_steps;
                tpo2[b] = temp_sum_tp2/n_small_steps;
                temp_sum_xr2 = 0;
                temp_sum_tr2 = 0;
                temp_sum_xp2 = 0;
                temp_sum_tp2 = 0;
            }
            m = 0;
        }
        
        }
    if (m>0) {
            b++;
            xro[b] = temp_sum_xr/m;
            tro[b] = temp_sum_tr/m;
            xpo[b] = temp_sum_xp/m;
            tpo[b] = temp_sum_tp/m;

            time[b] = step*v;
            if (sec_chr==1) {
                xro2[b] = temp_sum_xr2/m;
                tro2[b] = temp_sum_tr2/m;
                xpo2[b] = temp_sum_xp2/m;
                tpo2[b] = temp_sum_tp2/m;
            }
    }
   	}
