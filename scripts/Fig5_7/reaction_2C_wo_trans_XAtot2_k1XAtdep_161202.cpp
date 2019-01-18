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

int collisions (int &dimx_xp, deque<int> &xist_gene, deque<int> &tsix_gene, 
        int &sum_xp, int &sum_tp, boost::random::uniform_int_distribution<> &dist)
{
    int q;
    for (q=0; q<dimx_xp; q++)
            {
            if ((xist_gene[q] + tsix_gene[q])>1) {
                if (dist(gen)==0) {
                    xist_gene[q]=0;
                    sum_xp--;
                } else {
                    tsix_gene[q]=0;
                    sum_tp--;
                }
            }
        }
}

int elongation (int &dimx_xp, deque<int> &xist_gene, deque<int> &tsix_gene, int &sum_xp, int &sum_tp, 
        bool &xpol_ini, bool &tpol_ini, int &xist_rna, int &tsix_rna, 
        vector<int> &all_xist, int &temp_sum_xr, int &temp_sum_tr, int &temp_sum_xp, int &temp_sum_tp,
        boost::random::uniform_int_distribution<> &dist, int &step, double &p8)
{
    xist_gene.push_front(xpol_ini);
    sum_xp = sum_xp - xist_gene.back() + xist_gene.front();
    xist_rna = xist_rna + (int)xist_gene.back();
    xist_gene.pop_back();
    if (p8>0) {
        collisions (dimx_xp, xist_gene, tsix_gene, sum_xp, sum_tp, dist);
    }   
    tsix_gene.push_back(tpol_ini);
    sum_tp = sum_tp - tsix_gene.front() + tsix_gene.back();
    tsix_rna = tsix_rna + (int)tsix_gene.front();
    tsix_gene.pop_front();
    if (p8>0) {
        collisions (dimx_xp, xist_gene, tsix_gene, sum_xp, sum_tp, dist);
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
            *xist_pol_out2, *tsix_pol_out2, *xist_rna_out2, *tsix_rna_out2, *test, 
            *switch_off_out, *stable_switch_on_out; //, *xist_prom_out;
	const mwSize *dims;
	double *xp, *tp, *xr, *tr, *xpr, *xp2, *tp2, *xr2, *tr2, *xpr2, *p, *cp, 
            t, sil_threshold, t_diff, t_before;
    double *xpo, *tpo, *xro, *tro, *xpro, *xpo2, *tpo2, *xro2, *tro2, *xpro2, 
            *time, v = 1.0/1440.0, out_step, *to, *switch_off, *stable_sw_on, k_adv_sil;
	int dimx_xp, dimx_tp, xist_rna, tsix_rna, xist_prom, xist_rna2, tsix_rna2, xist_prom2;
	int n_steps,n_steps_before, n_out_steps, n_small_steps, sec_chr, step, sel_rx;
	
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
    t_before =  cp[1];
    t_diff =  cp[2];
    sil_threshold = cp[3];
    out_step = cp[4];
    sec_chr = cp[5];
    k_adv_sil = cp[6];
    
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
    
    n_out_steps = 1+ceil((t_diff+t_before)/out_step);
    n_small_steps = round(out_step/v);
    n_steps = (n_out_steps-1)*n_small_steps;
    n_steps_before = (ceil(t_before/out_step))*n_small_steps;
	
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
    switch_off_out = plhs[10] = mxCreateDoubleMatrix(6, 1, mxREAL);
    stable_switch_on_out = plhs[11] = mxCreateDoubleMatrix(1, 1, mxREAL);
        
    xpo2 = mxGetPr(xist_pol_out2);
    tpo2= mxGetPr(tsix_pol_out2);
    xro2 = mxGetPr(xist_rna_out2);
    tro2 = mxGetPr(tsix_rna_out2);
    to = mxGetPr(test);
    switch_off = mxGetPr(switch_off_out);
    stable_sw_on = mxGetPr(stable_switch_on_out);
    
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
    vector<double> rx(18), cum_rx(18);
    vector<int> all_xist(n_steps), all_xist2(n_steps), swo(6,0), counter1(6,0), counter2(6,0); // 6 integers with value 0
    
    int sil_xa = 0, sil_xa2 = 0, xa_sil_steps = round(p[11]/v), 
            sil_tsix = 0, sil_tsix2 = 0, 
	    sil_state_xa = 0, sil_state_xa2 = 0, sil_state_tsix = 0, sil_state_tsix2 = 0, stst = 0; 
    double xa = 1, xa2 = 0;
    if (sec_chr==1) {xa2 = 1;};
    double temp_sum_xa = xa + xa2;
    int m = 0, b = 0, temp_sum_xr = 0, temp_sum_tr = 0, temp_sum_xp = 0, temp_sum_tp = 0,
            temp_sum_xr2 = 0, temp_sum_tr2 = 0, temp_sum_xp2 = 0, temp_sum_tp2 = 0;
    all_xist[0] = xist_rna;
    xro[b] = xist_rna;
    tro[b] = tsix_rna;
    xpo[b] = sum_xp;
    tpo[b] = sum_tp;
    to[b] = temp_sum_xa;
    if (sec_chr==1) {
        all_xist2[0] = xist_rna2;
        xro2[b] = xist_rna2;
        tro2[b] = tsix_rna2;
        xpo2[b] = sum_xp2;
        tpo2[b] = sum_tp2;
    }
    //while (act_time<t_max)
    for (step=0; step<n_steps; step++)
        { 
			
		// Identify first timepoint at which stable state with 1 Xa and 1 Xi is reached
		if (xist_rna>=sil_threshold & xa ==0 & sil_tsix ==1 & xist_rna2<sil_threshold & xa2==1 & sil_tsix2==0 & stst==0) {stable_sw_on[0] = step*v; stst=1;}
		else { if (xist_rna2>=sil_threshold & xa2 ==0 & sil_tsix2 ==1 & xist_rna<sil_threshold & xa==1 & sil_tsix==0 & stst==0) {stable_sw_on[0] = step*v; stst=1;}
	}
	
		//counter how long Xist is already above threshold for switch_off 1: 2*XA, Tsix ON
        if (xist_rna>=sil_threshold & xa+xa2==2 & sil_tsix==0) {counter1[0] = counter1[0]+1;}
        else {counter1[0] = 0;}
        // X2
        if (xist_rna2>=sil_threshold & xa+xa2==2 & sil_tsix2==0) {counter2[0] = counter2[0]+1;}
        else {counter2[0] = 0;}
        //counter switch_off 2: 1*XA, Tsix ON
        if (xist_rna>=sil_threshold & xa+xa2==1 & sil_tsix==0) {counter1[1] = counter1[1]+1;}
        else {counter1[1] = 0;}
        if (xist_rna2>=sil_threshold & xa+xa2==1 & sil_tsix2==0) {counter2[1] = counter2[1]+1;}
        else {counter2[1] = 0;}
        //counter switch_off 3: 2*XA, Tsix OFF
        if (xist_rna>=sil_threshold & xa+xa2==2 & sil_tsix==1) {counter1[2] = counter1[2]+1;}
        else {counter1[2] = 0;}
        if (xist_rna2>=sil_threshold & xa+xa2==2 & sil_tsix2==1) {counter2[2] = counter2[2]+1;}
        else {counter2[2] = 0;}
        //counter switch_off 4: 1*XA, Tsix OFF
        if (xist_rna>=sil_threshold & xa+xa2==1 & sil_tsix==1) {counter1[3] = counter1[3]+1;}
        else {counter1[3] = 0;}
        if (xist_rna2>=sil_threshold & xa+xa2==1 & sil_tsix2==1) {counter2[3] = counter2[3]+1;}
        else {counter2[3] = 0;}
        //counter switch_off 5: 0*XA, Tsix ON
        if (xist_rna>=sil_threshold & xa+xa2==0 & sil_tsix==0) {counter1[4] = counter1[4]+1;}
        else {counter1[4] = 0;}
        if (xist_rna2>=sil_threshold & xa+xa2==0 & sil_tsix2==0) {counter2[4] = counter2[4]+1;}
        else {counter2[4] = 0;}
        //counter switch_off 6: 0*XA, Tsix OFF
        if (xist_rna>=sil_threshold & xa+xa2==0 & sil_tsix==1) {counter1[5] = counter1[5]+1;}
        else {counter1[5] = 0;}
        if (xist_rna2>=sil_threshold & xa+xa2==0 & sil_tsix2==1) {counter2[5] = counter2[5]+1;}
        else {counter2[5] = 0;}
        
        
        
        
        
	//If Tsix silencing has reached final state, Tsix is silenced //p[12]=parameter that determines # of intermediate states until silencing of tsix occurs
     	if (sil_state_tsix>=p[12]&xist_rna>=sil_threshold) {sil_tsix=1;} 
	//If silencing has not reached the final state but Xist rna level is below threshold system falls back into initial tsix silencing state
	else {
	if (xist_rna<sil_threshold&sil_state_tsix<p[12]) {sil_state_tsix=0; sil_tsix=0;}
	}
	//If XA silencing has reached final state, XA is silenced //p[11]=parameter that determines # of intermediate states until silencing of xa occurs
	if (sil_state_xa>=p[11]&xist_rna>=sil_threshold) {xa=0;} 
	//If silencing has not reached the final state but Xist rna level is below threshold system falls back into initial xa silencing state
	else {
	if (xist_rna<sil_threshold&sil_state_xa<p[11]) {sil_state_xa=0; xa=1;}	
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
	if (sil_state_xa2>=p[11]&xist_rna2>=sil_threshold) {xa2=0;} 
	else {
	if (xist_rna2<sil_threshold&sil_state_xa2<p[11]) {sil_state_xa2=0; xa2=1;}	
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
			
            //before induction of diff k1 10 times reduced
            if (step<n_steps_before) {rx[0] = p[0]*0.1*(xa+xa2+0.001)*(xist_prom==1)*(tsix_gene.front()==0);}
            else {rx[0] = p[0]*(xa+xa2+0.001)*(xist_prom==1)*(tsix_gene.front()==0);} //Xist initiation
            rx[1] = p[1]*(sil_tsix<1); //Tsix initiation
            rx[2] = p[7]*(xist_prom==2)*(tsix_gene.front()==0); // turning on the Xist promoter: p[7]=percentage of p[2] with which transition OFF state 2 -> ON state 1 happens (p[2]=rate with which transition OFF state 0 -> ON state 1 happens)
            rx[3] = p[3]*xist_rna; //Xist degradation
            rx[4] = p[4]*tsix_rna; //Tsix degradation
            /////////NEW!!!!!
            //rx[5] = p[9]*(xist_prom==1); //turning off the Xist promoter - basal rate
            
            //second chr
            if (step<n_steps_before) {rx[5] = (sec_chr==1)*p[0]*0.1*(xa+xa2+0.001)*(xist_prom2==1)*(tsix_gene2.front()==0);}
            else {rx[5] = (sec_chr==1)*p[0]*(xa+xa2+0.001)*(xist_prom2==1)*(tsix_gene2.front()==0);} //Xist initiation
            rx[6] = (sec_chr==1)*p[1]*(sil_tsix2<1);   //Tsix initiation
            rx[7] = (sec_chr==1)*(p[7]*(xist_prom2==2))*(tsix_gene2.front()==0); // turning on the Xist promoter
            rx[8] = (sec_chr==1)*p[3]*xist_rna2; //Xist degradation
            rx[9] = (sec_chr==1)*p[4]*tsix_rna2; //Tsix degradation
            //rx[11] = (sec_chr==1)*p[9]*(xist_prom2==1); //turning off the Xist promoter

	    //New!!! Transition between silencing states as Gillespie reaction
	    ////chromosome 1: now k_adv_sil= rate for transition between silencing states
	    rx[10] = k_adv_sil*(xist_rna>=sil_threshold)*(sil_state_xa<p[11]|sil_state_tsix<p[12]); //XA and Tsix silencing advance one state
	    ////if system is in last silencing state of the xa but xist rna fall below the threshold the silencing state can go back to the initial state with rate k[13]
	    rx[11] = p[13]*(xist_rna<sil_threshold)*(sil_state_xa>=p[11]);
	    //if system is in last silencing state of tsix but xist rna fall below the threshold the silencing state can go back to the initial state with rate k[14]
	    rx[12] = p[14]*(xist_rna<sil_threshold)*(sil_state_tsix>=p[12]);
	    ////second chromosome
	    rx[13] = (sec_chr==1)*k_adv_sil*(xist_rna2>=sil_threshold)*(sil_state_xa2<p[11]|sil_state_tsix2<p[12]);
	    rx[14] = (sec_chr==1)*p[13]*(xist_rna2<sil_threshold)*(sil_state_xa2>=p[11]);
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
                        xist_rna--;
                        //switch_off1:first timepoint at which Xist RNA falls back below threshold with Tsix on in cis at double XA dose
                        if (xist_rna==(sil_threshold-1) & swo[0]==0 & xa+xa2==2 & sil_tsix==0) {swo[0]=1; switch_off[0] = counter1[0]*v;} 
                        //switch_off2:first timepoint at which Xist RNA falls back below threshold with Tsix on in cis at single XA dose
                        if (xist_rna==(sil_threshold-1) & swo[1]==0 & xa+xa2==1 & sil_tsix==0) {swo[1]=1; switch_off[1] = counter1[1]*v;} 
                        //switch_off3:first timepoint at which Xist RNA falls back below threshold with Tsix off in cis at double XA dose
                        if (xist_rna==(sil_threshold-1) & swo[2]==0 & xa+xa2==2 & sil_tsix==1) {swo[2]=1; switch_off[2] = counter1[2]*v;} 
                        //switch off 4
                        if (xist_rna==(sil_threshold-1) & swo[3]==0 & xa+xa2==1 & sil_tsix==1) {swo[3]=1; switch_off[3] = counter1[3]*v;}
                        //switch off 5
                        if (xist_rna==(sil_threshold-1) & swo[4]==0 & xa+xa2==0 & sil_tsix==0) {swo[4]=1; switch_off[4] = counter1[4]*v;}
                        //switch off 6
                        if (xist_rna==(sil_threshold-1) & swo[5]==0 & xa+xa2==0 & sil_tsix==1) {swo[5]=1; switch_off[5] = counter1[5]*v;}
                         break;
                    case 4:
                        tsix_rna--; break;
                    //case 5:
                        //xist_prom = 0; break;
                                       
                    //second chr                        
                      
                   case 5:
                        xpol_ini2 = 1; break;
                    case 6:
                        tpol_ini2 = 1; break;
                    case 7:
                        xist_prom2 = 1; break;
                    case 8:
                        xist_rna2--; 
                        if (xist_rna2==(sil_threshold-1) & swo[0]==0 & xa+xa2==2 & sil_tsix2==0) {swo[0]=1; switch_off[0] = counter2[0]*v;} 
                        if (xist_rna2==(sil_threshold-1) & swo[1]==0 & xa+xa2==1 & sil_tsix2==0) {swo[1]=1; switch_off[1] = counter2[1]*v;} 
                        if (xist_rna2==(sil_threshold-1) & swo[2]==0 & xa+xa2==2 & sil_tsix2==1) {swo[2]=1; switch_off[2] = counter2[2]*v;} 
                        if (xist_rna2==(sil_threshold-1) & swo[3]==0 & xa+xa2==1 & sil_tsix2==1) {swo[3]=1; switch_off[3] = counter2[3]*v;}
                        if (xist_rna2==(sil_threshold-1) & swo[4]==0 & xa+xa2==0 & sil_tsix2==0) {swo[4]=1; switch_off[4] = counter2[4]*v;}
                        if (xist_rna2==(sil_threshold-1) & swo[5]==0 & xa+xa2==0 & sil_tsix2==1) {swo[5]=1; switch_off[5] = counter2[5]*v;}
                        break;
                    case 9:
                        tsix_rna2--; break;
                    //case 11:
                        //xist_prom2 = 0; break;
		    case 10:
			if (sil_state_xa<p[11]) {sil_state_xa++;}
			if (sil_state_tsix<p[12]) {sil_state_tsix++;}
			break;
		    case 11:
			sil_state_xa=0; xa=1; break;
		    case 12:
			sil_state_tsix=0; sil_tsix=0; break;
		    case 13:
			if (sil_state_xa2<p[11]) {sil_state_xa2++;}
			if (sil_state_tsix2<p[12]) {sil_state_tsix2++;}
			break;
		    case 14:
			sil_state_xa2=0; xa2=1; break;
		    case 15:
			sil_state_tsix2=0; sil_tsix2=0; break;
                    default:
                        break;
                }
            }
        }
        temp_sum_xa = temp_sum_xa + xa + xa2;
        //elongation and collisions       
        elongation (dimx_xp, xist_gene, tsix_gene, sum_xp, sum_tp, 
            xpol_ini, tpol_ini, xist_rna, tsix_rna, 
            all_xist, temp_sum_xr, temp_sum_tr, temp_sum_xp, temp_sum_tp, dist, step, p[8]);
        
               
        //second chromosome
        if (sec_chr==1) {
            elongation (dimx_xp, xist_gene2, tsix_gene2, sum_xp2, sum_tp2, 
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
            to[b] = temp_sum_xa/n_small_steps;
            time[b] = step*v;
            temp_sum_xr = 0;
            temp_sum_tr = 0;
            temp_sum_xp = 0;
            temp_sum_tp = 0;
            temp_sum_xa = 0;
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
            to[b] = temp_sum_xa/n_small_steps;
            time[b] = step*v;
            if (sec_chr==1) {
                xro2[b] = temp_sum_xr2/m;
                tro2[b] = temp_sum_tr2/m;
                xpo2[b] = temp_sum_xp2/m;
                tpo2[b] = temp_sum_tp2/m;
            }
    }
   	}
