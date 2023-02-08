//
//  main.cpp
//  NonSpatial-ABM
//
//  Created by Andrew Yates on 2020-11-02.
//
// Simple ABM of cell dynamics:
// -Division and death rates can depend on...
//      poolsize, time, cell age, time since last division, number of divisions
// - Tracks Ki67 and GFP expression (residual from Rag)
//   Latter decays exponentially with time, and halves when cell divides
//   Ki67 expression also decays exponentially/deterministically... cells can also be restricted
//   from dividing again too soon
//   We ALSO track a binary Ki67 hi/lo toggle which flips at rate beta 1/3.5d - analogous to ODE model of loss of Ki67
//   (So, we have Ki67_det and Ki67_stoch)
// - Also tracks clonal structure (not needed here)
// - Inputs new clones from the thymus.
//    All of that is handled in the function "update" which takes one (fixed, small) timestep and handles all the death
//    division and influx.
//
//    The basic principle is that you pre-allocate a huge fixed array of objects of class "cell" (cellstore).
//    And two lists that contain pointers to all the slots in this array that contain live cells
//    (cell_location_list_1, and cell_location_list_2); and another list of pointers to empty slots (spacelist).
//    Each timestep you update you take the current list of live cells, record all the deaths and divisions and
//    immigrations, and build a new  list of pointers to live cells from that, using and updating the list of free spaces
//    as you go. Then, toggle the lists so you update the old one on the next timestep.
//    There are prob more efficient ways to do it but this seemed reasonable. My macbook can handle about 100,000 cells.

#include <iostream>
#include <random>
#include <string>
//#include <iomanip>
#include <fstream>
//#include <sstream>

//#include <cstdlib>
#include "custom_functions.cpp"

using namespace std;

# define FILEPATH "/Users/sanketrane/Dropbox/Conference_Room_Files/Sanket/T_cell_ABM/Tregs/"
# define PARAMETER_FILENAME "parameters_ASM8rhovar_Ki67gap_90.txt"

# define OUTPUT_FILENAME "allout_Incumbent.csv"
# define GFP_DIST_FILENAME "GFP_dist_ASM8rhovar.csv"
# define GFP_DUMP_FILENAME "GFP_dump_ASM8rhovar.csv"

//# define GFP_KI67_FREQS_FILENAME "GFPKi67_freq.csv"

# define SCALE_OUTPUT 0.005f // Scale N0 and influx by this number to get actual simulated cell numbers
# define CLONEUNIVERSESIZE 100 // number of possible clonotypes/TCR sequences
# define TSTEP 0.04f // time step. needs to be small, bcs event probabilities are estimated as (process rate)*TSTEP
# define STORAGELISTLENGTH 100000 // max number of cells
# define LOG_KPOS_THRESHOLD -1.f // Define threshold of Ki67 positivity. Set to exp(-1) consistent with Sanket's ASM
// max value of deterministic Ki67 is 1 - > log equals zero. So Ki67+ cells lie between -1 and 0
# define KI67_LOSS_RATE 0.2857143f  // 1/3.5d
# define EXP_MINUS_1 0.3678794412f // handy constant
# define KPOS_THRESHOLD 0.3678794412f // handy constant
# define LOG2 0.69314718f // another handy constant

# define sep "  ,  " // column separator for output files (e.g. comma or tab)

// CD4 neutral model;  DIV  0.00456f  LOSS 0.07f, N0 = 907,000 N0, THYMIC_EXPORT_RATE 0.54
// CD8 rhovar ASM; N0 = 70983.6881, delta = 0.04533878, rho=0.061905344 * exp(0.164375558 a)
// ---------------------------------------------------------------------
// Cell object definition
class cell
{
private:
    int cloneID; // 1...CLONEUNIVERSESIZE
    float ki67_intens_norm; // normalized Ki67 expression. 
    bool location_thymus; // whether its a thymic or peripheral Treg
    int ndivs; // number of divisions since differentiation into Treg phenotype in the thymus
    float last_division_time; // time at which this cell last divided
    float export_time;// time at which its ancestor differentiatiated into Treg

public:
    cell(); // default constructor fn.
    cell(int p1, float p2, bool p3, bool p4, int p5, float p6, float p7)); // general constructor fn.
    void set_cloneID(int x){cloneID=x;}
    void set_ki67_intens_norm(float x){ki67_intens_norm=x;}
    void set_donor_derived(bool x){donor_derived=x;}
    void set_location_thymus(bool x){location_thymus=x;}
    void set_ndivs(int i){ndivs=i;}
    void set_last_division_time (float x) {last_division_time =x;}
    void set_export_time (float x) {export_time=x;}

    int get_cloneID() {return(cloneID);}
    float ki67_intens_norm() {return(logKi67_det);}
    bool get_donor_derived(){return (donor_derived);}
    bool get_location_thymus(){return (location_thymus);}
    int get_ndivs() {return(ndivs);}
    float get_last_division_time() {return(last_division_time);}
    float get_export_time() {return(export_time);}
};

cell::cell(){ //default constructor
    cloneID=0;
    ki67_intens_norm=runif<double>(0, KPOS_THRESHOLD);    // for a given cell ki67 exoression is chosen from a uniform distribution between 0 and exp(-1)
    donor_derived=false;
    location_thymus=true;
    ndivs=0;
    last_division_time=-99999.; // last division time not defined for undivided cells
    export_time=0.;
}

cell::cell(int p1, float p2, bool p3, bool p4, int p5, float p6, float p7){ //general constructor
    cloneID=p1;
    ki67_intens_norm=p2;
    donor_derived=p3;
    location_thymus=p4;
    ndivs=p5;
    last_division_time=p6;
    export_time=p7;
}

// parameters (all floats)
// 0: 4 or 8 (CD4 or CD8)
// 1: N0 (to be scaled by 0.005 to get actual cell numbers)
// 2: T0 (start time for simulations)
// 3: TMAX
// 4: g0  pop density at age zero at time T0. # thymic export rate at time t=1 = C* SP(1)  =  N(t=1)g(0)
// 5: N_densitydependence (scale for density dep div or loss - in "true" physiological numbers
// 6: delta0 (base death rate for cells of age zero)
// 7: r_delta (death rate goes with cell age as exp(-r_delta * age))
// 8: rho_0 (base death rate for cells of age zero)
// 9: r_rho (div rate goes with cell age as exp(-r_rho * age))

float lossrate(int poolsize, float current_time, float cell_age, float time_since_last_division, int ndivs, float params[]){
    float delta0=params[6];
    float r_delta= params[7];
    float host_age_at_export = current_time - cell_age;
    float N_densitydependence = params[5]; // work with "true" physiological numbers
    if (host_age_at_export < 0. ) host_age_at_export = 0.;
    // insert function etc. here to calculate host age dependent rho_0
    //return (delta0 * exp(-1.* cell_age * r_delta));
    return delta0;
}

float divrate(int poolsize, float current_time, float cell_age, float time_since_last_division, int ndivs, float params[]){
    float rho0=params[8];
    float r_rho= params[9];
    float host_age_at_export = current_time - cell_age;
    float N_densitydependence = params[5]; // work with "true" physiological numbers
    if (host_age_at_export < 0. ) host_age_at_export = 0.;
    // if (host_age_at_export <30.) rho0=rho0*0.02; else rho0=rho0*1.2;
    // insert function etc. here to calculate host age dependent rho_0
    //return (rho0 * exp(-1.* cell_age * r_rho));
    return rho0;
}

float sp_numbers(float t, float params[]){ // mature SP thymocyte numbers - spline defs
    // these are true numbers - need to scale down for simulation (done elsewhere)
    float  dpt0, T0=params[2], sp_numbers, theta0, nu;

    // define the splines for SP4 and SP8
    theta0  = 6.4;  nu = 0.0024;}
    dpt0 = t - T0;     // days post t0
    if (dpt0<0.) dpt0=0.;
    sp_numbers = pow(10, theta0) * exp(-1 * nu * dpt0); 
    return sp_numbers;
}

void update(cell *fromlist[], cell *tolist[], cell cellstore[], int clonerecord[], int *poolsize, cell *spacelist[],
            int *spacelistlength, float current_time, float params[], float THYMIC_EXPORT_RATE_CONSTANT)
{
    int i, j=0, cells_in, cloneID, ndivs, updated_poolsize=*poolsize; // j tracks number of live cells in this cycle
    float matureSP_Ki67pos_frac, runif, age, donor_derived, ki67_intens_norm, location_thymus, last_division_time, time_since_last_division, export_time, p_div, p_loss;
    cell *this_cell;
    cell *new_cell_location;
    bool dead, divide;

    // stock random number libraries
    std::random_device rd;
    std::mt19937 generator(rd()); // use the mersenne twister as the underlying RN generator
    // or ... fixed seed
    // std::mt19937 generator(0);

    // uniform distribution between 0 and 1; call with uniformrandom(generator)
    std::uniform_real_distribution<> uniformrandom(0, 1);

    // call these with Normalized_ki67_dist_SP(generator)
    std::uniform_real_distribution<> Ki67_dist_SP(0, exp(-1)); 

    for(i=0;i<*poolsize;i++){ // loop through all cells alive at the start of this timestep
        // fromlist is a list of pointers to live cells... it should have no gaps in it.
        this_cell=*(fromlist+i);

        // get this cell's attributes
        cloneID=(*this_cell).get_cloneID(); //similar to python syntax?
        ki67_intens_norm= (*this_cell).get_ki67_intens_norm();
        donor_derived= (*this_cell).get_donor_derived();
        location_thymus= (*this_cell).get_location_thymus();
        ndivs=(*this_cell).get_ndivs();
        last_division_time=(*this_cell).get_last_division_time();
        export_time=(*this_cell).get_export_time();
        age=current_time-export_time;
        time_since_last_division = current_time - last_division_time;

        // update loggfp and logki67 values - linear drop on log scale
        logKi67_det  = logKi67_det - TSTEP*KI67_LOSS_RATE; // Ki67 drops with time
        (*this_cell).set_logKi67_det(logKi67_det);

        // Now decide what the cell does this timestep
        // calculate its probs of division and loss ... = rate * timestep (if timestep small!)
        p_div  = divrate(*poolsize, current_time, age, time_since_last_division, ndivs, params)*TSTEP;
        p_loss = lossrate(*poolsize, current_time, age, time_since_last_division, ndivs, params)*TSTEP;
        // Carve the unit interval into (0, p_loss, (1-p_div), 1) and see where random # lands.
        // p_loss and p_div should both be small (short timestep) so shouldn't get into trouble here
        dead=false;
        divide=false;
        runif=uniformrandom(generator);
        if(runif<p_loss) dead=true;
        // hack - don't let cells re-divide before 0.2
        //if(runif>(1-p_div) && time_since_last_division>0.2) divide=true;
        if(runif>(1-p_div)) divide=true;

        if(dead){
            divide=false; // just in case of shenanigans with the probabilities above
            updated_poolsize--; // one less cell
            (*(clonerecord+cloneID))--; // and one less of this clone
            if(updated_poolsize<=0) {
                cout << "Poolsize equals zero" << endl;
                //exit(1);
            }
            // Add this dead cell's space to the free list
            (*spacelistlength)++;
            // add the address of this (dead) cell to the end of the spacelist
            spacelist[*spacelistlength - 1] = this_cell; // indexing of arrays starts at zero
        }

        if(divide){
            updated_poolsize++;
            (*(clonerecord+cloneID))++;
            // first update the parent cell
            (*this_cell).set_ndivs(ndivs+1);
            (*this_cell).set_last_division_time(current_time);
            (*this_cell).set_logGFP(logGFP - logGFP_halfstep); // GFP halves
            (*this_cell).set_logKi67_det(0.); // max Ki67 expression = 1
            (*this_cell).set_Ki67_stoch(true); // make it Ki67+ for the ODE contingent
            tolist[j]=this_cell; // put this cell in the updated list of live cells
            // now make a new cell in the first free slot, and set up a pointer
            // to it so we can find it next timestep
            if(*spacelistlength==0){ // there are no gaps in the array cellstore, so just add onto the end
                new_cell_location=cellstore + updated_poolsize -1;
            }
            else
            {
                new_cell_location=spacelist[*spacelistlength-1];// use the space on the top of the stack
                (*spacelistlength)--; //remove this space from the stack
            }
            tolist[j+1]=new_cell_location; // pointer to where this new cell will live
            (*new_cell_location).set_cloneID(cloneID);
            (*new_cell_location).set_ndivs(ndivs+1);
            (*new_cell_location).set_export_time(export_time); // time that ancestor left the thymus
            (*new_cell_location).set_last_division_time(current_time);
            (*new_cell_location).set_logGFP(logGFP - logGFP_halfstep); // same as sibling
            (*new_cell_location).set_logKi67_det(0.); // same as sibling
            (*new_cell_location).set_Ki67_stoch(true); // same as sibling
            j=j+2; // count 2 cells here
        }
        if(!dead && !divide) { //  survives and does nothing
            // if it's Ki67_stoch positive, does it become Ki67_stoch neg?
            if((*this_cell).get_Ki67_stoch()){
                if(uniformrandom(generator)<TSTEP*KI67_LOSS_RATE)  (*this_cell).set_Ki67_stoch(false);
            }
            //  add it to the list
            *(tolist+j)=*(fromlist+i);
            j++; // count one more cell
        }
    }
    //////////////////////////////////////////
    // now add new cells from the thymus (integer)
    // we work in rescaled numbers (i.e. scaled down from "true" numbers)
    //////////////////////////////////////////

    cells_in = lrint( THYMIC_EXPORT_RATE_CONSTANT * sp_numbers(current_time, params) * SCALE_OUTPUT * TSTEP);

    for(i=0;i<cells_in;i++){
        updated_poolsize++;
        cloneID=rand() % CLONEUNIVERSESIZE; // choose random clonotype
        if(*spacelistlength==0){
            // there are no gaps in the array cellstore, so just add onto the end
            new_cell_location=cellstore + updated_poolsize -1;
        }
        else
        {
            new_cell_location=spacelist[(*spacelistlength) -1];// space on the top of the stack
            (*spacelistlength)--; //remove this space from the stack
        }
        tolist[j]=new_cell_location;
        j++;
        (*(clonerecord+cloneID))++; // record a new cell of this clone

        //make random number to determine whether new cell is Ki67 high or low,
        // depending on Ki67 expression in SP cells
        runif=uniformrandom(generator);
        matureSP_Ki67pos_frac = sp_ki67(current_time, params);
        if(runif < matureSP_Ki67pos_frac)  { // It's Ki67+;
            // if(runif<0.5) { // It's Ki67+;
            // now choose a value uniform on log scale between 0 and -0.4 (Ki67+ threshold is -1)
            // Ki67 on absolute expression scale is in (0,1). So (-infinity, 0) on log scale
            // Sanket found that Ki67+ cells from the thymus need to have absolute expression in the range
            // 0.9 - 1 (otherwise they become Ki67- too quickly)
            // ln(0.9) = -0.10536
            (*new_cell_location).set_logKi67_det(-0.10536 * uniformrandom(generator));
            (*new_cell_location).set_Ki67_stoch(true);
            // set GFP expression from distn of Ki67+ Mature SPs
            (*new_cell_location).set_logGFP(LOG_GFPdist_Ki67HI_RTE(generator));
        }
        else // it's Ki67 negative and value doesn't really matter unless you're trying to model the
            // whole Ki67 distribution (and we're not - we just want to model Ki67+ cells well)
        {
            (*new_cell_location).set_logKi67_det(-2.);
            (*new_cell_location).set_Ki67_stoch(false);
            // set GFP expression from distn of Ki67- Mature SPs
            (*new_cell_location).set_logGFP(LOG_GFPdist_Ki67LO_RTE(generator));
        }
        (*new_cell_location).set_cloneID(cloneID);
        (*new_cell_location).set_ndivs(0);
        (*new_cell_location).set_export_time(current_time);
        (*new_cell_location).set_last_division_time(-999999.);
    }
    *poolsize=updated_poolsize;
}

void gfp_histogram(cell *loclist[], int poolsize,  float current_time, float params[]){
    // scan through the cells and make a histogram of GFP expression (bins on log scale)
    int i,j,ndivs,  khi_numbers=0, klo_numbers=0;
    int binned_khi[NUM_GFP_BINS], binned_klo[NUM_GFP_BINS];
    float ticks [NUM_GFP_BINS];
    float logGFP;
    cell *this_cell;
    bool Ki67hi, first_time = false;
    string Ki67_status;
    ifstream GFPfile_test;
    ofstream GFP_dist_file, GFP_dump_file;

    // Set up file file to append this info to
    // first check if it exists ('main' deletes it every time it runs)
    // if it doesn't exist, write the header
    // try to open file to read //

    const string filepath =  FILEPATH;
    const string GFP_dist_filename = GFP_DIST_FILENAME;
    const string GFP_dump_filename = GFP_DUMP_FILENAME;
    const string GFP_dist_fullpath = filepath + GFP_dist_filename;
    const string GFP_dump_fullpath = filepath + GFP_dump_filename;
    GFPfile_test.open(GFP_dist_fullpath);
    if(!GFPfile_test) first_time = true; else GFPfile_test.close(); // it doesn't exist yet
    GFP_dist_file.open(GFP_dist_fullpath, ios::out | ios::app);
    GFP_dump_file.open(GFP_dump_fullpath, ios::out | ios::app);
    if(first_time){
        GFP_dist_file << "Time" << sep << "Ki67" << sep << "Lower" << sep << " Upper" << sep <<"Mid "<< sep << "Count" << sep <<"Frequency" << endl;
        GFP_dump_file << "Time" << sep << "Ki67" << sep << "logGFP" << endl;
    }
    // make the tick marks for the bin boundaries. LogGFP is maximum zero. Here, 20 bins from -5 to zero.
    for(i=0;i<NUM_GFP_BINS;i++){
        ticks[i]=LOG_GFP_MIN + float(i) * (LOG_GFP_MAX-LOG_GFP_MIN)/float(NUM_GFP_BINS);
        binned_khi[i]=0; // number of cells in the i-th bin in Ki67 hi
        binned_klo[i]=0; // ditto for Ki67lo
    }
    // loop through all cells, get their Ki67 and GFP expression, and work out what bin it goes in.
    // bin zero  = less than zero
    // bin 1  = zero to x (first bin with). These two will be combined.

    for(i=0; i<poolsize; i++){
        this_cell=*(loclist+i);
        logGFP= (*this_cell).get_logGFP();
        Ki67hi=(*this_cell).get_Ki67_stoch();
        if(Ki67hi) Ki67_status="hi"; else Ki67_status="lo";
        //  dump raw GFP values to file
        if(logGFP<0.) logGFP=0.;
        GFP_dump_file << (int) floor(current_time) << sep << Ki67hi << sep << logGFP << endl;
        ndivs=(*this_cell).get_ndivs();
        if(Ki67hi) khi_numbers++; else klo_numbers++;
        if(logGFP<ticks[0] && Ki67hi)  binned_khi[0]++; // GFP below lower limit
        if(logGFP<ticks[0] && !Ki67hi) binned_klo[0]++;

        if(logGFP>ticks[NUM_GFP_BINS-1] && Ki67hi)  binned_khi[NUM_GFP_BINS-1]++; // greater than last tick
        if(logGFP>ticks[NUM_GFP_BINS-1] && !Ki67hi) binned_klo[NUM_GFP_BINS-1]++; // greater than last tick

        for(j=0; j<(NUM_GFP_BINS-1); j++){ // all the others (bins 1 ... NUM_GFP_BINS-2)
            if(logGFP >= ticks[j] && logGFP<ticks[j+1] && Ki67hi)  {binned_khi[j+1]++; break;}
            if(logGFP >= ticks[j] && logGFP<ticks[j+1] && !Ki67hi) {binned_klo[j+1]++; break;}
        }
    }

    GFP_dump_file.close();
    // add the GFP negative cells to the lowest bin
    binned_khi[1]=binned_khi[1] + binned_khi[0];
    binned_klo[1]=binned_klo[1] + binned_klo[0];

    // write both histograms to file
    for(j=0; j<(NUM_GFP_BINS-2); j++){
        // cells between ticks j and j+1  = binned [j+1]
        GFP_dist_file << (int) floor(current_time) << "," <<  "hi" << "," << ticks[j] << sep << ticks[j+1] << sep << (ticks[j]+ ticks[j+1])/2 << sep <<  binned_khi[j+1]
                      << sep << float(binned_khi[j+1])/float(khi_numbers) << endl;
        GFP_dist_file << (int) floor(current_time) << "," <<  "lo" << "," << ticks[j] << sep << ticks[j+1] << sep << (ticks[j]+ ticks[j+1])/2 << sep <<  binned_klo[j+1]
                      << sep << float(binned_klo[j+1])/float(klo_numbers) << endl;
    }

    GFP_dist_file << (int) floor(current_time) << "," <<  "hi" << "," << ticks[NUM_GFP_BINS-1] << sep << LOG_GFP_MAX << sep <<
                  (ticks[NUM_GFP_BINS-1] + LOG_GFP_MAX)/2 << sep << binned_khi[NUM_GFP_BINS-1] <<  sep <<
                  float(binned_khi[NUM_GFP_BINS-1])/float(khi_numbers) << endl;
    GFP_dist_file << (int) floor(current_time) << "," <<  "lo" << "," << ticks[NUM_GFP_BINS-1] << sep << LOG_GFP_MAX << sep <<
                  (ticks[NUM_GFP_BINS-1] + LOG_GFP_MAX)/2 << sep << binned_klo[NUM_GFP_BINS-1] <<  sep <<
                  float(binned_klo[NUM_GFP_BINS-1])/float(klo_numbers) << endl;
    GFP_dist_file.close();

    return;
}

// calculate Ki67 fraction (not used currently)
float ki67_frac(cell *loclist[], int poolsize,  float current_time, float params[]){
    int i, kpos=0;
    float logKi67_det;
    cell this_cell;
    // loop through all cells and look at ki67 expression
    for(i=0;i<poolsize;i++){
        this_cell=**(loclist+i);
        logKi67_det=this_cell.get_logKi67_det();
        if(logKi67_det>LOG_KPOS_THRESHOLD) kpos++;  // log Ki67 starts at log(1) =0;  this decays to -1 after 3.5d
    }
    return float(kpos)/float(poolsize); // fraction of cells that are ki67 pos
}

// calculate Ki67/GFP fractions and store in vector (passed in here as "y")
// using the deterministic model for Ki67 decay, as in the PDE model.
void ki67_GFP_fracs(cell *loclist[], int poolsize, float y[], float params[]){
    int i;
    int kpos_det_gpos=0, kpos_det_gneg=0, kneg_det_gneg=0;
    int kpos_stoch_gpos=0, kpos_stoch_gneg=0, kneg_stoch_gneg=0;
    float logKi67_det,logGFP;
    bool Ki67_stoch;
    cell this_cell;

    // loop through all cells and look at ki67 and GFP expression
    for(i=0;i<poolsize;i++){
        this_cell=**(loclist+i);
        logKi67_det=this_cell.get_logKi67_det();
        Ki67_stoch=this_cell.get_Ki67_stoch();
        logGFP=this_cell.get_logGFP();
        // quadrants based on deterministic loss of Ki67
        if(logKi67_det>-1.  && logGFP>LOG_GFPPOS_THRESHOLD) kpos_det_gpos++;
        if(logKi67_det>-1.  && logGFP<=LOG_GFPPOS_THRESHOLD) kpos_det_gneg++;
        if(logKi67_det<=-1. && logGFP<=LOG_GFPPOS_THRESHOLD) kneg_det_gneg++;
        // quadrants based on stochastic loss of Ki67
        if(Ki67_stoch  && logGFP>LOG_GFPPOS_THRESHOLD) kpos_stoch_gpos++;
        if(Ki67_stoch  && logGFP<=LOG_GFPPOS_THRESHOLD) kpos_stoch_gneg++;
        if(!Ki67_stoch && logGFP<=LOG_GFPPOS_THRESHOLD) kneg_stoch_gneg++;
    }
    y[0]= float(kpos_det_gpos + kpos_det_gneg)/float(poolsize);  // Ki67+ fraction (deterministic)
    y[1]= float(kpos_stoch_gpos + kpos_stoch_gneg)/float(poolsize);  // Ki67+ fraction (stochastic)
    y[2]= 1. - float(kpos_det_gneg +kneg_det_gneg)/float(poolsize);  // GFP+ fraction
    y[3]=float(kpos_det_gpos)/float(poolsize); //  Ki67+ GFP+ (deterministic)
    y[4]=float(kpos_det_gneg)/float(poolsize); // Ki67+ GFP- (deterministic)
    y[5]=float(kneg_det_gneg)/float(poolsize); // Ki67- GFP- (deterministic)
    y[6]=1. - (y[3]+y[4]+y[5]); // Ki67- GFP+ (deterministic)

    y[7]=float(kpos_stoch_gpos)/float(poolsize); //  Ki67+ GFP+ (stochastic)
    y[8]=float(kpos_stoch_gneg)/float(poolsize); // Ki67+ GFP- (stochastic)
    y[9]=float(kneg_stoch_gneg)/float(poolsize); // Ki67- GFP- (stochastic)
    y[10]=1. - (y[7]+y[8]+y[9]); // Ki67- GFP+ (stochastic)

    return;
}

// calculate Simpson's index  -  sum (1/freq^2) over clone
float simpsons(int poolsize, const int clonerecord[], float params[])
{
    int i;
    float temp=0., temp2;
    for (i=0;i<CLONEUNIVERSESIZE;i++){
        temp2=(float) clonerecord[i];
        temp+=(temp2*temp2);
    }
    temp2= ((float) poolsize)*((float) poolsize)/temp;
    return temp2;
}

// ---------------------------------------------------------------------------------------------------------------

int main (int argc, char * const argv[]) {

    int i=0, j=0, k=0;
    int STARTCELLS;
    int toggle, poolsize,spacelistlength=0,spec, cells_in;
    int clonerecord[CLONEUNIVERSESIZE]; // Track the number of clones of each type, because we can

    int num_GFP_output_times, GFP_output_times[10]; // output GFP expression profiles (up to 10, usually less)
    int num_GFPKi67_output_times, GFPKi67_output_times[100]; // output GFP expression profiles (up to 100)
    float T0, TMAX, sp_numbers_at_T0, THYMIC_EXPORT_RATE_CONSTANT, g0, current_time;
    float simp;
    float fraction_results[11]; // storage for Ki67+/- GFP+/- fracs, for the two defs of Ki67 status
    float params[12];
    cell cellstore[STORAGELISTLENGTH],this_cell; // This is where the cells actually are!
    // now define two list of pointers to these cell storage spots, and one to where the empty slots are
    cell *cell_location_list_1[STORAGELISTLENGTH],*cell_location_list_2[STORAGELISTLENGTH], *spacelist[STORAGELISTLENGTH];

    // bare filenames
    const string filepath =  FILEPATH;
    const string output_filename = OUTPUT_FILENAME; // for cell counts and GFP/Ki67 quadrant freqs
    const string parameter_filename = PARAMETER_FILENAME;
    const string GFP_dist_filename = GFP_DIST_FILENAME; // for binned distributions of GFP expression
    const string GFP_dump_filename = GFP_DUMP_FILENAME; // raw dump of GFP expression per cell, in Ki67 hi and lo

    // construct full paths to these files
    const string output_fullpath = filepath + output_filename;
    const string parameter_fullpath = filepath + parameter_filename;
    const string GFP_dist_fullpath = filepath + GFP_dist_filename;
    const string GFP_dump_fullpath = filepath + GFP_dump_filename;

    // file objects for reading in parameters, and output
    ifstream PAR_FILE;
    ofstream OUTPUT_FILE;
    ofstream GFP_DIST_FILE;
    ofstream GFP_DUMP_FILE;

    // remove any old output files
    remove(output_fullpath.c_str());
    remove(GFP_dist_fullpath.c_str());
    remove(GFP_dump_fullpath.c_str());

    //  random number libraries
    std::random_device rd;
    std::mt19937 generator(rd()); // random seed for mersenne twister
    //std::mt19937 generator(0); // fixed seed

    //uniform distribution between 0 and 1; call with uniformrandom(generator)
    std::uniform_real_distribution<> uniformrandom(0, 1);
    // GFP dists on RTE (ki67 hi and lo);
    std::normal_distribution<> LOG_GFPdist_Ki67HI_RTE(LOG_GFP_MEAN_KHI, LOG_GFP_SD_KHI);
    std::normal_distribution<> LOG_GFPdist_Ki67LO_RTE(LOG_GFP_MEAN_KLO, LOG_GFP_SD_KLO);
    // can use other dists - see e.g. http://www.cplusplus.com/reference/random/

    // read in parameters (all floats)
    // 0: 4 or 8 (CD4 or CD8)
    // 1: N0 (to be scaled by 0.005 to get actual cell numbers)
    // 2: T0 (start time for simulations)
    // 3: TMAX
    // 4: g0  pop density at age zero at time T0. # thymic export rate at time T0 = C* SP(T0)  =  N(t=T0)g(0)
    // 5: N_densitydependence (scale for density dep div or loss - in "true" physiological numbers
    // 6: delta0 (base death rate for cells of age zero)
    // 7: r_delta (death rate goes with cell age as exp(-r_delta * age))
    // 8: rho_0 (base death rate for cells of age zero)
    // 9: r_rho (div rate goes with cell age as exp(-r_rho * age))

    PAR_FILE.open(parameter_fullpath);
    for(i=0; i<10; i++)  PAR_FILE >> params[i];
    PAR_FILE.close();

    STARTCELLS=(int) params[1];
    poolsize=STARTCELLS;
    T0=params[2];
    TMAX=params[3];
    g0=params[4];
    sp_numbers_at_T0 = sp_numbers(T0, params);
    THYMIC_EXPORT_RATE_CONSTANT = ((float) poolsize)*g0/(sp_numbers_at_T0 * SCALE_OUTPUT);
    //THYMIC_EXPORT_RATE_CONSTANT = 0.;

    // set times at which to output GFP histograms and Ki67/GFP quadrant fractions
    num_GFP_output_times = 2;
    GFP_output_times[0]=5.;
    GFP_output_times[1]=15.;

    //num_GFPKi67_output_times= 20;
    //for(i=0; i<num_GFPKi67_output_times; i++){
    //    GFPKi67_output_times[i] = 5.+ (i*5.) ;
    //}

    // Set up pre-existing cells at time T0
    // initialise the list of pointers to each cell, and do some initialising of Ki67 and GFP and age at T0.
    for(i=0;i<poolsize;i++){ // loop through all cells
        spec=rand() % CLONEUNIVERSESIZE; // vaguely random number between 0 and CLONEUNIVERSESIZE-1
        (*(cellstore+i)).set_cloneID(spec);// give this cell this random TCR from the universe of TCRs
        // note there that "cellstore+i" is a pointer to the i-th cell in the array; and *(cellstore+i) is that cell itself.
        clonerecord[spec]++; // make a note that we've got one of these clonotypes
        cell_location_list_1[i]=cellstore+i;  // the i-th entry in cell_location_list_1 is now a pointer to this cell.
        spacelist[i]=nullptr; // space "i" is occupied
        // set up cell properties at T0. It could be more sophisticated than this.
        (*(cellstore+i)).set_export_time(T0 * uniformrandom(generator));// cells randomly aged between 0 and T0 days, at T0
        // Work out if it's Ki67 high or not - this impacts its GFP expression
        // make the Ki67 (binary = stochastic) fraction equal to 0.80737 = fitted Ki67+ fraction at T0=5 for neutral ODE model
        // start by assuming it's Ki67-
        (*(cellstore+i)).set_Ki67_stoch(false);
        (*(cellstore+i)).set_logKi67_det(-2); // less than -1 means Ki67 negative. So this can be anything
        (*(cellstore+i)).set_logGFP(LOG_GFPdist_Ki67LO_RTE(generator));

        if(uniformrandom(generator)>0.5) { // if it's Ki67 high, recalculate its GFP
            (*(cellstore + i)).set_Ki67_stoch(true);
            // deterministic model of Ki67 expression - positive cells are in Ki67+ interval (-0.1,0)
            (*(cellstore + i)).set_logKi67_det(-0.1 * uniformrandom(generator));
            (*(cellstore + i)).set_logGFP(LOG_GFPdist_Ki67HI_RTE(generator));
        }
        // cout << (*(cellstore + i)).get_Ki67_stoch() << sep << (*(cellstore + i)).get_logGFP() << endl;
    }

    //--------------------------------------------
    //main loop

    // we keep two lists of cell locations; we use one to update the other, and then flip-flop with this variable
    toggle=0;

    // Calculate Ki67/GFP fracs in population; store in vector 'fraction_results_stoch', with elements:
    // 0 = Ki67+ fraction, 1 = GFP+ fraction
    // 2 =Ki67+ GFP+, 3 = Ki67+ GFP-, 4 = Ki67- GFP-, 5 = Ki67- GFP+
    ki67_GFP_fracs(cell_location_list_1, poolsize, fraction_results, params);

    // Write main output file headers and initial conditions
    OUTPUT_FILE.open(output_fullpath);
    OUTPUT_FILE << "time , time.int, sim_counts , physiol_counts, sp.numbers , sp.ki67 , new.RTE , Ki67.det.pos , Ki67.stoch.pos ,GFP.pos , ";
    OUTPUT_FILE << "Ki67.det.pos.GFP.pos , Ki67.det.pos.GFP.neg , Ki67.det.neg.GFP.neg , Ki67.det.neg.GFP.pos , ";
    OUTPUT_FILE << "Ki67.stoch.pos.GFP.pos , Ki67.stoch.pos.GFP.neg , Ki67.stoch.neg.GFP.neg , Ki67.stoch.neg.GFP.pos" << endl;

    OUTPUT_FILE << T0 << sep  << (int) floor(T0) << sep << poolsize << sep  << ((float) poolsize)/SCALE_OUTPUT << sep << sp_numbers(T0,params) << sep << sp_ki67(T0,params) << sep << "NA" ;
    for (i=0; i<11; i++) OUTPUT_FILE << sep << fraction_results[i];
    OUTPUT_FILE << endl;

    //gfp_histogram(cell_location_list_1, poolsize, T0, params);
    // ------ now loop over timesteps
    for(current_time=T0+TSTEP;current_time<=TMAX; current_time+=TSTEP) {
        cout << current_time << endl;
        if (toggle == 0) {
            update(cell_location_list_1, cell_location_list_2, cellstore, clonerecord,
                   &poolsize, spacelist, &spacelistlength, current_time, params, THYMIC_EXPORT_RATE_CONSTANT);
            ki67_GFP_fracs(cell_location_list_2, poolsize, fraction_results, params);

        } else {
            update(cell_location_list_2, cell_location_list_1, cellstore, clonerecord,
                   &poolsize, spacelist, &spacelistlength, current_time, params, THYMIC_EXPORT_RATE_CONSTANT);
            ki67_GFP_fracs(cell_location_list_1, poolsize, fraction_results, params);
        }
        toggle = 1 - toggle;

        cells_in = lrint( THYMIC_EXPORT_RATE_CONSTANT * sp_numbers(current_time, params) * SCALE_OUTPUT * TSTEP);

        OUTPUT_FILE << current_time << sep << (int) floor(current_time) << sep << poolsize << sep <<  ((float) poolsize)/SCALE_OUTPUT << sep << sp_numbers(current_time, params)
                    << sep << sp_ki67(current_time, params) << sep << cells_in;
        for (i = 0; i < 11; i++) OUTPUT_FILE << sep << fraction_results[i];
        OUTPUT_FILE << endl;
        // Do we output gfp histgram to std out?
        if (j < num_GFP_output_times & GFP_output_times[j] == (int) floor(current_time)) {
            j++;
            // Make sure we use the most up-to-date list of cells
            if (toggle == 1) gfp_histogram(cell_location_list_2, poolsize, current_time, params);
            else gfp_histogram(cell_location_list_1, poolsize, current_time, params);
        }

    }
    cout << "... done!" << endl;
    return 0;
}
