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
namespace fs = std::filesystem;

//setting current WD and filepaths for input and output
auto wdir = fs::current_path();
fs::path parfile ("mytest_parameters.txt");
fs::path fullparfile_path = wdir / parfile;

// defining fixed variables and parameters
# define SCALE_OUTPUT 0.005f // Scale N0 and influx by this number to get actual simulated cell numbers
# define CLONEUNIVERSESIZE 100 // number of possible clonotypes/TCR sequences
//# define TSTEP 0.04f // time step. needs to be small, bcs event probabilities are estimated as (process rate)*TSTEP
# define TSTEP 0.4f // time step. needs to be small, bcs event probabilities are estimated as (process rate)*TSTEP
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
    float ki67_intens_log; // normalized Ki67 expression. 
    bool donor_derived; // whether its a donor or host BM derived
    bool location_thymus; // whether its a thymic or peripheral Treg
    int ndivs; // number of divisions since differentiation into Treg phenotype in the thymus
    float last_division_time; // time at which this cell last divided
    float export_time;// time at which its ancestor differentiatiated into Treg

public:
    cell(); // default constructor fn.
    cell(int p1, float p2, bool p3, bool p4, int p5, float p6, float p7); // general constructor fn.
    void set_cloneID(int x){cloneID=x;}
    void set_ki67_intens_log(float x){ki67_intens_log=x;}
    void set_donor_derived(bool x){donor_derived=x;}
    void set_location_thymus(bool x){location_thymus=x;}
    void set_ndivs(int i){ndivs=i;}
    void set_last_division_time (float x) {last_division_time =x;}
    void set_export_time (float x) {export_time=x;}

    int get_cloneID() {return(cloneID);}
    float get_ki67_intens_log() {return(ki67_intens_log);}
    bool get_donor_derived(){return (donor_derived);}
    bool get_location_thymus(){return (location_thymus);}
    int get_ndivs() {return(ndivs);}
    float get_last_division_time() {return(last_division_time);}
    float get_export_time() {return(export_time);}
};

cell::cell(){ //default constructor
    cloneID=0;
    ki67_intens_log=randunif<double>(0, LOG_KPOS_THRESHOLD);    // for a given cell ki67 exoression is chosen from a uniform distribution between 0 and exp(-1)
    donor_derived=false;
    location_thymus=true;
    ndivs=0;
    last_division_time=-99999.; // last division time not defined for undivided cells
    export_time=0.;
}

cell::cell(int p1, float p2, bool p3, bool p4, int p5, float p6, float p7){ //general constructor
    cloneID=p1;
    ki67_intens_log=p2;
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


float sp_numbers(float t, float params[]){ // mature SP thymocyte numbers - spline defs
    // these are true numbers - need to scale down for simulation (done elsewhere)
    float dpt0, T0=params[1], sp_numbers, theta0, nu;

    // define the splines for SP4 and SP8
    theta0  = 6.4;  nu = 0.0024;
    dpt0 = t - T0;     // days post t0
    if (dpt0<0.) dpt0=0.;
    sp_numbers = pow(10, theta0) * exp(-1 * nu * dpt0); 
    return sp_numbers;
}

float lossrate(int poolsize, float current_time, float cell_age, float time_since_last_division, int ndivs, float params[]){
    float delta0=params[5];
    float r_delta= params[6];
    float host_age_at_export = current_time - cell_age;
    float N_densitydependence = params[4]; // work with "true" physiological numbers
    if (host_age_at_export < 0. ) host_age_at_export = 0.;
    // insert function etc. here to calculate host age dependent rho_0
    //return (delta0 * exp(-1.* cell_age * r_delta));
    return delta0;
}

float divrate(int poolsize, float current_time, float cell_age, float time_since_last_division, int ndivs, float params[]){
    float rho0=params[7];
    float r_rho= params[8];
    float host_age_at_export = current_time - cell_age;
    float N_densitydependence = params[4]; // work with "true" physiological numbers
    if (host_age_at_export < 0. ) host_age_at_export = 0.;
    // if (host_age_at_export <30.) rho0=rho0*0.02; else rho0=rho0*1.2;
    // insert function etc. here to calculate host age dependent rho_0
    //return (rho0 * exp(-1.* cell_age * r_rho));
    return rho0;
}


int update(cell *fromlist[], cell *tolist[], cell cellstore[], int clonerecord[], int *poolsize, cell *spacelist[],
            int *spacelistlength, float current_time, float params[], float THYMIC_EXPORT_RATE_CONSTANT)
{
    int i, j=0, cells_in, cloneID, ndivs, updated_poolsize=*poolsize, bb; // j tracks number of live cells in this cycle
    float runif, age, ki67_intens_log, last_division_time, time_since_last_division, export_time, p_div, p_loss;
    cell *this_cell;
    cell *new_cell_location;
    bool dead, divide, location_thymus, donor_derived;

    std::cout << std::boolalpha;

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
        ki67_intens_log= (*this_cell).get_ki67_intens_log();
        donor_derived= (*this_cell).get_donor_derived();
        location_thymus= (*this_cell).get_location_thymus();
        ndivs=(*this_cell).get_ndivs();
        last_division_time=(*this_cell).get_last_division_time();
        export_time=(*this_cell).get_export_time();
        age=current_time-export_time;
        time_since_last_division = current_time - last_division_time;

        // update loggfp and logki67 values - linear drop on log scale
        ki67_intens_log  = ki67_intens_log - TSTEP*KI67_LOSS_RATE; // Ki67 drops with time
        (*this_cell).set_ki67_intens_log(ki67_intens_log);

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



        cout << "clone# " << cloneID << "div prop: " << 1-p_div << " divided: " << divide <<  " death prop: " << p_loss << " dead: "<< dead << " dice: " << runif << endl;

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
            (*this_cell).set_ki67_intens_log(0.); // max Ki67 expression = 1
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
            (*new_cell_location).set_ki67_intens_log(0.); // same as sibling
            j=j+2; // count 2 cells here
        }
        if(!dead && !divide) { //  survives and does nothing
            // if it's Ki67_stoch positive, does it become Ki67_stoch neg?
            //if((*this_cell).get_Ki67_stoch()){
            //    if(uniformrandom(generator)<TSTEP*KI67_LOSS_RATE)  (*this_cell).set_Ki67_stoch(false);
            //}
            //  add it to the list
            *(tolist+j)=*(fromlist+i);
            j++; // count one more cell
        }

        cout << "I am clone number " << cloneID << " and I come from donor is " << donor_derived 
        << ". I am " << age << " days old. I divided " << ndivs << " many times, most recently " 
        << time_since_last_division << " days ago. My Ki67 intensity is " << exp(ki67_intens_log) << endl;
    }

    return 0;
}

int main (int argc, char * const argv[]) {

    int i=0, j=0, k=0;
    int STARTCELLS;
    int toggle, poolsize,spacelistlength=0,spec, cells_in;
    int clonerecord[CLONEUNIVERSESIZE]; // Track the number of clones of each type, because we can

    float T0, TMAX, sp_numbers_at_T0, THYMIC_EXPORT_RATE_CONSTANT, g0, current_time;
    float params[4];
    cell cellstore[STORAGELISTLENGTH],this_cell; // This is where the cells actually are!
    // now define two list of pointers to these cell storage spots, and one to where the empty slots are
    cell *cell_location_list_1[STORAGELISTLENGTH],*cell_location_list_2[STORAGELISTLENGTH], *spacelist[STORAGELISTLENGTH];

    // file objects for reading in parameters, and output
    ifstream PAR_FILE;
    ofstream OUTPUT_FILE;

    //  random number libraries
    std::random_device rd;
    std::mt19937 generator(rd()); // random seed for mersenne twister
    //std::mt19937 generator(0); // fixed seed

    //uniform distribution between 0 and 1; call with uniformrandom(generator)
    std::uniform_real_distribution<> uniformrandom(0, 1);

    

    PAR_FILE.open(fullparfile_path);
    for(i=0; i<10; i++)  PAR_FILE >> params[i];
    PAR_FILE.close();

    STARTCELLS=(int) params[0];
    poolsize=STARTCELLS;
    T0=params[1];
    TMAX=params[2];
    g0=params[3];
    sp_numbers_at_T0 = sp_numbers(T0, params);
    THYMIC_EXPORT_RATE_CONSTANT = ((float) poolsize)*g0/(sp_numbers_at_T0 * SCALE_OUTPUT);

    cout << " I started at " << T0 << " I will run till " << TMAX << endl;

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
    }

    for(current_time=T0+TSTEP;current_time<=TMAX; current_time+=TSTEP) {
        cout << current_time << endl;
        update(cell_location_list_1, cell_location_list_2, cellstore, clonerecord,
        &poolsize, spacelist, &spacelistlength, current_time, params, THYMIC_EXPORT_RATE_CONSTANT);
    }
    cout << "... done!" << endl;
    return 0;
}