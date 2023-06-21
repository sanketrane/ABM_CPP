
//
//  NonSpatial-ABM
//
//  Written by Sanket Rane.
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


#include "custom_functions.cpp"

// defining fixed variables and parameters
# define SCALE_OUTPUT 0.005f // Scale N0 and influx by this number to get actual simulated cell numbers
# define CLONEUNIVERSESIZE 100 // number of possible clonotypes/TCR sequences
# define TSTEP 0.04f // time step. needs to be small, bcs event probabilities are estimated as (process rate)*TSTEP
# define STORAGELISTLENGTH 100000 // max number of cells
# define LOG_KPOS_THRESHOLD -1.f // Define threshold of Ki67 positivity. Set to exp(-1) consistent with Sanket's ASM
// max value of deterministic Ki67 is 1 - > log equals zero. So Ki67+ cells lie between -1 and 0
# define KI67_LOSS_RATE 0.2857143f  // 1/3.5d
# define Host_SP_Ki67 0.33f // the mean ki67 postive fraction in SP4 post 8 weeks of age.
# define EXP_MINUS_1 0.3678794412f // handy constant
# define KPOS_THRESHOLD 0.3678794412f // handy constant
# define LOG2 0.69314718f // another handy constant

# define sep "  ,  " // column separator for output files (e.g. comma or tab)


// Cell object definition
class cell
{
private:
    int cloneID; // 1...CLONEUNIVERSESIZE
    float ki67_intens_norm; // normalized Ki67 expression. 
    bool donor_derived; // whether its a donor or host BM derived
    bool location_thymus; // whether its a thymic or peripheral Treg
    int ndivs; // number of divisions since differentiation into Treg phenotype in the thymus
    float last_division_time; // time at which this cell last divided
    float export_time;// time at which its ancestor differentiatiated into Treg

public:
    cell(); // default constructor fn.
    cell(int p1, float p2, bool p3, bool p4, int p5, float p6, float p7); // general constructor fn.
    void set_cloneID(int x){cloneID=x;}
    void set_ki67_intens_norm(float x){ki67_intens_norm=x;}
    void set_donor_derived(bool x){donor_derived=x;}
    void set_location_thymus(bool x){location_thymus=x;}
    void set_ndivs(int i){ndivs=i;}
    void set_last_division_time (float x) {last_division_time =x;}
    void set_export_time (float x) {export_time=x;}

    int get_cloneID() {return(cloneID);}
    float get_ki67_intens_norm() {return(ki67_intens_norm);}
    bool get_donor_derived(){return (donor_derived);}
    bool get_location_thymus(){return (location_thymus);}
    int get_ndivs() {return(ndivs);}
    float get_last_division_time() {return(last_division_time);}
    float get_export_time() {return(export_time);}
};

cell::cell(){ //default constructor
    cloneID=0;
    ki67_intens_norm=randunif<float>(0, 1);    // for a given cell ki67 exoression is chosen from a uniform distribution -- this needs redfinition!!!
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
// 0: N0 (to be scaled by 0.005 to get actual cell numbers)
// 1: T0 (start time for simulations)
// 2: TMAX
// 3: g0  pop density at age zero at time T0. # thymic export rate at time t=1 = C* SP(1)  =  N(t=1)g(0)
// 4: psi -- perc capita rate of influx
// 5: delta0 (base death rate for cells of age zero)
// 6: r_delta (death rate goes with cell age as exp(-r_delta * age))
// 7: N_densitydependence (scale for density dep div or loss - in "true" physiological numbers
// 8: rho_0 (base death rate for cells of age zero)
// 9: r_rho (div rate goes with cell age as exp(-r_rho * age))
// 10: TBMT -- host age at bmt
// 11: rate of loss/division of an incumbent cell

 // spline1 --
float sp_numbers(float t, float params[]){ 
    // these are true numbers - need to scale down for simulation (done elsewhere)
    float dpt0, T0=params[1], sp_numbers, theta0, nu, psi=params[4];

    // thymic FoxP3 negative SP4 T cell numbers - spline defs
    theta0  = 6.4;  nu = 0.0024;
    dpt0 = t - T0;     // days post t0
    if (dpt0<0.) dpt0=0.;
    
    sp_numbers = psi * pow(10., theta0) * exp(-1. * nu * dpt0); 
    return sp_numbers;
}

 // spline2 --
float chi_spline(float t, float params[]) {
   // gives fraction of donor in the thymic FoxP3 negative SP4 T cells 
   float dpBMT, TBMT=params[10], chi_value, chiEst, qEst;

   // spline def for chimerism
   chiEst = 0.816; qEst = 0.063;
   dpBMT = t - TBMT;     // days post Bone Marrow Transplant
   if (dpBMT<0.) dpBMT=0.;

   chi_value = chiEst * (1. - exp(-qEst * (dpBMT - 0.)));
   if (dpBMT - 0. < 0.) chi_value = 0.;           // assumption: no donor cells seen in FoxP3neg SP4 compartment for ~ 10 days
     
   return chi_value;
 }

 float donor_eps_spline(float t, float params[]){
  // gives fraction of ki67 hi cells in the donor thymic FoxP3 negative SP4 T cells 
   float dpt0, T0=params[1],  eps_value, eps_0, eps_f;

  // spline def for ki67 fraction 
  eps_0 = 0.378;  eps_f = 0.068;
  dpt0 = t - T0;     // days post t0
  if (dpt0<0.) dpt0=0.;

  eps_value = exp(- eps_f * dpt0) + eps_0;
  return eps_value;
}

float lossrate(int poolsize, float current_time, float cell_age, float time_since_last_division, int ndivs, float params[]){
    // parame defs
    float delta0=params[5]; float r_delta= params[6];  
    float N_densitydependence = params[7]; // work with "true" physiological numbers
    float host_age_at_export = current_time - cell_age;
    
    if (host_age_at_export < 0. ) host_age_at_export = 0.;
    //return (delta0 * exp(-1.* cell_age * r_delta));
    return delta0;
}

float divrate(int poolsize, float current_time, float cell_age, float time_since_last_division, int ndivs, float params[]){
    // parame defs
    float rho0=params[8]; float r_rho= params[9];
    float N_densitydependence = params[7]; // work with "true" physiological numbers
    float host_age_at_export = current_time - cell_age;

    if (host_age_at_export < 0. ) host_age_at_export = 0.;
    //return (rho0 * exp(-1.* cell_age * r_rho));
    return rho0;
}

void update(cell *fromlist[], cell *tolist[], cell cellstore[], int clonerecord[], int *poolsize, cell *spacelist[],
            int *spacelistlength, float current_time, float params[], float THYMIC_EXPORT_RATE_CONSTANT){
    int i, j=0, cells_in, cloneID, ndivs, updated_poolsize; // j tracks number of live cells in this cycle
    float sp_chimerism, matureSP_Ki67pos_frac, runif, age, updated_age, ki67_intens_norm, last_division_time, time_since_last_division, export_time, p_div, p_loss;
    bool dead, divide, location_thymus, donor_derived;
    cell *this_cell; cell *new_cell_location;    // definitions for new and prexisting cells using the class cell

    updated_poolsize = *poolsize; 

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
        ki67_intens_norm  = ki67_intens_norm * exp(- TSTEP*KI67_LOSS_RATE); // Ki67 drops with time
        (*this_cell).set_ki67_intens_norm(ki67_intens_norm);

        // Now decide what the cell does this timestep
        // calculate its probs of division and loss ... = rate * timestep (if timestep small!)
        p_div  = divrate(*poolsize, current_time, age, time_since_last_division, ndivs, params)*TSTEP;
        p_loss = lossrate(*poolsize, current_time, age, time_since_last_division, ndivs, params)*TSTEP;
        // Carve the unit interval into (0, p_loss, (1-p_div), 1) and see where random # lands.
        // p_loss and p_div should both be small (short timestep) so shouldn't get into trouble here
        dead=false;
        divide=false;
        runif=randunif<float>(0, 1);
        if(runif<p_loss) dead=true;
        // hack - don't let cells re-divide before 0.2
        //if(runif>(1-p_div) && time_since_last_division>0.2) divide=true;
        if(runif>(1-p_div)) divide=true;
        
        if(dead){
            divide=false; // just in case of shenanigans with the probabilities above
            updated_poolsize--; // one less cell
            (*(clonerecord+cloneID))--; // and one less of this clone
            if(updated_poolsize<=0) {
                std::cout << "Poolsize equals zero" << '\n';
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
            (*this_cell).set_ki67_intens_norm(1.); // max Ki67 expression = 1
            (*this_cell).set_location_thymus(false); // progeny leaves thymus after division
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
            (*new_cell_location).set_ki67_intens_norm(1.); // max Ki67 expression = 1 
            (*new_cell_location).set_location_thymus(false); // progeny leaves thymus after division
            (*new_cell_location).set_donor_derived((*this_cell).get_donor_derived()); // same as parent and sibling
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
        updated_age=current_time -  (*this_cell).get_export_time();
    }

    ////////////////////////////////////////////
    // now add new cells from the thymus (integer)
    // we work in rescaled numbers (i.e. scaled down from "true" numbers)
    ////////////////////////////////////////////

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

        //make random number to determine whether new cell is donor or host,
        // depending on the chimerism in SP cells
        runif= randunif<float>(0., 1.);
        sp_chimerism = chi_spline(current_time, params);
        
        if(runif < sp_chimerism)  { // It's donor derived;
            (*new_cell_location).set_donor_derived(true);
            
            //make random number to determine whether new cell is Ki67 high or low,
            // depending on Ki67 expression in donor SP cells
            runif= randunif<float>(0., 1.);
            matureSP_Ki67pos_frac = donor_eps_spline(current_time, params);
            
            if(runif < matureSP_Ki67pos_frac)  { 
                // It's Ki67 positive;
                // now choose a value from a uniform distribution between and Ki67+ threshold i.e. exp(-1) and 1.
                // Ki67 on absolute expression scale is in (0,1). So (-infinity, 0) on log scale
                (*new_cell_location).set_ki67_intens_norm(randunif<float>(KPOS_THRESHOLD, 1));
            } else {
                // it's Ki67 negative and chose a value from a uniform distribution between  0 and Ki67+ threshold i.e. exp(-1).
                (*new_cell_location).set_ki67_intens_norm(randunif<float>(0, KPOS_THRESHOLD));
            }
        } else {
            // it's host derived
            (*new_cell_location).set_donor_derived(false);

            //make random number to determine whether new cell is Ki67 high or low,
            // depending on Ki67 expression in donor SP cells
            runif= randunif<float>(0., 1.);
            matureSP_Ki67pos_frac = Host_SP_Ki67;
            
            if(runif < matureSP_Ki67pos_frac)  { 
                // It's Ki67 positive;
                // now choose a value from a uniform distribution between and Ki67+ threshold i.e. exp(-1) and 1.
                // Ki67 on absolute expression scale is in (0,1). So (-infinity, 0) on log scale
                (*new_cell_location).set_ki67_intens_norm(randunif<float>(KPOS_THRESHOLD, 1));
            } else {
                // it's Ki67 negative and chose a value from a uniform distribution between  0 and Ki67+ threshold i.e. exp(-1).
                (*new_cell_location).set_ki67_intens_norm(randunif<float>(0, KPOS_THRESHOLD));
            }
        }

        // Other attributes of the new cell
        (*new_cell_location).set_cloneID(cloneID);
        (*new_cell_location).set_ndivs(0);
        (*new_cell_location).set_export_time(current_time);
        (*new_cell_location).set_last_division_time(-999999.);
    }
    // reassign the poolsize pointer to the unpdated number
    *poolsize=updated_poolsize;
}

// calculate Donor fraction
void donor_Ki67_frac(cell *loclist[], float current_time, int poolsize,  float y[], float params[]){
    int i, donorcount=0, donor_kpos=0, host_kpos=0;
    float donor_derived, ki67_intens_norm;
    cell this_cell;
    
    // loop through all cells and look at ki67 expression
    for(i=0;i<poolsize;i++){
        this_cell=**(loclist+i);
        donor_derived=this_cell.get_donor_derived();
        if(donor_derived)  {// bool variable
            donorcount++;  
            // check ki67 frac
            ki67_intens_norm=this_cell.get_ki67_intens_norm();
            if(ki67_intens_norm>KPOS_THRESHOLD) donor_kpos++;  
        } else {
            // its a host derived cell
            // check ki67 frac
            ki67_intens_norm=this_cell.get_ki67_intens_norm();
            if(ki67_intens_norm>KPOS_THRESHOLD) host_kpos++;  
        }
    }
    y[0] = float(donorcount)/float(poolsize); // fraction of cells that are donor derived
    y[1] = float(donorcount)/(float(poolsize) * chi_spline(current_time, params)); if (donorcount <= 0.0) y[1] = 0.0; // donor fraction normalised to chimerism in SP cells
    y[2] = float(donor_kpos)/float(donorcount); if (donorcount <= 0.0) y[2] = 0.0; // fraction of donor cells that are Ki67+
    y[3] = float(host_kpos)/(float(poolsize) - float(donorcount)); // fraction of host cells that are Ki67+
}

// calculate Ki67 fraction
void ki67_frac(cell *loclist[], int poolsize,  float y[], float params[]){
    int i, kpos=0;
    float ki67_intens_norm;
    cell this_cell;
    // loop through all cells and look at ki67 expression
    for(i=0;i<poolsize;i++){
        this_cell=**(loclist+i);
        ki67_intens_norm=this_cell.get_ki67_intens_norm();
        if(ki67_intens_norm>KPOS_THRESHOLD) kpos++;  // log Ki67 starts at log(1) =0;  this decays to -1 after 3.5d
    }
    y[0] = float(kpos)/float(poolsize); // fraction of cells that are ki67 pos
}

int main (int argc, char * const argv[]) {

    int i=0, j=0, k=0;
    int STARTCELLS;
    int toggle, poolsize,spacelistlength=0,spec, cells_in;
    int clonerecord[CLONEUNIVERSESIZE]; // Track the number of clones of each type, because we can

    float T0, TMAX, TBMT, sp_numbers_at_T0, THYMIC_EXPORT_RATE_CONSTANT, g0, current_time;
    float params[11];
    float fraction_results[4]; 
    cell cellstore[STORAGELISTLENGTH],this_cell; // This is where the cells actually are!
    // now define two list of pointers to these cell storage spots, and one to where the empty slots are
    cell *cell_location_list_1[STORAGELISTLENGTH],*cell_location_list_2[STORAGELISTLENGTH], *spacelist[STORAGELISTLENGTH];

    //setting current WD and filepaths for input and output
    std::string parfile ("mytest_parameters.txt");

    std::string outname = "output_csv/Neutral/outfile_";
    std::string const& arrayid = argv[1];
    std::string outfile = outname+arrayid+".csv";


    // file objects for reading in parameters, and output
    std::ifstream PAR_FILE;
    std::ofstream OUTPUT_FILE;

    // remove any old output files
    remove(outfile.c_str());

    PAR_FILE.open(parfile);
    for(i=0; i<12; i++)  PAR_FILE >> params[i];
    PAR_FILE.close();

    STARTCELLS=(int) params[0];
    poolsize=STARTCELLS;
    T0=params[1];
    TMAX=params[2];
    g0=params[3];
    TBMT=params[10];
    sp_numbers_at_T0 = sp_numbers(T0, params);
    THYMIC_EXPORT_RATE_CONSTANT = ((float) poolsize)*g0/(sp_numbers_at_T0 * SCALE_OUTPUT);

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
        (*(cellstore+i)).set_export_time(T0 * randunif<double>(0, 1));// cells randomly aged between 0 and T0 days, at T0
    }

    // we keep two lists of cell locations; we use one to update the other, and then flip-flop with this variable
    toggle=0;

    // Write main output file headers and initial conditions
    
    
    OUTPUT_FILE.open(outfile);
    OUTPUT_FILE << "time , time.int, sim_counts , physiol_counts, sp.numbers , new_RTE , Donor_fraction , Normalized_fd , Donor_Ki67_pos , Host_Ki67_pos " << '\n';

    OUTPUT_FILE << T0 << sep  << (int) floor(T0) << sep << poolsize << sep  << ((float) poolsize)/SCALE_OUTPUT << sep << sp_numbers(T0,params) << sep << 0.0 ;
    for (i=0; i<4; i++) OUTPUT_FILE << sep << fraction_results[i];
    OUTPUT_FILE << '\n';

    for(current_time=T0+TSTEP;current_time<=TMAX; current_time+=TSTEP) {
        //std::cout << current_time << '\n';
        if(fmod(current_time, 25) == 0.0){
            std::cout << current_time << '\n';
        }

        if (toggle == 0) {
            update(cell_location_list_1, cell_location_list_2, cellstore, clonerecord,
                   &poolsize, spacelist, &spacelistlength, current_time, params, THYMIC_EXPORT_RATE_CONSTANT);
            donor_Ki67_frac(cell_location_list_2, current_time, poolsize, fraction_results, params);

        } else {
            update(cell_location_list_2, cell_location_list_1, cellstore, clonerecord,
                   &poolsize, spacelist, &spacelistlength, current_time, params, THYMIC_EXPORT_RATE_CONSTANT);
            donor_Ki67_frac(cell_location_list_1, current_time, poolsize, fraction_results, params);
        }
        toggle = 1 - toggle;

        cells_in = lrint( THYMIC_EXPORT_RATE_CONSTANT * sp_numbers(current_time, params) * SCALE_OUTPUT * TSTEP);

        OUTPUT_FILE << current_time << sep << (int) floor(current_time) << sep << poolsize << sep <<  ((float) poolsize)/SCALE_OUTPUT << sep << sp_numbers(current_time, params)
                    << sep << cells_in;
        for (i = 0; i < 4; i++) OUTPUT_FILE << sep << fraction_results[i];
        OUTPUT_FILE << '\n';
    }
    

    std::cout << "... done!" << '\n';

    return 0;
}
