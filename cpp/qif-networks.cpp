//Libraries
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <random>
#include <vector>
#include <string>
#include "CNetwork/source/CNet.cpp" //Route local to this file
#include <chrono>

using namespace std;
using namespace std::chrono;

// --- PRECOMPILER OPTIONS --- //

//Auxiliary
#define FALSE 0
#define TRUE 1

//Deterministc Lorentzian or noisy Gaussian
#ifndef DETERMINISTIC
#define DETERMINISTIC FALSE
#endif

//Still unused, will be useful for parallelizing MF code later
#ifndef PARALLEL
#define PARALLEL 1
#endif

//Network model to use in this simulation
#define ERDOS 0
#define CLUSTERS 1
#define HMODULAR 2
#define HMCORE 3
#define CORE 4

//Network selection
#ifndef NETWORK 
#define NETWORK HMCORE
#endif


//Inhibitory percentage
#ifndef INH
const static double alpha = 0.0;
#define INH FALSE
#else
const static double alpha = INH;
#undef INH
#define INH TRUE
#endif

#ifndef INHBIAS
#define INHBIAS FALSE
const static double inhbias = -1.0;
#else
const static double inhbias = INHBIAS;
#undef INHBIAS
#define INHBIAS TRUE
#endif


//Modes of use 
#define TSERIES 0 //Obtain time series
#define DIAGRAM 1 //Compute phase diagram

#ifndef MODE
#define MODE TSERIES
#endif


// --- Global parameters --- //

unsigned int N = 1000;           //Total number of neurons
int Ne, Ni;         //Number of E/I neurons, respectively
double J = 10;          //Synaptic weight
double eta = -5.0;       //Intrinsic frequency
double sigma = 2.5;    //Noise intensity
string filename = "";   //Base filename for output


// --- Network parameters --- //

vector<double> network_params = vector<double>();   //Parameters that define the graph

// --- Neurons parameters --- //

const double Vp = 100.0;    //Threshold voltage
const double Vr = -Vp;      //Rest voltage

const double tau = 1e-3;    //Time that spikes are able to influence other neurons

const double refractory = 1.0 / Vp; //Neuron refractory period

//Structure that represents a neuron, so we can build a network of these
struct neuron
{
    double v;   //Membrane potential
    double tk = -1.0;  //Time of last spike, setting a default value
    double eta;
    double refractory = 0.0; //1.0 / Vp; //Refractory period 
    int inh = +1;      //Inhibitory neurons have this negative

};

// --- Simulation --- //

unsigned int N_raster = 400;   //Number of neurons to print in raster
vector<int> raster_indices; //Indices of the neurons that will be written to raster

const double h = 1e-4;      //Timestep
const double sqh = sqrt(h); //Stochastic timestep
const double therm_time = 10.0;     //Time to thermalize 10
const double measure_time = 30.0;   //Time to make measures 30
const double trelax = therm_time;    //Relaxation time 10
const double tinput = measure_time + trelax;    //Time of input 40
const double tmax (measure_time + 30 +  tinput);    //Simulation end


const double dt = 2e-2;     //Firing rate time binning
const int write_rate = dt / h; //To write data in the file

// --- Auxiliary stuff --- //

mt19937 gen(45323543);      //RNG 
uniform_real_distribution<double> ran_u(0.0, 1.0);  //Typical [0,1) distribution (antes neuronas_raster)
normal_distribution<double> ran_g(0.0, 1.0);        //Typical normal distribution (antes llamada ruido)

// --- Function declarations --- //

void initial_conditions(CNetwork<neuron> &net);
void step(CNetwork<neuron> &net, const double iext, double &t,  int &t_rate, double &avg_v, double &rate, double &re, double &ri, ofstream &file_potencial, ofstream &file_raster, ofstream &file_rate);
void step(vector<neuron> &net, const double iext, double &t, int &t_rate, double &rate, double &mean_v, double &mean_v2, double &mean_r, double &mean_r2);
void step_relax(CNetwork<neuron> &net, double &t,  const double iext);

void euler_maruyama(double &V, const double I);

//Computation modes
void compute_time_series();
void compute_phase_diagram(const bool use_intensity, const int sim_id, const double sigma0, const double sigmaf, const double dsigma);

//All network topologies (different implementation depending on NETWORK, see below)
void create_network(CNetwork<neuron> &net);
void set_indices_for_rasterplot();
string show_net_args(); //Helper function that shows you the network arguments when you don't remember them :)



// --- Main --- ///

int main(int argc, char* argv[])
{

    int i;
    
    #if MODE==TSERIES

    int nargs = 6;
    if (argc >= nargs)
    {
        N = stoi(argv[1]);
        N_raster = min(N_raster, N);
        J = stod(argv[2]);
        eta = stod(argv[3]);
        sigma = stod(argv[4]);
        filename = string(argv[5]);
        //Network parameters. Look specifically for them and load all in the parameter vector
        if (argc > nargs)
        {
            if (string(argv[nargs]) == "-net") 
            {
                for (i=nargs + 1; i < argc; i++) network_params.push_back(stod(argv[i])); 
            }   
        }
    }
    else 
    {
        cout << "Expected arguments: N, J, eta, sigma, filename, -net " + show_net_args()  << endl;
        return EXIT_SUCCESS;
    }

    auto start = high_resolution_clock::now(); 
    try 
    { 
        compute_time_series();
    }
    catch  (const exception &e) 
    {
        cout << "ERROR DURING EXECUTION. EXCEPTION DETECTED: " << endl << e.what() << endl;
        return EXIT_FAILURE;
    }
    auto end = high_resolution_clock::now(); 

    auto sim_time = duration_cast<milliseconds>(end - start);
    cout << endl << "Simulation time: " << sim_time.count() / 1000.0 << "s" << endl;

    #elif MODE==DIAGRAM

    double sigma0 = 0.0;
    double sigmaf = 1.0;
    int nsigma = 10;
    int id = 0;
    double dsigma;
    bool use_iext = false;

    int nargs = 10;
    if (argc >= nargs)
    {

        N = stoi(argv[1]);
        N_raster = min(N_raster, N);
        J = stod(argv[2]);
        eta = stod(argv[3]);
        sigma0 = stod(argv[4]);
        sigmaf = stod(argv[5]);

        nsigma = stoi(argv[6]);

        use_iext = bool(stoi(argv[7]));

        filename = string(argv[8]); 
        id = stoi(argv[9]);

        //Network parameters. Look specifically for them and load all in the parameter vector
        if (argc > nargs)
        {
            if (string(argv[nargs]) == "-net") for (i=nargs+1; i < argc; i++) network_params.push_back(stod(argv[i]));
        }
    }
    else 
    {
        cout << "Expected arguments: N, J, eta, sigma0, sigmaf, nsigma, use_iext, filename, id, -net " + show_net_args()  << endl;
        return EXIT_SUCCESS;
    }

    dsigma = (sigmaf - sigma0) / (1.0 * nsigma);

    auto start = high_resolution_clock::now(); 
    try
    {
        compute_phase_diagram(use_iext, id, sigma0, sigmaf, dsigma);
    }
    catch(const exception& e)
    {
        cout << "ERROR DURING EXECUTION. EXCEPTION DETECTED: " << endl << e.what() << endl;
        return EXIT_FAILURE;
    }
    
    auto end = high_resolution_clock::now(); 

    auto sim_time = duration_cast<milliseconds>(end - start);
    cout << endl << "Simulation time: " << sim_time.count() / 1000.0 << "s" << endl;
    #endif

    return EXIT_SUCCESS;
}



// --- Functions --- //

// Core of the simulation program, these function make the hard computations


//Initalize the network with 
void initial_conditions(CNetwork<neuron> &net)
{
    int i,j;
    double diff = (Vp - Vr); //Useful later to generate random initial conditions

    int is_inh; //Aux variable to determine if neuron is inhibitory

    double k_z;


    //Lorentzian distribution, taking in account that the term sigma in a Langevin equation corresponds
    //to a noise sqrt(2sigma) noise intensity
    #if DETERMINISTIC==TRUE
    cauchy_distribution<double> cauchy(eta, sigma);
    #endif

    #if INHBIAS==TRUE
    vector<double> k_weight(N);
    k_z = 0.0;
    for (i=0; i < N; i++)
    {
        k_weight[i] = pow(net.degree(i), inhbias);
        k_z += k_weight[i];
    }
    k_z = 1.0 / k_z;
    #endif

    for (i=0; i < N; i++) 
    {
        net[i].v = ran_u(gen) * diff + Vr; //Set all potentials
        net[i].tk = -1; //Reset spikes


        //Heterogeneity in distribution of frequencies or not
        #if DETERMINISTIC==TRUE
        net[i].eta = cauchy(gen);
        #else
        net[i].eta = eta;
        #endif

        //If we have some inhibitory neurons, initialize them
        #if INH==TRUE

        //Use random choice OR weight neurons according to its connectivity
        #if INHBIAS==FALSE
        is_inh = 2*(ran_u(gen) > alpha) - 1; 
        #else
        is_inh = ran_u(gen) <= k_weight[i] * k_z;
        #endif

        net[i].inh = is_inh;
        Ne += is_inh > 0;
        Ni += is_inh < 0;
        #endif
    }

    return;
}

//Integration step
void step(CNetwork<neuron> &net, const double iext, double &t, int &t_rate, double &avg_v, double &rate, double &re, double &ri, ofstream &file_potencial, ofstream &file_raster, ofstream &file_rate)
{
    static int j,k;               //Counter

    static double tdiff;        //Difference between times
    static double tdiff2;      //The same, but this time we need 2 ;)
    static int non_refractory;  //Number of active neurons
    static double s;            //Mean field spike actuvity
    
    static bool is_in_raster;   //True if the neuron is inside the raster_indices vector
    static int current_raster;  //Next index for the raster_index vector that will be written to file

    static int degree;          //Degree of a given node
    static int neigh;           //Index of the network neigh

    #if INH==TRUE
    bool is_inh;
    static double avg_ve, avg_vi;
    static int non_ref_e, non_ref_i;

    static double se, si;

    avg_ve = avg_vi = 0.0;
    non_ref_e = non_ref_i = 0;
    #endif


    //Then integrate all these equations using adequate input
    avg_v = 0.0;
    non_refractory = 0;
    current_raster = 0;
    for (j=0;j < N; j++)
    {
        //Use that raster_indices is sorted in order to know if index j in contained 
        //is_in_raster = (j == raster_indices[current_raster]);

        //Refractoriness for a time 2*refractory = 2 / Vp
        tdiff = t - net[j].tk;
        if (tdiff >= refractory)
        {
            //Compute input to this neuron
            s = 0.0;
            
            #if INH==TRUE
            se = si = 0.0;
            #endif
            
            degree = net.degree(j);
            for (k=0; k < degree; k++)
            {
                neigh = net.get_in(j, k); 
                tdiff2 = t - net[neigh].tk; 

                #if INH==FALSE
                s += net[neigh].inh *((tdiff2 > 0 && tdiff2 <= tau));
                #else
                se += (net[neigh].inh > 0) *((tdiff2 > 0 && tdiff2 <= tau));
                si += (net[neigh].inh < 0) *((tdiff2 > 0 && tdiff2 <= tau));
                #endif
            }

            #if INH==FALSE
            s *= J / (degree * tau);
            #else
            s = J * (se - si) / (degree * tau);
            #endif

            //Count non-refractory neurons
            non_refractory++;
            avg_v += net[j].v; 

            #if INH==TRUE
            is_inh = net[j].inh < 0;
            non_ref_e += not is_inh;
            non_ref_i += is_inh;

            avg_ve += (not is_inh) * net[j].v; 
            avg_vi += is_inh * net[j].v; 
            #endif

            //Integrate
            euler_maruyama(net[j].v, net[j].eta + iext + s);

            //If spike, set the spike time and reset. There is no problem in computing the rate now
            //due to the discretization of time
            if (net[j].v >= Vp)
            {
                net[j].tk = t + refractory;
                net[j].v = Vr;
                rate++;

                #if INH==TRUE
                re += not is_inh;
                ri += is_inh;
                #endif

                //This neuron made a spike and is inside raster, hence put it to raster plot!
                //We use current raster to keep the file always in interval [0,N_raster[
                //if (is_in_raster) file_raster << net[j].tk << " " << current_raster << endl;
                file_raster << net[j].tk << " " << j << endl;
            }

            //If the neuron was included, advance the raster_indices marker to next one
            //current_raster += is_in_raster;
        }

    }   //End of j<N for

    //Finish average potential over active neurons
    avg_v /= 1.0 * non_refractory;

    //Write firing rate and potential every 100 steps
    if (t_rate % write_rate == 0)
    {
        rate = 1.0*rate/(N*dt);

        #if INH==FALSE
        file_rate << t << " " << rate << endl;
        file_potencial << t << " " << avg_v << endl;
        #else
        avg_ve /= 1.0 * non_ref_e;
        avg_vi /= 1.0 * non_ref_i;

        re = 1.0 * re / (Ne * dt);
        ri = 1.0 * ri / (Ni * dt); 

        file_rate << t << " " << rate << " " << re << " " << ri << endl;
        file_potencial << t << " " << avg_v << " " << avg_ve << " " << avg_vi << endl;

        re = ri = 0;
        #endif
        
        rate = 0;
    }
    t_rate++;
} 


//Integration step
void step(CNetwork<neuron> &net, const double iext, double &t, int &t_rate, double &rate, double &mean_v, double &mean_v2, double &mean_r, double &mean_r2)
{
    static int j,k;               //Counter

    static double tdiff;        //Difference between times
    static double tdiff2;      //The same, but this time we need 2 ;)
    static int non_refractory;  //Number of active neurons
    static double s;            //Mean field spike actuvity
    static double avg_v;        //Average membrane potential

    static int degree;          //Degree of a given node
    static int neigh;           //Index of the network neigh

    //Then integrate all these equations using adequate input
    avg_v = 0.0;
    non_refractory = 0;
    for (j=0;j < N; j++)
    {
        //Refractoriness for a time 2*refractory = 2 / Vp
        tdiff = t - net[j].tk;
        if (tdiff >= refractory)
        {
            //Compute input to this neuron
            s = 0.0;
            degree = net.degree(j);
            for (k=0; k < degree; k++)
            {
                neigh = net.get_in(j, k); 
                tdiff2 = t - net[neigh].tk; 
                s += net[neigh].inh *((tdiff2 > 0 && tdiff2 <= tau));
            }
            s *= J / (degree * tau);

            //Count non-refractory neurons
            non_refractory++;
            avg_v += net[j].v; 

            //Integrate
            euler_maruyama(net[j].v, net[j].eta + iext + s);
            
            //If spike, set the spike time and reset. There is no problem in computing the rate now
            //due to the discretization of time
            if (net[j].v >= Vp)
            {
                net[j].tk = t + refractory;
                net[j].v = Vr;
                rate++;
            }
        }

    }   //End of j<N for

    //Finish average potential over active neurons
    avg_v /= 1.0 * non_refractory;

    //Write firing rate and potential every 100 steps
    if (t_rate % write_rate == 0)
    {
        rate *= 1.0/(N*dt);

        //Make measures
        mean_v += avg_v;
        mean_v2 += avg_v * avg_v;
        mean_r += rate;
        mean_r2 += rate*rate;

        rate = 0;
    }
    t_rate++;
} 




//Integration step used to relax the system. It is possible to excite it with an external current
//This one does not compute compute ANY average.
void step_relax(CNetwork<neuron> &net, double &t, const double iext)
{
    static int j,k;          //Counter

    static double tdiff, tdiff2; //Difference between times
    static int non_refractory;   //Number of active neurons
    static double s;             //Mean field spike actuvity

    static int degree;           //Degree of a given node
    static int neigh;           //Index of the network neigh


    //Then integrate all these equations using adequate input
    for (j=0;j < N; j++)
    {
        tdiff = t - net[j].tk;
        if (tdiff >= refractory)
        {
            //Compute input to this neuron
            s = 0.0;
            degree = net.degree(j);
            for (k=0; k < degree; k++)
            {
                neigh = net.get_in(j, k); 
                tdiff2 = t - net[neigh].tk; 
                s += net[neigh].inh * ((tdiff2 > 0 && tdiff2 <= tau));
            }
            s *= J / (degree * tau);


            euler_maruyama(net[j].v, net[j].eta + iext + s);
            
            if (net[j].v >= Vp)
            {
                net[j].tk = t + refractory;
                net[j].v = Vr;
            }
        }
    }

} 


//Integrate using the Euler-Maruyama stochastic integrator
#if DETERMINISTIC==TRUE
inline void euler_maruyama(double &V, const double I)
{
    V += h*(V*V+I);
}
#else
inline void euler_maruyama(double &V, const double I)
{
    V += h*(V*V+I) + sqh * sigma * ran_g(gen);
}
#endif











// --- PROGRAM MODES --- //

// Functions below do the complete simulations depending on the specified MODE. Be sure to check what
//does each function. They should have no arguments whatsoever, receiving only the global parameters that
//were set by the code or by the terminal

//This mode compute time series, including mean rate and potential, as well as a raster plot. 
void compute_time_series()
{
    //Auxiliary
    int i,j;        
    const double iext = 3.0;

    //Time counters
    double t;       //Time
    int t_rate;     //Firing rate counter

    //Neuronal Network
    CNetwork<neuron> network(N);   //Network of N neurons

    //Observables
    double avg_v, rate; //Average potential and rate
    double re, ri;

    //Output
    ofstream file_rate; //Firing rate output
    ofstream file_raster; //Raster plot output
    ofstream file_potencial; //Average potential output

    
    //Create a network of N neurons.
    cout << "Setting up system" << endl;
    create_network(network); //Can throw exceptions
    cout << "Network created." << endl;
    initial_conditions(network);

    //Select raster plot neurons (choose first N_raster, other options are possible)
    set_indices_for_rasterplot();

    //System thermalization
    cout << "Relaxing system..." << endl;
    for (t=0; t < trelax; t +=h) step_relax(network, t, 0.0);   //Then let the system relax without any external input

    //Create output files
    file_rate.open(filename + "_fr");
    file_rate << "#Parameters: N=" << N << ", eta=" << eta << ", J=" << J << ", sigma=" << sigma << endl;
    file_potencial.open(filename + "_v");
    file_potencial << "#Parameters: N=" << N << ", eta=" << eta << ", J=" << J << ", sigma=" << sigma << endl;
    file_raster.open(filename + "_raster");
    file_raster << "#Parameters: N=" << N << ", eta=" << eta << ", J=" << J << ", sigma=" << sigma << endl;

    //Finally do the simulation storing the measures
    cout << "Start simulation" << endl;
    avg_v = 0.0;//Start V average
    rate = 0;   //Start FR counter
    t_rate = 0; //Start time for FR counter

    //Allows to include external input if needed
    for (t=trelax; t <= tinput; t += h) step(network, 3.0, t, t_rate, avg_v, rate, re, ri, file_potencial, file_raster, file_rate);
    for (t=tinput; t <= tmax; t += h) step(network, 0.0, t, t_rate, avg_v, rate, re, ri, file_potencial, file_raster, file_rate);
    
    //Close files
    file_rate.close();
    file_raster.close();
    file_potencial.close();

    cout << "Simulation finished" << endl;
    return;
}


void compute_phase_diagram(const bool use_intensity, const int sim_id, const double sigma0, const double sigmaf, const double dsigma)
{
    //Auxiliary
    int i,j;        
    const double iext = 3.0;

    double rate; //Current count of spikes
    const double n_measures = measure_time / (h * write_rate); //Number of measures we did

    //Time counters
    double t;       //Time
    int t_rate;     //Firing rate counter

    //Observables
    double mean_v, mean_v2;
    double mean_r, mean_r2;

    //Neuronal Network
    CNetwork<neuron> network(N);   //Network of N neurons

    //Output
    ofstream results;

    //Create output files
    results.open(filename + "_" + to_string(sim_id));

    //Do all the simulatins over the same network
    create_network(network); //Can throw exceptions

    //Create N neurons
    for (sigma = sigma0; sigma <= sigmaf; sigma += dsigma)
    {
        cout << sigma << endl;
        initial_conditions(network); //Init the system

        //System thermalization
        for (t=0; t < trelax; t +=h) step_relax(network, t, 0.0);   //Then let the system relax without any external input

        //If there is external intensity, apply it. 
        if (use_intensity) for (t=trelax; t <= tinput; t += h)  step_relax(network, t, iext);
        else for (t=trelax; t <= tinput; t += h)  step_relax(network, t, 0.0);

        //Last step is let the system decay with no intensity, measuring it.
        rate = 0;               //Start FR counter
        t_rate = 0;             //Start time for FR counter
        mean_v = mean_v2 = 0.0; //Observable init
        mean_r = mean_r2 = 0.0;
        for (t=tinput; t <= tmax; t += h)  step(network, 0.0, t, t_rate, rate, mean_v, mean_v2, mean_r, mean_r2);

        //Finish averages and write mean and variance to file
        mean_v /= n_measures;
        mean_v2 /= n_measures;
        mean_r /= n_measures;
        mean_r2 /= n_measures;

        results << sigma << " " << mean_v << " " << mean_v2 - mean_v*mean_v << " " << mean_r << " " << mean_r2 - mean_r*mean_r << endl;
    }
    results.close();

    cout << "Simulation finished" << endl;
    return;
}
















// --- NETWORK TOPOLOGIES --- //

//The code below creates different network topologies dependent on the NETWORK specified.
//The second vectors, params, contain the desired parameters for each network. Be sure to check each function
//to know how it works.
//If an incorrect number of params is fed to the function, it just throws an exception and exits.


#if NETWORK==ERDOS

void create_network(CNetwork<neuron> &net)
{
    if (network_params.size() != 2) throw runtime_error("Critical error: incorrect number of network parameters. Expected 2");
    else
    {
        double k = network_params[0];
        unsigned int seed = int(network_params[1]);

        net.create_erdos_renyi(N, k, seed);
    }
}


//In Erdos-Renyi, it does not matter which indices you choose, so take the first ones
void set_indices_for_rasterplot()
{
    int i;
    raster_indices = vector<int>(N_raster);
    for (i=0; i < N_raster; i++) raster_indices[i] = i;
}

string show_net_args()
{
    return "k, seed";
}


#elif NETWORK==CLUSTERS

void create_network(CNetwork<neuron> &net)
{
    if (network_params.size() != 4) throw runtime_error("Critical error: incorrect number of network parameters. Expected 4");
    else
    {

        //Network parameters
        const double kc = network_params[0]; //Cluster connectivity
        const unsigned int seed = int(network_params[1]);
        const unsigned int n_clusters = int(network_params[2]);
        const double ki = network_params[3]; //Cluster-to-outside ratio. 0 perfect clusters, more random as r grows...

        int i,j; //Counters
        
        //Create the nodes
        net.add_nodes(N);
        
        //Cluster variables
        unsigned int neurons_per_cluster = int(N / n_clusters);
        int c_index_i, c_index_j; //Cluster index per neurons i,j
        double p, pcluster, pout; //Connection probabilities

        //Set connection probabilities
        pcluster = kc / neurons_per_cluster;
        pout = ki / (N - neurons_per_cluster);

        //Erdos-Renyi construction with the adequate prob.
        for (i=0; i < N; i++)
        {
            c_index_i = int (i / neurons_per_cluster);

            for (j=i+1; j < N; j++)
            {
                c_index_j = int (j / neurons_per_cluster);

                p = c_index_i == c_index_j ? pcluster : pout; 
                if (ran_u(gen) <= p) net.add_link(i,j);
            }
        }

    } //endelse
}

//For clusters, set the same number of neurons from each cluster
void set_indices_for_rasterplot()
{
    const unsigned int n_clusters = int(network_params[2]);
    unsigned int neurons_per_cluster = int(N / n_clusters);

    const unsigned int n_cluster_raster =  int(N_raster / n_clusters);
    int i,j;

    raster_indices = vector<int>(N_raster);


    for (i=0; i < n_cluster_raster; i++)
    {
        for (j=0; j < n_clusters; j++)
        {
            raster_indices[i + j*n_cluster_raster] = i + j*neurons_per_cluster;
        }
        
    }
}

string show_net_args()
{
    return "k_cluster, seed, num_clusters, ratio_kc_2_ki";
}



#elif NETWORK==CORE
void create_network(CNetwork<neuron> &net)
{
    if (network_params.size() != 4) throw runtime_error("Critical error: incorrect number of network parameters. Expected 4");
    else
    {
        const double kc = network_params[0];                              //Cluster connectivity
        const unsigned int seed = int(network_params[1]);
        const unsigned int n_clusters = int(network_params[2]);           //Number of clusters
        const unsigned int neurons_core_cluster = int(network_params[3]); //Neurons from each cluster belonging to the core
        
        int i,j; //Counters

        net.add_nodes(N);
        //Cluster variables
        unsigned int neurons_per_cluster = int(N / n_clusters);
        int c_index_i, c_index_j; //Cluster index per neurons i,j
        double p, pcluster; //Connection probabilities
        bool core_i, core_j; //Indicates if the node belongs to the core

        //Set connection probabilities
        pcluster = kc / neurons_per_cluster;
        
        //Erdos-Renyi construction with the adequate prob.
        for (i=0; i < N; i++)
        {
            c_index_i = int(i / neurons_per_cluster);
            core_i = (i % neurons_per_cluster) < neurons_core_cluster;
            
            for (j=i+1; j < N; j++)
            {
                core_j = (j % neurons_per_cluster) < neurons_core_cluster;
                if (core_i && core_j) 
                {
                    net.add_link(i,j);
                }
                else
                {
                    c_index_j = int(j / neurons_per_cluster);
                    p = c_index_i == c_index_j ? pcluster : 0; 

                    if (ran_u(gen) <= p) net.add_link(i,j);
                }
            }   
        }
    }//end else
}

//It is exactly the same as in the CLUSTER case
void set_indices_for_rasterplot()
{
    const unsigned int n_clusters = int(network_params[2]);
    unsigned int neurons_per_cluster = int(N / n_clusters);

    const unsigned int n_cluster_raster =  int(N_raster / n_clusters);
    int i,j;

    raster_indices = vector<int>(N_raster);

    for (i=0; i < n_cluster_raster; i++)
    {
        for (j=0; j < n_clusters; j++)
        {
            raster_indices[i + j*n_cluster_raster] = i + j*neurons_per_cluster;
        }
        
    }
}

string show_net_args()
{
    return "k_cluster, seed, num_clusters, neruons_in_core";
}



#elif NETWORK==HMODULAR

//Modular network has H levels. Level H=0 is clusters of neurons. H=1 is connectivity among clusters, and so on.
//Connectivity: the vector {k0, k1, ..., kH} means that each neuron has in average k0 connections with level H=0, k1
//conections with level H=1, and so on. Ideally it should look like {20, 2, 1}, or even lower values of higher level clusters.
//The number of elements in each levels goes as {n0, n1, ..., nH}, and means that level j as nj macro-components. Therefore,
//vector {300, 10, 2} means 2 modules, where each module has 10 clusters, and each cluster 300 neurons.
//This is, a total number of 2 modules, 20 clusters, and 6000 neurons.


//Auxiliary function. Given the indices of two neurons and information about the hierarchical level,
//return the level at which they are connected .
//Take in account a each neuron belongs to levels (for example) {3,0,0} means  cluster 3, module 0. 
//Last index is always 0 since the last level is the whole network, so in worst-case nlevels-1 is returned
int same_level(const int index_i, const int index_j, const int nlevels, const vector<int> &neurons_level)
{
    int m;
    int level_index_i, level_index_j;

    m=0;
    do 
    {
        level_index_i = int(index_i / neurons_level[m]);
        level_index_j = int(index_j / neurons_level[m]);

        m++;
    } while (m < nlevels && level_index_i != level_index_j);

    return m-1;
}


void create_network(CNetwork<neuron> &net)
{
    int i,j;
    int nlevels;

    //First of all get number of hierarchical levels
    if (network_params.size() == 0) throw runtime_error("Critical error: incorrect number of network parameters. Expected 4");
    
    nlevels = network_params[0];      

    vector<int> n_elements_level(nlevels);      //Number of modules in each level(lower level modules are neurons)
    vector<int> n_neurons_level(nlevels);       //Number of neurons in each level (lower level modules are neurons). 
    vector<double> k_level(nlevels);            //Connectivity for each level. 
    vector<double> p_level(nlevels);            //Probability of connection for each level. 

    //Aux variable to store at which level neurons coincide
    int coincidence_level;

    //This should have, as arguments: nlevels {list_n_elements} {list_k}, so 1+2*nlevels arguments. If not, throw error
    if (1 + 2*nlevels != network_params.size()) throw runtime_error("Critical error: incorrect number of network parameters. Expected 4");

    //Fill the information provided 
    for (i=0; i < nlevels; i++)
    {
        //network_params uses i+1 since 0 is nlevels
        n_elements_level[i] = int(network_params[i+1]);

        //Get total number of neurons until this level. First level is just the number of neuron, so make distinc cases 
        n_neurons_level[i] = i > 0 ? n_elements_level[i] * n_neurons_level[i-1] : n_elements_level[i]; 

        k_level[i] = network_params[i+1+nlevels];
    }

    //Last level contains all neurons, so make a consistency check
    if (N != n_neurons_level[nlevels-1]) throw runtime_error("Critical error: the number of neurons indicated does not coincide with hierarchical level computation, that was " + to_string(n_neurons_level[nlevels-1]));
    
    net.add_nodes(N);

    //Once we have the total number of neurons and the connectivities get the connection probabilities
    for (i=0; i < nlevels; i++) p_level[i] = k_level[i] / (1.0 * n_neurons_level[i]);

    //Erdos-Renyi-like construction with the adequate prob.
    for (i=0; i < N; i++)
    {
        for (j=i+1; j < N; j++)
        {
            //Get at which level both neurons coincide. Higher level is the whole network, with index nlevels-1
            coincidence_level = same_level(i, j, nlevels, n_neurons_level);

            //Add a link with the probability corresponding to that level
            if (ran_u(gen) < p_level[coincidence_level])
            {
                net.add_link(i, j);
            }
        }   
    }

    net.write_graphml("generated-network");
}

//STILL NOT IMPLEMENTED
void set_indices_for_rasterplot()
{
    const unsigned int n_clusters = int(network_params[2]);
    unsigned int neurons_per_cluster = int(N / n_clusters);

    const unsigned int n_cluster_raster =  int(N_raster / n_clusters);
    int i,j;

    raster_indices = vector<int>(N_raster);

    for (i=0; i < n_cluster_raster; i++)
    {
        for (j=0; j < n_clusters; j++)
        {
            raster_indices[i + j*n_cluster_raster] = i + j*neurons_per_cluster;
        }
        
    }
}

string show_net_args()
{
    return "nlevels {list_n_elements} {list_k}, where both lists have nlevels elements";
}

#elif NETWORK==HMCORE

//Read the comments of HMODULAR to read all information about how modules are constructed. Here it is the same, but all modules are constructed in a scale-free way.


//Auxiliary function. Given the indices of two neurons and information about the hierarchical level,
//return the level at which they are connected .
//Take in account a each neuron belongs to levels (for example) {3,0,0} means  cluster 3, module 0. 
//Last index is always 0 since the last level is the whole network, so in worst-case nlevels-1 is returned
int same_level(const int index_i, const int index_j, const int nlevels, const vector<int> &neurons_level)
{
    int m;
    int level_index_i, level_index_j;

    m=0;
    do 
    {
        level_index_i = int(index_i / neurons_level[m]);
        level_index_j = int(index_j / neurons_level[m]);

        m++;
    } while (m < nlevels && level_index_i != level_index_j);

    return m-1;
}

//Auxiliary function to get the cumulative probabilities for each neuron at each hierarchical level
void get_node_cumulative(const int neurons_cluster, const double beta, vector<double> &cumuweight)
{
    int i;

    double hna;
    double power_beta;

    if (beta < 0.0)
    {
        vector<double> weights(neurons_cluster);
        cumuweight = vector<double>(neurons_cluster, 0.0);

        //Store weights and compute the complete sum in order to compute later the probabilities
        hna = 1.0;
        weights[0] = 1.0;
        for (i=2; i <= neurons_cluster; i++) 
        {
            power_beta = pow(i, beta);
            hna += power_beta;
            weights[i-1] = power_beta;
        }

        //Then directly get the cumulatives
        cumuweight[0] = weights[0] / hna;
        for (i=1; i < neurons_cluster; i++) cumuweight[i] = cumuweight[i-1] + weights[i] / hna;
    }
    else
    {
        vector<double> weights(neurons_cluster, 1.0 / neurons_cluster);
        cumuweight = vector<double>(neurons_cluster, 0.0);

        cumuweight[0] = weights[0];
        for (i=1; i < neurons_cluster; i++) cumuweight[i] = cumuweight[i-1] + weights[i];
    }

}

//Return the index of a random node inside the cluster of specified index
int get_random_node(const int cluster, const int neurons_cluster, const vector<double> &cumuweight)
{
    //Get a random_number
    double r = ran_u(gen);

    //lower_bound gets the index where I would insert it in the (sorted) vector of cumulative probabilities, 
    //hence getting the (iterator) index of the selected one.
    auto lower = lower_bound(cumuweight.begin(), cumuweight.end(), r);

    //Return it as integer and not as iterator
    return distance(cumuweight.begin(), lower) + cluster * neurons_cluster;
}

//Create the hierarchical network.
//Input arguments: nlevels {list_n_elements} {list_k}
void create_network(CNetwork<neuron> &net)
{
    int i,j,level,link,cluster; //Counters
    int nclusters;

    int nlevels; //Number of hierarchical levels

    //Auxiliary variables
    int nlinks;             //Number of links in this level
    int comm1, comm2;       //Indices of two randomly selected communities
    double beta;            //Exponent of this level

    //First of all get number of hierarchical levels
    if (network_params.size() == 0) throw runtime_error("Critical error: incorrect number of network parameters. Expected 4");
    
    nlevels = network_params[0];      

    vector<int> n_elements_level(nlevels);      //Number of modules in each level(lower level modules are neurons)
    vector<int> n_neurons_level(nlevels);       //Number of neurons in each level (lower level modules are neurons). 
    vector<double> k_level(nlevels);            //Connectivity for each level. 
    vector<double> gamma_level(nlevels);        //Power-law exponents
    vector<double> p_level(nlevels);            //Probability of connection for each level. 


    //This should have, as arguments: nlevels {list_n_elements} {list_k}, so 1+2*nlevels arguments. If not, throw error
    if (1 + 3*nlevels != network_params.size()) throw runtime_error("Critical error: incorrect number of network parameters. Expected 4");

    //Fill the information provided 
    for (i=0; i < nlevels; i++)
    {
        //network_params uses i+1 since 0 is nlevels
        n_elements_level[i] = int(network_params[i+1]);

        //Get total number of neurons until this level. First level is just the number of neuron, so make distinc cases 
        n_neurons_level[i] = i > 0 ? n_elements_level[i] * n_neurons_level[i-1] : n_elements_level[i]; 

        //Set connectivity and power law exponents.
        k_level[i] = network_params[i+1+nlevels] ;
        gamma_level[i] = network_params[i+1+2*nlevels];
    }


    //Last level contains all neurons, so make a consistency check
    if (N != n_neurons_level[nlevels-1]) throw runtime_error("Critical error: the number of neurons indicated does not coincide with hierarchical level computation, that was " + to_string(n_neurons_level[nlevels-1]));
    
    //Once the consistency check is done, add nodes and select weights for scale free network
    net.add_nodes(N);
    vector<double> cumulative;

    //Now construct the network. We do first the first level, connecting the individual neurons
    //Number of links to generate in this level (number of links in each cluster)
    //We use the number of links of the random graph as an estimator
    nlinks = 0.5 * k_level[0] * n_neurons_level[0]; 

    //Beta to use in this level. Values gamma < 1 lead to beta=0.0, i.e. random connectivity
    beta = gamma_level[0] > 1.0 ? 1.0 / (gamma_level[0] - 1.0) : 0.0; 

    //Fill the node weights for this level
    get_node_cumulative(n_neurons_level[0], -beta, cumulative); 

    //Compute the total number of clusters at the first level
    nclusters = 1;
    for(i=1; i < nlevels; i++) nclusters *= n_elements_level[i];

    //Then set the links. This part will be different in higher levels. 
    //Each cluster will have nlinks, in order to have the selected k in this cluster!
    for (cluster=0; cluster < nclusters; cluster++)
    {
        link = 0;
        while (link < nlinks)
        {
            //Select two neurons from two different communities. Ensure they are different
            do
            {
                i = get_random_node(cluster, n_neurons_level[0], cumulative);
                j = get_random_node(cluster, n_neurons_level[0], cumulative);
            } while (i == j);

            //If there is no link between these two, add it 
            if (net.get_link_index(i,j) == -1) 
            {
                net.add_link(i,j);
                link++;
            }
        }
    }


    //Finally do the same for higher levels
    for (level=1; level < nlevels; level++)
    {
        //These lines are the same as in level 0
        nlinks = 0.5 * k_level[level] * n_neurons_level[level]; 
        beta = gamma_level[level] > 1.0 ? 1.0 / (gamma_level[level] - 1.0) : 1.0; 

        //Set the cumulative probability of selecting each node
        get_node_cumulative(n_neurons_level[level-1], -beta, cumulative); 

        //Compute the number of modules in this level. Now each 'cluster' (module) contains a
        //a set of clusters of the previous level. Last level only has 1 cluster.
        //Use that it was previosly computed by multipliying everybody, to 'substract' the number of elements in this level
        nclusters /= n_elements_level[level];

        //Finally, for each module, set the links
        for (cluster=0; cluster < nclusters; cluster++)
        {
            link = 0;
            while (link < nlinks)
            {                
                //Select two nodes from two different communities, from the current cluster. 
                //Ensure they are different communities
                do
                {
                    comm1 = int(ran_u(gen) * n_elements_level[level]) + cluster * n_elements_level[level];
                    comm2 = int(ran_u(gen) * n_elements_level[level]) + cluster * n_elements_level[level];
                } while (comm1 == comm2);

                //And now obtain a neuron from that community. Check that the link is yet non-existent
                i = get_random_node(comm1, n_neurons_level[level-1], cumulative);
                j = get_random_node(comm2, n_neurons_level[level-1], cumulative);

                //If everything is correct, then make a link
                if (net.get_link_index(i,j) == -1) 
                {
                    net.add_link(i,j);
                    link++;
                }
            }
        }
    }

    //Record the used graph
    net.write_graphml("generated_network", vector<string>());
}

//STILL NOT IMPLEMENTED
void set_indices_for_rasterplot()
{
    const unsigned int n_clusters = int(network_params[2]);
    unsigned int neurons_per_cluster = int(N / n_clusters);

    const unsigned int n_cluster_raster =  int(N_raster / n_clusters);
    int i,j;

    raster_indices = vector<int>(N_raster);

    for (i=0; i < n_cluster_raster; i++)
    {
        for (j=0; j < n_clusters; j++)
        {
            raster_indices[i + j*n_cluster_raster] = i + j*neurons_per_cluster;
        }
    }
}

string show_net_args()
{
    return "nlevels {list_n_elements} {list_k} {gamma_k}, where all lists have nlevels elements";
}

#endif