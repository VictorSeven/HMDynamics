/* -------------------------------------------------------- *
*  This code generates hierarchical-modular networks with a *
*  core-periphery structure, based on a work by             *
*  Zamora-Lopez et al., 2016, SciRep.                       *
*  Implementation: Victor Buendía, 2021                     *
* --------------------------------------------------------- */

#define MAXNETSIZE 10000

#define SINGLE 0 
#define DIAGRAM 1 
#define TIMETRACE 2
#define TIMESERIES 3
#define CHIMERA 4 

#ifndef MODE
#define MODE CHIMERA 
#endif 

//Include libraries
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <fstream>
#include <random>
#include <vector>
#include "CNetwork/source/CNet.cpp" //Route local to this file. Find it on github.com/VictorSeven/CNetwork

using namespace std;


// --- Declare functions --- //
bool set_up_network(CNetwork<double> &net, const string filename);
void initial_conditions(CNetwork<double> &net);

void step_no_measures(CNetwork<double> &net);
void step(CNetwork<double> &net);
void step_chimera(CNetwork<double> &net, const int n_moduli, const int osc_per_modulus);

void simulate_single(CNetwork<double> &net, ofstream &output, const double control);
void simulate_single_chimera(CNetwork<double> &net, ofstream &output, const int n_moduli, const int osc_per_modulus);

void simulate_diagram(CNetwork<double> &net, double &selected_var, const double var0, const double varf, const int nvar, ofstream &output);
void time_traces(CNetwork<double> &net, ofstream &output, const int ntraces, const double duration);
void time_traces_chimera(CNetwork<double> &net, ofstream &output, const int n_moduli);

// --- Random number generation --- //
mt19937 gen(85345385434);
uniform_real_distribution<double> ran_u(0.0, 1.0);
normal_distribution<double> ran_g(0.0, 1.0);
cauchy_distribution<double> lorentzian(0.0, 1.0);
normal_distribution<double> gaussian;

// --- Global constants --- //
int N;
int osc_per_modulus = 100;

vector<double> w, a;
double w0,delta,q,a0,s;
string networkname = "hmrandom";
string filename = "kuramoto";

const double dt = 0.01;
const double sqdt = sqrt(dt);

double tf = 4000.0;
double trelax = 600.0;
const double tmeasure = 5.0;
const int nitswindow = 50;

const complex<double> I = complex<double>(0.0, 1.0);

complex<double> z;
double r, psi;

vector<complex<double>> z_modulus;

// --- Main --- ///
// This is used to construct and export the desired network to a file

int main(int argc, char* argv[])
{
    //Counters
    int i,j;

    //Define the network and the path to it
    CNetwork<double> net(MAXNETSIZE);

    //Input/output stuff
    bool correct_setup;

    ofstream output;    //File to write stuff



    #if MODE==SINGLE

        if (argc == 8)
        {
            w0    = stod(argv[1]);
            delta = stod(argv[2]);
            s     = stod(argv[3]);
            a0     = stod(argv[4]);
            q     = stod(argv[5]);
            networkname = string(argv[6]);
            filename = string(argv[7]); 
        }

        //Set up the network, including initial conditions
        correct_setup = set_up_network(net, networkname);
        if (!correct_setup) 
        {
            cout << "[HMDYNAMICS]: Network not loaded, execution aborted" << endl;
            return EXIT_FAILURE;
        }

        output.open(filename);
        simulate_single(net, output);
        output.close();
    #elif MODE==DIAGRAM
        
        double si, sf;
        double ai, af;
        int ns, na;
        bool variable_a;

        if (argc == 11)
        {
            w0    = stod(argv[1]);
            delta = stod(argv[2]);
            q     = stod(argv[3]);
            variable_a = stoi(argv[4]);
            if (!variable_a)
            {
                a0     = stod(argv[5]);
                si    = stod(argv[6]);
                sf    = stod(argv[7]);
                ns    = stoi(argv[8]);
            }
            else 
            {
                s     = stod(argv[5]);
                ai    = stod(argv[6]);
                af    = stod(argv[7]);
                na    = stoi(argv[8]);
            }
            networkname = string(argv[9]);
            filename = string(argv[10]); 
        }
        else
        {
            cout << "[HMDYNAMICS]: incorrect number of arguments" << endl;
            return EXIT_SUCCESS;        
        }

        //Set up the network, including initial conditions
        correct_setup = set_up_network(net, networkname);

        if (!correct_setup) 
        {
            cout << "[HMDYNAMICS]: Network not loaded, execution aborted" << endl;
            return EXIT_SUCCESS;
        }
        
        if (variable_a) simulate_diagram(net, a0, ai, af, na, output);
        else simulate_diagram(net, s, si, sf, ns, output);
//        simulate_diagram(net, s0, sf, ns, output);
    #elif MODE==TIMESERIES
        double t;
        if (argc == 10)
        {
            w0    = stod(argv[1]);
            delta = stod(argv[2]);
            s     = stod(argv[3]);
            a0    = stod(argv[4]);
            q     = stod(argv[5]);
            tf    = stod(argv[6]);
            trelax= stod(argv[7]);
            networkname = string(argv[8]);
            filename = string(argv[9]); 
        }
        else
        {
            cout << "[HMDYNAMICS]: incorrect number of arguments" << endl;
            return EXIT_SUCCESS;        
        }

        //Set up the network, including initial conditions
        correct_setup = set_up_network(net, networkname);

        if (!correct_setup) 
        {
            cout << "[HMDYNAMICS]: Network not loaded, execution aborted" << endl;
            return EXIT_SUCCESS;
        }

        cout << "a" << endl;

        initial_conditions(net);

        cout << "a" << endl;
        for (t=0; t < trelax; t += dt) step_no_measures(net);
        output.open(filename);
        for (t=0; t < tf; t += dt)
        {
            step(net);
            output << r << " " << 1.0 + sin(psi) << endl;
        }
        output.close();

    #elif MODE==TIMETRACE
        int ntraces, trace_duration;

        if (argc == 11)
        {
            w0    = stod(argv[1]);
            delta = stod(argv[2]);
            s     = stod(argv[3]);
            a0    = stod(argv[4]);
            q     = stod(argv[5]);
            ntraces = stoi(argv[6]);
            trace_duration = stod(argv[7]);
            tf    = stod(argv[8]);
            networkname = string(argv[9]);
            filename = string(argv[10]); 
        }
        else
        {
            cout << "[HMDYNAMICS]: incorrect number of arguments" << endl;
            return EXIT_SUCCESS;        
        }

        //Set up the network, including initial conditions
        correct_setup = set_up_network(net, networkname);

        if (!correct_setup) 
        {
            cout << "[HMDYNAMICS]: Network not loaded, execution aborted" << endl;
            return EXIT_SUCCESS;
        }

        time_traces(net, output, ntraces, trace_duration);
    #elif MODE==CHIMERA
        int n_moduli;

        if (argc == 9)
        {
            w0    = stod(argv[1]);
            delta = stod(argv[2]);
            s     = stod(argv[3]);
            a0    = stod(argv[4]);
            q     = stod(argv[5]);
            osc_per_modulus = stoi(argv[6]);
            networkname = string(argv[7]);
            filename = string(argv[8]); 
        }
        else
        {
            cout << "[HMDYNAMICS]: incorrect number of arguments" << endl;
            return EXIT_SUCCESS;        
        }

        //Set up the network, including initial conditions
        correct_setup = set_up_network(net, networkname);
        n_moduli = net.get_node_count() / osc_per_modulus;

        if (!correct_setup) 
        {
            cout << "[HMDYNAMICS]: Network not loaded, execution aborted" << endl;
            return EXIT_SUCCESS;
        }

        //simulate_single_chimera(net, output, n_moduli, osc_per_modulus);
        time_traces_chimera(net, output, n_moduli);

    #endif
    return EXIT_SUCCESS;
}


// --- Core functions --- //
// These function do all the core hard work

bool set_up_network(CNetwork<double> &net, const string path_to_network)
{
    bool network_loaded_ok = false;


    //Read network from file
    network_loaded_ok = net.read_mtx(path_to_network);
    if (!network_loaded_ok) return false;

    N = net.get_node_count();

    return network_loaded_ok;
}

void initial_conditions(CNetwork<double> &net)
{
    int i;
    double x,y;
    int n_moduli;
    
    //gaussian = normal_distribution<double>(w0, delta);
    normal_distribution<double> gaussian1 = normal_distribution<double>(w0, delta);
    normal_distribution<double> gaussian2 = normal_distribution<double>(-w0, delta);

    n_moduli = N / osc_per_modulus;
    vector<double> wmoduli = vector<double>(n_moduli);  

    for (i=0; i < n_moduli; i++)
    {
        //wmoduli[i] = gaussian(gen);
        wmoduli[i] = ran_u(gen) < 0.5 ? gaussian1(gen) : gaussian2(gen);
    }

    //Set initial conditions
    x = y = 0.0;
    w = vector<double>(N);
    a = vector<double>(N);
    for (i=0; i < N; i++) 
    {
        net[i] = 2.0 * M_2_PI * ran_u(gen);
        w[i] = wmoduli[i / osc_per_modulus]; 
        a[i] = a0*w[i]/w0;

        x += cos(net[i]);
        y += sin(net[i]);
    }

    z = complex<double>(x,y) / (1.0*N);
    r = abs(z);
    psi = arg(z);
    
    return;
}

void step(CNetwork<double> &net)
{
    int i,j,neigh; 

    double x,y;
    
    double coupling;
    z = complex<double>(0.0, 0.0);
    x = y = 0.0;

    //Copy the state to avoid overwriting
    vector<double> old_state = net.get_values();

    for (i=0; i < N; i++)
    {
        coupling = 0.0;
        for (j=0; j < net.degree(i); j++)
        {
            neigh = net.get_in(i, j);
            coupling += sin(old_state[neigh] - old_state[i]);
        }

        coupling *= q / net.degree(i);
        net[i] += dt * (w[i] + a[i]*sin(old_state[i]) + coupling) + sqdt * ran_g(gen) * s;

        x += cos(net[i]);
        y += sin(net[i]);
    }

    z = complex<double>(x,y) / (1.0*N);
    r = abs(z);
    psi = arg(z);
}

void step_no_measures(CNetwork<double> &net)
{
    int i,j,neigh; 

    double coupling;

    //Copy the state to avoid overwriting
    vector<double> old_state = net.get_values();

    for (i=0; i < N; i++)
    {
        coupling = 0.0;
        for (j=0; j < net.degree(i); j++)
        {
            neigh = net.get_in(i, j);
            coupling += sin(old_state[neigh] - old_state[i]);
        }

        coupling *= q / net.degree(i);
        net[i] += dt * (w[i] + a[i]*sin(old_state[i]) + coupling) + sqdt * ran_g(gen) * s;

    }
}

void step_chimera(CNetwork<double> &net, const int n_moduli, const int osc_per_modulus)
{
    int i,j,neigh; 
    double coupling;

    z_modulus = vector<complex<double>>(n_moduli, complex<double>(0.0,0.0));


    //Copy the state to avoid overwriting
    vector<double> old_state = net.get_values();

    for (i=0; i < N; i++)
    {
        coupling = 0.0;
        for (j=0; j < net.degree(i); j++)
        {
            neigh = net.get_in(i, j);
            coupling += sin(old_state[neigh] - old_state[i]);
        }

        coupling *= q / net.degree(i);
        net[i] += dt * (w[i] + a[i]*sin(old_state[i]) + coupling) + sqdt * ran_g(gen) * s;

        z_modulus[i / osc_per_modulus] += complex<double>(cos(net[i]),sin(net[i])); 
    }


}



// --- Functions for phase diagrams and so on

void simulate_single(CNetwork<double> &net, ofstream &output, const double control)
{
    const int measure_its = tmeasure / dt;   //Number of iterations to do between measurements
    int nmeasures;                           //Number of measures we did
    double t;

    int i;

    //Initial conditions
    initial_conditions(net);

    //Variables to make averages
    double av_r  = 0.0;    
    double av_r2 = 0.0;

    //Then make simulations. First relaxation, then measurement.
    for (t = 0.0; t <= trelax; t += dt) step_no_measures(net);
    t = 0.0;
    nmeasures = 0;
    while(t < tf)
    {
        //Allow a bit of time between measures
        for (i = 0; i < measure_its - 1; i++)
        {
            step_no_measures(net);
            t += dt;
        }

        //Last step and measure
        step(net);
        av_r  += r;
        av_r2 += r*r;
        nmeasures++;
        t += dt;
    }

    av_r /= 1.0 * nmeasures;
    av_r2 /= 1.0 * nmeasures;
    output << control << " " << av_r << " " << av_r2 - av_r*av_r << endl;
    
    return;
}

//Simulate a single output, storing the Kuramoto order parameter of each moduli directly
void simulate_single_chimera(CNetwork<double> &net, ofstream &output, const int n_moduli, const int osc_per_modulus)
{
    const int measure_its = tmeasure / dt;   //Number of iterations to do between measurements
    double t;

    int i,j;

    //Initial conditions
    initial_conditions(net);

    //Then make simulations. First relaxation, then measurement.
    for (t = 0.0; t <= trelax; t += dt) step_no_measures(net);

    t = 0.0;
    output.open(filename);
    while(t < tf)
    {
        //Allow a bit of time between measures
        for (i = 0; i < measure_its - 1; i++)
        {
            step_no_measures(net);
            t += dt;
        }

        //Measure for a short time window
        for (i=0; i < nitswindow; i++)
        {
            step_chimera(net, n_moduli, osc_per_modulus);
            t += dt;
            output << t << " ";
            for (j=0; j < n_moduli; j++)
            {
                //output << real(z_modulus[j])/osc_per_modulus << " " << imag(z_modulus[j]) / osc_per_modulus << " ";
                output << abs(z_modulus[j]) / osc_per_modulus << " ";
            }
            output << endl;
        }
    }
    output.close();

    return;
}

//Generate a phase diagram for the selected parameter
void simulate_diagram(CNetwork<double> &net, double &selected_var, const double var0, const double varf, const int nvar, ofstream &output)
{
    const double dvar = (varf - var0) / (1.0 * nvar);
    output.open(filename);
    for (selected_var=var0; selected_var < varf; selected_var += dvar)
    {
        cout << selected_var << endl;
        simulate_single(net, output, selected_var);
    }
    output.close();
}

//Save time traces
void time_traces(CNetwork<double> &net, ofstream &output, const int ntraces, const double duration)
{
    int i,trace;
    double t;

    //Initial conditions
    initial_conditions(net);

    //First relaxation, then measurement.
    for (t = 0.0; t < trelax; t += dt) step(net);

    for (trace=0; trace < ntraces; trace++)
    {
        //Write to file part of the system evolution
        output.open(filename + "-trace" + to_string(trace));
        for (t = 0.0; t < duration; t += dt)
        {
            step(net);
            for (i=0; i < N; i++) output << net[i] << " ";
            output << endl;
        }
        output.close();

        //Wait a bit between traces
        for (t = 0.0; t < tf; t += dt) step(net);
    } 
}

void time_traces_chimera(CNetwork<double> &net, ofstream &output, const int n_moduli)
{
    int i,j;
    double t;

    //Initial conditions
    initial_conditions(net);

    //First relaxation, then measurement.
    for (t = 0.0; t < trelax; t += dt) step(net);

    t = 0.0;
    output.open(filename);
    while(t < tf)
    {
        step_chimera(net, n_moduli, osc_per_modulus);
        t += dt;
        for (j=0; j < n_moduli; j++)
        {
            output << real(z_modulus[j])/osc_per_modulus << " " << imag(z_modulus[j]) / osc_per_modulus << " ";
        }
        output << endl;
    }
    output.close();
}