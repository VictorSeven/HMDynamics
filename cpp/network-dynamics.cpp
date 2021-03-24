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

#ifndef MODE
#define MODE DIAGRAM 
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

void step(CNetwork<double> &net);

void simulate_single(CNetwork<double> &net, ofstream &output);
void simulate_diagram(CNetwork<double> &net, const double s0, const double sf, const int ns, ofstream &output);
void simulate_diagram(CNetwork<double> &net, double &selected_var, const double var0, const double varf, const int nvar, ofstream &output);
void time_traces(CNetwork<double> &net, ofstream &output, const int ntraces, const double duration);

// --- Random number generation --- //
mt19937 gen(85345385434);
uniform_real_distribution<double> ran_u(0.0, 1.0);
normal_distribution<double> ran_g(0.0, 1.0);
cauchy_distribution<double> lorentzian(0.0, 1.0);

// --- Global constants --- //
int N;

vector<double> w;
double w0,delta,q,a,s;
string networkname = "hmrandom";
string filename = "kuramoto";


const double dt = 0.01;
const double sqdt = sqrt(dt);

double tf = 1000.0;
double trelax = 100.0;
const double tmeasure = 1.0;

const complex<double> I = complex<double>(0.0, 1.0);

complex<double> z;
double r, psi;

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
            a     = stod(argv[4]);
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
        
        double s0, sf;
        double a0, af;
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
                a     = stod(argv[5]);
                s0    = stod(argv[6]);
                sf    = stod(argv[7]);
                ns    = stoi(argv[8]);
            }
            else 
            {
                s     = stod(argv[5]);
                cout << s << endl;
                a0    = stod(argv[6]);
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
        
        if (variable_a) simulate_diagram(net, a, a0, af, na, output);
        else simulate_diagram(net, s, s0, sf, ns, output);
//        simulate_diagram(net, s0, sf, ns, output);
    #elif MODE==TIMETRACE
        int ntraces, trace_duration;

        if (argc == 11)
        {
            w0    = stod(argv[1]);
            delta = stod(argv[2]);
            s     = stod(argv[3]);
            a     = stod(argv[4]);
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
        
    #endif
    return EXIT_SUCCESS;
}


// --- Core functions --- //
// These function do all the core hard work

bool set_up_network(CNetwork<double> &net, const string path_to_network)
{
    int i;
    bool network_loaded_ok = false;

    double x,y;

    //Read network from file
    network_loaded_ok = net.read_mtx(path_to_network);
    if (!network_loaded_ok) return false;
    
    N = net.get_node_count();
    
    //Link stuff
    //L = net.get_link_count();
    //for (i=0; i < L; i++) net.get_link(i) = q; 

    //Distribution of frequencies
    lorentzian = cauchy_distribution<double>(w0, delta);

    //Set initial conditions
    x = y = 0.0;
    w = vector<double>(N);
    for (i=0; i < N; i++) 
    {
        net[i] = 2.0 * M_2_PI * ran_u(gen);
        w[i] = lorentzian(gen);

        x += cos(net[i]);
        y += sin(net[i]);
    }

    z = complex<double>(x,y) / (1.0*N);
    r = abs(z);
    psi = arg(z);
    
    return true;
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
        net[i] += dt * (w[i] + a*sin(old_state[i]) + coupling) + sqdt * ran_g(gen) * s;

        x += cos(net[i]);
        y += sin(net[i]);
    }

    z = complex<double>(x,y) / (1.0*N);
    r = abs(z);
    psi = arg(z);
}

// --- Functions for phase diagrams and so on

void simulate_single(CNetwork<double> &net, ofstream &output)
{
    cout << a << " " << s << endl;
    const int measure_its = tmeasure / dt;   //Number of iterations to do between measurements
    const int nmeasures = tf / tmeasure;     //Number of measures we did
    double t;

    int i;

    //Variables to make averages
    double av_r  = 0.0;    
    double av_r2 = 0.0;

    //Then make simulations. First relaxation, then measurement.
    for (t = 0.0; t <= trelax; t += dt) step(net);

    t = 0.0;
    while(t < tf)
    {
        //Allow a bit of time between measures
        for (i = 0; i < measure_its - 1; i++)
        {
            step(net);
            t += dt;
        }

        //Last step and measure
        step(net);
        av_r  += r;
        av_r2 += r*r;
        t += dt;
    }

    av_r /= 1.0 * nmeasures;
    av_r2 /= 1.0 * nmeasures;
    output << q << " " << av_r << " " << av_r2 - av_r*av_r << endl;
    
    return;
}

void simulate_diagram(CNetwork<double> &net, const double s0, const double sf, const int ns, ofstream &output)
{
    const double ds = (sf - s0) / (1.0 * ns);
    cout << ds << " " << filename << endl;
    output.open(filename);
    for (s=s0; s < sf; s += ds)
    {
        simulate_single(net, output);
    }
    output.close();
}

void simulate_diagram(CNetwork<double> &net, double &selected_var, const double var0, const double varf, const int nvar, ofstream &output)
{
    const double dvar = (varf - var0) / (1.0 * nvar);
    output.open(filename);
    for (selected_var=var0; selected_var < varf; selected_var += dvar)
    {
        simulate_single(net, output);
    }
    output.close();
}

void time_traces(CNetwork<double> &net, ofstream &output, const int ntraces, const double duration)
{
    int i,trace;
    double t;

    //First relaxation, then measurement.
    for (t = 0.0; t < trelax; t += dt) step(net);

    for (trace=0; trace < ntraces; trace++)
    {
        //Write to file part of the system evolution
        output.open(filename + "_trace" + to_string(trace));
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