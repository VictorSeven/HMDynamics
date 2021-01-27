/* -------------------------------------------------------- *
*  This code generates hierarchical-modular networks with a *
*  core-periphery structure, based on a work by             *
*  Zamora-Lopez et al., 2016, SciRep.                       *
*  Implementation: Victor Buendía, 2021                     *
* --------------------------------------------------------- */

//TODO: no medir el promedio a cada paso temporal, tal vez cada 5 o 10. Asegurar el criterio de Nyquist
//Usar desenrrollado de bucles para la funcion step

#define MAXNETSIZE 10000

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

void step_relaxation(CNetwork<double> &net);
void step(CNetwork<double> &net);


// --- Random number generation --- //
mt19937 gen(85345385434);
uniform_real_distribution<double> ran_u(0.0, 1.0);

// --- Global constants --- //
int N;

double w,q,s;

const double dt = 0.01;
const double sqdt = sqrt(0.01);


const complex<double> I = complex<double>(0.0, 1.0);

complex<double> z;
double r, psi;

// --- Main --- ///
// This is used to construct and export the desired network to a file

int main(int argc, char* argv[])
{
    //Counters
    int i,j;

    //Time-related variables
    double t;

    const double trelax = 100.0;            //Relaxation time
    const double tf = 100.0;                //Simulation time
    const double tmeasure = 1.0;            //Time between measurements
    const int measure_its = tmeasure / dt;   //Number of iterations to do between measurements

    //Variables to make averages
    double av_r  = 0.0;    
    double av_r2 = 0.0;

    //Define the network and the path to it
    CNetwork<double> net(MAXNETSIZE);
    string filename = "hmrandom";

    //Input/output stuff
    bool correct_setup;

    ofstream output;

    //Set up the network, including initial conditions
    correct_setup = set_up_network(net, filename);
    if (!correct_setup) 
    {
        cout << "[HMDYNAMICS]: Network not loaded, execution aborted" << endl;
        return EXIT_FAILURE;
    }

    //Then make simulations. First relaxation, then measurement.
    for (t = 0.0; t <= trelax; t += dt) step_relaxation(net);

    t = 0.0;
    while(t < tf)
    {
        //Allow a bit of time between measures
        for (i = 0; i < measure_its - 1; i++)
        {
            step_relaxation(net);
            t += dt;
        }

        //Last step and measure
        step(net);
        av_r  += r;
        av_r2 += r*r;
        t += dt;
    }

    output.open(filename);
    output << q << " " << av_r << " " << av_r2 << endl;
    output.close();

    return EXIT_SUCCESS;
}

bool set_up_network(CNetwork<double> &net, const string filename)
{
    int i;
    bool network_loaded_ok = false;

    double x,y;

    //int L;

    //Read network from file
    network_loaded_ok = net.read_mtx(filename);
    if (!network_loaded_ok) return false;
    
    N = net.get_node_count();
    //L = net.get_link_count();

    //for (i=0; i < L; i++) net.get_link(i) = q; 

    //Set initial conditions
    z = complex<double>(0.0, 0.0);
    x = y = 0.0;
    for (i=0; i < N; i++) 
    {
        net[i] = 2.0 * M_2_PI * ran_u(gen);

        x += cos(net[i]);
        y += sin(net[i]);
    }

    z = complex<double>(x,y) / N;
    r = abs(z);
    psi = arg(z);
    
    return true;
}


void step_relaxation(CNetwork<double> &net)
{
    int i; 

    double x,y;
    
    z = complex<double>(0.0, 0.0);
    x = y = 0.0;

    for (i=0; i < N; i++)
    {
        net[i] = dt * (w + q * r * sin(psi - net[i])) + sqdt * ran_g(gen) * s;

        x += cos(net[i]);
        y += sin(net[i]);
    }

    z = complex<double>(x,y) / N;
    r = abs(z);
    psi = arg(z);
}


void step(CNetwork<double> &net)
{
    int i; 

    double x,y;
    
    z = complex<double>(0.0, 0.0);
    x = y = 0.0;

    for (i=0; i < N; i++)
    {
        net[i] = dt * (w + q * r * sin(psi - net[i])) + sqdt * ran_g(gen) * s;

        x += cos(net[i]);
        y += sin(net[i]);
    }

    z = complex<double>(x,y) / N;
    r = abs(z);
    psi = arg(z);
}