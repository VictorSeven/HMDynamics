/* -------------------------------------------------------- *
*  This code generates hierarchical-modular networks with a *
*  core-periphery structure, based on a work by             *
*  Zamora-Lopez et al., 2016, SciRep.                       *
*  Implementation: Victor Buendía, 2021                     *
* --------------------------------------------------------- */

#define MAXNETSIZE 10000

//Include libraries
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <random>
#include <vector>
#include "CNetwork/source/CNet.cpp" //Route local to this file. Find it on github.com/VictorSeven/CNetwork

using namespace std;

// --- Declare functions --- //
void set_up_network(CNetwork<double> &net, const string filename)


// --- Random number generation --- //
mt19937 gen(85345385434);
uniform_real_distribution<double> ran_u(0.0, 1.0);

// --- Global constants --- //
int N;


// --- Main --- ///
// This is used to construct and export the desired network to a file

int main(int argc, char* argv[])
{
    CNetwork<double> &net(MAXNETSIZE);
    string filename = "prueba";


    set_up_network(net, filename);

    return EXIT_SUCCESS;
}

void set_up_network(CNetwork<double> &net, const string filename)
{
    int i;

    //Read network from file
    net.read_mtx(filename);
    N = net.get_node_count();

    //Set initial conditions
    for (i=0; i < N; i++) net[i] = 2.0 * M_2_PI * ran_u(gen);
}
