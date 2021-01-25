/* -------------------------------------------------------- *
*  This code generates hierarchical-modular networks with a *
*  core-periphery structure, based on a work by             *
*  Zamora-Lopez et al., 2016, SciRep.                       *
*  Implementation: Victor Buendía, 2021                     *
* --------------------------------------------------------- */

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
void show_help();

int same_level(const int index_i, const int index_j, const int nlevels, const vector<int> &neurons_level);
void hm_random(CNetwork<> &net, const int nlevels, const vector<int> &levels, const vector<double> &k_level);

void get_node_cumulative(const int neurons_cluster, const double beta, vector<double> &cumuweight);
int get_random_node(const int cluster, const int neurons_cluster, const vector<double> &cumuweight);
void hm_core(CNetwork<> &net, const int nlevels, const vector<int> &levels, const vector<double> &k_level, const vector<double> &gamma_level);

// --- Random number generation --- //
mt19937 gen(85345385434);
uniform_real_distribution<double> ran_u(0.0, 1.0);

// --- Global constants --- //
int N;


// --- Main --- ///
// This is used to construct and export the desired network to a file

int main(int argc, char* argv[])
{
    int i;

    int nlevels = 0;
    bool generate_random = true;
    int nargs;

    CNetwork<> net;

    vector<int> levels;
    vector<double> klevel;
    vector<double> glevel;

    if (argv[1] == "-help" || argv[1] == "--help" || argv[1] == "-h" || argv[1] == "--h") 
    {
        show_help();
        return EXIT_SUCCESS;
    }

    generate_random = !(string(argv[1]) == "--core");
    nlevels = stoi(argv[2]);

    if (generate_random)
    {
        nargs = 3 + nlevels * 2;
        if (nargs != argc)
        {
            cout << nargs << " ! " << argc << endl;
            cout << "Wrong number of arguments. Correct format:" << endl;
            show_help();
            return EXIT_SUCCESS;
        }

        levels = vector<int>(nlevels);
        klevel = vector<double>(nlevels);

        N=1;
        for (i=0; i < nlevels; i++)
        {
            levels[i] = stoi(argv[3 + i]);
            N *= levels[i];
            klevel[i] = stod(argv[3 + nlevels + i]);
        }

        net = CNetwork<>(N);
        net.add_nodes(N);

        hm_random(net, nlevels, levels, klevel);
    }
    else
    {
        nargs = 3 + nlevels * 3;
        if (nargs != argc)
        {
            cout << "Wrong number of arguments. Correct format:" << endl;
            show_help();
            return EXIT_SUCCESS;
        }

        levels = vector<int>(nlevels);
        klevel = vector<double>(nlevels);
        glevel = vector<double>(nlevels);

        N=1;
        for (i=0; i < nlevels; i++)
        {
            levels[i] = stoi(argv[3 + i]);
            N *= levels[i];
            klevel[i] = stod(argv[3 + nlevels + i]);
            glevel[i] = stod(argv[3 + nlevels*2 + i]);
        }
        net = CNetwork<>(N);
        net.add_nodes(N);

        hm_core(net, nlevels, levels, klevel, glevel);
    }

    return EXIT_SUCCESS;
}

void show_help()
{
    cout << "create-network [mode] [M] [n1 n2 ... nM] [k1 k2 ... kM] [g1 g2 ... gM]" << endl << endl;
    cout << "[mode]: either --random or --core. Selects the network to construct." << endl;
    cout << "[M]: number of hierarchical levels. Integer." << endl;
    cout << "[n1 n2 ... nM]: number of elements in each hierarchical level. Integers." << endl;
    cout << "[k1 k2 ... kM]: connectivity in each hierarchical level. Floats." << endl << endl;
    cout << "[g1 g2 ... gM]: Only in core mode. Scale-free exponent in each hierarchical level. Floats." << endl << endl;
    cout << "Examples:" << endl;
    cout  << "create-network --random 3   40 20 5   10.0 2.0 1.5" << endl;
    cout <<  "create-network --core 2  40 20   10.0 2.0   0.0 1.5" << endl << endl;
    cout << "Please read documentation for more info on parameters and its behaviour." << endl;
}



// --- Functions --- ///
// Implementation of the network creation itself. Please read documentation notes for implementation details


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

//Random modular graph
void hm_random(CNetwork<> &net, const int nlevels, const vector<int> &levels, const vector<double> &k_level)
{
    int i,j;

    //First of all get number of hierarchical levels
    vector<int> n_neurons_level(nlevels);       //Number of neurons in each level (lower level modules are neurons). 
    vector<double> p_level(nlevels);            //Probability of connection for each level. 

    //Aux variable to store at which level neurons coincide
    int coincidence_level;

    //Store total number of neurons up to level k and connection probability 
    //First level is just the number of neurons inside the cluster
    n_neurons_level[0] = levels[0];
    p_level[0] = k_level[0] / (1.0 * n_neurons_level[i]);
    for (i=1; i < nlevels; i++)
    {
        n_neurons_level[i] = levels[i] * n_neurons_level[i-1];
        p_level[i] = k_level[i] / (1.0 * n_neurons_level[i]); //Connection probability
    }
    
    //Erdos-Renyi-like construction with the adequate prob.
    for (i=0; i < N; i++)
    {
        for (j=i+1; j < N; j++)
        {
            //Get at which level both neurons coincide. 
            //Higher level is the whole network, with index nlevels-1
            coincidence_level = same_level(i, j, nlevels, n_neurons_level);

            //Add a link with the probability corresponding to that level
            if (ran_u(gen) < p_level[coincidence_level])
            {
                net.add_link(i, j);
            }
        }   
    }

    //Record the used graph
    net.write_graphml("hmrandom", vector<string>());
    net.write_mtx("hmrandom");

}

//Core-periphery

//Auxiliary function to get the cumulative probabilities for each neuron at each hierarchical level
void get_node_cumulative(const int neurons_cluster, const double beta, vector<double> &cumuweight)
{
    int i;

    double normalization;

    if (beta < 0.0)
    {
        vector<double> weights(neurons_cluster);
        cumuweight = vector<double>(neurons_cluster + 1, 0.0); //+1 to start the cumulative at 0

        normalization = 0.0;

        //Store weights and compute the complete sum in order to compute later the probabilities
        for (i=0; i < neurons_cluster; i++) 
        {
            weights[i] = pow(i+1, beta);
            normalization += weights[i];
        }

        //Then directly get the cumulatives
        cumuweight[0] = 0.0;
        cumuweight[1] = weights[0] / normalization;
        for (i=2; i < neurons_cluster; i++) cumuweight[i] = cumuweight[i-1] + weights[i-1] / normalization;

    }
    else
    {
        //Random connectivity, all same weight
        double weight =  1.0 / neurons_cluster;
        cumuweight = vector<double>(neurons_cluster, 0.0);
        for (i=0; i < neurons_cluster; i++) cumuweight[i] = i * weight;
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
    return distance(cumuweight.begin(), lower) - 1 + cluster * neurons_cluster;
}

//Create the hierarchical network.
//Input arguments: nlevels {list_n_elements} {list_k}
void hm_core(CNetwork<> &net, const int nlevels, const vector<int> &levels, const vector<double> &k_level, const vector<double> &gamma_level)
{
    int i,j,level,link,cluster; //Counters
    int nclusters;

    //Auxiliary variables
    int nlinks;             //Number of links in this level
    int comm1, comm2;       //Indices of two randomly selected communities
    double beta;            //Exponent of this level

    vector<int> n_neurons_level(nlevels);       //Number of neurons in each level (lower level modules are neurons). 

    vector<double> cumulative;                  //Cumulative weight for obtaining random numbers from power-law distrib.

    //Fill the information provided 
    n_neurons_level[0] = levels[0];
    for (i=1; i < nlevels; i++)
    {
        n_neurons_level[i] = levels[i] * n_neurons_level[i-1]; 
    }
    
    //Now do the connections level by level
    for (level=0; level < nlevels; level++)
    {
        //Get the number of links adequate to this level, expected randomly, and set the power-law exponents
        nlinks = 0.5 * k_level[level] * n_neurons_level[level]; 
        beta = gamma_level[level] > 1.0 ? 1.0 / (gamma_level[level] - 1.0) : 0.0;

        //Compute cumulative weight of all the neurons inside each cluster
        //Cumulative is a vector of size "numbers of neurons in this level"
        get_node_cumulative(n_neurons_level[level], -beta, cumulative); 
 
        //Get the number of clusters at this hierarchical level
        nclusters = N / n_neurons_level[level];

        //Make connections inside those clusters
        for (cluster=0; cluster < nclusters; cluster++)
        {
            link = 0;
            //Set all necessary links...
            while (link <= nlinks)
            {
                //Select two neurons from the same community using the scale-free cumulative weight
                do
                {
                    i = get_random_node(cluster, n_neurons_level[level], cumulative);
                    j = get_random_node(cluster, n_neurons_level[level], cumulative);
                } while (i == j);

                if (i >= N || j >= N || i < 0 || j < 0)
                {
                    cout << i << " " << j << " " << cluster << " " << level << " " << n_neurons_level[level] << endl;
                    return;
                }

                //Connect them if they are not already
                if (net.get_link_index(i,j) == -1) 
                {
                    net.add_link(i,j);
                    link++;
                }
            }
        }

        //Now we have connected the neurons inside this level, connect the clusters between them in next level!
    }

    //Record the used graph
    net.write_graphml("hmcore", vector<string>());
    net.write_mtx("hmcore");



/*

    //Now construct the network. We do first the first level, connecting the individual neurons
    //Number of links to generate in this level (number of links in each cluster)
    //We use the number of links of the random graph as an estimator
    nlinks = 0.5 * k_level[0] * n_neurons_level[0]; 

    //Beta to use in this level. Values gamma < 1 lead to beta=0.0, i.e., random connectivity
    beta = gamma_level[0] > 1.0 ? 1.0 / (gamma_level[0] - 1.0) : 0.0; 

    //Fill the node weights for this level
    get_node_cumulative(n_neurons_level[0], -beta, cumulative); 

    //At the first level, each neuron works as a cluster itself
    //nclusters = N;
    //for(i=1; i < nlevels; i++) nclusters *= n_elements_level[i];
    //TODO CHECK: creo que este nclusters se puede obtener de forma mas simple
*/
    //Then set the links. This part will be different in higher levels. 
    //Each cluster will have nlinks, in order to have the selected k in this cluster!
    /*for (cluster=0; cluster < nclusters; cluster++)
    {
        link = 0;
        while (link < nlinks)
        {
            //Select two neurons from the same community using the scale-free cumulative weight
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
        get_node_cumulative(n_neurons_level[level], -beta, cumulative); //TODO CHECK first it was level-1

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
                if (n_elements_level[level] > 2)
                {
                    do
                    {
                        comm1 = int(ran_u(gen) * n_elements_level[level]) + cluster * n_elements_level[level];
                        comm2 = int(ran_u(gen) * n_elements_level[level]) + cluster * n_elements_level[level];
                    } while (comm1 == comm2);
                }
                else if (n_elements_level[level] == 2)
                {
                    comm1 = 0;
                    comm2 = 1;
                }
                else
                {
                    comm1 == comm2 = 0;
                }
                
                //And now obtain a neuron from those communities. 
                i = get_random_node(comm1, n_neurons_level[level-1], cumulative);
                j = get_random_node(comm2, n_neurons_level[level-1], cumulative);

                //Check that the link is yet non-existent. If everything is correct, then make a link
                if (net.get_link_index(i,j) == -1) 
                {
                    net.add_link(i,j);
                    link++;
                }
            }
        }
    }
*/
    //Record the used graph
    //net.write_graphml("hmcore", vector<string>());
    //net.write_mtx("hmcore");

}