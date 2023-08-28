#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;
//physical parameters
double m = 1.;
double omega= 1.;
double r = 0.0;
// Yoshida coefficients
double denominatore = 2.-cbrt(2.);
double w_0 = -(cbrt(2.))/denominatore;
double w_1 = 1./denominatore;
double c_1 = w_1/2.;
double c_2 = (w_0 + w_1)/2.;
double c_3 = c_2;
double c_4 = c_1;
double d_1 = w_1;
double d_2 = w_0;
double d_3 = w_1;

//simulation parameters
double h = 0.001;
int steps = 100000;
int num_samples = 10000;


double forza(double q){
    return -m*omega*omega*q-r*q*q*q;
}

double momento(double p){
    return p/m;
}

int main()
{
    // ensemble averaged trajectories
    vector<double> q_ensemble(steps, 0.0);
    vector<double> p_ensemble(steps, 0.0);
    // file to read initial conditions
    ifstream file_condizioni_iniziali_q("condizioni_iniziali_q.txt"); // Apri il file in modalità di lettura
    ifstream file_condizioni_iniziali_p("condizioni_iniziali_p.txt"); // Apri il file in modalità di lettura

    // initial conditions
    vector<double> condizione_iniziale_q;
    vector<double> condizione_iniziale_p;
    ///////////////////////////////////////////
    if (file_condizioni_iniziali_q.is_open()) {
        double valore;
        while (file_condizioni_iniziali_q >> valore) {
            condizione_iniziale_q.push_back(valore);
        }
        file_condizioni_iniziali_q.close();
    }
    else {
        cerr << "Impossibile aprire il file delle condizioni iniziali q." << endl;
    }
    ///////////////////////////////////////////
    ///////////////////////////////////////////
    if (file_condizioni_iniziali_p.is_open()) {
        double valore;
        while (file_condizioni_iniziali_p >> valore) {
            condizione_iniziale_p.push_back(valore);
        }
        file_condizioni_iniziali_p.close();
    }
    else {
        cerr << "Impossibile aprire il file delle condizioni iniziali p." << endl;
    }
    //////////////////////////////////
    // file to memorize data
    ofstream file_q("q_trajectory_ensemble.txt");
    ofstream file_p("p_trajectory_ensemble.txt");

    //N stands for the number of samples in the ensemble
    for (int n = 0; n < num_samples; n++){
        // variables where I will put initial conditions

        // vectors containing trajectories
        vector<double> q(steps, 0.0);
        vector<double> p(steps, 0.0);
        // reading initial conditions
        q[0] = condizione_iniziale_q[n];
        p[0] = condizione_iniziale_p[n];
        // Yoshida-4 algorithm
        for(int i = 0; i < steps-1; i++){
            // temporary variables
            // step-1
            double q_temp = q[i] + c_1*momento(p[i])*h;
            double p_temp = p[i] + d_1*forza(q_temp)*h;
            // step-2
            q_temp = q_temp + c_2*momento(p_temp)*h;
            p_temp = p_temp + d_2*forza(q_temp)*h;
            // step-3
            q_temp = q_temp + c_3*momento(p_temp)*h;
            p_temp = p_temp + d_3*forza(q_temp)*h;
            // step-4
            q[i+1] = q_temp + c_4*momento(p_temp)*h;
            p[i+1] = p_temp;
        }
        // update ensemble vectors
        for (int i = 0; i < steps; i++){
            q_ensemble[i] += q[i];
            p_ensemble[i] += p[i];
        }
    }



    // da riadattare per metterci dentro le medie di ensemble
    for(int i = 0; i < steps; i++){
        q_ensemble[i] /= num_samples;
        p_ensemble[i] /= num_samples;
    }


    for(int i = 0; i < steps; i++){

        if (file_q.is_open()) {
            file_q << q_ensemble[i] << endl; // Scrive la variabile nel file e va a capo
        }
        else {
            cerr << "Impossibile aprire il file delle q." << endl;
        }
        if (file_p.is_open()) {
            file_p << p_ensemble[i] << endl; // Scrive la variabile nel file e va a capo
        }
        else {
            cerr << "Impossibile aprire il file delle p." << endl;
        }
    }

    //closing files
    file_q.close();
    file_p.close();


    return 0;
}
