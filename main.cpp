#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

//physical parameters
double m = 1.;
double omega = 1.;
double r = 0.0001; // to vary depending on the degree of nonlinearity
double k_B = 1.;
double T = 50.;
double beta = 1/(k_B*T);
double tau = 1.0;


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
double h = 0.01;
int steps = 50000;
int num_samples = 100000; // this parameter has to be the same of the associated jupiter notebook


double forza(double q, double zeta, double p){
    //return -m*omega*omega*q-r*q*q*q;
    return -m*omega*omega*q - zeta*p;
}

double momento(double p){
    return p/m;
}

double controllo(double p){
    return ((p*p)/(m*k_B*T) - 1)*(1/tau);
}

int main()
{
    // 1. define vectors that will contain ensemble trajectories
    // ensemble averaged trajectories

    // 1.1 dynamical variables
    vector<double> q_ensemble(steps, 0.0);
    vector<double> p_ensemble(steps, 0.0);
    vector<double> zeta_ensemble(steps, 0.0);

    // 1.2 observables
    vector<double> q_square_ensemble(steps, 0.0);
    vector<double> p_square_ensemble(steps, 0.0);
    vector<double> omega_ensemble(steps, 0.0);

    // 1.3 correlations functions
    vector<double> corr_q_square(steps, 0.0);
    vector<double> corr_p_square(steps, 0.0);
    vector<double> corr_omega(steps, 0.0);

    // 1.4 integral omega
    vector <double> integral_omega_ensemble(steps, 0.0);


    // 2. define file to read initial conditions from
    // file to read initial conditions
    ifstream file_condizioni_iniziali_q("condizioni_iniziali_q.txt"); // Apri il file in modalità di lettura
    ifstream file_condizioni_iniziali_p("condizioni_iniziali_p.txt"); // Apri il file in modalità di lettura

    // 2.1 define variable to store the i-th sample of the initial conditions
    // initial conditions
    vector<double> condizione_iniziale_q;
    vector<double> condizione_iniziale_p;

    // 2.2 fill the array of q initial conditions
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

    // 2.3 fill the array of p initial conditions
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

    cout << "Condizioni iniziali lette." << endl;

    // 2.4 check whether initial conditions are properly initialized
    //for(int i = 0; i < num_samples; i++){
        //cout << i << "-th position initial condition: " << condizione_iniziale_q[i] << endl;
        //cout << i << "-th momentum initial condition: " << condizione_iniziale_p[i] << endl;
        //cout << endl;
        //}


    // 3. define the files to store ensemble averaged data
    // file to memorize data

    // 3.1 file to memorize dynamical variables
    ofstream file_q("q_trajectory_ensemble.txt");
    ofstream file_p("p_trajectory_ensemble.txt");
    ofstream file_zeta("zeta_trajectory_ensemble.txt");

    // 3.2 file to memorize observables
    ofstream file_q_square("q_square_trajectory_ensemble.txt");
    ofstream file_p_square("p_square_trajectory_ensemble.txt");
    ofstream file_omega("omega_trajectory_ensemble.txt");

    // 3.3 file to memorize correlation functions
    ofstream file_corr_q_square("corr_q_square_ensemble.txt");
    ofstream file_corr_p_square("corr_p_square_ensemble.txt");
    ofstream file_corr_omega("corr_omega_ensemble.txt");

    // 3.4 file to memorize the integral omega
    ofstream file_integral_omega("integral_omega_ensemble.txt");


    // 4. start the ensemble simulation
    cout << "Iniziata ensemble simulation" << endl;
    for (int n = 0; n < num_samples; n++){

        // 4.1 for each cycle I define a fresh array to store the n-th realization trajectory
        // vectors containing trajectories
        vector<double> q(steps, 0.0);
        vector<double> p(steps, 0.0);
        vector<double> zeta(steps, 0.0);

        // vectors containing the elements to compute the integral
        vector<double> rettangoli(steps, 0.0);
        vector<double> integral_omega(steps, 0.0);

        // 4.2 store the n-th initial condition in the first element of the array
        // reading initial conditions
        q[0] = condizione_iniziale_q[n];
        p[0] = condizione_iniziale_p[n];
        zeta[0] = 1.;
        double omega_zero = zeta[0]*(1-(p[0]*p[0])/(m*k_B*T));

        // 4.3 Runge-Kutta-4 algorithm
        for(int i = 0; i < steps-1; i++){
            // define temporary k-variables to advance in time
            double kappa[3][4]; // 3 variabili dinamiche - 4 step RK

            // 4.3.1 step-1
            kappa[0][0] = momento(p[i]);
            kappa[1][0] = forza(q[i], zeta[i], p[i]);
            kappa[2][0] = controllo(p[i]);

            // 4.3.2 step-2
            kappa[0][1] = momento(p[i]+(h/2)*kappa[0][0]);
            kappa[1][1] = forza(q[i]+(h/2)*kappa[1][0], zeta[i]+(h/2)*kappa[1][0], p[i]+(h/2)*kappa[1][0]);
            kappa[2][1] = controllo(p[i]+(h/2)*kappa[2][0]);

            // 4.3.3 step-3
            kappa[0][2] = momento(p[i] + (h/2)*kappa[0][1]);
            kappa[1][2] = forza(q[i]+(h/2)*kappa[1][1], zeta[i]+(h/2)*kappa[1][1], p[i]+(h/2)*kappa[1][1]);
            kappa[2][2]= controllo(p[i]+(h/2)*kappa[2][1]);

            // 4.3.4 step-4
            kappa[0][3] = momento(p[i]+ h*kappa[0][2]);
            kappa[1][3] = forza(q[i]+h*kappa[1][2], zeta[i]+h*kappa[1][2], p[i]+h*kappa[1][2]);
            kappa[2][3] = controllo(p[i]+h*kappa[2][2]);

            // 4.3.5 new trajectory point
            q[i+1] = q[i]+(h/6)*(kappa[0][0] + 2*kappa[0][1] + 2*kappa[0][2] + kappa[0][3]);
            p[i+1] = p[i]+(h/6)*(kappa[1][0] + 2*kappa[1][1] + 2*kappa[1][2] + kappa[1][3]);
            zeta[i+1] = zeta[i]+(h/6)*(kappa[2][0] + 2*kappa[2][1] + 2*kappa[2][2] + kappa[2][3]);


        }


        // 4.4 update ensemble vectors after a whole evolution has occurred
        for (int i = 0; i < steps; i++){

            // 4.4.1 dynamical variables
            q_ensemble[i] += q[i];
            p_ensemble[i] += p[i];
            zeta_ensemble[i] += zeta[i];

            // 4.4.2 observables
            q_square_ensemble[i] += q[i]*q[i];
            p_square_ensemble[i] += p[i]*p[i];
            omega_ensemble[i] += zeta[i]*(1-(p[i]*p[i])/(m*k_B*T));

            // 4.4.3 correlation functions
            corr_q_square[i] += omega_zero*q[i]*q[i];
            corr_p_square[i] += omega_zero*p[i]*p[i];
            corr_omega[i] += omega_zero*zeta[i]*(1-(p[i]*p[i])/(m*k_B*T));


        }

        // 4.5 calculate the integral omega
        // 4.5.1. calculate the rectangles
        for(int i = 0; i < steps; i++){
            rettangoli[i] = h*zeta[i]*(1-(p[i]*p[i])/(m*k_B*T));
        }
        //4.5.2 calculate the integral function
        integral_omega[0] = rettangoli[0];
        for(int i = 0; i < steps-1; i++){
            integral_omega[i+1] = integral_omega[i] + rettangoli[i];
        }
        // 4.5.3 update the ensemble integral function
        for(int i = 0; i < steps; i++){
            integral_omega_ensemble[i] += integral_omega[i];
        }
    }
    // 4.6 close the ensemble simulation
    cout << "Finita ensemble simulation" << endl;



    // 5. Normalize the ensemble averaged trajectories
    for(int i = 0; i < steps; i++){
        // 5.1 dynamical variables
        q_ensemble[i] /= num_samples;
        p_ensemble[i] /= num_samples;
        zeta_ensemble[i] /= num_samples;

        // 5.2 observables
        q_square_ensemble[i] /= num_samples;
        p_square_ensemble[i] /= num_samples;
        omega_ensemble[i] /= num_samples;

        // 5.3 correlation functions
        corr_q_square[i] /= num_samples;
        corr_p_square[i] /= num_samples;
        corr_omega[i] /= num_samples;

        // 5.4 integral omega function
        integral_omega_ensemble[i] /= num_samples;
    }

    cout << "Inizio memorizzazione dati" << endl;
    // 6. Memorize the ensemble averaged trajectories
    for(int i = 0; i < steps; i++){


        // 6.1 storing dynamical variables
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

        if (file_zeta.is_open()){
            file_zeta << zeta_ensemble[i] << endl;
        }
        else {
            cerr << "Impossibile aprire il file delle zeta." << endl;
        }
        ////////////////////////////////////////////

        // 6.2 storing observables
        if (file_q_square.is_open()){
            file_q_square << q_square_ensemble[i] << endl;
        }
        else {
            cerr << "Impossibile aprire il file delle q quadrate." << endl;
        }
        if (file_p_square.is_open()){
            file_p_square << p_square_ensemble[i] << endl;
        }
        else {
            cerr << "Impossibile aprire il file delle p quadrate." << endl;
        }
        if (file_omega.is_open()){
            file_omega << omega_ensemble[i] << endl;
        }
        else {
            cerr << "Impossibile aprire il file delle omega." << endl;
        }
        ///////////////////////////////////////////////////////////////

        // 6.3 storing correlation functions
        if (file_corr_q_square.is_open()){
            file_corr_q_square << corr_q_square[i] << endl;
        }
        else {
            cerr << "Impossibile aprire il file delle correlazioni delle q quadrate." << endl;
        }
        if (file_corr_p_square.is_open()){
            file_corr_p_square << corr_p_square[i] << endl;
        }
        else {
            cerr << "Impossibile aprire il file delle correlazioni delle p quadrate." << endl;
        }
        if (file_corr_omega.is_open()){
            file_corr_omega << corr_omega[i] << endl;
        }
        else {
            cerr << "Impossibile aprire il file delle correlazioni delle omega." << endl;
        }
        //////////////////////////////////////

        //6.4 storing the ensemble averaged integral omega
        if (file_integral_omega.is_open()){
            file_integral_omega << integral_omega_ensemble[i] << endl;
        }
        else{
            cerr << "Impossible aprire il file della omega integrale" << endl;
        }




    }

    // 7. Close the files of ensemble averaged trajectories.
    cout << "Fine memorizzazione dati" << endl;

    // 7.1 dynamical variables
    file_q.close();
    file_p.close();
    file_zeta.close();

    // 7.2 observables
    file_q_square.close();
    file_p_square.close();
    file_omega.close();

    // 7.3 correlations
    file_corr_q_square.close();
    file_corr_p_square.close();
    file_corr_omega.close();

    // 7.4 integrals
    file_integral_omega.close();


    return 0;
}
