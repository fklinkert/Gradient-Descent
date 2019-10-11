/*
* Group 2 HW2
* Federico Klinkert(fgk2106)
* Mateo Gomez(mg4010)
* Haoxuan Huang(hh2773)
* Jiaxin Lin(jl5304)
*/



#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <ctime>
using namespace std;
using namespace std::chrono;


void split(vector<double> &vec, string line, char delimeter){
    // split the whole line by delimeter and store them into vector of doubles
    int i = 0;
    int start = 0;
    while(i < line.length()){
        if(line[i] == delimeter){
            string str = line.substr(start, i-start);
            start = i+1;
            vec.push_back(stod(str));
        }
        i ++;
    }
}

void getData(vector<vector<double> > &datalist, string filename){
    //get data from <filename.csv> and store them as 2d vector of doubles
    ifstream file(filename);

    string line = "";

    while (getline(file, line)){
        vector<double> vec;
        split(vec, line, ',');
        datalist.push_back(vec);
    }

    file.close();
}

void getDelta(vector<vector<double> > datalist, vector<vector<double> > &delta){
    for(int r = 0; r < datalist.size(); r++){
        vector<double> line;
        for(int c = 1; c < datalist[0].size(); c++){
            line.push_back(datalist[r][c] - datalist[r][c-1]);
        }
        delta.push_back(line);
    }
}

double g(vector<double> alpha, vector<double> delta, int t){
    double res = delta[t];
    for(int i = t-1, j = 0; i >= t-10; i--, j++){
        res -= alpha[j]*delta[i];
    }

    return res;
}

void getGradient(vector<double> &gradient, vector<double> delta, vector<double> alpha){
    for(int i = 0; i < 10; i++){
        gradient[i] = 0;
    }
    for(int t = 10; t < delta.size(); t++){
        double const_term = 2 * g(alpha, delta, t);
        for(int i = 0; i < 10; i++){
            gradient[i] += const_term * (-1) * delta[t-i-1];
        }
    }
}
double f(vector<double> alpha, vector<double> delta){
    double res = 0;
    for(int t = 10; t < delta.size(); t++){
        res += pow(g(alpha, delta, t),2);
    }
    return res;
}

double dot(vector<double> vec1, vector<double> vec2){
    double res = 0;
    for(int i = 0; i < vec1.size(); i++){
        res += vec1[i]*vec2[i];
    }
    return res;
}

double getNorm(vector<double> v){
    double res = 0;
    for(int i = 0; i < v.size(); i++){
        res += v[i]*v[i];
    }

    return sqrt(res);
}

void backtracking(vector<double> &alpha, const vector<double> &old_alpha, const vector<double> &delta_alpha,
    const vector <double> &delta, const vector<double> &gradient){
    double a = 0.25, b = 0.5;
    double step = 1;
    while(1){
        for(int k = 0; k < 10; k++){
            alpha[k] = old_alpha[k] + step * delta_alpha[k];
        }
        double lhs = f(alpha, delta);
        double rhs = f(old_alpha, delta) + a * step * dot(gradient, delta_alpha);

        if(lhs <= rhs){
            break;
        }
        else{
            step *= b;
        }
    }
}

void accelerated(vector<double> &alpha, vector<double> &old_alpha,
    const vector<double> &delta, double e) {
    vector<double> gradient(10, 1);
    vector<double> old_y(10, 0);
    double gradient_norm = 100;
    int k = 1;
    while(gradient_norm > e){
        // obtain gradient
        getGradient(gradient, delta, alpha);

        vector<double> delta_alpha;
        // assign new alpha and store the old one
        for(int k = 0; k < 10; k++){
            delta_alpha.push_back(-gradient[k]);
            old_alpha[k] = alpha[k];
        }

        //accelerated part
        vector<double> y(10, 0.0);
        double s = k/(k+3.0);
        double t = 0.00085;

        for(int k = 0; k < 10; k++){
            y[k] = old_alpha[k] + t * delta_alpha[k];
        }

        for(int i = 0; i < 10; i++){
            alpha[i] = y[i] + s * (y[i] - old_y[i]);
            old_y[i] = y[i];
        }


        gradient_norm = getNorm(gradient);
        k ++;
    }
}
void accelerate( vector<double> &alpha, const vector<double> &old_alpha, const vector<double> &delta_alpha,
    const vector<vector <double> > &delta, int i, const vector<double> &gradient, vector<double> &old_y, int k){
    
}

void stochastic(vector<double> &alpha, vector<double> &old_alpha,
    const vector<double> &delta, double e, double N) {
    double gradient_norm = 100;
    double t = 0.00000005*N;
    while(gradient_norm > e){
        //generate a random index k
        int k = rand() % (delta.size() - 10 - 1) + 10;
        vector<double> gradient(10, 0.0);

        // calculate the constant part for all 10 alphas
        double const_term = 2 * g(alpha, delta, k);

        // getting gradient and update alpha
        for(int i = 0; i < 10; i++){
            gradient[i] = -const_term * delta[k-i-1];
            alpha[i] -= t*gradient[i];
        }
        getGradient(gradient, delta, alpha);

        gradient_norm = getNorm(gradient);
    }
    cout << e << "; t: " << t <<endl;
}

void gradientDescent(vector<double> &alpha, vector<double> &old_alpha,
    const vector<double> &delta, double e) {
    vector<double> gradient(10, 1);
    double gradient_norm = 100;
    while(gradient_norm > e){
        // obtain gradient
        getGradient(gradient, delta, alpha);

        vector<double> delta_alpha;
        // assign new alpha and store the old one
        for(int k = 0; k < 10; k++){
            delta_alpha.push_back(-gradient[k]);
            old_alpha[k] = alpha[k];
        }

        //backtracking to get the correct step
        backtracking(alpha, old_alpha, delta_alpha, delta, gradient);

        gradient_norm = getNorm(gradient);
    }
}



int main(int argc, char ** argv){
    // Get starting timepoint
    auto start = high_resolution_clock::now();

    if (argc != 2){
        cerr << "usage: ./a.out <filename>\n";
        exit(1);
    }
    cout << fixed;
    cout << setprecision(8);

    // create a datalist and fetch data from the required file
    vector<vector<double> > datalist;
    vector<vector<double> > delta;

    getData(datalist, argv[1]);

    getDelta(datalist, delta);

    srand(time(NULL));
    double N = double(delta[0].size()-10);
    vector<vector<double> > alphas;
    for(int i = 0; i < datalist.size(); i++){
        // 1000 rows of assets
        double e = 0.01;
        vector<double> alpha(10, 0);
        vector<double> gradient(10, 1);
        vector<double> old_alpha(10, 0);
        

        // gradient descent
        gradientDescent(alpha, old_alpha, delta[i], e);


        //accelerated
        accelerated(alpha, old_alpha, delta[i], e);
        

        // stochastic gradient descent
        stochastic(alpha, old_alpha, delta[i], e, N);

    
        alphas.push_back(alpha);

    }



    for(int i = 0; i < 10; i++){
        cout << setprecision(10) << alphas[1][i] << endl;
    }
    cout << endl;


    // Get ending timepoint
    auto stop = high_resolution_clock::now();
    // Get duration. Substart timepoints to
    // get durarion. To cast it to proper unit
    // use duration cast method
    auto duration = duration_cast<microseconds>(stop - start);
  
    cout << "Time taken by function: "
         << duration.count() /1000000.0 << " seconds" << endl;

    return 0;
}
