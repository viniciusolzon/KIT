#ifndef MLP_H
#define MLP_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <chrono>

#include <random>

using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::setprecision;
using std::fixed;
using std::numeric_limits;

struct Data{
    double **matrizAdj;
    int vertices;
};

struct Solucao{
    vector<int> sequence;
    double cost;
};

struct InsertionInfo{
    int inserir;
    int remover;
    double cost = 0.0;
};

struct Subsequence{
    double T = 0, C = 0;
    int W = 0;
    int first = 0, last = 0;

    inline static Subsequence Concatenate(Subsequence &sigma_1, Subsequence &sigma_2, Data *d){
        Subsequence sigma;
        double temp = d->matrizAdj[sigma_1.last][sigma_2.first];
        sigma.W = sigma_1.W + sigma_2.W;
        sigma.T = sigma_1.T + temp + sigma_2.T;
        sigma.C = sigma_1.C + sigma_2.W * (sigma_1.T + temp) + sigma_2.C;
        sigma.first = sigma_1.first;
        sigma.last = sigma_2.last;

    return sigma;
    }
};

class MLP{
private:
    int maxIter;
    int maxiterILS;

    Data dist;
    Solucao final_solution;
    vector<vector<Subsequence>> subseq_matrix;

public:
    MLP(Data d); // CONSTRUTOR
    ~MLP(); // DESTRUTOR

    double epsilon(double a, double b);
    bool improve(double value_1, double value_2);
    void calcularcost(Solucao *s, Data *d);

    Solucao Construcao(Data *d);
    vector<InsertionInfo> calcularCustoInsercao(Solucao *s, Data *d, vector<int> CL);

    void BuscaLocal(Solucao *s, Data *d);
    double calculateSwapCost(Solucao *s, Data *d, int first, int second);
    bool bestImprovementSwap(Solucao *s, Data *d);
    double calculate2OptCost(Solucao *s, Data *d, int first, int second);
    bool bestImprovement2Opt(Solucao *s, Data *d);
    double calculateOrOptCost(Solucao *s, Data *d, int first, int second, int amount);        
    bool bestImprovementOrOpt(Solucao *s, Data *d, int amount);

    Solucao Pertubacao(Solucao *s, Data *d);

    void solve();
    void show_solution();
};

bool compare(const InsertionInfo &a, const InsertionInfo &b);

void UpdateAllSubseq(Solucao *s, vector<vector<Subsequence>> &subseq_matrix, Data *d);

#endif
