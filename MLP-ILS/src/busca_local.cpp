#include "MLP.h"

double MLP :: calculateSwapCost(Solucao *s, Data *d, int first, int second){
    Subsequence sigma1, sigma2, sigma3, major_sigma;

    if(second == first + 1){
        sigma1 = Subsequence :: Concatenate(subseq_matrix[0][first - 1], subseq_matrix[second][second], d);
        sigma2  = Subsequence :: Concatenate(sigma1, subseq_matrix[first][first], d);

        major_sigma = Subsequence :: Concatenate(sigma2, subseq_matrix[second + 1][d->vertices], d);
    }
    else{
        sigma1 = Subsequence :: Concatenate(subseq_matrix[0][first - 1], subseq_matrix[second][second], d);
        sigma2  = Subsequence :: Concatenate(sigma1, subseq_matrix[first + 1][second - 1], d);
        sigma3 = Subsequence :: Concatenate(sigma2, subseq_matrix[first][first], d); 

        major_sigma = Subsequence :: Concatenate(sigma3, subseq_matrix[second + 1][d->vertices], d);
    }
	return major_sigma.C;
}

bool MLP :: bestImprovementSwap(Solucao *s, Data *d){
    double bestcost = s->cost;
    double cost = 0.0;
    int best_i, best_j;
    for(int i = 1; i < s->sequence.size() - 2; i++){
        for(int j = i + 1; j < s->sequence.size() - 1; j++){
            cost = calculateSwapCost(s, d, i, j);
            if(improve(bestcost, cost)){
                bestcost = cost;
                best_i = i;
                best_j = j;
            }
        }
    }
    if(bestcost < s->cost){
        int aux = s->sequence[best_i];
        s->sequence[best_i] = s->sequence[best_j];
        s->sequence[best_j] = aux;
        s->cost = bestcost;
        return true;
    }

    return false;
}

double MLP :: calculate2OptCost(Solucao *s, Data *d, int first, int second){
   Subsequence sigma_1, major_sigma;

   sigma_1 = Subsequence :: Concatenate(subseq_matrix[0][first - 1], subseq_matrix[second][first], d);
   major_sigma = Subsequence :: Concatenate(sigma_1, subseq_matrix[second + 1][d->vertices], d);

   return major_sigma.C;
}

bool MLP :: bestImprovement2Opt(Solucao *s, Data *d){
    double bestcost = s->cost, cost = 0.0;
    int best_i, best_j;             

    for(int i = 1; i < s->sequence.size() - 2; i++){
        for(int j = i + 1; j < s->sequence.size() - 1; j++){
            cost = calculate2OptCost(s, d, i, j);
            if(improve(bestcost, cost)){
                best_i = i;
                best_j = j;
                bestcost = cost;
            }
        }
    }
    if(bestcost < s->cost){
        reverse(s->sequence.begin() + best_i, s->sequence.begin() + best_j + 1);
        s->cost = bestcost;
        return true;
    }

    return false;
}

double MLP :: calculateOrOptCost(Solucao *s, Data *d, int first, int second, int amount){
    Subsequence sigma_1, sigma_2, complete_sub;

    if (second > first){
        sigma_1 = Subsequence :: Concatenate(subseq_matrix[0][first - 1], subseq_matrix[first + amount][second + amount - 1], d);
        sigma_2 = Subsequence :: Concatenate(sigma_1, subseq_matrix[first][first + amount - 1], d);
        complete_sub = Subsequence :: Concatenate(sigma_2, subseq_matrix[second + amount][d->vertices], d);
    }
    else if (first > second){
        sigma_1 = Subsequence ::  Concatenate(subseq_matrix[0][second - 1], subseq_matrix[first][first + amount - 1], d);
        sigma_2 = Subsequence :: Concatenate(sigma_1, subseq_matrix[second][first - 1], d);
        complete_sub = Subsequence :: Concatenate(sigma_2, subseq_matrix[first + amount][d->vertices], d);
    }

   return complete_sub.C;
}

bool MLP :: bestImprovementOrOpt(Solucao *s, Data *d, int amount){
    double bestcost = s->cost, cost = 0.0;
    int best_i, best_j;
    for(int i = 1; i < s->sequence.size() - amount; i++){
        for(int j = 1; j < s->sequence.size() - amount; j++){
            if(i == j)
                continue;
            
            cost = calculateOrOptCost(s, d, i, j, amount);
            
            if(improve(bestcost, cost)){
                best_i = i;
                best_j = j;
                bestcost = cost;
            }
        }
    }

    if(bestcost < s->cost){
        vector<int> section(s->sequence.begin() + best_i, s->sequence.begin() + best_i  + amount);
        s->sequence.erase(s->sequence.begin() + best_i, s->sequence.begin() + best_i + amount);
        s->sequence.insert(s->sequence.begin() + best_j, section.begin(), section.end());

        s->cost = bestcost;
        
        return true;
    }

    return false;
}

void MLP :: BuscaLocal(Solucao *s, Data *d){
    vector<int> NL = {1, 2, 3, 4, 5};
    bool improved = false;
    while(NL.empty() == false){
        int n = rand() % NL.size(); // Gera aleatório
        switch(NL[n]){
            case 1:
                improved = bestImprovementSwap(s, d);
                break;
            case 2:
                improved = bestImprovement2Opt(s, d);
                break;
            case 3:
                improved = bestImprovementOrOpt(s, d, 1);
                break;
            case 4:
                improved = bestImprovementOrOpt(s, d, 2);
                break;
            case 5:
                improved = bestImprovementOrOpt(s, d, 3);
                break;
        }
        if(improved){
            NL = {1, 2, 3, 4, 5}; // Reinicia
            UpdateAllSubseq(s, subseq_matrix, d);
        }
        else
            NL.erase(NL.begin() + n);
    }
}
