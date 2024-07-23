#include "MLP.h"

MLP :: MLP(Data d){ // Construtor
    dist = d;
    maxIter = 10;
    maxiterILS = d.vertices >= 100 ? d.vertices : d.vertices;
    subseq_matrix = vector<vector<Subsequence>> (dist.vertices + 1, vector<Subsequence>(dist.vertices + 1));
}

MLP :: ~MLP(){} // Destrutor

bool compare(const InsertionInfo &a, const InsertionInfo &b){
    return a.cost < b.cost;
}

double MLP :: epsilon(double a, double b){
    return fabs(a + b) * numeric_limits<double>::epsilon() * 15;
}

bool MLP :: improve(double value_1, double value_2){
    return (value_1 - value_2) > epsilon(value_1, value_2);
}

void MLP :: calcularcost(Solucao *s, Data *d){
    s->cost = 0.0;
    for(int i = 0; i < d->vertices; i++)
        s->cost += (d->vertices - i) * d->matrizAdj[s->sequence[i]][s->sequence[i+1]];
}

void UpdateAllSubseq(Solucao *s, vector<vector<Subsequence>> &subseq_matrix, Data *d){
    int n = s->sequence.size();
    // subsequencias de um unico no
    for (int i = 0; i < n; i++){
        // int v = s->sequence[i];
        subseq_matrix[i][i].W = (i > 0);
        subseq_matrix[i][i].C = 0;
        subseq_matrix[i][i].T = 0;
        subseq_matrix[i][i].first = s->sequence[i];
        subseq_matrix[i][i].last = s->sequence[i];
    }

    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            subseq_matrix[i][j] = Subsequence :: Concatenate(subseq_matrix[i][j-1], subseq_matrix[j][j], d);

    
    for (int i = n - 1; i >= 0; i--)
        for (int j = i - 1; j >= 0; j--)
            subseq_matrix[i][j] = Subsequence :: Concatenate(subseq_matrix[i][j+1], subseq_matrix[j][j], d);
}

void MLP :: solve(){
    auto start = std::chrono::high_resolution_clock::now(); // Inicia o cronômetro
    srand(time(NULL)); // Para conseguir gerar números aleatórios
    Solucao bestOfAll; // Cria a solução que guardará a melhor solução possível

    for(int i = 0; i < maxIter; i++){
        Solucao current = Construcao(&dist); // Cria a solução atual ( Muita gula e pouca aleatoriedade )

        Solucao best = current; // Cria a solução que guardará a melhor solução dessa iteração
        if(i == 0)
            bestOfAll.cost = current.cost; // Se for a primeira solução criada ela já é atribuída para a melhor solução possível

        int iterILS = 0;
        while(iterILS <= maxiterILS){ // Enquanto houver melhoras nessa solução, reinicia-se o processo de melhora
            BuscaLocal(&current, &dist); // Tenta melhorar a solução atual
            if(improve(best.cost, current.cost)){ // Se a solução atual foi melhorada a melhor solução atual é atualizada pela atual
                best = current;
                iterILS = 0;
            }
            current = Pertubacao(&best, &dist); // Pertuba de forma aleatória a solução atual pra ver se ao alterá-la é possível melhorá-la
            
            if(current.cost < 0){
                show_solution();
                break;
            }

            iterILS++;
        }

        if(improve(bestOfAll.cost, best.cost)) // Se a melhor solução atual for melhor que a melhor solução possível, ela se torna a melhor solução possível
            bestOfAll = best;
    }

    final_solution = bestOfAll;

    auto end = std::chrono::high_resolution_clock::now(); // Para o cronômetro
    std::chrono::duration<double, std::milli> float_ms = end - start; // Calcula o tempo do cronômetro
    cout << "\nExecution time:  " << float_ms.count() / 1000.0000000000000 << " seconds" << "\n";

}

void MLP :: show_solution(){
    cout << "Sequence:  ";
    for(int i = 0; i < final_solution.sequence.size() - 1; i++)
        cout << final_solution.sequence[i] << " -> ";
    cout << final_solution.sequence.back() << endl;

    cout << "Cost:  " << setprecision(2) << fixed << final_solution.cost << "\n";
}
