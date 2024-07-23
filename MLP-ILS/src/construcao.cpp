#include "MLP.h"

Solucao MLP :: Construcao(Data *d){
    Solucao current;
    vector<int> CL;

    for (int i = 2; i <= d->vertices; i++){
        CL.push_back(i);
    }
    // Sequencia começa com 1
    current.sequence.push_back(1);

    random_shuffle(CL.begin(), CL.end());
    for(int i = 0; i < 2; i++){
        current.sequence.push_back(CL[i]);
        CL.erase(remove(CL.begin(), CL.end(), CL[i]), CL.end());
    }
    // E termina com 1
    current.sequence.push_back(1);

    vector<InsertionInfo> insertioncost;

    while(!CL.empty()){
        insertioncost = calcularCustoInsercao(&current, d, CL);
        sort(insertioncost.begin(), insertioncost.end(), compare);

        double alpha = (double) (rand() % 25) / 100;
        int selecionado = (alpha == 0) ? 0 : rand() % ((int) ceil(alpha * insertioncost.size())) ;

        current.sequence.insert(current.sequence.begin() + insertioncost[selecionado].remover + 1, insertioncost[selecionado].inserir);

        CL.erase(remove(CL.begin(), CL.end(), insertioncost[selecionado].inserir), CL.end());
    }

    UpdateAllSubseq(&current, subseq_matrix, d);
    current.cost = subseq_matrix[0][d->vertices].C;

    return current;
}

vector<InsertionInfo> MLP :: calcularCustoInsercao(Solucao *s, Data *d, vector<int> CL){
    InsertionInfo info;
    vector<InsertionInfo> custoInsercao;

    int l = 0;
    for(int a = 0, b = 1; l < s->sequence.size() - 1; a++, b++){
        int i = s->sequence[a];
        int j = s->sequence[b];
        
        for(auto k : CL){
            info.inserir = k;
            info.remover = a;
            info.cost = d->matrizAdj[i][k] + d->matrizAdj[j][k] - d->matrizAdj[i][j];

            custoInsercao.push_back(info);
        }
        l++;
    }

    return custoInsercao;
}
