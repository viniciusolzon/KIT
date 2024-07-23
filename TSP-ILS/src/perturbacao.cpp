#include "ILS.h"

Solucao ILS :: Pertubacao(Solucao *s, Data *d){
    Solucao s_copy;
    s_copy.sequence = s->sequence;
    s_copy.cost = s->cost;

    int lim = ceil(s_copy.sequence.size() / 10.0);

    int size1, start1, end1, size2, start2, end2;
    size1 = std::max(2, rand() % lim);
    start1 = rand()%(s_copy.sequence.size() - size1 - 1) + 1;
    end1 = start1 + size1 - 1;
    size2 = std::max(2, rand() % lim);

    if(s_copy.sequence.size() - 1 - end1 > size2 && start1 > size2){
        int coin = rand() % 2;
        if(coin == 0){
            int min = end1 + 1;
            int max = s_copy.sequence.size() - size2 - 1;
            start2 = rand()%(max - min + 1) + min;
            end2 = start2 + size2 - 1;
        }
        else{
            int max = start1 - size2;
            int min = 1;
            start2 = rand()%(max - min + 1) + min;
            end2 = start2 + size2 - 1;
        }
    }
    else if(s_copy.sequence.size() - 1 - end1 > size2){
            int min = end1 + 1;
            int max = s_copy.sequence.size() - size2 - 1;
            start2 = rand()%(max - min + 1) + min;
            end2 = start2 + size2 - 1;
    }   else{
            int max = start1 - size2 - 1;
            int min = 1;
            start2 = rand()%(max - min + 1) + min;
            end2 = start2 + size2 - 1;
        }
    
    vector<int> segment1(s_copy.sequence.begin() + start1, s_copy.sequence.begin() + start1 + size1);
    vector<int> segment2(s_copy.sequence.begin() + start2, s_copy.sequence.begin() + start2 + size2);
    
    double origin, after, cost;

    if(start1 + size1 == start2){
        origin = d->matrizAdj[s_copy.sequence[start1 - 1]][s_copy.sequence[start1]]
        + d->matrizAdj[s_copy.sequence[start1 + size1 - 1]][s_copy.sequence[start2]]
        + d->matrizAdj[s_copy.sequence[start2 + size2 - 1]][s_copy.sequence[start2 + size2]];

        after = d->matrizAdj[s_copy.sequence[start1 - 1]][s_copy.sequence[start2]]
        + d->matrizAdj[s_copy.sequence[start2 + size2 - 1]][s_copy.sequence[start1]]
        + d->matrizAdj[s_copy.sequence[start1 + size1 - 1]][s_copy.sequence[start2 + size2]];
    }
    else if(start2 + size2 == start1){
        origin = d->matrizAdj[s_copy.sequence[start2 - 1]][s_copy.sequence[start2]]
        + d->matrizAdj[s_copy.sequence[start2 + size2 - 1]][s_copy.sequence[start1]]
        + d->matrizAdj[s_copy.sequence[start1 + size1 - 1]][s_copy.sequence[start1 + size1]];
         
        after = d->matrizAdj[s_copy.sequence[start2 - 1]][s_copy.sequence[start1]]
        + d->matrizAdj[s_copy.sequence[start1 + size1 - 1]][s_copy.sequence[start2]]
        + d->matrizAdj[s_copy.sequence[start2 + size2 - 1]][s_copy.sequence[start1 + size1]];
    }
        else{
            origin = d->matrizAdj[s_copy.sequence[start1 - 1]][s_copy.sequence[start1]]
            + d->matrizAdj[s_copy.sequence[start1 + size1 - 1]][s_copy.sequence[start1 + size1]]
            + d->matrizAdj[s_copy.sequence[start2 - 1]][s_copy.sequence[start2]]
            + d->matrizAdj[s_copy.sequence[start2 + size2 - 1]][s_copy.sequence[start2 + size2]];

            after = d->matrizAdj[s_copy.sequence[start1 - 1]][s_copy.sequence[start2]]
            + d->matrizAdj[s_copy.sequence[start2 + size2 - 1]][s_copy.sequence[start1 + size1]]
            + d->matrizAdj[s_copy.sequence[start2 - 1]][s_copy.sequence[start1]]
            + d->matrizAdj[s_copy.sequence[start1 + size1 - 1]][s_copy.sequence[start2 + size2]];
    }

    if(start1 > start2){
        s_copy.sequence.erase(s_copy.sequence.begin() + start1, s_copy.sequence.begin() + start1 + size1);
        s_copy.sequence.insert(s_copy.sequence.begin() + start1, segment2.begin(), segment2.end());
        s_copy.sequence.erase(s_copy.sequence.begin() + start2, s_copy.sequence.begin() + start2 + size2);
        s_copy.sequence.insert(s_copy.sequence.begin() + start2, segment1.begin(), segment1.end());
    }
    else{
        s_copy.sequence.erase(s_copy.sequence.begin() + start2, s_copy.sequence.begin() + start2 + size2);
        s_copy.sequence.insert(s_copy.sequence.begin() + start2, segment1.begin(), segment1.end());
        s_copy.sequence.erase(s_copy.sequence.begin() + start1, s_copy.sequence.begin() + start1 + size1);
        s_copy.sequence.insert(s_copy.sequence.begin() + start1, segment2.begin(), segment2.end());
    }

    cost = after - origin;
    s_copy.cost += cost;

    return s_copy;
}
