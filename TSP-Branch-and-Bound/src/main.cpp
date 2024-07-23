#include <vector>
#include <list>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <iterator>
#include <chrono>
#include <queue>

using std::cout;
using std::endl;
using std::vector;
using std::list;
using std::pair;
using std::setprecision;
using std::fixed;
using std::numeric_limits;

#include "data.h"
#include "hungarian.h"

typedef struct Node{ // Estrutura de cada nó da árvore a ser utilizada no algoritmo Branch and Bound Combinatório
    vector<pair<int, int>> forbidden_arcs; // Vector de arcos proibidos, terão o custo 9999999 para não serem escolhidos
    vector<vector<int>> subtours; // Conjunto de subtours da solução gerada pelo Algoritmo Húngaro
    
    double lower_bound = 0; // Lower Bound do nó (ou custo total da solução do Algoritmo Húngaro)
    int pick = 0; // Indicador de qual subtour foi escolhido, quase sempre vai ser o menor aqui nesse código
    bool feasible; // Se a solução do gerada pelo Assignment Problem é viável para o TSP ou não
    double cost; // Custo total da solução 

	// bool operator == (const Node& rhs){
	// 	this->subtours == rhs.subtours;
	// }
}Node;

// Função auxiliar
bool is_feasible(Node node){ // indica se a solução é viável ou não
	return node.feasible;
}

// Função auxiliar
void show_subtours(Node node){ // exibe apenas os subtours da solução
    cout << "Subtours count: " << node.subtours.size() << "\n";
    cout << "Subtours:\n";
    for(int i = 0; i < node.subtours.size(); i++){
        for(int j = 0; j < node.subtours[i].size(); j++){
            cout << node.subtours[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}

// Função auxiliar
void show_arcs(Node node){ // exibe apenas os arcos proibidos da solução
	for(int i = 0; i < node.forbidden_arcs.size(); i++){  
		cout << "Forbidden arc" << i << "]: ";
		cout << node.forbidden_arcs[i].first  << " " << node.forbidden_arcs[i].second << "\n";
	}

	cout << "\n";
}

// Função auxiliar
void show_lowerbound(Node node){ // exibe apenas o lower bound da solução
	cout << "Lower Bound: " << 	node.lower_bound << "\n" << "\n";
}

// Função auxiliar
void show_info(Node node){ // exibe todas as informações da solução
    cout << "\nA new solution has been generated;\n";
    show_subtours(node);
    cout << "Cost: " << setprecision(2) << fixed << node.cost << "\n";
    if(node.feasible)
        cout << "Solution is feasible\n";
    else
        cout << "Solution isn't feasible\n";
}

// Função auxiliar
int subtour_pick(const vector<vector<int>> &subtour_vector){ // escolhe o subtour de menor tamanho para continuar o algoritmo
    int smallest;
    for(int i = 1; i < subtour_vector.size(); i++){
        if(subtour_vector[i].size() < subtour_vector[0].size()){ // se o subtour atual for menor que o primeiro subtour, ele é "escolhido"
            smallest = i; // salva o índice de qual é o menor subtour
        }
    }
    return smallest;
}

// Função core
vector<vector<int>> generate_subtours(hungarian_problem_t p){ // gera o vector que guarda os subtours
	vector<vector<int>> subtour_vector;
    vector<int> available_nodes(p.num_rows);

    // Cria vector dos nós ainda disponíveis
    for(int i = 0; i < p.num_rows; i++) // preenche vector dos nós dispońiveis
        available_nodes[i] = i + 1;

    while(!available_nodes.empty()){ // Vai executando até analisar todos os nós
        int last_node = available_nodes[0]; // Primeiro nó do tour
        int node_aux = last_node;

        vector<int> tour; // tour que será gerado
        tour.push_back(last_node); // insere o primeiro nó dos disponíveis no tour que está sendo formado

        do{
            for(int j = 0; j < p.num_rows; j++){
                if(p.assignment[node_aux - 1][j] == 1){
                    node_aux = j + 1;
					 // infelizmente n tem simples função p remover um elemento específico, tem q utilizar esse método
                    available_nodes.erase(remove(available_nodes.begin(), available_nodes.end(), node_aux), available_nodes.end());
                        break;
                }
            }
            tour.push_back(node_aux); // insere nó no tour

        }while(node_aux != last_node); // até o "último" nó encontrar com o primeiro
        	subtour_vector.push_back(tour); // insere tour no vector de tours
    }
    
    // for(int i = 0; i < subtour_vector.size(); i++){ // pra debugar
    //     for(int j = 0; j < subtour_vector[i].size(); j++){
    //         cout << subtour_vector[i][j] << " -> ";
    //     }
    //     cout << "\n";
    // }

    return subtour_vector; // retorna o vector que cnotém todos os subtours gerados para a solução em questão
}

// Função auxiliar
// Node find_lowest_bound(list<Node> tree){
// 	Node lowest_node = tree.begin();
// 	double menor_lb = numeric_limits<double>::infinity(); // o primeiro lower_bound encontrado já será atribuído a essa variável devido ao seu custo "infinito"

// 	for(int i = 0; i < tree.size(); i++){ // varre a lista procurando o nó com o menor lower bound
// 		if(tree[i].cost - menor_lb < numeric_limits<double>::epsilon()){
// 			menor_lb = tree[i].cost;
// 			lowest_node = tree[i];
// 		}
// 	}

// 	return lowest_node;
// }

// Função core
void getSolutionHungarian(Node &node, Data *data){
	// Precisamos criar arcos proibidos, e p/ isso a gente vai ter que alterar os custos da matriz de custos, todavida não podemos alterar a matriz
	// original pra não dar problema no algoritmo, por isso a gente vai criar uma matriz de custos "provisória" aqui, de forma que a gente consiga
	// gerar uma criar a matriz alterada, gerar uma solução a partir dessa matriz com os arcos proibidos e dps a gente só libera ela
	double **cost = new double *[data->getDimension()]; // Aloca espaço na memória para a matriz de custo
	for(int i = 0; i < data->getDimension(); i++){ // Preenche a matriz de custo
		cost[i] = new double[data->getDimension()];
		for(int j = 0; j < data->getDimension(); j++){
			cost[i][j] = data->getDistance(i, j);
		}
	}

	for(int i = 0; i < node.forbidden_arcs.size(); i++){ //primeiro limita a escolha de certos arcos, o custo desses são 999999..., por isso eles não são escolhidos
		cost[node.forbidden_arcs[i].first-1][node.forbidden_arcs[i].second-1] = 99999999;
	}

	hungarian_problem_t p;
	int mode = HUNGARIAN_MODE_MINIMIZE_COST;
	hungarian_init(&p, cost, data->getDimension(), data->getDimension(), mode); // carregando o problema

	double obj_value = hungarian_solve(&p); // depois guarda o custo da solução gerada pelo algoritmo húngaro
    node.cost = obj_value; // e atribui esse custo à solução sendo gerada
    
    node.subtours = generate_subtours(p); // atribui o vector com os subtours gerados à solução atual
	node.feasible = node.subtours.size() == 1 ? true : false; // se tiver só um subtour a solução é declarada com válida
	
	for (int i = 0; i < data->getDimension(); i++){
		delete [] cost[i];
	}
	delete [] cost;
	hungarian_free(&p);


	int smallest_subtour = 0;
    for(int i = 1; i < node.subtours.size(); i++){
        if(node.subtours[i].size() < node.subtours[0].size()){ // se o subtour atual for menor que o primeiro subtour, ele é "escolhido"
            smallest_subtour = i; // salva o índice de qual é o menor subtour
        }
    }
	node.pick = smallest_subtour;
}

Node BranchBound_BFS(Data *data, double **cost){
	Node root; // cria o nó raíz da árvore
    Node best_solution; // cria o nó em que será guardada a melhor solução
	best_solution.cost = numeric_limits<double>::infinity(); // seta o custo da melhor solução inicialmente para "infinito"
	list<Node> tree;
	tree.push_back(root); // cria a árvore que será preenchida e feita a busca nela
    // a "árvore" será na verdade uma lista em que iremos inserindo e removendo no começo e/ou no final de acordo com o algortimo

	double upper_bound = 99999999; // "qualquer", ou seja, o primeiro upper_bound encontrado já será atribuído a essa variável devido ao seu custo

	while(!tree.empty()){ // equanto todos os nós da árvore não forem analisados
    	// Node node = tree.front(); // vai pelo primeiro nó -> BFS strategy
		
    	auto node = tree.begin(); // vai pelo primeiro nó -> BFS strategy

		getSolutionHungarian(*node, data); // calcula a solução do subtour atual

        // se for debugar descomenta a linha de baixo, mas o melhor pra debuggar é logo aqui embaixo na linha 231
        // show_info(*node);

		if(node->cost > upper_bound){ // solução gerada pelo algoritmo húngaro do assignment problema é descartada pois o custo é maior que o upper bound
			tree.erase(node);
			continue;
		}
		if(node->feasible == true){ // se a solução for viável
			if(node->cost < best_solution.cost){// e se o custo da solução for menor que o o custo da melhor solução encontrada até agora
                best_solution = *node; // atualizamos a melhor solução
				upper_bound = node->cost; // upper_bound é atualizado com o custo dessa nova solução
                // Solução viável foi gerada, caso queira olhar as soluções geradas na execução do algoritmo descomenta "show_info()"
                // show_info(best_solution);
			}

			tree.erase(node); // esse nó já foi analisado então tira ele da árvore
			continue;
		}
        // Solução inviável foi gerada, caso queira olhar as soluções geradas na execução do algoritmo descomenta "show_info()"
        // else{
                // show_info(*node);
        // }
        

		for(int i = 0; i < node->subtours[node->pick].size() - 1; i++){ // Adiciona folhas/filhos na árvore
			Node new_node;
			new_node.forbidden_arcs = node->forbidden_arcs; // atribui os arcos proibidos do nó da árvore ao nó atual que será inserido no final da árvore

			pair<int, int> new_forbidden_arc;
			new_forbidden_arc.first = node->subtours[node->pick][i];
			new_forbidden_arc.second = node->subtours[node->pick][i + 1];

			new_node.forbidden_arcs.push_back(new_forbidden_arc);
			tree.insert(tree.end(), new_node); // insere novos nós na árvore de busca
		}
		tree.erase(node);
	}

    return best_solution; // returns the best solution
}

Node BranchBound_DFS(Data *data, double **cost){ // executa o algoritmo utilizando busca por profundidade DFS(Depth First Search)
	Node root; // cria o nó raíz da árvore
    Node best_solution; // cria o nó em que será guardada a melhor solução
	best_solution.cost = numeric_limits<double>::infinity(); // seta o custo da melhor solução inicialmente para "infinito"
	list<Node> tree;
	tree.push_back(root); // cria a árvore que será preenchida e feita a busca nela
    // a "árvore" será na verdade uma lista em que iremos inserindo e removendo no começo e/ou no final de acordo com o algortimo

	double upper_bound = 99999999; // "qualquer", ou seja, o primeiro upper_bound encontrado já será atribuído a essa variável devido ao seu custo

	while(!tree.empty()){ // equanto todos os nós da árvore não forem analisados
    	// auto node = tree.end(); // assim não funciona, precisa ser assim aí na linha de baixo
    	auto node = prev(tree.end()); // vai pelo último nó -> DFS strategy

		getSolutionHungarian(*node, data); // calcula a solução do subtour atual

        // se for debugar descomenta a linha de baixo, mas o melhor pra debuggar é logo aqui embaixo na linha 231
        // show_info(*node);

		if(node->cost > upper_bound){ // solução gerada pelo algoritmo húngaro do assignment problema é descartada pois o custo é maior que o upper bound
			tree.erase(node);
			continue;
		}
		if(node->feasible == true){ // se a solução for viável
			if(node->cost < best_solution.cost){// e se o custo da solução for menor que o o custo da melhor solução encontrada até agora
                best_solution = *node; // atualizamos a melhor solução
				upper_bound = node->cost; // upper_bound é atualizado com o custo dessa nova solução
                // Solução viável foi gerada, caso queira olhar as soluções geradas na execução do algoritmo descomenta "show_info()"
                // show_info(best_solution);
			}

			tree.erase(node); // esse nó já foi analisado então tira ele da árvore
			continue;
		}
        // Solução inviável foi gerada, caso queira olhar as soluções geradas na execução do algoritmo descomenta "show_info()"
        // else{
        //         // show_info(node);
        // }
        
		// aqui tava dando mt problema pq o (node.subtours[i].size()) nunca era < que (node.subtours[0].size()) no final da função getSolutionHungarian() então
		// o node->pick não era escolhido, fazendo com que ele ficasse com o valor inicializado, e eu não setava ele inicialmente pra 0, então ele pegava lixo
		// de memória. Logo, esse loop aqui de baixo nunca era inicializado e bugava o código. Pra arrumar só seto o pick inicial pra zero, ou seja, se ele não
		// achar nenhum tour menor pra escolher, (todos são do mesmo tamanho), ele escolhe o primeiro.
		for(int i = 0; i < node->subtours[node->pick].size() - 1; i++){ // Adiciona folhas/filhos na árvore
			Node new_node;
			new_node.forbidden_arcs = node->forbidden_arcs; // atribui os arcos proibidos do nó da árvore ao nó atual que será inserido no final da árvore

			pair<int, int> new_forbidden_arc;
			new_forbidden_arc.first = node->subtours[node->pick][i];
			new_forbidden_arc.second = node->subtours[node->pick][i + 1];

			new_node.forbidden_arcs.push_back(new_forbidden_arc);
			tree.insert(tree.end(), new_node); // insere novos nós na árvore de busca
		}

		tree.erase(node);
	}

    return best_solution; // returns the best solution
}

Node BranchBound_LB(Data *data, double **cost){ // executa o algoritmo utilizando a estratégia de busca pelo menor Lower Bound
	Node root; // cria o nó raíz da árvore
    Node best_solution; // cria o nó em que será guardada a melhor solução
	best_solution.cost = numeric_limits<double>::infinity(); // seta o custo da melhor solução inicialmente para "infinito"
	list<Node> tree;
	tree.push_back(root); // cria a árvore que será preenchida e feita a busca nela
    // a "árvore" será na verdade uma lista em que iremos inserindo e removendo no começo e/ou no final de acordo com o algortimo

	double upper_bound = 99999999; // "qualquer", ou seja, o primeiro upper_bound encontrado já será atribuído a essa variável devido ao seu custo

	while(!tree.empty()){ // equanto todos os nós da árvore não forem analisados
    	
		// Node *node; // tentei fazer isso funcionar mas dá um erro esquisito
		// double menor_lb = numeric_limits<double>::infinity();
		// for(int i = 0; i < tree.size(); i++){
		// 	if(tree[i].cost - menor_lb < numeric_limits<double>::epsilon()){
		// 		menor_lb = tree[i].cost;
		// 		*node = tree[i];
		// 	}
		// }

		auto node = tree.begin(); // iterador que vai ser rutilizado na busca pelo menor lower bound na lista dos nós(árvore)
		double menor_lb = numeric_limits<double>::infinity();
		for(auto it = tree.begin(); it != tree.end(); it++){ // varre a árvore procurando o nó com o menor lower bound para continuar a busca 
			if(it->cost - menor_lb < numeric_limits<double>::epsilon()){
				menor_lb = it->cost;
				node = it;
			}
		}

		getSolutionHungarian(*node, data); // calcula a solução do subtour atual

        // se for debugar descomenta a linha de baixo, mas o melhor pra debuggar é logo aqui embaixo na linha 231
        // show_info(*node);

		if(node->cost > upper_bound){ // solução gerada pelo algoritmo húngaro do assignment problema é descartada pois o custo é maior que o upper bound
			tree.erase(node); // esse nó já foi analisado então tira ele da árvore
			continue;
		}
		if(node->feasible == true){ // se a solução for viável
			if(node->cost < best_solution.cost){// e se o custo da solução for menor que o o custo da melhor solução encontrada até agora
                best_solution = *node; // atualizamos a melhor solução
				upper_bound = node->cost; // upper_bound é atualizado com o custo dessa nova solução
                // Solução viável foi gerada, caso queira olhar as soluções geradas na execução do algoritmo descomenta "show_info()"
                // show_info(best_solution);
			}

			tree.erase(node); // esse nó já foi analisado então tira ele da árvore
			continue;
		}
        // Solução inviável foi gerada, caso queira olhar as soluções geradas na execução do algoritmo descomenta "show_info()"
        // else{
        //         // show_info(node);
        // }
        
		for(int i = 0; i < node->subtours[node->pick].size() - 1; i++){ // Adiciona folhas/filhos na árvore
			Node new_node;
			new_node.forbidden_arcs = node->forbidden_arcs; // atribui os arcos proibidos do nó da árvore ao nó atual que será inserido no final da árvore

			pair<int, int> new_forbidden_arc;
			new_forbidden_arc.first = node->subtours[node->pick][i];
			new_forbidden_arc.second = node->subtours[node->pick][i + 1];

			new_node.forbidden_arcs.push_back(new_forbidden_arc);
			tree.insert(tree.end(), new_node); // insere novos nós na árvore de busca
		}

		tree.erase(node);
	}

    return best_solution; // returns the best solution
}

void BranchBound_show_solution(Node best_solution){ // exibe informações sobre a solução final encontrada
    cout << "Final Solution: ";
    for(int i = 0; i < best_solution.subtours.size(); i++){
        for(int j = 0; j < best_solution.subtours[i].size(); j++){
            cout << best_solution.subtours[i][j] << " -> ";
        }
        cout << "\n";
    }
    cout << "Cost: "  << setprecision(2) << fixed << best_solution.cost << "\n";
}

int main(int argc, char **argv){
    // O usuário define qual tipo de estrateǵia será utilizada na exedcução do algortimo
    cout << "Choose what strategy to use in the Branch and Bound algorithm:\n";
    cout << "\t|1|- Breath First Search\n\t|2|- Depth First Search\n\t|3|- Lower Bound Search\n-> ";
    int choice;
    cin >> choice;
    getchar();
    if(choice == 1)
        cout << "Using Breath First Search;\n";
    else if(choice == 2)
            cout << "Using Depth First Search;\n";
        else
            cout << "Using Lower Bound Search;\n";

	Data *data = new Data(argc, argv[1]);
	data->readData();

    auto start = std::chrono::high_resolution_clock::now(); // Inicia o cronômetro

	double **cost = new double *[data->getDimension()]; // Aloca espaço na memória para a matriz de custo
	for(int i = 0; i < data->getDimension(); i++){ // Preenche a matriz de custo
		cost[i] = new double[data->getDimension()];
		for(int j = 0; j < data->getDimension(); j++){
			cost[i][j] = data->getDistance(i, j);
		}
	}

    // dependendo da esolha do usuário, a técnica de execução do algoritmo varia
    if(choice == 1)
        BranchBound_show_solution(BranchBound_BFS(data, cost));
    else if(choice == 2)
            BranchBound_show_solution(BranchBound_DFS(data, cost));
        else
            BranchBound_show_solution(BranchBound_LB(data, cost));

    auto end = std::chrono::high_resolution_clock::now(); // Para o cronômetro
    std::chrono::duration<double, std::milli> float_ms = end - start; // Calcula o tempo do cronômetro
    cout << "Execution time: " << float_ms.count() / 1000.0000000000000 << " seconds" << "\n";

	for (int i = 0; i < data->getDimension(); i++) // Libera a memória previamente alocada
		delete[] cost[i];
	delete[] cost;
	delete data;
	return 0;
}
