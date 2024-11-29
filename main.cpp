#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <limits>
#include <fstream>

using namespace std;

const int NUM_VERTICES = 100;      // Кількість вершин
const int MIN_WEIGHT = 5;          // Мінімальна вага ребра
const int MAX_WEIGHT = 150;        // Максимальна вага ребра
const int POPULATION_SIZE = 100;   // Розмір популяції
const int MAX_GENERATIONS = 100;  // Максимальна кількість поколінь
const double MUTATION_RATE = 0.2;  // Ймовірність мутації
const int ELITE_COUNT = 10;        // Кількість елітних рішень у популяції

struct Edge {
    int to;
    int weight;
};

struct Individual {
    vector<int> path;
    int fitness;
};

vector<vector<Edge>> graph(NUM_VERTICES);

void generateGraph() {
    srand(time(0));
    for (int i = 0; i < NUM_VERTICES; ++i) {
        for (int j = i + 1; j < NUM_VERTICES; ++j) {
            int weight = rand() % (MAX_WEIGHT - MIN_WEIGHT + 1) + MIN_WEIGHT;
            graph[i].push_back({j, weight});
            graph[j].push_back({i, weight});
        }
    }
}

int evaluatePath(const vector<int>& path) {
    int totalWeight = 0;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        bool found = false;
        for (const auto& edge : graph[path[i]]) {
            if (edge.to == path[i + 1]) {
                totalWeight += edge.weight;
                found = true;
                break;
            }
        }
        if (!found) return numeric_limits<int>::max();
    }
    return totalWeight;
}

vector<int> generateRandomPath(int start, int end) {
    vector<int> path = {start};
    while (path.back() != end) {
        int current = path.back();
        if (graph[current].empty()) break;
        int next = graph[current][rand() % graph[current].size()].to;
        if (find(path.begin(), path.end(), next) == path.end()) {
            path.push_back(next);
        }
    }
    if (path.back() != end) path.clear();
    return path;
}

vector<Individual> initializePopulation(int start, int end) {
    vector<Individual> population;
    for (int i = 0; i < POPULATION_SIZE; ++i) {
        vector<int> path = generateRandomPath(start, end);
        if (!path.empty()) {
            population.push_back({path, evaluatePath(path)});
        }
    }
    return population;
}

Individual onePointCrossover(const Individual& parent1, const Individual& parent2, int start, int end) {
    vector<int> childPath = {start};
    size_t cutPoint = rand() % parent1.path.size();
    for (size_t i = 1; i < cutPoint; ++i) {
        childPath.push_back(parent1.path[i]);
    }
    for (int gene : parent2.path) {
        if (find(childPath.begin(), childPath.end(), gene) == childPath.end()) {
            childPath.push_back(gene);
        }
        if (gene == end) break;
    }
    return {childPath, evaluatePath(childPath)};
}

Individual twoPointCrossover(const Individual& parent1, const Individual& parent2, int start, int end) {
    vector<int> childPath = {start};
    int cut1 = rand() % parent1.path.size();
    int cut2 = rand() % parent1.path.size();
    if (cut1 > cut2) swap(cut1, cut2);

    for (int i = cut1; i < cut2; ++i) {
        childPath.push_back(parent1.path[i]);
    }
    for (int gene : parent2.path) {
        if (find(childPath.begin(), childPath.end(), gene) == childPath.end()) {
            childPath.push_back(gene);
        }
        if (gene == end) break;
    }
    return {childPath, evaluatePath(childPath)};
}

Individual orderedCrossover(const Individual& parent1, const Individual& parent2, int start, int end) {
    vector<int> childPath(parent1.path.size(), -1);
    size_t cut1 = rand() % parent1.path.size();
    size_t cut2 = rand() % parent1.path.size();
    if (cut1 > cut2) swap(cut1, cut2);

    for (size_t i = cut1; i <= cut2; ++i) {
        childPath[i] = parent1.path[i];
    }
    size_t index = (cut2 + 1) % parent1.path.size();
    for (size_t i = 0; i < parent2.path.size(); ++i) {
        int gene = parent2.path[(cut2 + 1 + i) % parent2.path.size()];
        if (find(childPath.begin(), childPath.end(), gene) == childPath.end()) {
            childPath[index] = gene;
            index = (index + 1) % childPath.size();
        }
    }
    return {childPath, evaluatePath(childPath)};
}

void swapMutation(Individual& individual) {
    if ((double)rand() / RAND_MAX < MUTATION_RATE) {
        int i = rand() % individual.path.size();
        int j = rand() % individual.path.size();
        swap(individual.path[i], individual.path[j]);
        individual.fitness = evaluatePath(individual.path);
    }
}

void shiftMutation(Individual& individual) {
    if ((double)rand() / RAND_MAX < MUTATION_RATE) {
        int start = rand() % individual.path.size();
        int end = rand() % individual.path.size();
        if (start > end) swap(start, end);
        rotate(individual.path.begin() + start, individual.path.begin() + end, individual.path.end());
        individual.fitness = evaluatePath(individual.path);
    }
}

void twoOptImprovement(Individual& individual) {
    for (size_t i = 0; i < individual.path.size() - 1; ++i) {
        for (size_t j = i + 1; j < individual.path.size(); ++j) {
            vector<int> newPath = individual.path;
            reverse(newPath.begin() + i, newPath.begin() + j + 1);
            int newFitness = evaluatePath(newPath);
            if (newFitness < individual.fitness) {
                individual.path = newPath;
                individual.fitness = newFitness;
            }
        }
    }
}

void replaceVertexImprovement(Individual& individual) {
    for (size_t i = 1; i < individual.path.size() - 1; ++i) {
        for (const auto& edge : graph[individual.path[i - 1]]) {
            vector<int> newPath = individual.path;
            newPath[i] = edge.to;
            int newFitness = evaluatePath(newPath);
            if (newFitness < individual.fitness) {
                individual.path = newPath;
                individual.fitness = newFitness;
            }
        }
    }
}

void geneticAlgorithm(int start, int end, const string& outputFilename) {
    vector<Individual> population = initializePopulation(start, end);
    Individual bestIndividual = *min_element(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
        return a.fitness < b.fitness;
    });

    ofstream outFile(outputFilename);
    outFile << "Generation,BestFitness" << endl;

    for (int generation = 0; generation < MAX_GENERATIONS; ++generation) {
        vector<Individual> newPopulation;

        sort(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
            return a.fitness < b.fitness;
        });

        for (int i = 0; i < ELITE_COUNT; ++i) {
            newPopulation.push_back(population[i]);
        }

        while (newPopulation.size() < POPULATION_SIZE) {
            int parent1 = rand() % population.size();
            int parent2 = rand() % population.size();

            // Схрещування
            Individual child = orderedCrossover(population[parent1], population[parent2], start, end);

            // Мутація
            swapMutation(child);

            // Локальне покращення
            replaceVertexImprovement(child);

            newPopulation.push_back(child);
        }

        population = newPopulation;

        Individual currentBest = *min_element(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
            return a.fitness < b.fitness;
        });
        if (currentBest.fitness < bestIndividual.fitness) {
            bestIndividual = currentBest;
        }

        outFile << generation + 1 << "," << bestIndividual.fitness << endl;
        if (generation % 10 == 0) {
            cout << "Generation " << generation + 1 << ": Best fitness = " << bestIndividual.fitness << endl;
        }
    }

    outFile.close();
    cout << "Best path fitness: " << bestIndividual.fitness << endl;
}

int main() {
    generateGraph();
    int start = 0;
    int end = NUM_VERTICES - 1;
    geneticAlgorithm(start, end, "results.csv");
    return 0;
}