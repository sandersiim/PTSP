#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <queue>
#include <cmath>
#include <random>
#include <functional>
#include <algorithm>
#include <queue>
#include "vector.h"
#include "util.h"

using namespace std;

mt19937 ranEng;

// Globals here

const size_t populationSize = 50;
const size_t numberOfParents = 3;
const size_t numberOfNew = 80;

struct Trajectory {
    double fitness;
    int size;
    vector<Action> acts;
};

vector<Vector> cities;
vector<Trajectory> population;

// distributions used
normal_distribution<double> normalDist(0,1);
uniform_int_distribution<> geneDist(0,4);
uniform_int_distribution<> disMating(0, numberOfParents - 1);
uniform_int_distribution<> disPopulation(0, populationSize - 1);

void mutateSwap(Trajectory & t) {
    uniform_int_distribution<> d(2, t.size);
    auto index = d(ranEng);
    auto tmp = t.acts[index-2];
    t.acts[index-2] = t.acts[index-1];
    t.acts[index-1] = tmp;
}

void mutateFlip(Trajectory & t) {
    uniform_int_distribution<> d(1, t.size);
    t.acts[d(ranEng) - 1] = static_cast<Action>(geneDist(ranEng));
}

void initPopulation() {
    random_device d;
    ranEng.seed(d());
    Trajectory t;
    ifstream f;
    f.open("best_pid.txt");
    int n;
    f >> n;
    int a;
    for (size_t i = 0; i < n; ++i) {
        f >> a;
        t.acts.emplace_back(static_cast<Action>(a));
    }
    f.close();

    t.size = t.acts.size();

    population.resize(populationSize);
    for (auto & p : population) {
        p = t;
        mutateSwap(p);
    }
}

void fitness(Trajectory & t) {
    Vector pos(160, 120);
    Vector speed;

    int cityIndex = 0;
    int steps = 0;

    t.size = t.acts.size();
    for (auto & a : t.acts) {
        if (touched(pos, cities[cityIndex])) {
            ++cityIndex;
            if (cityIndex == cities.size()) {
                t.size = steps + 1;
                break;
            }
        }
        ++steps;
        applyAction(pos, speed, a);
    }

    double e = 0;

    if (cityIndex != cities.size()) {
        e = 1.0 / sqrt(distSqr(pos, cities[cityIndex]));
    }
    double fitness = cityIndex + e + ( 1.0 / t.size);
    t.fitness = fitness;
}

void evaluateAll() {
    for (auto & p : population) {
        fitness(p);
    }
}

Trajectory & tournamentSelect() {
    Trajectory & t1 = population[disPopulation(ranEng)];
    Trajectory & t2 = population[disPopulation(ranEng)];
    Trajectory & t3 = population[disPopulation(ranEng)];
    double f1 = t1.fitness;
    double f2 = t2.fitness;
    double f3 = t3.fitness;
    if (f1 > f2) {
        return f1 > f3 ? t1 : t3;
    }
    else {
        return f2 > f3 ? t2 : t3;
    }
}

void naturalSelection(vector<Trajectory> & parents)
{
    parents.clear();
    while (parents.size() < numberOfParents) {
        parents.push_back(tournamentSelect());
    }
}

void reproduce(vector<Trajectory> & matingPool) {
    for (int i = 0; i < populationSize; ++i) {
        population[i] = matingPool[disMating(ranEng)];
        mutateSwap(population[i]);
        //mutateFlip(population[i]);
    }
}

void readInput(const char *file) {
    ifstream f;
    f.open(file);
    int n;
    f >> n;

    cities.resize(n);
    for (size_t i = 0; i < n; ++i) {
        f >> cities[i].x >> cities[i].y;
    }

    f.close();
}

void printResult(vector<Action> & actions) {
    ofstream out;
    out.open("result.txt");
    out << actions.size() << endl;

    for (auto a: actions) {
        out << a << endl;
    }

    out.close();
}

Trajectory & findBestSoltion() {
    double bestfit = 0;
    Trajectory * best = nullptr;
    double fit;
    for (auto & p : population) {
        if (p.fitness > bestfit) {
            bestfit = p.fitness;
            best = &p;
        }
    }
    return *best;
}

int main(int argc, char**argv) {

    readInput(argv[1]);

    initPopulation();
    vector<Trajectory> parents;
    parents.resize(numberOfParents);
    for (size_t i = 0; i < 1000; ++i) {
        evaluateAll();
        naturalSelection(parents);
        reproduce(parents);
    }

    Trajectory & p = findBestSoltion();
    cout << p.fitness << endl;
    cout << p.size << endl;
    printResult(p.acts);

    return 0;
}
