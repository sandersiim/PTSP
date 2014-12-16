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
#include <cassert>
#include "vector.h"
#include "util.h"

using namespace std;

mt19937 ranEng;

// Globals here

const size_t populationSize = 25;
const size_t numberOfParents = 5;

struct Trajectory {
    double fitness;
    int size;
    vector<Action> acts;
};


double globalMax;
Trajectory globalBest;

vector<Vector> cities;
vector<Trajectory> population;

// distributions used
normal_distribution<double> normalDist(0,1);
uniform_int_distribution<> geneDist(0,4);
uniform_int_distribution<> disMating(0, numberOfParents - 1);
uniform_int_distribution<> disPopulation(0, populationSize - 1);
binomial_distribution<> swapDist(40, 0.3);
bernoulli_distribution flipDist(0.05);

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
        applyAction(pos, speed, a);
        if (touched(pos, cities[cityIndex])) {
            ++cityIndex;
            if (cityIndex == cities.size()) {
                t.size = steps + 1;
                break;
            }
        }
        ++steps;
    }

    double e = 0;

    if (cityIndex != cities.size()) {
        e = 1.0 / sqrt(distSqr(pos, cities[cityIndex]));
    }
    double fitness = (cityIndex/1.0) + e + ( 1.0 / t.size);
    assert(fitness > 0 && fitness < cityIndex + 1);
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
    for (int i = 0; i < numberOfParents; ++i) {
        parents[i] = tournamentSelect();
    }
}

void reproduce(vector<Trajectory> & matingPool) {
    for (int i = 0; i < populationSize; ++i) {
        population[i] = matingPool[disMating(ranEng)];
        int swaps = swapDist(ranEng);
        while (swaps--) {
            mutateSwap(population[i]);
        }
        if (flipDist(ranEng))
            mutateFlip(population[i]);
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
    double bestfit = globalBest.fitness;
    Trajectory * best = &globalBest;
    double fit;
    for (auto & p : population) {
        if (p.fitness > bestfit) {
            bestfit = p.fitness;
            best = &p;
        }
    }
    return *best;
}

inline void calcStatistics(const vector<Trajectory> & pop, double & avg,
                           double & max)
{
    max = 0;
    double sum = 0;
    for (const auto & p : pop) {
        sum += p.fitness;
        if (p.fitness > max) {
            max = p.fitness;
            if (max > globalMax) {
                globalBest = p;
                globalMax = max;
            }
        }
    }
    avg = sum / pop.size();
}

int main(int argc, char**argv) {

    readInput(argv[1]);

    initPopulation();
    evaluateAll();
    vector<Trajectory> parents;
    parents.resize(numberOfParents);
    double avg, max;
    bool printStat = false;
    for (size_t i = 0; i < 200000; ++i) {
        printStat = !(i % 0xff);
        if (printStat) {
            calcStatistics(population, avg, max);
            cout << "Pop: avg: " << avg << "\tmax: " << max << "\t";
        }
        naturalSelection(parents);
        if (printStat) {
            calcStatistics(parents, avg, max);
            cout << "Par: avg: " << avg << "\tmax: " << max << endl;
        }
        reproduce(parents);
        evaluateAll();
    }

    Trajectory & p = findBestSoltion();
    cout << p.fitness << endl;
    cout << p.size << endl;
    p.acts.resize(p.size);
    printResult(p.acts);

    return 0;
}
