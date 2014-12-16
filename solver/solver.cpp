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
#include <iterator>
#include "vector.h"
#include "util.h"

using namespace std;

mt19937 ranEng;

// Globals here

struct Trajectory {
    double fitness;
    int size;
    vector<Action> acts;
    Vector pos;
    Vector speed;
    int cityIndex;
};


double globalMax;

vector<Action> globalActs;

int windowStart;
int windowEnd;
const size_t populationSize = 10;
const size_t numberOfParents = 4;

const int timeWindow = 50;
const int totalSteps = 1000;
const int timeIncresement = 15;
const int generations = 1000;


vector<Vector> cities;
vector<Trajectory> population;

// distributions used
normal_distribution<double> normalDist(0,1);
uniform_int_distribution<> geneDist(0,4);
uniform_int_distribution<> disMating(0, numberOfParents - 1);
uniform_int_distribution<> disPopulation(0, populationSize - 1);
binomial_distribution<> swapDist(80, 0.5);
bernoulli_distribution flipDist(0.07);
uniform_int_distribution<> swapChooseDist(0, timeWindow-2);
uniform_int_distribution<> flipChooseDist(0, timeWindow-1);

void mutateSwap(Trajectory & t) {
    auto index = swapChooseDist(ranEng) + windowStart;
    if (index > windowEnd - 2)
        return;
    auto tmp = t.acts[index+1];
    t.acts[index+1] = t.acts[index];
    t.acts[index] = tmp;
}

void mutateFlip(Trajectory & t) {
    auto index = flipChooseDist(ranEng) + windowStart;
    if (index > windowEnd - 1)
        return;
    t.acts[index] = static_cast<Action>(geneDist(ranEng));
}

void initPopulation() {
    random_device d;
    ranEng.seed(d());
    Trajectory t;
    t.pos = Vector(160,120);
    t.speed = 0;
    t.cityIndex = 0;
#if 0
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
#endif

    t.size = t.acts.size();

    t.size = totalSteps;

    population.resize(populationSize);
    for (auto & p : population) {
        p = t;
        p.acts.resize(p.size);
        for (auto & a : p.acts) {
            a = static_cast<Action>(geneDist(ranEng));
        }
        /*
        for (int i = 0; i < 5; ++i)
            mutateFlip(p);
            */
    }
}

inline int simulate(vector<Action> & acts, Vector & pos, Vector & speed,
              int & cityIndex)
{
    if (cityIndex == cities.size())
        return windowStart;
    for (int i = windowStart; i < windowEnd; ++i) {
        applyAction(pos, speed, acts[i]);
        if (touched(pos, cities[cityIndex])) {
            ++cityIndex;
            if (cityIndex == cities.size()) {
                return i;
                break;
            }
        }
    }
    return windowEnd;
}

Trajectory & findBestSoltion();

void moveTime(vector<Trajectory> & pop) {
    int timeend = windowEnd;
    windowEnd = windowStart + timeIncresement;
    if (windowEnd > totalSteps)
        windowEnd = totalSteps;
    for (auto & p : pop) {
        p.size = simulate(p.acts, p.pos, p.speed, p.cityIndex);

        double e = 0;
        if (p.cityIndex < cities.size()) {
            e = 1.0 / sqrt(distSqr(p.pos, cities[p.cityIndex]));
        }
        double fitness = (p.cityIndex/1.0) + e + ( 1.0 / p.size);
        assert(fitness > 0 && fitness < p.cityIndex + 1);
        p.fitness = fitness;
    }
    windowEnd = timeend;

    Trajectory & best = findBestSoltion();
    population[0] = best;
    copy_n(best.acts.begin() + windowStart, timeIncresement,
           back_inserter(globalActs));

    for (auto & p : pop) {
        p.pos = best.pos;
        p.speed = best.speed;
        p.cityIndex = best.cityIndex;
        p.size = best.size;
    }
}

void fitness(Trajectory & t) {
    Vector pos = t.pos;
    Vector speed = t.speed;

    int cityIndex = t.cityIndex;

    t.size = simulate(t.acts, pos, speed, cityIndex);

    double e = 0;

    if (cityIndex < cities.size()) {
        e = 1.0 / sqrt(distSqr(pos, cities[cityIndex]));
    }
    double fitness = (cityIndex/1.0) + e + ( 1.0 / t.size);
    assert(fitness > 0 && fitness < cityIndex + 1);
    t.fitness = fitness;
}

void evaluateAll() {
    for (auto & p : population) {
        fitness(p);
        if (p.fitness > globalMax) {
            globalMax = p.fitness;
            population[0] = p;
        }
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
    for (int i = 1; i < populationSize; ++i) {
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

inline void calcStatistics(const vector<Trajectory> & pop, double & avg,
                           double & max)
{
    max = 0;
    double sum = 0;
    for (const auto & p : pop) {
        sum += p.fitness;
        if (p.fitness > max) {
            max = p.fitness;
#if 0
            if (max > globalMax) {
                globalBest = p;
                globalMax = max;
            }
#endif
        }
    }
    avg = sum / pop.size();
}

int main(int argc, char**argv) {

    readInput(argv[1]);

    initPopulation();
    vector<Trajectory> parents;
    parents.resize(numberOfParents);
    double avg, max;
    bool printStat = false;
    for (size_t t = 0; t + timeIncresement < totalSteps; t += timeIncresement) {
        if (t != 0) {
            cout << "Increase timewindow!" << endl;
            moveTime(population);
        }
        windowStart = t;
        windowEnd = t + timeWindow;
        if (windowEnd > totalSteps)
            windowEnd = totalSteps;

        for (size_t i = 0; i < generations; ++i) {
            evaluateAll();
            printStat = !(i % 0xfff);
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
        }
    }

    evaluateAll();

    Trajectory & p = findBestSoltion();
    copy_n(p.acts.begin() + windowStart, windowEnd - windowStart,
           back_inserter(globalActs));

    Trajectory fin;
    fin.pos = Vector(160,120);
    fin.speed = 0;
    fin.cityIndex = 0;
    fin.acts = globalActs;
    windowStart = 0;
    windowEnd = totalSteps;
    fitness(fin);
    cout << fin.fitness << endl;
    globalActs.resize(fin.size + 4);
    cout << globalActs.size() << endl;
    printResult(globalActs);

    return 0;
}
