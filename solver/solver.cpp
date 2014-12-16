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

const size_t populationSize = 50;
const size_t numberOfParents = 5;

struct Trajectory {
    double fitness;
    int size;
    vector<Action> acts;
    Vector pos;
    Vector speed;
    int cityIndex;
    int steps;
};


double globalMax;
Trajectory globalBest;

vector<Action> globalActs;

int windowStart;
int windowEnd;
const int timeWindow = 100;
const int totalSteps = 1000;
const int timeIncresement = 10;

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
    uniform_int_distribution<> d(windowStart, windowEnd);
    auto index = d(ranEng);
    auto tmp = t.acts[index+1];
    t.acts[index+1] = t.acts[index];
    t.acts[index] = tmp;
}

void mutateFlip(Trajectory & t) {
    uniform_int_distribution<> d(windowStart, windowEnd);
    t.acts[d(ranEng)] = static_cast<Action>(geneDist(ranEng));
}

void initPopulation() {
    random_device d;
    ranEng.seed(d());
    Trajectory t;
    t.pos = Vector(160,120);
    t.speed = 0;
    t.cityIndex = 0;
    t.steps = 0;
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

inline void simulate(vector<Action> & acts, Vector & pos, Vector & speed,
              int & cityIndex, int & steps)
{
    for (int i = steps; i < windowEnd; ++i) {
        applyAction(pos, speed, acts[i]);
        if (touched(pos, cities[cityIndex])) {
            ++cityIndex;
            if (cityIndex == cities.size()) {
                break;
            }
        }
        ++steps;
    }
}

Trajectory & findBestSoltion();

void moveTime(vector<Trajectory> & pop) {
    int timeend = windowEnd;
    windowEnd = windowStart + timeIncresement;
    for (auto & p : pop) {
        simulate(p.acts, p.pos, p.speed, p.cityIndex, p.steps);
        p.size = p.steps;

        double e = 0;

        if (p.cityIndex != cities.size()) {
            e = 1.0 / sqrt(distSqr(p.pos, cities[p.cityIndex]));
        }
        double fitness = (p.cityIndex/1.0) + e + ( 1.0 / p.size);
        assert(fitness > 0 && fitness < p.cityIndex + 1);
        p.fitness = fitness;
    }
    windowEnd = timeend;

    Trajectory & best = findBestSoltion();
    copy_n(best.acts.begin() + windowStart, timeIncresement,
           back_inserter(globalActs));
#if 0
    cout << "Local: ";
    copy_n(globalActs.begin(), windowStart,
           ostream_iterator<Action>(cout, " "));
    copy_n(best.acts.begin() + windowStart, timeWindow,
           ostream_iterator<Action>(cout, " "));
    cout << endl;
    cout << "Global: ";
    copy(globalActs.begin(), globalActs.end(),
           ostream_iterator<Action>(cout, " "));
    cout << endl;
#endif
    for (auto & p : pop) {
        //simulate(p.acts, p.pos, p.speed, p.cityIndex, p.steps);
        //p.size = p.steps;
        p.pos = best.pos;
        p.speed = best.speed;
        p.cityIndex = best.cityIndex;
        p.steps = best.steps;
        p.size = best.size;
    }
}

void fitness(Trajectory & t) {
    Vector pos = t.pos;
    Vector speed = t.speed;

    int cityIndex = t.cityIndex;
    int steps = t.steps;

    simulate(t.acts, pos, speed, cityIndex, steps);

    t.size = steps;

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
    for (size_t t = 0; t + timeWindow < totalSteps; t += timeIncresement) {
        if (t != 0) {
            cout << "Increase timewindow!" << endl;
            moveTime(population);
        }
        windowStart = t;
        windowEnd = t + timeWindow;

        for (size_t i = 0; i < 20000; ++i) {
            evaluateAll();
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
        }
    }

    evaluateAll();

    Trajectory & p = findBestSoltion();
    copy_n(p.acts.begin() + windowStart, timeWindow,
           back_inserter(globalActs));

    cout << p.fitness << endl;
    globalActs.resize(p.size);
    cout << globalActs.size() << endl;
    printResult(globalActs);

    return 0;
}
