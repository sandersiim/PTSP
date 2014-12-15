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
uniform_real_distribution<> distP(-8, 80);
uniform_real_distribution<> distI(-2, 2);
uniform_real_distribution<> distD(-5, 50);
uniform_real_distribution<> distPre(-2, 4);

const size_t populationSize = 20000;
const size_t numberOfParents = 1400;
const size_t numberOfNew = 8000;
const size_t cutOff = 78;

int optMax = 3;
int optMin = 0;

vector<Vector> cities;

struct Chromosome {
    double P;
    double I;
    double D;
    double predict;
};

std::ostream& operator<<(std::ostream& os, const Chromosome& c)
{
    // write obj to stream
    os << "{" << c.P << ", " << c.I << ", " << c.D << ", " << c.predict << "}";
    return os;
}

Vector integral;
Vector prevError;

void getFout(Vector & pos, Vector & speed, const Vector & target,
             const Chromosome & chromo, vector<Action> *actions)
{
    Vector predictedPos = pos + (speed * chromo.predict);

    Vector posError = target - predictedPos;
    Vector derivate = (posError - prevError) / dt;
    integral += posError * dt;

    Vector output = chromo.P * posError + chromo.I * integral + chromo.D * derivate;

    prevError = posError;

    double ax = abs(output.y) / abs(output.x);
    double ay = abs(output.x) / abs(output.y);

    Action act = NOP;
    if (ax > 0.9) {
        act = output.y > 0 ? DOWN : UP;
    }
    if (ay > 0.9) {
        act = output.x > 0 ? RIGHT : LEFT;
    }
    applyAction(pos, speed, act);
    if (actions != nullptr) {
        actions->push_back(act);
    }
}

using Pilot = Chromosome;

Pilot bestPilot;
vector<Pilot> population;
normal_distribution<double> normalDist(0,1);
uniform_int_distribution<> geneDist(0,3);

void mutateChromo(Chromosome & c) {
    int gene = geneDist(ranEng);
    double value = normalDist(ranEng);

    switch (gene) {
        case 0:
            c.P += value;
            break;
        case 1:
            c.I += value;
            break;
        case 2:
            c.D += value;
            break;
        case 3:
            c.predict += value;
            break;
    }
}

void initChromo(Chromosome & c) {
    c.P = distP(ranEng);
    c.I = distI(ranEng);
    c.D = distD(ranEng);
    c.predict = distPre(ranEng);
}

void initPopulation() {
    random_device d;
    ranEng.seed(d());
    population.resize(populationSize);
    for (auto & p : population) {
        initChromo(p);
    }
}

int countFail=0;

double evaluate(const Pilot & p, vector<Action> *path) {
    Vector pos(160, 120);
    Vector speed;

    double time;
    int cutMult = 0;

    int visited = 0;
    bool cityReached;

    int actions = 0;
    for (size_t i = 0; i < cities.size(); ++i) {
        integral = 0;
        prevError = 0;
        ++cutMult;
        cityReached = false;
        do {
            getFout(pos, speed, cities[i], p, path);
            ++actions;
            if (touched(pos, cities[i])) {
                cityReached = true;
                ++visited;
            }
        } while (!cityReached && actions < cutOff * cutMult);
        if (actions < cutOff * cutMult) {
            time = actions;
        }
        else {
            break;
        }
    }
    if (visited < cities.size()) {
        countFail++;
        double t = (1.0 / (time + cutOff)) + visited / 1000.0;
        //cout << "fail: " << t << endl;
        return t;
    }
    return (1.0 / (time + cutOff)) + visited / 1000.0;
}

struct FitPilot {
    FitPilot () { }
    FitPilot (Pilot & p, double f) : fitness(f), pilot(&p) { }
    double fitness;
    Pilot * pilot;
    bool operator< (const FitPilot & r) const {
        return fitness < r.fitness;
    }
};

uniform_int_distribution<> crossOverDist(0, 3);
void crossover(const Pilot & p1, const Pilot & p2, Pilot & c1, Pilot & c2) {
}

void evaluateIndividuals(vector<double> & fitness) {
    countFail = 0;
    for (size_t i = 0; i < populationSize; ++i) {
       fitness[i] = evaluate(population[i], nullptr);
    }
}

void naturalSelection(vector<double> & fitness, vector<Pilot> & parents,
                      vector<FitPilot> & parentFitness)
{
    discrete_distribution<int> selection(fitness.begin(), fitness.end());
    for (size_t i = 0; i < numberOfParents; ++i) {
        int s = selection(ranEng);
        parents[i] = population[s];
        parentFitness[i] = FitPilot(parents[i], fitness[i]);
    }
}

uniform_int_distribution<> disParents(0, numberOfParents - 1);
uniform_int_distribution<> disPart(optMin, optMax);
void reproduce(vector<FitPilot> & parents) {
    for (int i = 0; i < populationSize; i += 2) {
        if (i < numberOfNew) {
            auto & p = *parents[disParents(ranEng)].pilot;
            mutateChromo(p);
            population[i] = p;
        }
        else {
            population[i] = *parents[disParents(ranEng)].pilot;
        }
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

Pilot & findBestPilot() {
    double bestfit = 0;
    Pilot * best = nullptr;
    double fit;
    for (auto & p : population) {
        fit = evaluate(p, nullptr);
        if (fit > bestfit) {
            bestfit = fit;
            best = &p;
        }
    }
    return *best;
}

int main(int argc, char**argv) {

    readInput(argv[1]);

    initPopulation();
    vector<double> fitness;
    vector<Pilot> parents;
    vector<FitPilot> parentFitness;
    fitness.resize(populationSize);
    parents.resize(numberOfParents);
    parentFitness.resize(numberOfParents);
    for (size_t i = 0; i < 10; ++i) {
        evaluateIndividuals(fitness);
        naturalSelection(fitness, parents, parentFitness);
        reproduce(parentFitness);
        auto bf = max_element(fitness.begin(), fitness.end());
        double bestFitness = *bf;
        double avgFitness = accumulate(fitness.begin(), fitness.end(), 0.0)
                / populationSize;
        cout << bestFitness << "\t" << avgFitness << "\tfailes: "
            << countFail << endl;
#if 1
        optMax = min<int>(cities.size() - 1, bestFitness + 1);
        optMin = max<int>(avgFitness - 2, 0);
#else
        optMax = cities.size() - 1;
        optMin = 0;
#endif
        crossOverDist = uniform_int_distribution<> (optMin, optMax);
        disPart = uniform_int_distribution<> (optMin, optMax);
    }

    Pilot & p = findBestPilot();
    vector<Action> path;
    cout << evaluate(p, &path) << endl;
    cout << path.size() << endl;
    cout << p << endl;
    printResult(path);

    return 0;
}
