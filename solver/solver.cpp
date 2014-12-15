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

using namespace std;

const double dt = 0.1; // sqrt(0.1)
const double rad = 5;

mt19937 ranEng;
uniform_real_distribution<> distP(0, 50);
uniform_real_distribution<> distI(0, 40);
uniform_real_distribution<> distD(0, 50);
uniform_real_distribution<> distPre(0, 20);

const size_t populationSize = 200;
const size_t numberOfParents = 160;
const size_t numberOfNew = 80;
const size_t cutOff = 78;

int optMax = 3;
int optMin = 0;

vector<Vector> cities;

inline double distSqr(const Vector & l, const Vector & r) {
    double xd = l.x - r.x;
    double yd = l.y - r.y;
    return xd*xd + yd*yd;
}

enum Action {
    NOP,
    UP,
    RIGHT,
    DOWN,
    LEFT
};

inline Vector toVect(Action & a) {
    switch (a) {
        case UP:
            return Vector(0,-1);
        case DOWN:
            return Vector(0,1);
        case LEFT:
            return Vector(-1,0);
        case RIGHT:
            return Vector(1,0);
    }
    return Vector(0,0);
}

inline bool touched(const Vector & p, size_t index) {
    double xd = cities[index].x - p.x;
    double yd = cities[index].y - p.y;
    return xd * xd + yd * yd < rad * rad;
}

class State {
public:
    State () : time(0) { }
    State (const State & prev, Action act) : action(act), prevState(&prev) {
        Vector a = toVect(act);
        speed = prev.speed + a;
        pos = prev.pos + prev.speed * dt + 0.5 * a * dt * dt;
        time = prev.time + 1;
    }
    Vector pos;
    Vector speed;
    const State * prevState = nullptr;
    Action action;
    int time;
};

vector<Action> actions;

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

void getFout(const State & st, State & next, const Vector & target,
             const Chromosome & chromo)
{
    Vector predictedPos = st.pos + (st.speed * chromo.predict);

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
    next = State(st, act);
    actions.push_back(act);
}

using Pilot = vector<Chromosome>;

void printPilot(const Pilot & pl) {
    for (auto & c : pl) {
        cout << c;
    }
    cout << endl;
}

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

void initPilot(Pilot & p) {
    for (auto & c : p) {
        initChromo(c);
    }
}

void initPopulation() {
    random_device d;
    ranEng.seed(d());
    population.resize(populationSize);
    size_t nrOfChromos = cities.size();
    for (auto & p : population) {
        p.resize(nrOfChromos);
        initPilot(p);
    }
}

int countFail=0;

double evaluate(const Pilot & p) {
    actions.clear();
    State start;
    start.pos = Vector(160, 120);
    State next;
    State cur = start;

    double time;
    int cutMult = 0;

    int visited = 0;
    bool cityReached;

    for (size_t i = 0; i < cities.size(); ++i) {
        integral = 0;
        prevError = 0;
        ++cutMult;
        cityReached = false;
        do {
            getFout(cur, next, cities[i], p[i]);
            cur = next;
            if (touched(cur.pos, i)) {
                cityReached = true;
                ++visited;
            }
        } while (!cityReached && cur.time < cutOff * cutMult);
        if (cur.time < cutOff * cutMult) {
            time = next.time;
        }
        else {
            break;
        }
    }
    if (visited < cities.size()) {
        countFail++;
        double t = (1.0 / (time + cutOff)) + visited;
        //cout << "fail: " << t << endl;
        return t;
    }
    return (1.0 / (time + cutOff)) + visited;
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

uniform_int_distribution<> crossOverDist(0, optMax);
void crossover(const Pilot & p1, const Pilot & p2, Pilot & c1, Pilot & c2) {
    int i1 = crossOverDist(ranEng);
    int i2 = crossOverDist(ranEng);
    int from = min(i1,i2);
    int to = max(i1,i2);
    for (int i = 0; i < cities.size(); ++i) {
        if (i >= to && i <= from) {
            c1[i] = p2[i];
            c2[i] = p1[i];
        }
        else {
            c1[i] = p1[i];
            c2[i] = p2[i];
        }
    }
}

void evaluateIndividuals(vector<double> & fitness) {
    countFail = 0;
    for (size_t i = 0; i < populationSize; ++i) {
       fitness[i] = evaluate(population[i]);
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
            mutateChromo(p[disPart(ranEng)]);
            population[i] = p;
        }
        else {
            crossover(*parents[disParents(ranEng)].pilot,
                    *parents[disParents(ranEng)].pilot,
                    population[i], population[i+1]);
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

void printResult() {
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
        fit = evaluate(p);
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
    for (size_t i = 0; i < 1000; ++i) {
        evaluateIndividuals(fitness);
        naturalSelection(fitness, parents, parentFitness);
        reproduce(parentFitness);
        double bestFitness = *max_element(fitness.begin(), fitness.end());
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
    cout << evaluate(p) << endl;
    cout << actions.size() << endl;

    printResult();

    return 0;
}
