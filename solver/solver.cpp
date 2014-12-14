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

const double dt = 0.31622776601683794; // sqrt(0.1)
const double rad = 5;

mt19937 ranEng;
uniform_real_distribution<> dist(0, 20);

const size_t populationSize = 300;

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
        speed = prev.speed + a * dt;
        pos = prev.pos + prev.speed * dt + 0.5 * a * dt * dt;
        visited.resize(prev.visited.size());
        for (size_t i = 0; i < visited.size(); ++i) {
            if (prev.visited[i]) {
                visited[i] = true;
            }
            else if (touched(pos, i)) {
                visited[i] = true;
            }
        }
        time = prev.time + 1;
    }
    Vector pos;
    Vector speed;
    vector<bool> visited;
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
    //double dist = distSqr(target, st.pos);
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

void initPilot(Pilot & p) {
    for (auto & c : p) {
        c.P = dist(ranEng);
        c.I = dist(ranEng);
        c.D = dist(ranEng);
        c.predict = dist(ranEng);
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

double evaluate(const Pilot & p) {
    actions.clear();
    State start;
    start.pos = Vector(160, 120);
    start.visited.resize(cities.size());
    State next = start;
    State tmp = start;

    for (size_t i = 0; i < cities.size(); ++i) {
        integral = 0;
        prevError = 0;
        do {
            getFout(tmp, next, cities[i], p[i]);
            tmp = next;
        } while (!next.visited[i] && next.time < 2000);
    }
    if (next.time == 2000)
        return 0.01;
    double visited = accumulate(next.visited.begin(), next.visited.end(), 0.0);
    return 1.0/actions.size() + visited;
}

struct FitPilot {
    FitPilot (Pilot & p, double f) : fitness(f), pilot(&p) { }
    double fitness;
    Pilot * pilot;
    bool operator< (const FitPilot & r) const {
        return fitness < r.fitness;
    }
};

int findSelect(const vector<bool> & selected, int & index, bool p) {
    while (index < populationSize && selected[index] == p) {
        ++index;
    }
    if (index == populationSize) {
        return -1;
    }
    return index++;
}

void crossover(const Pilot & p1, const Pilot & p2, Pilot & c1, Pilot & c2) {
    uniform_int_distribution<> dis(0, cities.size() - 1);
    int i1 = dis(ranEng);
    int i2 = dis(ranEng);
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

void evaluateIndividuals() {
    vector<double> fitness;
    fitness.resize(populationSize);
    for (size_t i = 0; i < populationSize; ++i) {
       fitness[i] = evaluate(population[i]);
    }
    cout << *max_element(fitness.begin(), fitness.end()) << endl;
    discrete_distribution<int> selection(fitness.begin(), fitness.end());
    vector<bool> selected;
    selected.resize(populationSize);
    cout << "selection:" ;
    for (size_t i = 0; i < populationSize / 2; ++i) {
        int s;
        do {
            s = selection(ranEng);
        } while (selected[s]);
        cout << s << " ";
        selected[s] = true;
    }
    cout << endl;
    int firstRemaining = 0;
    int firstFree = 0;
    for(;;) {
        int p1 = findSelect(selected, firstRemaining, true);
        int p2 = findSelect(selected, firstRemaining, true);
        int c1 = findSelect(selected, firstFree, false);
        int c2 = findSelect(selected, firstFree, false);
        if (max(fitness[c1], fitness[c2]) > min(fitness[p1], fitness[2]))
            continue;
        if (p2 == -1 || c2 == -1)
            break;
        cout << "fit: " << fitness[p1] << " " << fitness[p2] << " " << fitness[c1]
            << " " << fitness[c2] << endl;
        crossover(population[p1], population[p2], population[c1], population[c2]);
    }
    int free;
    while ((free = findSelect(selected, firstFree, false)) != -1) {
        initPilot(population[free]);
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
    for (size_t i = 0; i < 100; ++i) {
        evaluateIndividuals();
    }

    Pilot & p = findBestPilot();
    evaluate(p);

    printResult();

    return 0;
}
