#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <queue>
#include <cmath>
#include <random>
#include <functional>
#include <queue>
#include "vector.h"

using namespace std;

const double dt = 0.31622776601683794; // sqrt(0.1)
const double rad = 5;

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
    cout << act << "\t" << output << "\t" << next.pos <<  endl;
    actions.push_back(act);
}

using Pilot = vector<Chromosome>;

mt19937 ranEng;

const size_t populationSize = 500;
vector<Pilot> population;

void initPopulation() {
    random_device d;
    ranEng.seed(d());
    uniform_real_distribution<> dist(0, 20);
    auto g = bind(dist, ranEng);
    population.resize(populationSize);
    size_t nrOfChromos = cities.size();
    for (auto & p : population) {
        p.resize(nrOfChromos);
        for (auto & c : p) {
            c.P = g();
            c.I = g();
            c.D = g();
            c.predict = g();
        }
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
    return 1.0/actions.size();
}

struct FitPilot {
    FitPilot (Pilot & p, double f) : fitness(f), pilot(&p) { }
    double fitness;
    Pilot * pilot;
    bool operator< (const FitPilot & r) const {
        return fitness < r.fitness;
    }
};

priority_queue<FitPilot> testedPilots;

void evaluateIndividuals() {
    testedPilots=priority_queue<FitPilot>();
    for (auto & p : population) {
        double fitness = evaluate(p);
        testedPilots.emplace(FitPilot(p,fitness));
    }
}

void naturalSelection() {
}

void recombine() {
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

int main(int argc, char**argv) {

    readInput(argv[1]);

    initPopulation();
    for (size_t i = 0; i < 1; ++i) {
        evaluateIndividuals();
        naturalSelection();
        recombine();
    }

    evaluateIndividuals();
    Pilot & p = *testedPilots.top().pilot;
    evaluate(p);

    printResult();

    return 0;
}
