#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <queue>
#include <cmath>

using namespace std;

const double dt = 0.31622776601683794; // sqrt(0.1)
const double rad = 5;

struct Vector {
    Vector () : x(0), y(0) {}
    Vector(double x, double y) : x(x), y(y) {}
    double x;
    double y;

    Vector & operator+= (const Vector & rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }
    friend Vector operator+ (const Vector & lhs, const Vector & rhs) {
        return Vector(lhs.x + rhs.x, lhs.y + rhs.y);
    }
    friend Vector operator- (const Vector & lhs, const Vector & rhs) {
        return Vector(lhs.x - rhs.x, lhs.y - rhs.y);
    }
    friend Vector operator* (const Vector & lhs, double s) {
        return Vector(lhs.x * s, lhs.y * s);
    }
    friend Vector operator/ (const Vector & lhs, double s) {
        return Vector(lhs.x / s, lhs.y / s);
    }
    friend Vector operator* (double s, const Vector & rhs) {
        return Vector(rhs.x * s, rhs.y * s);
    }
};

std::ostream& operator<<(std::ostream& os, const Vector& v)
{
    // write obj to stream
    os << v.x << " " << v.y;
    return os;
}

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

struct chromosome {
    double P;
    double I;
    double D;
    double predict;
};

chromosome chromo = { 10, 0, 3, 2};
Vector integral;
Vector prevError;

void getFout(const State & st, State & next, size_t target) {
    Vector t = cities[target];
    double dist = distSqr(t, st.pos);
    Vector predictedPos = st.pos + (st.speed * chromo.predict);

    Vector posError = t - predictedPos;
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

int main(int argc, char**argv) {

    ifstream f;
    f.open(string(argv[1]));
    int n;
    f >> n;

    cities.resize(n);
    for (size_t i = 0; i < n; ++i) {
        f >> cities[i].x >> cities[i].y;
    }

    f.close();

    State start;
    start.pos = Vector(160, 120);
    start.visited.resize(cities.size());
    State next = start;
    State tmp = start;

    for (size_t i = 0; i < cities.size(); ++i) {
        do {
        getFout(tmp, next, i);
        tmp = next;
        } while (!next.visited[i] && next.time < 2000);
        cout << "Visited city " << i << endl;
    }

    ofstream out;
    out.open("result.txt");
    out << actions.size() << endl;

    for (auto a: actions) {
        out << a << endl;
    }

    out.close();

    return 0;
}
