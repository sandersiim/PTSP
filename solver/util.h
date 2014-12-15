#ifndef UTIL_H
#define UTIL_H

#include "vector.h"

const double dt = 0.1;
const double rad = 5;

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

inline bool touched(const Vector & p, const Vector & c) {
    double xd = c.x - p.x;
    double yd = c.y - p.y;
    return xd * xd + yd * yd < rad * rad;
}

inline void applyAction(Vector & pos, Vector & speed, Action act) {
        Vector a = toVect(act);
        pos += speed * dt + 0.5 * a * dt * dt;
        speed += a;
}

#endif
