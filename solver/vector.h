#ifndef VECTOR_H
#define VECTOR_H

struct Vector {
    Vector () : x(0), y(0) {}
    Vector(double x, double y) : x(x), y(y) {}
    double x;
    double y;

    Vector & operator= (const double & d) {
        x = d;
        y = d;
        return *this;
    }
    Vector & operator+= (const Vector & rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }
    friend inline Vector operator+ (const Vector & lhs, const Vector & rhs) {
        return Vector(lhs.x + rhs.x, lhs.y + rhs.y);
    }
    friend inline Vector operator- (const Vector & lhs, const Vector & rhs) {
        return Vector(lhs.x - rhs.x, lhs.y - rhs.y);
    }
    friend inline Vector operator* (const Vector & lhs, double s) {
        return Vector(lhs.x * s, lhs.y * s);
    }
    friend inline Vector operator/ (const Vector & lhs, double s) {
        return Vector(lhs.x / s, lhs.y / s);
    }
    friend inline Vector operator* (double s, const Vector & rhs) {
        return Vector(rhs.x * s, rhs.y * s);
    }
};

std::ostream& operator<<(std::ostream& os, const Vector& v)
{
    // write obj to stream
    os << v.x << " " << v.y;
    return os;
}

#endif
