#include <assert.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>

#include "vector.h"

using namespace std;

vector<Vector> cities;
using CityOrder = std::vector<Vector>;
Vector startPos(160, 120);

double distance(Vector& a, Vector& b) {
    return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y));
}

double dotProduct(Vector& a, Vector& b) {
    return a.x*b.x + a.y*b.y;
}

double magnitude(Vector& a) {
    return sqrt(a.x*a.x + a.y*a.y);
}

//angle between vectors. always returns value <= M_PI
//http://math.stackexchange.com/questions/361412/finding-the-angle-between-three-points
double angle(Vector& a, Vector& b) {
    return acos(dotProduct(a, b)/(magnitude(a)*magnitude(b)));
}

// from http://mathforum.org/library/drmath/view/54323.html
Vector circleCenterPoint(Vector& a, Vector& b, Vector& c) {
    float temp = b.x*b.x+b.y*b.y;
    double ab = (a.x*a.x + a.y*a.y - temp)/2.0d;
    double bc = (temp - c.x*c.x - c.y*c.y)/2.0d;
    double det = (a.x-b.x)*(b.y-c.y)-(b.x-c.x)*(a.y-b.y);
    if (abs(det) < 1.0e-6) {
        assert(false);
    }
    det = 1/det;
    double center_x = (ab*(b.y-c.y)-bc*(a.y-b.y))*det;
    double center_y = ((a.x-b.x)*bc-(b.x-c.x)*ab)*det;
    return Vector(center_x, center_y);
}

double evaluateCityOrder(CityOrder& order) {
    double totalTime = 0;
    Vector currentPos = startPos;
    Vector pos1, pos2, circleCenter;
    Vector vec_a,vec_b,vec_c,vec_up,vec_down;
    double radius, trajectory_angle;

    Vector newSpeed;
    Vector previousSpeed = Vector(0, 0);
    for (size_t i = 0; i < order.size()-1; i++) {
        pos1 = order[i];
        pos2 = order[i+1];

        circleCenter = circleCenterPoint(currentPos, pos1, pos2);
        radius = distance(currentPos, circleCenter);
        //speed = sqrt(radius);

        vec_a = Vector(currentPos.x - circleCenter.x, currentPos.y - circleCenter.y);
        vec_b = Vector(pos1.x - circleCenter.x, pos1.y - circleCenter.y);
        vec_c = Vector(pos2.x - circleCenter.x, pos2.y - circleCenter.y);
        vec_up = Vector(0, radius);
        vec_down = Vector(0, -radius);

        trajectory_angle = angle(vec_a, vec_b);

        bool clockWise, pos2IsInBetween;
        if (currentPos.x >= circleCenter.x && pos1.x >= circleCenter.x) {
            clockWise = currentPos.y < pos1.y;
            pos2IsInBetween = pos2.x >= circleCenter.x &&
                              pos2.y > min(currentPos.y, pos1.y) &&
                              pos2.y < max(currentPos.y, pos1.y);
        } else if (currentPos.x < circleCenter.x && pos1.x < circleCenter.x) {
            clockWise = currentPos.y > pos1.y;
            pos2IsInBetween = pos2.x < circleCenter.x &&
                              pos2.y > min(currentPos.y, pos1.y) &&
                              pos2.y < max(currentPos.y, pos1.y);
        } else if (currentPos.x >= circleCenter.x && pos1.x < circleCenter.x) {
            clockWise = angle(vec_a, vec_down)+angle(vec_b, vec_down) <= M_PI;
            if (clockWise) {
                if (pos2.x >= circleCenter.x) {
                    pos2IsInBetween = angle(vec_a,vec_down) > angle(vec_c,vec_down);
                } else {
                    pos2IsInBetween = angle(vec_b,vec_down) > angle(vec_c,vec_down);
                }
            } else {
                if (pos2.x >= circleCenter.x) {
                    pos2IsInBetween = angle(vec_a,vec_up) > angle(vec_c,vec_up);
                } else {
                    pos2IsInBetween = angle(vec_b,vec_up) > angle(vec_c,vec_up);
                }
            }
        } else {//if (currentPos.x < circleCenter.x && pos1.x >= circleCenter.x) {
            clockWise = angle(vec_a, vec_up)+angle(vec_b, vec_up) <= M_PI;
            if (clockWise) {
                if (pos2.x < circleCenter.x) {
                    pos2IsInBetween = angle(vec_a,vec_up) > angle(vec_c,vec_up);
                } else {
                    pos2IsInBetween = angle(vec_b,vec_up) > angle(vec_c,vec_up);
                }
            } else {
                if (pos2.x < circleCenter.x) {
                    pos2IsInBetween = angle(vec_a,vec_down) > angle(vec_c,vec_down);
                } else {
                    pos2IsInBetween = angle(vec_b,vec_down) > angle(vec_c,vec_down);
                }
            }
        }

        if (pos2IsInBetween) {
            trajectory_angle = 2*M_PI - trajectory_angle;
            clockWise = !clockWise;
        }

        //arc_length = trajectory_angle*radius;
        totalTime += trajectory_angle*sqrt(radius);

        if (clockWise && currentPos.y >= circleCenter.y) {
            newSpeed = Vector(1.0d, -(vec_a.x/vec_a.y));
        } else if (clockWise && currentPos.y < circleCenter.y) {
            newSpeed = Vector(-1.0d, (vec_a.x/vec_a.y));
        } else if (!clockWise && currentPos.y >= circleCenter.y) {
            newSpeed = Vector(-1.0d, (vec_a.x/vec_a.y));
        } else {
            newSpeed = Vector(1.0d, -(vec_a.x/vec_a.y));
        }

        newSpeed = newSpeed*(sqrt(radius)/magnitude(newSpeed));

        Vector differenceVector = Vector(newSpeed.x-previousSpeed.x, newSpeed.y-previousSpeed.y);
        totalTime += magnitude(differenceVector);

        if (clockWise && pos1.y >= circleCenter.y) {
            previousSpeed = Vector(1.0d, -(vec_b.x/vec_b.y));
        } else if (clockWise && pos1.y < circleCenter.y) {
            previousSpeed = Vector(-1.0d, (vec_b.x/vec_b.y));
        } else if (!clockWise && pos1.y >= circleCenter.y) {
            previousSpeed = Vector(-1.0d, (vec_b.x/vec_b.y));
        } else {
            previousSpeed = Vector(1.0d, -(vec_b.x/vec_b.y));
        }

        previousSpeed = previousSpeed*(sqrt(radius)/magnitude(previousSpeed));

        if (i == order.size()-1) {
            double bc_angle;
            if (clockWise) {
                if (pos1.x >= circleCenter.x && pos2.x >= circleCenter.x) {
                    assert(pos1.y < pos2.y);
                    bc_angle = angle(vec_b, vec_c);
                } else if(pos1.x < circleCenter.x && pos2.x < circleCenter.x) {
                    assert(pos1.y > pos2.y);
                    bc_angle = angle(vec_b, vec_c);
                } else if (pos1.x >= circleCenter.x && pos2.x < circleCenter.x) {
                    bc_angle = angle(vec_b, vec_down) + angle(vec_c, vec_down);
                } else {
                    bc_angle = angle(vec_b, vec_up) + angle(vec_c, vec_up);
                }
            } else {
                if (pos1.x >= circleCenter.x && pos2.x >= circleCenter.x) {
                    assert(pos1.y > pos2.y);
                    bc_angle = angle(vec_b, vec_c);
                } else if (pos1.x < circleCenter.x && pos2.x < circleCenter.x) {
                    assert(pos1.y < pos2.y);
                    bc_angle = angle(vec_b, vec_c);
                } else if (pos1.x >= circleCenter.x && pos2.x < circleCenter.x) {
                    bc_angle = angle(vec_b, vec_up) + angle(vec_c, vec_up);
                } else {
                    bc_angle = angle(vec_b, vec_down) + angle(vec_c, vec_down);
                }
            }

            totalTime += bc_angle*sqrt(radius);
        }
    }

    return totalTime;
}

const CityOrder findBestTrajectory() {
    //TODO teh genetic algo...
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

void printResult(CityOrder& order) {
    ofstream out;
    out.open("result.txt");
    out << order.size() << endl;

    for (auto a: order) {
        out << a.x << " " << a.y << endl;
    }

    out.close();
}

int main(int argc, char**argv) {

    //readInput(argv[1]);

    Vector A, B, C, D;

    cities.resize(4);
    A = Vector(100,140);
    B = Vector(185,145);
    C = Vector(230,131);
    D = Vector(190,40);

    cities = {A,B,C,D};

    cout << "A,B,C,D = " << evaluateCityOrder(cities) << endl;

    cities = {D,C,B,A};

    cout << "D,C,B,A = " << evaluateCityOrder(cities) << endl;

    cities = {A,B,D,C};

    cout << "A,B,D,C = " << evaluateCityOrder(cities) << endl;

    cities = {B,C,D,A};

    cout << "B,C,D,A = " << evaluateCityOrder(cities) << endl;

    cities = {A,D,C,B};

    cout << "A,D,C,B = " << evaluateCityOrder(cities) << endl;


}
