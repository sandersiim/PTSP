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

// from http://forums.techarena.in/software-development/1415154.htm
Vector circleCenterPoint(Vector& p, Vector& q, Vector& r) {
//    cout << p << ", " << q << ", " << r << endl;
    double p1, y11, yd1, pd1, p2, y2, yd2, pd2, ox, oy;
    p1 = (q.x + p.x) / 2;
    y11 = (q.y + p.y) / 2;
    yd1 = q.x - p.x;
    pd1 = -(q.y - p.y);
    //***********************
    p2 = (r.x + q.x) / 2;
    y2 = (r.y + q.y) / 2;
    yd2 = r.x - q.x;
    pd2 = -(r.y - q.y);
    //****************************
    ox = (y11 * pd1 * pd2 + p2 * pd1 * yd2 - p1 * yd1 * pd2 - y2 * pd1 * pd2)/ (pd1 * yd2 - yd1 * pd2);
    oy = (ox - p1) * yd1 / pd1 + y11;
    return Vector(ox, oy);
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
//        cout << "center: " << circleCenter << endl;
//        cout << "radius: " << radius << endl;
        //speed = sqrt(radius);

        //TODO: Handle situation where points are on the same line

        vec_a = Vector(currentPos.x - circleCenter.x, currentPos.y - circleCenter.y);
        vec_b = Vector(pos1.x - circleCenter.x, pos1.y - circleCenter.y);
        vec_c = Vector(pos2.x - circleCenter.x, pos2.y - circleCenter.y);
        vec_up = Vector(0, -radius);
        vec_down = Vector(0, radius);

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
//        cout << "Added trajectory time: " << trajectory_angle*sqrt(radius) << endl;

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
//        cout << "Added velocity change time: " << magnitude(differenceVector) << endl;

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

        if (i == order.size()-2) {
            double bc_angle;
            if (clockWise) {
                if (pos1.x >= circleCenter.x && pos2.x >= circleCenter.x) {
                    assert(pos1.y < pos2.y);
                    bc_angle = angle(vec_b, vec_c);
                } else if(pos1.x < circleCenter.x && pos2.x < circleCenter.x) {
                    assert(pos1.y > pos2.y);
                    bc_angle = angle(vec_b, vec_c);
                } else if (pos1.x >= circleCenter.x && pos2.x < circleCenter.x) {
//                    cout << "Hello" << endl;
//                    cout << vec_b << ", " << vec_down << endl;
//                    cout << angle(vec_b, vec_down) << endl;
//                    cout << angle(vec_c, vec_down) << endl;
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
//            cout << "Added last arc time: " <<  bc_angle*sqrt(radius) << endl;

        } else {
            currentPos = pos1;
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


    //Test 3 points on line
    /*A = Vector(100,140);
    B = Vector(185,141);
    C = Vector(230,141);
    D = Vector(240,141);

    cities = {A,B,D,C};

    cout << "A,B,D,C = " << evaluateCityOrder(cities) << endl;*/




}
