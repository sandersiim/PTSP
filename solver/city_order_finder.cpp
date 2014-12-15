#include <array>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <set>
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
Vector circleCenterPoint(Vector p, Vector q, Vector r) {
//    cout << p << ", " << q << ", " << r << endl;
    if ((int)q.y == (int)p.y) {
        Vector tmp = r;
        r = q;
        q = tmp;
    } else if ((int)q.y == (int)r.y) {
        Vector tmp = q;
        q = p;
        p = tmp;
    }
    double p1, y11, yd1, pd1, p2, y2, yd2, pd2, ox, oy;
    p1 = (q.x + p.x) / 2.0d;
    y11 = (q.y + p.y) / 2.0d;
    yd1 = q.x - p.x;
    pd1 = -(q.y - p.y);
    //***********************
    p2 = (r.x + q.x) / 2.0d;
    y2 = (r.y + q.y) / 2.0d;
    yd2 = r.x - q.x;
    pd2 = -(r.y - q.y);
    //****************************
    ox = (y11 * pd1 * pd2 + p2 * pd1 * yd2 - p1 * yd1 * pd2 - y2 * pd1 * pd2)/ (pd1 * yd2 - yd1 * pd2);
    oy = (ox - p1) * yd1 / pd1 + y11;
    return Vector(ox, oy);
}

bool pointsAreOnLine(Vector& a, Vector& b, Vector& c) {
    return abs((a.y - b.y) * (a.x - c.x) - (a.y - c.y) * (a.x - b.x)) <= 1e-9;
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

        //TODO: Handle situation where points are on the same line

        if (pointsAreOnLine(currentPos, pos1, pos2)) {
            if (i == 0) {
                double dist = distance(currentPos, pos1);
                totalTime += sqrt(2*dist);
                previousSpeed = Vector(pos1.x-currentPos.x, pos1.y-currentPos.y);
                previousSpeed = previousSpeed*(sqrt(2.0d)/sqrt(dist));
                if (i == order.size()-2) {
                    assert(false);
                }
                continue;
            }
            if (i == order.size()-2) {
                newSpeed = Vector(pos1.x-currentPos.x, pos1.y-currentPos.y);
                newSpeed = newSpeed*(magnitude(previousSpeed)/magnitude(newSpeed));
                Vector differenceVector = Vector(newSpeed.x - previousSpeed.x, newSpeed.y - previousSpeed.y);
                totalTime += magnitude(differenceVector);

                double dist = distance(currentPos, pos1) + distance(pos1, pos2);
                double v2 = (newSpeed.x*newSpeed.x + newSpeed.y*newSpeed.y);
                double time = -sqrt(v2) + sqrt(v2 + 2*dist);
                totalTime += time;
                return 1.0d/totalTime;
            }
            if (i > 0 && i < order.size()-2) {
                assert(false);
            }
            continue;
        }

        circleCenter = circleCenterPoint(currentPos, pos1, pos2);
        radius = distance(currentPos, circleCenter);
//        cout << "center: " << circleCenter << endl;
//        cout << "radius: " << radius << endl;
        //speed = sqrt(radius);


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
                    assert(pos1.y < pos2.y || pos2.y < currentPos.y);
                    if (pos1.y < pos2.y) {
                        bc_angle = angle(vec_b, vec_c);
                    } else {
                        bc_angle = M_PI + angle(vec_b, vec_down) + angle(vec_c, vec_up);
                    }
                } else if(pos1.x < circleCenter.x && pos2.x < circleCenter.x) {
                    assert(pos1.y > pos2.y || pos2.y > currentPos.y);
                    if (pos1.y > pos2.y) {
                        bc_angle = angle(vec_b, vec_c);
                    } else {
                        bc_angle = M_PI + angle(vec_b, vec_up) + angle(vec_c, vec_down);
                    }
                } else if (pos1.x >= circleCenter.x && pos2.x < circleCenter.x) {
                    bc_angle = angle(vec_b, vec_down) + angle(vec_c, vec_down);
                } else {
                    bc_angle = angle(vec_b, vec_up) + angle(vec_c, vec_up);
                }
            } else {
                if (pos1.x >= circleCenter.x && pos2.x >= circleCenter.x) {
                    assert(pos1.y > pos2.y || pos2.y > currentPos.y);
                    if (pos1.y > pos2.y) {
                        bc_angle = angle(vec_b, vec_c);
                    } else {
                        bc_angle = M_PI + angle(vec_b, vec_up) + angle(vec_c, vec_down);
                    }
                } else if (pos1.x < circleCenter.x && pos2.x < circleCenter.x) {
                    assert(pos1.y < pos2.y || pos2.y < currentPos.y);
                    if (pos1.y < pos2.y) {
                        bc_angle = angle(vec_b, vec_c);
                    } else {
                        bc_angle = M_PI + angle(vec_b, vec_down) + angle(vec_c,vec_up);
                    }
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

    return 1.0d/totalTime;
}

//double alternativeEvaluate(CityOrder& order) {
//    double totalTime = 0;
//    Vector currentPos = startPos;
//    Vector pos1;

//    Vector newSpeed;
//    Vector previousSpeed = Vector(0, 0);
//    Vector vec_ab;

//    double trajectoryTime, trajectorySpeed, dist;
//    for (size_t i = 0; i < order.size()-1; i++) {
//        pos1 = order[i];

//        dist = distance(currentPos, pos1);

//        trajectorySpeed = sqrt(2.0d*dist)/2.0d;

//        trajectoryTime = dist/trajectorySpeed;
//        totalTime += trajectoryTime;

//        vec_ab = Vector(pos1.x - currentPos.x, pos1.y - currentPos.y);

//        newSpeed = vec_ab/sqrt(2.0d*dist);

//        Vector differenceVector = Vector(newSpeed.x-previousSpeed.x, newSpeed.y-previousSpeed.y);
//        totalTime += magnitude(differenceVector);
////        cout << "Added velocity change time: " << magnitude(differenceVector) << endl;

//        previousSpeed = newSpeed;

//        currentPos = pos1;
//    }

//    pos1 = order.back();

//    newSpeed = Vector(pos1.x - currentPos.x, pos1.y - currentPos.y);
//    newSpeed = newSpeed*(magnitude(previousSpeed)/magnitude(newSpeed));

//    Vector differenceVector = Vector(newSpeed.x-previousSpeed.x, newSpeed.y-previousSpeed.y);
//    totalTime += magnitude(differenceVector);

//    dist = distance(currentPos, pos1);

//    double v2 = newSpeed.x*newSpeed.x + newSpeed.y*newSpeed.y;

//    totalTime += -sqrt(v2) + sqrt(v2 + 2*dist);

//    return totalTime;
//}

const CityOrder findBestCityOrder() {
    //TODO teh genetic algo...
    CityOrder bestOrder;
    CityOrder candidate;
    array<size_t, 3> firstThreeCities;
    candidate.resize(3);
    set<size_t> visitedCities;

    double bestFitness, fitness;
    bestFitness = 0;

    for (size_t i = 0; i < cities.size(); i++) {
        for (size_t j = 0; j < cities.size(); j++) {
            if (j == i) continue;
            for (size_t k = 0; k < cities.size(); k++) {
                if (k == i || k == j) continue;
                candidate[0] = cities[i];
                candidate[1] = cities[j];
                candidate[2] = cities[k];
                fitness = evaluateCityOrder(candidate);
                if (fitness > bestFitness) {
                    bestFitness = fitness;
                    bestOrder = candidate;
                    firstThreeCities = {i, j, k};
                }
            }
        }
    }

    visitedCities.insert(firstThreeCities[0]);
    visitedCities.insert(firstThreeCities[1]);
    visitedCities.insert(firstThreeCities[2]);


    CityOrder newBestOrder = bestOrder;
    size_t bestNextCity;
    while (visitedCities.size() < cities.size()) {
        bestFitness = 0;
        bestOrder.resize(bestOrder.size()+1);
        for (size_t i = 0; i < cities.size(); i++) {
            if (visitedCities.find(i) == visitedCities.end()) {
                bestOrder[bestOrder.size()-1] = cities[i];
                fitness = evaluateCityOrder(bestOrder);
                if (fitness > bestFitness) {
                    newBestOrder = bestOrder;
                    bestNextCity = i;
                }
            }
        }
        bestOrder = newBestOrder;
        visitedCities.insert(bestNextCity);
    }
    return bestOrder;
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
    out.open("best_order.txt");
    out << order.size() << endl;

    for (auto a: order) {
        out << a.x << " " << a.y << endl;
    }

    out.close();
}

int main(int argc, char**argv) {

    readInput(argv[1]);

    CityOrder bestOrder = findBestCityOrder();

    printResult(bestOrder);

    cout << evaluateCityOrder(bestOrder) << endl;
    cout << 1.0d/evaluateCityOrder(bestOrder) << endl;


    /*Vector A, B, C, D;

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


    cout << "Alternative evaluation" << endl;

    cities = {A,B,C,D};

    cout << "A,B,C,D = " << alternativeEvaluate(cities) << endl;

    cities = {D,C,B,A};

    cout << "D,C,B,A = " << alternativeEvaluate(cities) << endl;

    cities = {A,B,D,C};

    cout << "A,B,D,C = " << alternativeEvaluate(cities) << endl;

    cities = {B,C,D,A};

    cout << "B,C,D,A = " << alternativeEvaluate(cities) << endl;

    cities = {A,D,C,B};

    cout << "A,D,C,B = " << alternativeEvaluate(cities) << endl;*/


    //Test 3 points on line
    /*A = Vector(100,140);
    B = Vector(185,141);
    C = Vector(230,141);
    D = Vector(240,141);

    cities = {A,B,D,C};

    cout << "A,B,D,C = " << evaluateCityOrder(cities) << endl;*/




}
