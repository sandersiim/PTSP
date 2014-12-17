#include <array>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <set>
#include <vector>
#include <random>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>


#include "vector.h"
#define RANDOM_MUTATE
using namespace std;

mt19937 ranEng;

//struct City {
//    City () { }
//    City (Vector & _coords, size_t idx) : coords(_coords), index(idx) { }
//    Vector coords;
//    size_t index;
//};
//ostream& operator<<(ostream& os, const City* c)
//{
//    // write obj to stream
//    os << c->coords;
//    return os;
//}
using City = size_t;
vector<Vector> cities;
using CityOrder = vector<City>;
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
   //s if( b.x == 0 && b.y == 0) return 0;
    //return acos(dotProduct(a, b)/(magnitude(a)*magnitude(b)));
    return dotProduct(a, b)/(magnitude(a)*magnitude(b));
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
    //return abs((a.y - b.y) * (a.x - c.x) - (a.y - c.y) * (a.x - b.x)) <= 5.0d;

    double sum_x = a.x+b.x+c.x;
    double sum_y = a.y+b.y+c.y;
    double mean_x = (sum_x) / 3;
    double mean_y = (sum_y) / 3;

    double varx = (a.x*a.x + b.x*b.x + c.x*c.x) - sum_x * mean_x;
    double cov = (a.x*a.y + b.x*b.y + c.x*c.y) - sum_x * mean_y;

    // check for zero varx
    //y = kx + c
    double k = cov / varx;
    double d = mean_y - k * mean_x;
    double dista, distb, distc;
    dista = abs(k*a.x - a.y + d)/sqrt(k*k + 1);
    distb = abs(k*b.x - b.y + d)/sqrt(k*k + 1);
    distc = abs(k*c.x - c.y + d)/sqrt(k*k + 1);

    return dista < 4.5d && distb < 4.5d && distc < 4.5d;

}

template <typename T>
string vectorToString(vector<T>& vec) {
    ostringstream oss;

    if (!vec.empty())
    {
        // Convert all but the last element to avoid a trailing ","
        copy(vec.begin(), vec.end()-1,
            ostream_iterator<T>(oss, ", "));

        // Now add the last element with no delimiter
        oss << vec.back();
    }

    return oss.str();
}

struct Fitness {
    Fitness() { }
    Fitness ( double f) : fitness(f) { }
    Fitness (vector<double> & _times, double f) : fitness(f), times(_times) { }
    double fitness;
    vector<double> times;
    bool operator< (const Fitness & r) const {
        return fitness < r.fitness;
    }
};
double penaltyForTurning = 4.0d;

Fitness evaluateCityOrder2(CityOrder& order) {
    double totalTime = 0;
    Vector currentPos = startPos;
    Vector pos1, pos2, circleCenter;
    Vector vec_a,vec_b,vec_c,vec_up,vec_down;
    double radius, trajectory_angle;

    Vector newSpeed;
    Vector previousSpeed = Vector(0, 0);
    bool lastThreeOnLine = false;
    Vector lastCircleCenter;
    double lastRadius;
    double bc_angle, last_bc_angle;
    bool lastWasClockWise;
    vector<double> timeValues (order.size(), 0);
    double currentTime;
    for (size_t i = 0; i < order.size()-1; i++) {
        pos1 = cities[order[i]];
        pos2 = cities[order[i+1]];

        if (pointsAreOnLine(currentPos, pos1, pos2)) {
            //cout << i << endl;
            if (i == 0) {
                lastThreeOnLine = true;
                double dist = distance(currentPos, pos1);
                currentTime = sqrt(2*dist);
                totalTime += currentTime;
                timeValues[i] = currentTime;
                previousSpeed = Vector(pos1.x-currentPos.x, pos1.y-currentPos.y);
                previousSpeed = previousSpeed*(sqrt(2.0d)/sqrt(dist));
                if (i == order.size()-2) {
                    assert(false);
                }
                currentPos = pos1;
                continue;
            }
            if (i == order.size()-2) {
                if (lastThreeOnLine) {
                    double dist = distance(currentPos, pos1) + distance(pos1, pos2);
                    double v2 = (previousSpeed.x*previousSpeed.x + previousSpeed.y*previousSpeed.y);
                    double time = -sqrt(v2) + sqrt(v2 + 2*dist);
                    currentTime += time;
                    totalTime += currentTime;
                    timeValues[i] = currentTime/2.0;
                    timeValues[i+1] = currentTime/2.0; // Hack
                    return Fitness(timeValues, 1.0d/totalTime);
                }
                lastThreeOnLine = true;

                newSpeed = Vector(pos1.x-currentPos.x, pos1.y-currentPos.y);
                newSpeed = newSpeed*(magnitude(previousSpeed)/magnitude(newSpeed));
                Vector differenceVector = Vector(newSpeed.x - previousSpeed.x, newSpeed.y - previousSpeed.y);
                currentTime = penaltyForTurning*magnitude(differenceVector);

                double dist = distance(currentPos, pos1) + distance(pos1, pos2);
                double v2 = (newSpeed.x*newSpeed.x + newSpeed.y*newSpeed.y);
                double time = -sqrt(v2) + sqrt(v2 + 2*dist);
                currentTime += time;
                totalTime += currentTime;
                timeValues[i] = currentTime/2.0;
                timeValues[i+1] = currentTime/2.0; // Hack
                return Fitness(timeValues, 1.0d/totalTime);
            }
            if (i > 0 && i < order.size()-2) {
                if (lastThreeOnLine) {
                    newSpeed = Vector(pos1.x-currentPos.x, pos1.y-currentPos.y);
                    newSpeed = newSpeed*(magnitude(previousSpeed)/magnitude(newSpeed));
                    Vector differenceVector = Vector(newSpeed.x - previousSpeed.x, newSpeed.y - previousSpeed.y);
                    currentTime = penaltyForTurning*magnitude(differenceVector);

                    double dist = distance(currentPos, pos1);
                    double v2 = (newSpeed.x*newSpeed.x + newSpeed.y*newSpeed.y);
                    double time = -sqrt(v2) + sqrt(v2 + 2*dist);
                    currentTime += time;
                    totalTime += currentTime;
                    timeValues[i] = currentTime;
                    previousSpeed = newSpeed;
                    previousSpeed = previousSpeed*(sqrt(v2 + 2*dist)/magnitude(previousSpeed));

                } else {
                    lastThreeOnLine = true;
                    currentTime = last_bc_angle*sqrt(lastRadius);
                    totalTime += currentTime;
                    timeValues[i] = currentTime;
                    vec_b = Vector(pos1.x - circleCenter.x, pos1.y - circleCenter.y);

                    if (lastWasClockWise && currentPos.y >= lastCircleCenter.y) {
                        previousSpeed = Vector(1.0d, -(vec_b.x/vec_b.y));
                    } else if (lastWasClockWise && currentPos.y < lastCircleCenter.y) {
                        previousSpeed = Vector(-1.0d, (vec_b.x/vec_b.y));
                    } else if (!lastWasClockWise && currentPos.y >= lastCircleCenter.y) {
                        previousSpeed = Vector(-1.0d, (vec_b.x/vec_b.y));
                    } else {
                        previousSpeed = Vector(1.0d, -(vec_b.x/vec_b.y));
                    }

                    previousSpeed = previousSpeed*(sqrt(lastRadius)/magnitude(previousSpeed));
                }
                currentPos = pos1;
            }
            continue;
        }
        lastThreeOnLine = false;

//        cout << currentPos << "; " << pos1 <<"; " << pos2 << endl;

        if (distance(currentPos, pos1) < 12.0d) {

        }

        circleCenter = circleCenterPoint(currentPos, pos1, pos2);
        radius = distance(currentPos, circleCenter);
        lastCircleCenter = circleCenter;
        lastRadius = radius;
        assert(!(std::isnan(circleCenter.x) || std::isnan(circleCenter.y)));
//        cout << "center: " << circleCenter << endl;
//        cout << "radius: " << radius << endl;
        //speed = sqrt(radius);


        vec_a = Vector(currentPos.x - circleCenter.x, currentPos.y - circleCenter.y);
        vec_b = Vector(pos1.x - circleCenter.x, pos1.y - circleCenter.y);
        vec_c = Vector(pos2.x - circleCenter.x, pos2.y - circleCenter.y);
        vec_up = Vector(0, -radius);
        vec_down = Vector(0, radius);

//        cout << vec_a << endl;
//        cout << vec_b << endl;
//        cout << vec_c << endl;
//        cout << vec_up << endl;
//        cout << vec_down << endl;

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

        lastWasClockWise = clockWise;

        //arc_length = trajectory_angle*radius;
        currentTime = trajectory_angle*sqrt(radius);
//        cout << "Added trajectory time: " << trajectory_angle*sqrt(radius) << endl;
//        cout << "iter " << i << endl;
//        cout << "clockWise: " << clockWise << endl;
//        cout << "pos2IsInBetween: " << pos2IsInBetween << endl;

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
        currentTime += penaltyForTurning*magnitude(differenceVector);
        totalTime += currentTime;
        timeValues[i] = currentTime;
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
        last_bc_angle = bc_angle;

        if (i == order.size()-2) {
            currentTime = bc_angle*sqrt(radius);
            totalTime += currentTime;
            timeValues[i+1] = currentTime;
//            cout << "Added last arc time: " <<  bc_angle*sqrt(radius) << endl;

        } else {
            currentPos = pos1;
        }
    }

    return Fitness(timeValues, 1.0d/totalTime);
}
double k = 3.4d;
Fitness evaluateCityOrder(CityOrder& order) {
    double totalTime = 0;
    Vector currentPos = startPos;
    Vector pos1, pos2;

    Vector vec1, vec2;
    vector<double> timeValues (order.size(), 0);

    double currentTime, dist1, dist2, rotationAngle;
    for (size_t i = 0; i < order.size()-1; i++) {
        pos1 = cities[order[i]];
        pos2 = cities[order[i+1]];

        dist1 = distance(currentPos, pos1);
        dist2 = distance(pos1, pos2);

        vec1 = Vector(currentPos.x-pos1.x, currentPos.y-pos1.y);
        vec2 = Vector(pos2.x-pos1.x, pos2.y-pos1.y);

        rotationAngle = angle(vec1, vec2);

        currentTime = (sqrt(dist1) + sqrt(dist2))*(k + rotationAngle);
        //currentTime += 4*dist*(rotationAngle/M_PI);

        totalTime += currentTime;

        timeValues[i] = currentTime;

        currentPos = pos1;
    }

    return Fitness(timeValues, 1.0d/totalTime);
}

CityOrder nearestNeighbour(CityOrder& startingOrder) {
    set<size_t> visitedCities;
    set<size_t> remainingCities;

    Fitness bestFitness, fitness;

//    visitedCities.insert(city1);
//    visitedCities.insert(city2);
//    visitedCities.insert(city3);
    for (size_t i = 0; i < startingOrder.size(); i++) {
        visitedCities.insert(startingOrder[i]);
    }

    CityOrder bestOrder = startingOrder;
    CityOrder newBestOrder = bestOrder;
    size_t bestNextCity;
    while (visitedCities.size() < cities.size()) {
        bestFitness = Fitness(0);
        if (visitedCities.size() == cities.size()-3) {
            bestOrder.resize(bestOrder.size()+3);
            for (size_t i = 0; i < cities.size(); i++) {
                if (visitedCities.find(i) == visitedCities.end()) {
                    for (size_t j = 0; j < cities.size(); j++) {
                        if (i != j && visitedCities.find(j) == visitedCities.end()) {
                            for (size_t k = 0; k < cities.size(); k++) {
                                if (k != i && k != j && visitedCities.find(k) == visitedCities.end()) {
                                    bestOrder[bestOrder.size()-3] = i;
                                    bestOrder[bestOrder.size()-2] = j;
                                    bestOrder[bestOrder.size()-1] = k;
                                    fitness = evaluateCityOrder(bestOrder);
                                    if (bestFitness < fitness ) {
                                        bestFitness = fitness;
                                        newBestOrder = bestOrder;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            return newBestOrder;
        } else {
            bestOrder.resize(bestOrder.size()+4);
            for (size_t i = 0; i < cities.size(); i++) {
                if (visitedCities.find(i) == visitedCities.end()) {
                    for (size_t j = 0; j < cities.size(); j++) {
                        if (i != j && visitedCities.find(j) == visitedCities.end()) {
                            for (size_t k = 0; k < cities.size(); k++) {
                                if (k != i && k != j && visitedCities.find(k) == visitedCities.end()) {
                                    for (size_t l = 0; l < cities.size(); l++) {
                                        if (l != i && l != j && l != k && visitedCities.find(l) == visitedCities.end()) {
                                            bestOrder[bestOrder.size()-4] = i;
                                            bestOrder[bestOrder.size()-3] = j;
                                            bestOrder[bestOrder.size()-2] = k;
                                            bestOrder[bestOrder.size()-1] = l;
                                            fitness = evaluateCityOrder(bestOrder);
                                            if (bestFitness < fitness ) {
                                                bestFitness = fitness;
                                                bestNextCity = i;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

            }
            newBestOrder.resize(newBestOrder.size() + 1);
            newBestOrder[newBestOrder.size()-1] = bestNextCity;
            bestOrder = newBestOrder;
            visitedCities.insert(bestNextCity);
        }
    }
    return bestOrder;
}

//const CityOrder findBestOrderNearestNeighbour() {
//    //TODO teh genetic algo...
//    CityOrder bestOrder;
//    CityOrder candidate;
//    array<size_t, 3> firstThreeCities;
//    candidate.resize(3);

//    double bestFitness, fitness;
//    bestFitness = 0;

//    for (size_t i = 0; i < cities.size(); i++) {
//        for (size_t j = 0; j < cities.size(); j++) {
//            if (j == i) continue;
//            for (size_t k = 0; k < cities.size(); k++) {
//                if (k == i || k == j) continue;
//                candidate[0] = cities[i];
//                candidate[1] = cities[j];
//                candidate[2] = cities[k];
//                fitness = evaluateCityOrder(candidate);
//                if (fitness > bestFitness) {
//                    bestFitness = fitness;
//                    bestOrder = candidate;
//                    firstThreeCities = {i, j, k};
//                }
//            }
//        }
//    }

//    //return nearestNeighbour(firstThreeCities[0], firstThreeCities[1], firstThreeCities[2]);
//    return bestOrder;
//}


struct FitCityOrder {
    FitCityOrder () { }
    FitCityOrder (CityOrder & _order) : order(_order) { }
    FitCityOrder (CityOrder & _order, Fitness& f) : fitness(f), order(_order) { }
    Fitness fitness;
    CityOrder order;
    bool operator< (const FitCityOrder & r) const {
        return fitness < r.fitness;
    }
};


using Chromosome = FitCityOrder;
using Population = vector<Chromosome>;

Population population;

#ifdef RANDOM_MUTATE
size_t populationSize = 2000;
size_t numberOfParents = 1000;
size_t numberOfChildren = 1000;
#else
size_t populationSize = 200;
size_t numberOfParents = 200;
size_t numberOfChildren = 600;
#endif


void evaluateIndividuals(vector<double> & fitness) {
    for (size_t i = 0; i < populationSize; ++i) {
        fitness[i] = population[i].fitness.fitness;
    }
}

void naturalSelection(vector<double> & fitness, vector<FitCityOrder> & parents)
{
    discrete_distribution<int> selection(fitness.begin(), fitness.end());
    for (size_t i = 0; i < numberOfParents; ++i) {
        int s = selection(ranEng);
        parents[i] = population[s];
    }
}

uniform_int_distribution<> randomCity;


void mutateChromo3(Chromosome& chromo) {

    size_t worstIndex = distance(chromo.fitness.times.begin(),
                                 max_element(chromo.fitness.times.begin(), chromo.fitness.times.end()));
    size_t i = randomCity(ranEng);
    while(worstIndex == i) {
        i = randomCity(ranEng);
    }

    size_t movedElement = chromo.order[worstIndex];
    swap(chromo.order[i], chromo.order[worstIndex]);
    chromo.fitness = evaluateCityOrder(chromo.order);
}


void mutateChromo(Chromosome& chromo) {
    size_t i = randomCity(ranEng);
    size_t j = i;
    while(j == i) {
        j = randomCity(ranEng);
    }

    swap(chromo.order[i], chromo.order[j]);
    chromo.fitness = evaluateCityOrder(chromo.order);
}

uniform_int_distribution<> randomCity1;

void mutateChromo2(Chromosome& chromo) {
//    size_t i = randomCity1(ranEng);

//    swap(chromo.order[i], chromo.order[i+1]);
//    shuffle(chromo.order.begin(), chromo.order.end(), ranEng);

    mutateChromo(chromo);
//#ifdef RANDOM_MUTATE
    mutateChromo(chromo);
//#endif

    chromo.fitness = evaluateCityOrder(chromo.order);

}

uniform_real_distribution<> randomProbability(0, 1);
#ifdef RANDOM_MUTATE
double mutationProbability = 0.01;
#else
double mutationProbability = 0.12;
#endif
void mutatePopulation() {
    double pr;
    for (auto &c : population) {
        pr = randomProbability(ranEng);
        if (pr < mutationProbability) {
            mutateChromo2(c);
        }
    }

}

uniform_int_distribution<> distParents(0, numberOfParents - 1);
uniform_int_distribution<> distWindow(0, 8);
size_t windowSize = 22;
void reproduce(vector<FitCityOrder> & parents) {
    vector<FitCityOrder> selectionPool (numberOfParents+numberOfChildren);
    size_t parent1, parent2, windowStart;
    double fitness1, fitness2;
    for (int i = 0; i < selectionPool.size(); i ++) {
        if (i < numberOfChildren) {
            Chromosome child = parents[i%numberOfParents];
#ifdef RANDOM_MUTATE
            mutateChromo3(child);
            mutateChromo3(child);
#else
            mutateChromo(child);
#endif
            selectionPool[i] = child;

            /*parent1 = distParents(ranEng);
            parent2 = parent1;
            while (parent2 == parent1) {
                parent2 = distParents(ranEng);
            }
            windowStart = distWindow(ranEng);
            fitness1 = fitness2 = 0;
            vector<size_t> availableIndices;
            for (size_t i = windowStart; i < windowStart+windowSize; i++) {
                fitness1 += parents[parent1].fitness.times[i];
                fitness2 += parents[parent2].fitness.times[i];
            }*/

            /*if (fitness1 > fitness2) {
                for (size_t i = 0; i < cities.size(); i++) {
                    if (i < windowStart || i >= windowStart+windowSize) {
                        availableIndices.push_back(parents[parent1].order[i]->index);
                    }
                }
            } else {
                for (size_t i = 0; i < cities.size(); i++) {
                    if (i < windowStart || i >= windowStart+windowSize) {
                        availableIndices.push_back(parents[parent2].order[i]->index);
                    }
                }
            }

            shuffle(availableIndices.begin(), availableIndices.end(), ranEng);*/
            /*CityOrder childCityOrder (cities.size());
            for (size_t i = 0; i < cities.size(); i++) {
                if (i < windowStart) {
                    childCityOrder[i] = cities[availableIndices[i]];
                } else if (i >= windowStart && i < windowStart+windowSize) {
                    if (fitness1 > fitness2) {
                        childCityOrder[i] = parents[parent1].order[i];
                    } else {
                        childCityOrder[i] = parents[parent2].order[i];
                    }
                } else {
                    childCityOrder[i] = cities[availableIndices[i-windowSize]];
                }
            }*/
            /*CityOrder childCityOrder (windowSize);
            if (fitness1 > fitness2) {
                for (size_t i = 0; i < windowSize; i++) {
                    childCityOrder[i] = parents[parent1].order[windowStart+i];
                }
            } else {
                for (size_t i = 0; i < windowSize; i++) {
                    childCityOrder[i] = parents[parent2].order[windowStart+i];
                }
            }

            childCityOrder = nearestNeighbour(childCityOrder);

            Fitness childFitness = evaluateCityOrder(childCityOrder);
            selectionPool[i] = Chromosome(childCityOrder, childFitness);*/

            /*fitness1 = fitness2 = 0;
            for (size_t i = 0; i < windowSize; i++) {
                fitness1 += parents[parent1].fitness.fitness;
                fitness2 += parents[parent2].fitness.fitness;
            }
            CityOrder childCityOrder (windowSize);
            if (fitness1 > fitness2) {
                for (size_t i = 0; i < windowSize; i++) {
                    childCityOrder[i] = parents[parent1].order[i];
                }
            } else {
                for (size_t i = 0; i < windowSize; i++) {
                    childCityOrder[i] = parents[parent2].order[i];
                }
            }
            childCityOrder = nearestNeighbour(childCityOrder);
            Fitness childFitness = evaluateCityOrder(childCityOrder);
            selectionPool[i] = Chromosome(childCityOrder, childFitness);*/
        } else {
            selectionPool[i] = population[i-numberOfChildren];
        }
    }
    sort(selectionPool.rbegin(), selectionPool.rend());
    for (size_t i = 0; i < populationSize; i++) {
        population[i] = selectionPool[i];
    }
}



CityOrder canonicalOrder;

void initChromo(Chromosome & c) {
    c.order = canonicalOrder;
    shuffle(c.order.begin(), c.order.end(), ranEng);
    //c.order.resize(20);

    //c.order = nearestNeighbour(c.order);
    c.fitness = evaluateCityOrder(c.order);
    /*while (c.fitness < 0.0014d) {
        shuffle(c.order.begin(), c.order.end(), ranEng);
        c.fitness = evaluateCityOrder(c.order);
    }*/
}

void initPopulation() {
    random_device d;
    ranEng.seed(d());
    population.resize(populationSize);
    for (auto & p : population) {
        initChromo(p);
    }
    CityOrder order;
    //Seed population with nearestNeighbour result
#ifndef RANDOM_MUTATE
    population[0].order = nearestNeighbour(order);
    population[0].fitness = evaluateCityOrder(population[0].order);
#endif
}


void readInput(const char *file) {
    ifstream f;
    f.open(file);
    int n;
    f >> n;

    cities.resize(n);
    canonicalOrder.resize(n);
    for (size_t i = 0; i < n; ++i) {
        f >> cities[i].x >> cities[i].y;
        canonicalOrder[i] = i;
    }

    f.close();
    randomCity = uniform_int_distribution<> (0, n-1);
    randomCity1 = uniform_int_distribution<> (0, n-2);

}

CityOrder readOrder(const char *file) {
    ifstream f;
    f.open(file);
    int n;
    f >> n;

    vector<Vector> orderCities (n);

    for (size_t i = 0; i < n; ++i) {
        f >> orderCities[i].x >> orderCities[i].y;
    }

    CityOrder order;
    for (auto c : orderCities) {
        for (size_t i = 0; i < cities.size(); i++) {
            if (c.x == cities[i].x && c.y == cities[i].y) {
                order.push_back(i);
                break;
            }
        }
    }

    f.close();
    return order;

}

void printResult(CityOrder& order) {
    ofstream out;
    out.open("best_order.txt");
    out << order.size() << endl;

    for (auto a: order) {
        out << cities[a].x << " " << cities[a].y << endl;
    }

    out.close();
}

size_t numberOfGenerations = 10000;

int main(int argc, char**argv) {

    readInput(argv[1]);

    if (argc > 2) {
        CityOrder order = readOrder(argv[2]);
        cout << evaluateCityOrder(order).fitness << endl;
        cout << 1.0d/evaluateCityOrder(order).fitness << endl;
        return 0;
    }

    vector<double> fitness;
    vector<FitCityOrder> parents;
    fitness.resize(populationSize);
    parents.resize(numberOfParents);

    Fitness bestFitness(0);
    Chromosome bestChromosome;
    //while( true ) {
        initPopulation();
        evaluateIndividuals(fitness);

        auto bf = max_element(population.begin(), population.end());
        FitCityOrder currentGenerationBestChromosome = *bf;
        Fitness currentGenerationBestFitness = bestChromosome.fitness;
        if (bestFitness < currentGenerationBestFitness) {
            bestFitness = currentGenerationBestFitness;
            bestChromosome = currentGenerationBestChromosome;
            double avgFitness = accumulate(fitness.begin(), fitness.end(), 0.0)
                    / populationSize;
            cout << currentGenerationBestFitness.fitness << "\t" << avgFitness << endl;
            printResult(bestChromosome.order);
        }

        cout << "New Population" << endl;
        for (size_t i = 0; i < numberOfGenerations; ++i) {
            naturalSelection(fitness, parents);
            reproduce(parents);
            mutatePopulation();
            evaluateIndividuals(fitness);
            auto bf = max_element(population.begin(), population.end());
            currentGenerationBestChromosome = *bf;
            currentGenerationBestFitness = currentGenerationBestChromosome.fitness;
            double avgFitness = accumulate(fitness.begin(), fitness.end(), 0.0)
                    / populationSize;
            cout << currentGenerationBestFitness.fitness << "\t" << avgFitness << "\t" << i << endl;

            if (bestFitness < currentGenerationBestFitness) {
                bestFitness = currentGenerationBestFitness;
                bestChromosome = currentGenerationBestChromosome;

                printResult(bestChromosome.order);
            }

        }
    //}

    cout << "Best fitness result: " << bestChromosome.fitness.fitness << endl;

    printResult(bestChromosome.order);

    /*for (auto p : population) {
        cout << vectorToString(p) << endl;
    }*/

//    CityOrder bestOrder = {};
//    bestOrder = nearestNeighbour(bestOrder);

//    printResult(bestOrder);

//    cout << evaluateCityOrder(bestOrder).fitness << endl;
//    cout << 1.0d/evaluateCityOrder(bestOrder).fitness << endl;


//    Vector A, B, C, D;

//    cities.resize(4);
//    A = Vector(100,144);
//    B = Vector(185,145);
//    C = Vector(230,146);
//    D = Vector(190,40);
//    cities = {A,B,C,D};

//    CityOrder order;

//    order = {0,1,2,3};
//    Fitness f = evaluateCityOrder(order);
//    cout << "A,B,C,D = " << f.fitness << " ; " << vectorToString(f.times) << endl;

//    order = {3,2,1,0};
//    f = evaluateCityOrder(order);
//    cout << "D,C,B,A = " << f.fitness  << " ; " <<vectorToString(f.times) << endl;

//    order = {0,1,3,2};
//    f = evaluateCityOrder(order);
//    cout << "A,B,D,C = " << f.fitness << " ; " <<vectorToString(f.times)  << endl;

//    order = {1,2,3,0};
//    f = evaluateCityOrder(order);
//    cout << "B,C,D,A = " << f.fitness  << " ; " <<vectorToString(f.times)  << endl;

//    order = {0,3,2,1};
//    f = evaluateCityOrder(order);
//    cout << "A,D,C,B = " << f.fitness  <<" ; " <<vectorToString(f.times)  << endl;





    //Test 3 points on line
    /*A = Vector(100,140);
    B = Vector(185,141);
    C = Vector(230,141);
    D = Vector(240,141);

    cities = {A,B,D,C};

    cout << "A,B,D,C = " << evaluateCityOrder(cities) << endl;*/




}
