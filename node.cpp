//
// Created by vlad on 22.04.19.
//
#include "node.h"
#include <algorithm>

Node::Node(int numb):mu(0),n(numb),z(1e+06),t(n,0),m(0){}


int Node::getMu() const {
    return mu;
}


void Node::setMu(int var) {
    mu = var;
}


double Node::getZ() const{
    return z;
}


int Node::getM() const {
    return m;
}


void Node::setM(int var) {
    m = var;
}


void Node::findMin() {
    if(!t.empty()) {
        double min = t[0];
        for (int i = 1; i < t.size(); ++i) {
            if (min > t[i]) min = t[i];
        }
        z = min;
    }
    else z = 1e+06;
}


void Node::subtractValue(double val) {
    auto begin = t.begin();
    for(; begin != t.end(); ++begin){
        *begin -= val;
    }
    t.erase(std::remove(t.begin(),t.end(),0),t.end());
    findMin();
}


void Node::addValue(double val) {
    t.push_back(val);
    findMin();
}