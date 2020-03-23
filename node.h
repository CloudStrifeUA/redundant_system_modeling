//
// Created by vlad on 22.04.19.
//

#ifndef CMO_NODE_H
#define CMO_NODE_H

#include <vector>

class Node{

public:

explicit Node(int n);
int getMu() const ;
void setMu(int mu);
double getZ() const;
int getM() const ;
void setM(int m);
void subtractValue(double val);  // Віднімає задане значення від усіх елементів вектора
void addValue(double val);       // Додає нове значення у вектор

private:

int mu;                          // Кількість зайнятих каналів
int n;                           // Кількість каналів
std::vector<double> t;           // Вектор який зберігає час обслуговування в кожному каналі
double z;                        // Мінімальний час до завершення обслуговування
int m;                           // Число заявок в черзі
void findMin();                  // Знаходить мінімальний час до завершення обслуговування, якщо вектор порожній то
                                 // встановлює мінімальний час 10^6
};
#endif //CMO_NODE_H