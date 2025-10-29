#pragma once
#include<iostream>
#include<cmath>
#include "constants.h"

class vec
{
public:
    vec(double x, double y, double z) : x(x), y(y), z(z) {}
    vec() :x(0), y(0), z(0) {}
    vec(const vec& other);
    vec& operator = (const vec& other);
    double size() const;
    double x;
    double y;
    double z;
};

class M
{
public:
    M(int row, int column);

    ~M();

    M(const M& other);

    double* operator[] (int i) const;

    M& operator = (const M& other);

    M operator* (const M& other);

    void T();

    vec operator*(const vec& other);

    int Row() const;

    int Column() const;

protected:
    double** arr;
    int row;
    int column;
};

class MX : public M
{
public:
    MX(double fi);
};
class MY : public M
{
public:
    MY(double fi);
};
class MZ : public M
{
public:
    MZ(double fi);
};


void showM(const M& a);

void showV(const vec& a);
