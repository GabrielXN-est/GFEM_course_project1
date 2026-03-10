#ifndef BODY_FUNC_H
#define BODY_FUNC_H

class body_functions
{
    public:
    int grau {};
    body_functions(int d) : grau{d} {}
    virtual double operator()(double x) =0;
};

class body_function0 : public body_functions
{
    public:
    body_function0() : body_functions(0) {}
    double operator()(double x){return 0.0;}
};

class body_function1 : public body_functions
{
    public:
    body_function1() : body_functions(1) {}
    double operator()(double x){return x;}
};

class body_function3 : public body_functions
{
    public:
    body_function3() : body_functions(3) {}
    double operator()(double x)
    {return x*x*x-6.*x*x-x+12;}
};
class body_function10 : public body_functions
{
    public:

    double alpha {50.};
    double xb{0.5};

    body_function10(double a, double b) : alpha{a}, xb{b}, body_functions(10){}

    double operator()(double x)
    {
        return 2*alpha/(1+alpha*alpha*(x-xb)*(x-xb)) +
        2*(1-x)*alpha*alpha*alpha*(x-xb)/
        ((1+alpha*alpha*(x-xb)*(x-xb))*(1+alpha*alpha*(x-xb)*(x-xb)));
    }
};
#endif