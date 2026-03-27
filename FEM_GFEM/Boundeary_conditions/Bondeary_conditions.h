#ifndef BONDEARY_CONDITIONS_H
#define BONDEARY_CONDITIONS_H

#include <vector>
#include "node.h"

class BC_displacement
{
    public:
    int id {};
    Node* no;
    int dof {};
    double value {};

    BC_displacement (){}

    void assign_node (Node* node){no = node;}

    void clear()
    {
        id = dof = 0;
        value = 0.;
    }
};

class BC_load
{
    public:
    int id {};
    Node* no;
    int dof {};
    double value {};

    BC_load (){}
    
    void assign_node (Node* node){no = node;}

    void clear()
    {
        id = dof = 0;
        value = 0.;
    }
};
#endif