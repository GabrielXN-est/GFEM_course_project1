#include "Enrichment.h"
#include "mesh.h"

// Construtores
void Enrichment::assign_to_node(Node& node)
{node.enr.push_back(create_copy(node));}

Enrichment* Sukumar_enrichment_1D::create_copy(Node& node)
{return (new Sukumar_enrichment_1D(id, xGamma));}

Enrichment* polinomial_enrichment_1D::create_copy(Node& node)
{return (new polinomial_enrichment_1D(id, grau, shifted, scaled, &node));}

Enrichment* Moes_enrichment_1D::create_copy(Node& node)
{return (new Moes_enrichment_1D(id, xGamma, &node));}

Enrichment* Pair_enrichment_1D::create_copy(Node& node)
{return (new Pair_enrichment_1D(enr1->create_copy(node), enr2->create_copy(node)));}

polinomial_enrichment_1D::polinomial_enrichment_1D(int index, int g, bool sh, bool sc, Node* node) :
 Enrichment{index, g}, grau {g}, shifted {sh}, scaled {sc}, node_ptr {node}
{
    if (scaled)
    {
        for (Element* el: node_ptr->vicinal_elements)
            {node_ptr->biggest_vicinal_element_size = 
                max(node_ptr->biggest_vicinal_element_size, el->Nod_list[el->Nod_list.size()-1]->x - el->Nod_list[0]->x);}
    }
}

double Sukumar_derivate(double x, double xGamma)
{
    if (x < xGamma)
        return -1.;
    else if (x >= xGamma)
        return 1.;
}

double Sukumar_enrichment_1D::D(double x) {return Sukumar_derivate(x, xGamma);}

double Moes_enrichment_1D::operator()(double x)
{
    double result {0};
    std::vector<double> dist {};

    for (Element* el : node_ptr->vicinal_elements)
    {
        if ((el->Nod_list[0]->x < x && el->Nod_list[0]->x + el->el_size > x) || el->Nod_list[0]->x == x || el->Nod_list[0]->x + el->el_size == x)
        {
            shape_functions* sf {el->get_shape_func()};
            sf->operator()(el->mapping(x, el->Nod_list[0]->x, el->el_size));
            sf->mont_vector();
            dist.reserve(sf->size());
            for (Node* node: el->Nod_list)
                {dist.push_back(node->x-xGamma);}
            for (std::size_t i {0}; i < sf->size(); i++)
                {result += sf->vec[i]*dist[i];}
            result = -std::abs(result);
            for (std::size_t i {0}; i < sf->size(); i++)
                {result +=sf->vec[i]*std::abs(dist[i]);}
            return result;
        }
    }
}

double Moes_enrichment_1D::D(double x)
{
    std::vector<double> dist {};
    double result {0};
    double temp {0}, dxidx {0};

    for (Element* el : node_ptr->vicinal_elements)
    {
        if ((el->Nod_list[0]->x < x && el->Nod_list[0]->x + el->el_size > x) || el->Nod_list[0]->x == x || el->Nod_list[0]->x + el->el_size == x)
        {
            shape_functions* sf {el->get_shape_func()};
            shape_functions* Dsf {el->get_D_shape_func()};

            sf->operator()(el->mapping(x, el->Nod_list[0]->x, el->el_size));
            sf->mont_vector();
            Dsf->operator()(el->mapping(x, el->Nod_list[0]->x, el->el_size));
            Dsf->mont_vector();

            dxidx = 2/el->el_size;

            dist.reserve(sf->size());
            for (Node* node: el->Nod_list)
                {dist.push_back(node->x-xGamma);}

            for (std::size_t i {0}; i < sf->size(); i++)
                {temp += sf->vec[i]*dist[i];}

            if (temp >= 0)
            {
                for (std::size_t i {0}; i < sf->size(); i++)
                    {   result += std::abs(dist[i])*Dsf->vec[i] *dxidx;
                        result -= dist[i]*Dsf->vec[i]*dxidx;}
            }
            else
            {
                for (std::size_t i {0}; i < sf->size(); i++)
                    {   result += std::abs(dist[i])*Dsf->vec[i] * dxidx;
                        result += dist[i]*Dsf->vec[i]*dxidx;}
            }
            return result;
        }
    }
}


