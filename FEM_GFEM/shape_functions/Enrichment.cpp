#include "Enrichment.h"
#include "mesh.h"

void Sukumar_enrichment_1D::assign_to_node(Node& node)
{node.enr.push_back(new Sukumar_enrichment_1D(id, xGamma));}

void polinomial_enrichment_1D::assign_to_node(Node& node)
{node.enr.push_back(new polinomial_enrichment_1D(id, grau));}