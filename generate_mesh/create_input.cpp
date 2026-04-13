#include <fstream>
#include <string>
#include <vector>
#include <cmath>

// to apply a bc to more than one dof per node, repete the node in the _bcs and _bcs_pos
// Exlim is the maximum x for which the respective E is valid
void generate_input(std::string filename, int nel, int porder, std::string eltype,
    double L, const std::vector<double>& E, const std::vector<double>& Exlim, double A, double C,
    const std::vector<double>& d_bcs, const std::vector<int>& d_bcs_pos, const std::vector<int>& d_bcs_dofs, // dirichilet boundary conditions
    const std::vector<double>& f_bcs, const std::vector<int>& f_bcs_pos, const std::vector<int>& f_bcs_dofs, // Neumann Boundary conditions
    int bf_func_id, double alpha=0.0, double xb=0.0, double xi = 0.0, 
    double xgamma = 0.0, int porder_Enrichment = 1, std::vector<double> geom_enr = {})
{
    // flags
    bool enriched_interface {eltype == "pGFEMBar_WD_S" || eltype == "pGFEMBar_WD_M" || eltype == "pGFEMBar_2Proj"|| eltype == "pSGFEMBar_2Proj"};
    bool GFEM {enriched_interface || eltype == "pGFEMBar" || eltype == "pGFEMBar_sc"};
    
    // open file
    std::ofstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file");
    }

    // print header
    file << filename << "\n";

    //get n de nodes
    int nnodes {};
    int n_per_el {};

    if (eltype == "lBar" || GFEM)
    {
        nnodes = nel*porder + 1;
        n_per_el = porder+1;
    }
    else if (eltype == "pBar")
    {
        nnodes = nel*2 + 1;
        n_per_el = 3;
    }
    else
        throw std::invalid_argument("Unexpected element type (" + eltype + ")");

    // nodes description
    double xi_el {0}, xf_el {0};

    if (enriched_interface && geom_enr.size() != 0)
    {
        file << "nodes - nnodes ndim; nodeID x-coord enrID\n" << nnodes << " 1\n";
        for (int i {0}; i < nnodes; i++)
        {
            xi_el = i * L / (nnodes-1) + xi;
            if (xi_el >= geom_enr[0] && xi_el <= geom_enr[1])
                file << i+1 << " " << xi_el << " 1\n";
            else
                file << i+1 << " " << xi_el << "\n";
        }
    }
    else
    {
        file << "nodes - nnodes ndim; nodeID x-coord\n" << nnodes << " 1\n";
        for (int i {0}; i < nnodes; i++)
        {file << i+1 << " " << i * L / (nnodes-1) + xi << "\n";}
    }

    // elements description
    if (enriched_interface && geom_enr.size() == 0)
        file << "nelem; elemID Type propID x-Gamma nodes\n";
    else
        file << "nelem; elemID Type propID nodes\n";
    file << nel << "\n";

    for (int i {0}; i < nel; i++)
    {
        xi_el = i*(n_per_el-1)*L / (nnodes-1) + xi;
        xf_el = (i+1)*(n_per_el-1)*L / (nnodes-1) + xi;
        // verifica se a interface não pertence ao elemento
        // Enriquecimentos no projeto 2 validos para singularidades nos nós
         // Só vai colocar elementos enriquecidos na inteface sem atrapalhar as condições de contorno
         // if it have geometrical enrichment, non polinomial enrichments defined in the nodes 
        if ((enriched_interface) && (((geom_enr.size() == 0) && (((eltype == "pGFEMBar_2Proj" || eltype == "pSGFEMBar_2Proj") && !(xi_el <= xgamma && xf_el >= xgamma)) ||
                    ((eltype != "pGFEMBar_2Proj" && eltype != "pSGFEMBar_2Proj") && !(xi_el < xgamma && xf_el > xgamma)) ||
                    (eltype == "pGFEMBar_WD_S" && nel == 1))) || (geom_enr.size() != 0)))
            {file << i+1 << " " << "pGFEMBar" << porder+1;}
        else
            {file << i+1 << " " << eltype << porder+1;}
        if (GFEM)
            {file << "_" << porder_Enrichment;}
        //propriedades
        file << " " << 1 << " ";
        //xGamma
        if (enriched_interface && geom_enr.size() == 0)
            {file << xgamma << " ";}
        //Nodes
        for (int j {0}; j < n_per_el; j++)
            {file << i*(n_per_el-1) + j + 1 << " ";}
        file << "\n";
    }

    // nodal defined enrichments description
    if (enriched_interface && geom_enr.size() != 0)
    {
        if (eltype == "pGFEMBar_WD_S" || eltype == "pGFEMBar_WD_M")
            file << "nenrichments; enrID type xGamma\n";
        else
            file << "nenrichments; enrID type\n";

        if (eltype == "pGFEMBar_WD_S" || eltype == "pGFEMBar_WD_M" || eltype == "pGFEMBar_2Proj"|| eltype == "pSGFEMBar_2Proj")
            file << 1 << "\n" << 1;

        if (eltype == "pGFEMBar_WD_S")
            file << " ESuk " << xgamma << "\n";
        else if (eltype == "pGFEMBar_WD_M")
            file << " EMoes " << xgamma << "\n";
        else if (eltype == "pGFEMBar_2Proj")
            file << " EProject2\n";
        else if (eltype == "pSGFEMBar_2Proj")
            file << " E_SGFEM_Project2\n";
    }

    // properties description
    file << "properties - nprop; propID type ";
    for (int i {0}; i < E.size(); i++)
        {file << "E ";}
    for (int i {0}; i < Exlim.size(); i++)
        {file << "Exlim ";}
    if (bf_func_id == 10)
        file << "A C bf_fun alpha xb\n";
    else
        file << "A C bf_fun\n";

    file << 1 << "\n";

    file << 1 << " Mat" << "Bar" << " ";
    for (double Ei: E)
        {file << Ei << " ";}
    for (double Exlimi: Exlim)
        {file << Exlimi << " ";}

    file << A << " " << C << " " << bf_func_id << " ";
    if (bf_func_id == 10)
        {file << alpha << " " << xb;}
    file << "\n";

    // constraints description
    file << "constraints - nconstr;constrID nodeID dof value\n";
    file << d_bcs.size() << "\n";
    for (int i {0}; i < d_bcs.size();)
    {        
        file << i+1 << " ";
        if (d_bcs_pos[i] ==0)
                file << 1 << " ";
        else
            {file << nnodes << " ";}
        file << d_bcs_dofs[i] << " " <<d_bcs[i++] << "\n";
    }

    // loads description
    file << "loads - nload; loadID nodeID dof value\n";
    file << f_bcs.size() << "\n";
    for (int i {0}; i < f_bcs.size(); i++)
    {        
        file << i+1 << " ";
        if (f_bcs_pos[i] ==0)
            {file << 1 << " ";}
        else
            {file << nnodes << " ";}
        file << f_bcs_dofs[i] << " "<< f_bcs[i] << "\n";
    }
}