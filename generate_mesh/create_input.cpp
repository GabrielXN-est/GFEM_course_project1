#include <fstream>
#include <string>
#include <vector>

// to apply a bc to more than one dof per node, repete the node in the _bcs and _bcs_pos
void generate_input(std::string filename, int nel, int porder, std::string eltype,
    double L, std::vector<double> E, std::vector<double> Exlim, double A, double C,
    std::vector<double> d_bcs, std::vector<int> d_bcs_pos, std::vector<int> d_bcs_dofs, // dirichilet boundary conditions
    std::vector<double> f_bcs, std::vector<int> f_bcs_pos, std::vector<int> f_bcs_dofs, // Neumann Boundary conditions
    int bf_func_id, double alpha=0.0, double xb=0.0, double xi = 0.0, 
    double xgamma = 0.0, int porder_Enrichment = 1)
{
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

    if (eltype == "lBar" || eltype == "pGFEMBar" || eltype == "pGFEMBar_WD_S" || eltype == "pGFEMBar_WD_M")
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
        throw std::invalid_argument("Unexpected element type");

    // nodes description
    file << "nodes - nnodes ndim; nodeID x-coord\n";
    file << nnodes << " 1\n";
    for (int i {0}; i < nnodes; i++)
    {
        file << i+1 << " " << i * L / (nnodes-1) + xi << "\n";
    }

    // elements description
    if (eltype == "pGFEMBar_WD_S" || eltype == "pGFEMBar_WD_M")
        file << "nelem; elemID Type propID x-Gamma nodes\n";
    else
        file << "nelem; elemID Type propID nodes\n";
    file << nel << "\n";
    for (int i {0}; i < nel; i++)
    {
        //ID, tipo e ordem polinomial
        if ((eltype == "pGFEMBar_WD_S" || eltype == "pGFEMBar_WD_M") && 
            (i*n_per_el*L / (nnodes-1) + xi < xgamma && (i+1)*n_per_el*L / (nnodes-1) + xi > xgamma ))
            {file << i+1 << " " << "pGFEMBar" << porder+1;}
        else
            {file << i+1 << " " << eltype << porder+1;}
        if (eltype == "pGFEMBar" || eltype == "pGFEMBar_WD_S" || eltype == "pGFEMBar_WD_M")
            {file << "_" << porder_Enrichment;}
        //propriedades
        file << " " << 1 << " ";
        //xGamma
        if (eltype == "pGFEMBar_WD_S" || eltype == "pGFEMBar_WD_M")
            {file << xgamma;}
        //Nodes
        for (int j {0}; j < n_per_el; j++)
            {file << i*(n_per_el-1) + j + 1 << " ";}
        file << "\n";
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
    for (int i {0}; i < d_bcs.size(); i++)
    {        
        file << i+1 << " ";
        if (d_bcs_pos[i] ==0)
            {file << 1 << " ";}
        else
            {file << nnodes << " ";}
        file << d_bcs_dofs[i] << " " <<d_bcs[i] << "\n";
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