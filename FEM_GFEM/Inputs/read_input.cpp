// ajustar entrada de E

#include <fstream>
#include "mesh.h"
#include "read_input.h"
#include "bars.h"
#include "Bondeary_conditions.h"
#include "body_func.h"
#include "Enrichment.h"
#include "node.h"

body_functions* get_body_function(int id, double alpha = 0., double xb = 0.)
{
    switch (id)
    {
        case 0:
            return new body_function0();
        case 1:
            return new body_function1();
        case 3:
            return new body_function3();
        case 10:
            return new body_function10(alpha, xb);
        case 12:
            return new body_function12();
        default:
            throw std::invalid_argument("Invalid body function ID");
    }
}

bool is_num(char c)
{
    if (c == '0' || c == '1'|| c == '2' || c == '3'|| c == '4' || c == '5'|| c == '6' || c == '7'|| c == '8' || c == '9')
        return true;
    else 
        return false;
}

int get_number_of (std::string& line)
{
    std::string temp {};
    for (std::size_t i {0}; i < line.size(); i++)
    {
        if (!(line[i] == ' ' || line[i] == '\t'))
            temp += line[i];

        if (line[i] == ' ' || line[i] == '\t' || i == line.size()-1)
        {
            if (temp != "")
            {
                return std::stoi(temp);
            }
        }
    }
    throw std::invalid_argument("No number found in line");
}

Node* get_node_by_id (int id, std::vector<Node>& nodevec)
{
    for (Node& no: nodevec)
    {
        if (no.id == id)
            return &no;
    }
    throw std::invalid_argument("Node not found");
}

template <typename T>
void pointer_vector_clear(std::vector<T*>& vec)
{
    for (T* obj : vec)
        delete obj;
}

// input de dados
void read_input (const std::string& filename, Mesh& mesh)
{
    // get file
    std::ifstream in_file(filename);
    if (!in_file)
    {
        // Print an error and exit
        throw std::runtime_error("input file could not be opened for reading!\n");
    }

    // extraindo cabeçalho
    std::getline(in_file, mesh.name);

    std::string line {};
    std::string temp {};

    std::string type;
    std::vector<std::string> header_arg;
    std::vector<std::string> leg;

    int case_i {0};
    std::size_t marker {};

    int id;
    double x_coord;

    std::vector<Node>& node_vec{mesh.nodes};
    std::vector<Element*>& el_vec{mesh.c_bars};
    std::vector<BC_displacement>& bc_vec{mesh.bc_ds};
    std::vector<BC_load>& loads{mesh.bc_l};

    std::vector<Enrichment*> temp_enr_vec{};

    std::vector<Properties> pr_vec{};

    std::vector<int> nid {};
    std::vector<int> lnid {};

    while (std::getline(in_file, line))
    {
        // extrair tipo, argumentos do cabeçalho e legenda
        case_i = 0;
        for (std::size_t i {0}; i < line.size(); i++)
        {
            if (!(line[i] == ' ' || line[i] == '\t' || line[i] == ';' || line[i] == '-'))
            {
                temp += line[i];
            }
            
            if (line[i] == ' ' || line[i] == '\t' || i == line.size()-1)
            {
                if (temp != "")
                {
                    if (case_i == 0)
                        type = temp;
                    else if (case_i == 1)
                        header_arg.push_back(temp);
                    else
                        leg.push_back(temp);

                    temp = "";
                }
            }
            else if (line[i] == '-')
            {
                if (case_i == 0)
                {
                    case_i++;
                    if (temp != "")
                    {
                        type = temp;
                        temp = "";
                    }
                }
                else
                    temp += line[i];
            }
            else if (line[i] == ';')
            {
                if (temp != "")
                {
                    if (case_i == 0)                        
                        type = temp;
                    else if (case_i == 1)
                        header_arg.push_back(temp);
                    temp = "";
                }
                case_i =2;
            }
        }
        
        // extrair objetos de acordo com o tipo
        case_i = 0;
        if (type == "nodes")
        {
            int nnodes {};
            int ndim {};

            temp = "";

            std::getline(in_file, line);

            // extrair argumentos do cabeçalho
            for (std::size_t i {0}; i < line.size(); i++)
            {
                if (!(line[i] == ' ' || line[i] == '\t'))
                    temp += line[i];
                
                if (line[i] == ' ' || line[i] == '\t' || i == line.size()-1)
                {
                    if (temp != "")
                    {
                        if (header_arg[case_i] == "nnodes")
                            nnodes = std::stoi(temp);
                        else if (header_arg[case_i] == "ndim")
                        {
                            ndim = std::stoi(temp);
                            if (ndim != 1)
                                throw std::invalid_argument("Only 1D nodes implemented");
                        }
                        else
                            throw std::invalid_argument("Unexpected header argument for nodes section");

                        case_i++;
                        temp = "";
                    }
                }
            }

            // extrair ordem dos tipos dados
            std::vector<int> order {};
            for (std::string& arg : leg)
            {
                if (arg == "nodeID")
                    order.push_back(0);
                else if (arg == "x-coord")
                    order.push_back(1);
                else if (arg == "y-coord")
                    order.push_back(2);
                else if (arg == "enrID")
                    order.push_back(3);
                else
                    throw std::invalid_argument("Unexpected legend argument for nodes section");
            }

            std::vector<int> enr_ids {};
            // criando nós
            for (int j {0}; j < nnodes; j++)
            {
                std::getline(in_file, line);

                temp = "";
                case_i = 0;

                for (std::size_t i {0}; i < line.size(); i++)
                {
                    if (!(line[i] == ' ' || line[i] == '\t'))
                        temp += line[i];

                    if (line[i] == ' ' || line[i] == '\t' || i == line.size()-1)
                    {
                        if (temp != "")
                        {
                            if (order[case_i] == 0)
                                id = std::stoi(temp);
                            else if (order[case_i] == 1)
                                x_coord = std::stod(temp);
                            else if (order[case_i] == 3)
                                enr_ids.push_back(std::stoi(temp));
                            temp = "";
                            case_i++;
                        }
                    }
                }
                node_vec.emplace_back(id, x_coord, enr_ids);
                enr_ids.clear();
            }
        }
            
        else if (type == "nelem")
        {
            // extrair número de elementos
            std::getline(in_file, line);

            int nelem {get_number_of(line)};

            // extrair ordem dos tipos dados
            std::vector<int> order {};
            for (std::string& arg : leg)
            {
                if (arg == "elemID")
                    order.push_back(0);
                else if (arg == "Type")
                    order.push_back(1);
                else if (arg == "propID")
                    order.push_back(2);
                else if (arg == "nodes")
                    order.push_back(3);
                else if (arg == "x-Gamma")
                    order.push_back(4);
                else
                    throw std::invalid_argument("Unexpected legend argument for elements section");
            }

            // criar elementos
            int case_j {4}; // número de nós por elemento

            int id {};
            std::string type {};
            int propID {};
            std::vector<int> nodes {};
            int shape_func_order {};
            int Eshape_func_order {};
            double xgamma {};

            for (int k {0}; k < nelem; k++)
            {
                std::getline(in_file, line);

                temp = "";
                case_i = 0;

                //obtendo parametros da linha
                for (std::size_t i {0}; i < line.size(); i++)
                {
                    if (!(line[i] == ' ' || line[i] == '\t'))
                        temp += line[i];
                    if (line[i] == ' ' || line[i] == '\t' || i == line.size()-1)
                    {
                        if (temp != "")
                        {
                            if (order[case_i] == 0)
                            {
                                id = std::stoi(temp);
                                case_i++;
                            }
                            else if (order[case_i] == 1)
                            {
                                std::string temp2 {};
                                std::string temp3 {};
                                while (is_num(temp.back()))
                                {
                                    temp2 += temp.back();
                                    temp.pop_back();
                                    shape_func_order = std::stoi(temp2)-1;
                                }

                                if (temp.back() == '_')
                                {   
                                    temp.pop_back();
                                    while (is_num(temp.back()))
                                    {
                                        temp3 += temp.back();
                                        temp.pop_back();
                                        shape_func_order = std::stoi(temp3)-1;
                                        Eshape_func_order = std::stoi(temp2);
                                    }
                                }
                                if (temp == "pBar" || temp == "lBar" || temp == "pGFEMBar" || temp == "pGFEMBar_sc" || temp == "pGFEMBar_WD_S" || temp == "pGFEMBar_WD_M")
                                {
                                    type = temp;
                                    if (temp == "lBar"  || temp == "pGFEMBar" || temp == "pGFEMBar_sc"  || temp == "pGFEMBar_WD_S" || temp == "pGFEMBar_WD_M")
                                        case_j = shape_func_order+1;
                                    else
                                        case_j = 3;

                                    case_i++;
                                }
                            }
                            else if (order[case_i] == 2)
                            {
                                propID = std::stoi(temp);
                                case_i++;
                            }
                            else if (order[case_i] == 3)
                            {                                
                                nodes.push_back(std::stoi(temp));
                                case_j--;
                                if (case_j <= 0)
                                    case_i++;
                            }
                            else if (order[case_i] == 4)
                            {                                
                                xgamma = std::stod(temp);
                                case_i++;
                            }
                            temp = "";
                        }
                    }
                }

                // descrevendo elementos
                if (type == "lBar")
                    el_vec.push_back(new lagrangian_bar(id, nodes, propID, shape_func_order));
                else if (type == "pBar")
                    el_vec.push_back(new p_hier_bar(id, nodes, propID, shape_func_order));
                else if (type == "pGFEMBar")
                    el_vec.push_back(new p_GFEM_bar(id, nodes, propID, shape_func_order, Eshape_func_order));
                else if (type == "pGFEMBar_sc")
                    el_vec.push_back(new p_GFEM_bar(id, nodes, propID, shape_func_order, Eshape_func_order, true, true));
                else if (type == "pGFEMBar_WD_S")
                    el_vec.push_back(new p_GFEM_bar_weak_disc(id, nodes, propID, shape_func_order, Eshape_func_order, new Sukumar_enrichment_1D(0, xgamma)));
                else if (type == "pGFEMBar_WD_M")
                    el_vec.push_back(new p_GFEM_bar_weak_disc(id, nodes, propID, shape_func_order, Eshape_func_order, new Moes_enrichment_1D(0, xgamma)));
                else
                    throw std::invalid_argument("Unexpected element type (" + type + ")");
                nodes.clear();
            }
        }

        else if (type == "nenrichments")
        {
            // extrair número de elementos
            std::getline(in_file, line);

            int nenr {get_number_of(line)};

            // extrair ordem dos tipos dados
            std::vector<int> order {};
            for (std::string& arg : leg)
            {
                if (arg == "enrID")
                    order.push_back(0);
                else if (arg == "Type")
                    order.push_back(1);
                else if (arg == "xGamma")
                    order.push_back(2);
                else
                    throw std::invalid_argument("Unexpected legend argument for elements section");
            }

            // criar padrões de enriquecimento
            int id {};
            std::string type {};
            double xGamma {};

            for (int k {0}; k < nenr; k++)
            {
                std::getline(in_file, line);

                temp = "";
                case_i = 0;

                //obtendo parametros da linha
                for (std::size_t i {0}; i < line.size(); i++)
                {
                    if (!(line[i] == ' ' || line[i] == '\t'))
                        temp += line[i];
                    if (line[i] == ' ' || line[i] == '\t' || i == line.size()-1)
                    {
                        if (temp != "")
                        {
                            if (order[case_i] == 0)
                            {
                                id = std::stoi(temp);
                                case_i++;
                            }
                            else if (order[case_i] == 1)
                            {
                                type = temp;
                                case_i++;
                            }
                            else if (order[case_i] == 2)
                            {
                                xGamma = std::stod(temp);
                                case_i++;
                            }
                            temp = "";
                        }
                    }
                }

                // descrevendo elementos
                if (type == "ESuk")
                    temp_enr_vec.push_back(new Sukumar_enrichment_1D(id, xGamma));
                else
                    throw std::invalid_argument("Unexpected enrichment type (" + type + ")");
            }
        }

        else if (type == "properties")
        {
            temp = "";

            std::getline(in_file, line);
            
            int nprop {get_number_of(line)};

            // extrair ordem dos tipos dados
            std::vector<int> order {};
            for (std::string& arg : leg)
            {
                if (arg == "propID")
                    order.push_back(0);
                else if (arg == "type")
                    order.push_back(1);
                else if (arg == "E") //colocar mais de um E e Exlim
                    order.push_back(2);
                else if (arg == "Exlim")
                    order.push_back(3);
                else if (arg == "A")
                    order.push_back(4);
                else if (arg == "C")
                    order.push_back(5);
                else if (arg == "bf_fun")
                    order.push_back(6);
                else if (arg == "alpha")
                    order.push_back(7);
                else if (arg == "xb")
                    order.push_back(8);
                else
                    throw std::invalid_argument("Unexpected legend argument for elements section");
            }

            for (int j {0}; j < nprop; j++)
            {
                std::getline(in_file, line);

                temp = "";
                case_i = 0;

                pr_vec.emplace_back();

                int bf_func_id {};
                double alpha {0.};
                double xb {0.};

                //obtendo parametros da linha
                for (std::size_t i {0}; i < line.size(); i++)
                {
                    if (!(line[i] == ' ' || line[i] == '\t'))
                        temp += line[i];
                    
                    if (line[i] == ' ' || line[i] == '\t' || i == line.size()-1)
                    {
                        if (temp != "")
                        {
                            if (order[case_i] == 0)
                                pr_vec[j].id = std::stoi(temp);
                            else if (order[case_i] == 1)
                            {
                                if (temp == "MatBar")
                                    type = temp;
                                else
                                    throw std::invalid_argument("Unexpected property type (" + temp + ")");
                            }
                            else if (order[case_i] == 2)
                                pr_vec[j].E.push_back(std::stod(temp));
                            else if (order[case_i] == 3)
                                pr_vec[j].Exlim.push_back(std::stod(temp));
                            else if (order[case_i] == 4)                             
                                pr_vec[j].A = std::stod(temp);
                            else if (order[case_i] == 5)                              
                                pr_vec[j].C = std::stod(temp);
                            else if (order[case_i] == 6)                               
                                bf_func_id = std::stoi(temp);
                            else if (order[case_i] == 7)                              
                                alpha = std::stod(temp);
                            else if (order[case_i] == 8)                            
                                xb = std::stod(temp);
                            temp = "";
                            case_i++;
                        }
                    }
                }
                pr_vec[j].bf_func = get_body_function(bf_func_id, alpha, xb);
            }
        }

        else if (type == "constraints")
        {
            std::getline(in_file, line);
            
            int nconst {get_number_of(line)};

            // extrair ordem dos tipos dados
            std::vector<int> order {};
            for (std::string& arg : leg)
            {
                if (arg == "constrID")
                    order.push_back(0);
                else if (arg == "nodeID")
                    order.push_back(1);
                else if (arg == "dof")
                    order.push_back(2);
                else if (arg == "value")
                    order.push_back(3);
                else
                    throw std::invalid_argument("Unexpected legend argument for elements section");
            }

            for (int j {0}; j < nconst; j++)
            {
                std::getline(in_file, line);

                temp = "";
                case_i = 0;

                BC_displacement bc {};

                //obtendo parametros da linha
                for (std::size_t i {0}; i < line.size(); i++)
                {
                    if (!(line[i] == ' ' || line[i] == '\t'))
                        temp += line[i];
                    if (line[i] == ' ' || line[i] == '\t' || i == line.size()-1)
                    {
                        if (temp != "")
                        {
                            if (order[case_i] == 0)
                                bc.id = std::stoi(temp);
                            else if (order[case_i] == 1)
                                nid.push_back(std::stoi(temp));
                            else if (order[case_i] == 2)
                                bc.dof = std::stoi(temp) -1;
                            else if (order[case_i] == 3)                               
                                bc.value = std::stod(temp);
                            temp = "";
                            case_i++;
                        }
                    }
                }
                bc_vec.push_back(bc);
            }
        }

        else if (type == "loads")
        {
            std::getline(in_file, line);
            
            int nconst {get_number_of(line)};

            // extrair ordem dos tipos dados
            std::vector<int> order {};
            for (std::string& arg : leg)
            {
                if (arg == "loadID")
                    order.push_back(0);
                else if (arg == "nodeID")
                    order.push_back(1);
                else if (arg == "dof")
                    order.push_back(2);
                else if (arg == "value")
                    order.push_back(3);
                else
                    throw std::invalid_argument("Unexpected legend argument for elements section");
            }

            for (int j {0}; j < nconst; j++)
            {
                std::getline(in_file, line);

                temp = "";
                case_i = 0;

                BC_load bc {};

                //obtendo parametros da linha
                for (std::size_t i {0}; i < line.size(); i++)
                {
                    if (!(line[i] == ' ' || line[i] == '\t'))
                        temp += line[i];
                    
                    if (line[i] == ' ' || line[i] == '\t' || i == line.size()-1)
                    {
                        if (temp != "")
                        {
                            if (order[case_i] == 0)
                                bc.id = std::stoi(temp);
                            else if (order[case_i] == 1)
                                lnid.push_back(std::stoi(temp));
                            else if (order[case_i] == 2)
                                bc.dof = std::stoi(temp) -1;
                            else if (order[case_i] == 3)                    
                                bc.value = std::stod(temp);
                            temp = "";
                            case_i++;
                        }
                    }
                }
                loads.push_back(bc);
            }
        }
        
        header_arg.clear();
        leg.clear();
        type = "";
    }

    for (std::size_t i {0}; i<bc_vec.size();i++)
    {
        bc_vec[i].assign_node(get_node_by_id (nid[i], node_vec));
    }
    for (std::size_t i {0}; i<loads.size();i++)
    {
        loads[i].assign_node(get_node_by_id (lnid[i], node_vec));
    }
    
    sort_Enr_by_id(temp_enr_vec);
    // atribuir os enriquecimentos aos nós
    for (Node& no: node_vec)
    {
        for (int enr_id: no.enr_ids)
        {
            for (Enrichment* enr: temp_enr_vec)
            {
                if (enr->id == enr_id)
                {
                    enr->assign_to_node(no);
                    break;
                }
            }
        }
    }

    int dof0 {0};
    for (Element* el: el_vec)
        {el->start_el(node_vec, dof0, pr_vec);}
        
    // extrair número de dofs
    mesh.set_dofs(dof0);

    pointer_vector_clear(temp_enr_vec);
}