// ajustar entrada de E

#include <fstream>
#include "mesh.h"
#include "read_input.h"
#include "constrained_bar.h"
#include "Bondeary_conditions.h"
#include "body_func.h"

body_functions* get_body_function(int id)
{
    switch (id)
    {
        case 0:
            return new body_function0();
        case 1:
            return new body_function1();
        case 3:
            return new body_function3();
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

#include "node.h"

Node* get_node_by_id (int id, std::vector<Node>& nodevec)
{
    for (Node& no: nodevec)
    {
        if (no.id == id)
            return &no;
    }
    throw std::invalid_argument("Node not found");
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
    std::vector<Constrained_bar>& el_vec{mesh.c_bars};
    std::vector<BC_displacement>& bc_vec{mesh.bc_ds};
    std::vector<BC_load>& loads{mesh.bc_l};

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
                else
                    throw std::invalid_argument("Unexpected legend argument for nodes section");
            }

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
                            temp = "";
                            case_i++;
                        }
                    }
                }
                node_vec.emplace_back(id, x_coord);
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
                else
                    throw std::invalid_argument("Unexpected legend argument for elements section");
            }

            // criar elementos
            int case_j {4}; // número de nós por elemento

            int id {};
            std::string type {};
            int propID {};
            std::vector<int> nodes {};

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
                                while (is_num(temp.back()))
                                {
                                    temp2 += temp.back();
                                    temp.pop_back();
                                }
                                
                                if (temp == "pBar")
                                {
                                    type = temp;
                                    case_j = std::stoi(temp2);
                                    case_i++;
                                }
                                else
                                    throw std::invalid_argument("Unexpected element type");
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
                                if (case_j == 0)
                                    case_i++;
                            }
                            temp = "";
                        }
                    }
                }

                // descrevendo elementos
                if (type == "pBar")
                    el_vec.emplace_back(id, nodes, propID);
                nodes.clear();
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
                double alpha {};
                double xb {};

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
                                if (temp == "MatpBar")
                                    type = temp;
                                else
                                    throw std::invalid_argument("Unexpected element type");
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
                if (bf_func_id == 10)
                    pr_vec[j].bf_func = new body_function10(alpha, xb);
                else
                    pr_vec[j].bf_func = get_body_function(bf_func_id);
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
    for (Constrained_bar& bar: el_vec)
    {
        bar.get_nodes(node_vec);
        for (Properties& pr: pr_vec)
        {
            if (pr.id == bar.prop_id)
            {
                bar.get_properties(pr);
                break;
            }
        }
        bar.start_el();
    }
}