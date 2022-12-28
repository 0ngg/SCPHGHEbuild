#include"scphghe/temp_struct_cfd.h"
#include"windows.h"

using namespace cfd;

// util
template <typename U, class V>
std::vector<U> get_map_keys(V map)
{
    std::vector<U> get;
    for(auto const& entry : map)
    {
        get.push_back(entry.first);
    };
    return get;
};
double calc_fluid_prop(std::string what, double P, double T, double W)
{
    // load, export, and definition of function HAPropsSI from CoolProp.DLL
    // Define DLL functions
    typedef double (WINAPI *HAPropsSI)(const char* Output, const char* Name1, double Prop1, const char* Name2, double Prop2, const char* Name3, double Prop3);
    // addresses
    HAPropsSI HAPropsSIAddress;
    // load DLL; change this path as needed
    HINSTANCE CoolPropDll = LoadLibraryA("C:\\Users\\LENOVO\\Documents\\TA\\sandbox\\make_scheme_fix\\lib\\CoolProp.dll");
    // number of HAPropsSI args bytes = 8 (double) * 3 + 4 (char) * 4 = 40
    HAPropsSIAddress = (HAPropsSI) GetProcAddress(CoolPropDll, "_HAPropsSI@40"); // change address (check CoolPropLib.h)

    // rho, miu, cp, eps, alpha
    // W = absolute humidity
    if(what.compare(std::string("rho")) == 0)
    {
        double get = (*HAPropsSIAddress) ("Vha", "P", P, "T", T, "W", W);
        FreeLibrary(CoolPropDll);
        return 1 / get;
    }
    else if(what.compare(std::string("miu")) == 0)
    {
        double get = (*HAPropsSIAddress) ("mu", "P", P, "T", T, "W", W);
        FreeLibrary(CoolPropDll);
        return get;
    }
    else if(what.compare(std::string("cp")) == 0)
    {
        double get = (*HAPropsSIAddress) ("cp_ha", "P", P, "T", T, "W", W);
        FreeLibrary(CoolPropDll);
        return get;
    }
    else if(what.compare(std::string("k")) == 0)
    {
        double get = (*HAPropsSIAddress) ("k", "P", P, "T", T, "W", W);
        FreeLibrary(CoolPropDll);
        return get;
    }
    else if(what.compare(std::string("alpha")) == 0)
    {
        // k / (rho * cp)
        double k__ = (*HAPropsSIAddress) ("k", "P", P, "T", T, "W", W);
        double rho__ = 1 / (*HAPropsSIAddress) ("Vha", "P", P, "T", T, "W", W);
        double cp__ = (*HAPropsSIAddress) ("cp_ha", "P", P, "T", T, "W", W);
        FreeLibrary(CoolPropDll);
        return k__ / (rho__ * cp__);
    }
    else
    {
        FreeLibrary(CoolPropDll);
        return 0.0;
    };
};

// axes
double square_sum (double x, double y)
{
    return x + y * y;
};
double sqrt_sum(coor v)
{
    return pow(pow(v(0), 2) + pow(v(1), 2) + pow(v(2), 2), 0.5);
};
axes::axes(make<double>::sp_mat& in)
{
    this->x = in;
    this->y = in;
    this->z = in;
};
axes::axes(const axes& other)
{
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
};
axes& axes::operator()(int row, int col, coor& in)
{
    this->x.coeffRef(row, col) = in(0);
    this->y.coeffRef(row, col) = in(1);
    this->z.coeffRef(row, col) = in(2);
    return *this;
};
axes& axes::operator=(int x)
{
    return *this;
};
axes& axes::operator+=(const axes& other)
{
    this->x += other.x;
    this->y += other.y;
    this->z += other.z;
    return *this;
};
coor axes::axes_to_coor(int row, int col)
{
    coor ax_coor(this->x.coeffRef(row, col), this->y.coeffRef(row, col), this->z.coeffRef(row, col));
    return ax_coor;
};
double axes::axes_to_val(int row, int col)
{
    make<double>::vec ax_vec;
    ax_vec.push_back(this->x.coeffRef(row, col));
    ax_vec.push_back(this->y.coeffRef(row, col));
    ax_vec.push_back(this->z.coeffRef(row, col));
    coor ax_coor(ax_vec.data());
    double ax_val = sqrt_sum(ax_coor);
    return ax_val;
};

// user
make<make<std::string>::vec>::vec parse_csv(std::string filename)
{
    // output strings, revert to double later
    std::ifstream data(filename);
    std::string line;
    make<make<std::string>::vec>::vec parsed;
    while(std::getline(data, line))
    {
        std::stringstream lineStream(line);
        std::string element;
        make<std::string>::vec parsed_row;
        while(std::getline(lineStream, element, ','))
        {
            parsed_row.push_back(element);
        };
        parsed.push_back(parsed_row);
    };
    return parsed;
};
user::user(double P_init_in, double T_init_in, double W_init_in, std::string source_file,
        std::string prop_solid_file): P_init(P_init_in), T_init(T_init_in), W_init(W_init_in)
{
    this->read_source_csv(source_file);
    this->read_solid_prop_csv(prop_solid_file);
};
void user::read_source_csv(std::string source_file)
{
    make<make<std::string>::vec>::vec parsed = parse_csv(source_file);
    // make<make<make<double>::map_int>::map_int>::map_str face_source; // .... - time step - value
    // make<make<make<double>::map_int>::map_str>::map_str cell_source; // .... - time step - value
    for(auto i = parsed.begin(); i != parsed.end(); i++)
    {
        make<std::string>::vec temp = *i;
        char test_ifdigit = temp[1][0];
        if(std::isdigit(test_ifdigit))
        {
            // face_source
            make<std::string>::vec this_face_source_keys = get_map_keys<std::string, make<make<make<double>::map_int>::map_int>::map_str>(this->face_source);
            if(std::find(this_face_source_keys.begin(), this_face_source_keys.end(), temp[0]) == this_face_source_keys.end())
            {
                this->face_source.insert({temp[0], make<make<double>::map_int>::map_int()});
                make<double>::map_int value_vec;
                int ctd = 0;
                for(int j = 2; j < temp.size(); j++)
                {
                    value_vec.insert({ctd, std::stod(temp[j])});
                    ctd += 1;
                };
                this->face_source[temp[0]].insert({std::stoi(temp[1]), make<double>::map_int(value_vec)});
            }  
            else
            {
                make<double>::map_int value_vec;
                int ctd = 0;
                for(int j = 2; j < temp.size(); j++)
                {
                    value_vec.insert({ctd, std::stod(temp[j])});
                    ctd += 1;
                };
                this->face_source[temp[0]].insert({std::stoi(temp[1]), make<double>::map_int(value_vec)});
            };
        }
        else
        {
            // cell_source
            make<std::string>::vec this_cell_source_keys = get_map_keys<std::string, make<make<make<double>::map_int>::map_str>::map_str>(this->cell_source);
            if(std::find(this_cell_source_keys.begin(), this_cell_source_keys.end(), temp[0]) == this_cell_source_keys.end())
            {
                this->cell_source.insert({temp[0], make<make<double>::map_int>::map_str()});
                make<double>::map_int value_vec;
                int ctd = 0;
                for(int j = 2; j < temp.size(); j++)
                {
                    value_vec.insert({ctd, std::stod(temp[j])});
                    ctd += 1;
                };
                this->cell_source[temp[0]].insert({temp[1], make<double>::map_int(value_vec)});
            }
            else
            {
                make<double>::map_int value_vec;
                int ctd = 0;
                for(int j = 2; j < temp.size(); j++)
                {
                    value_vec.insert({ctd, std::stod(temp[j])});
                    ctd += 1;
                };
                this->cell_source[temp[0]].insert({temp[1], make<double>::map_int(value_vec)});  
            };
        };
    };
};
void user::read_solid_prop_csv(std::string prop_file)
{
    make<make<std::string>::vec>::vec parsed = parse_csv(prop_file);
    for(auto i = parsed.begin(); i != parsed.end(); i++)
    {
        make<std::string>::vec temp = *i;
        char test_ifdigit = temp[1][0];
        if(std::isdigit(test_ifdigit))
        {
            // s2s_eps
            this->s2s_eps.insert({std::stoi(temp[1]), std::stod(temp[2])});
        }
        else
        {
            // solid_k
            this->solid_k.insert({temp[1], std::stod(temp[2])});
        };
    };
};
void user::update_source(int current_time, scheme*& scheme_ref)
{
    binfo& source__ = *scheme_ref->source;
    for(std::pair<std::string, make<double>::map_int> entry_face : source__.face_value)
    {
        for(std::pair<int, double> entry_face_int : entry_face.second)
        {
            entry_face_int.second = this->face_source[entry_face.first][entry_face_int.first][current_time];
        };
    };
    for(std::pair<std::string, make<double>::map_str> entry_cell : source__.cell_value)
    {
        for(std::pair<std::string, double> entry_cell_str : entry_cell.second)
        {
            entry_cell_str.second = this->cell_source[entry_cell.first][entry_cell_str.first][current_time];
        };
    };
};

// scheme other functions
bool is_unique(std::string s)
{
    return !s.empty() && std::find_if(s.begin(), s.end(),
    [](char c) {return std::isdigit(c);}) != s.end(); 
};
int is_neighbor(make<int>::vec& v1, make<int>::vec& v2)
{
    for(auto i = v1.begin(); i != v1.end(); i++)
    {
        if(std::find(v2.begin(), v2.end(), *i) != v1.end())
        {
            return *i;
        };
    };
    return -1;
};
make<std::pair<std::string, int>>::vec is_boundary(int fid, make<make<make<int>::vec>::map_int>::map_str& fid_, make<std::string>::vec& reserve)
{
    make<std::pair<std::string, int>>::vec get = make<std::pair<std::string, int>>::vec();
    for(std::pair<std::string, make<make<int>::vec>::map_int> entry1 : fid_)
    {
        if(std::find(reserve.begin(), reserve.end(), entry1.first) != reserve.end())
        {
            for(std::pair<int, make<int>::vec> entry2 : fid_[entry1.first])
            {
                for(auto i = entry2.second.begin(); i != entry2.second.end(); i++)
                {
                    if(*i == fid)
                    {
                        std::pair<std::string, int> temp(entry1.first, entry2.first);
                        get.push_back(temp);
                    };
                };
            };
        };
    };
    return get;
};
struct quad
{
    make<coor>::vec centroid;
    make<double>::vec area;
};
void make_tri(finfo& face_, make<coor>::map_int& nodes_)
{
    const make<int>::vec& fnode_ = face_.fnode;
    coor centroid(0.0, 0.0, 0.0);
    coor v1 = nodes_[fnode_[1]] - nodes_[fnode_[0]];
    coor v2 = nodes_[fnode_[2]] - nodes_[fnode_[0]];
    coor vcross = v1.cross(v2);
    for(auto i = fnode_.begin(); i != fnode_.end(); i++)
    {
        centroid += nodes_[*i];
    };
    centroid = centroid/3;
    vcross(0) = pow(vcross(0), 2);
    vcross(1) = pow(vcross(1), 2);
    vcross(2) = pow(vcross(2), 2);
    double area = 0.5 * pow(vcross.block(0, 0, 3, 1).sum(), 0.5);
    coor vnorm = v1.cross(v2);
    vnorm.norm();
    face_.fcentroid = centroid;
    face_.fnormal = vnorm;
    face_.farea = area;
};
void make_quad(finfo& face_, make<coor>::map_int& nodes_)
{
    quad fquad;
    make<int>::vec& fnode_ = face_.fnode;
    std::pair<int, int> diag1;
    coor acc; double len; double rhs;
    for(int i = 0; i < 3; i++)
    {
        for(int j = i; j < 4; j++)
        {
            acc = nodes_[fnode_[i]] + nodes_[fnode_[j]];
            rhs = sqrt_sum(acc);
            if(rhs > len)
            {
                len = rhs;
                diag1.first = fnode_[i];
                diag1.second = fnode_[j];
            };
        };
    };
    double area; make<int>::vec diag2;
    for(auto i = fnode_.begin(); i != fnode_.end(); i++)
    {
        if(*i != diag1.first && *i != diag1.second)
        {
            coor centroid(0.0, 0.0, 0.0);
            centroid += nodes_[diag1.first] + nodes_[diag1.second] + nodes_[*i];
            centroid = centroid/3;
            coor v1 = nodes_[diag1.first] - nodes_[*i];
            coor v2 = nodes_[diag1.second] - nodes_[*i];
            coor vcross = v1.cross(v2);
            area = 0.5 * sqrt_sum(vcross);
            fquad.centroid.push_back(centroid);
            fquad.area.push_back(area);
        }
        else
        {
            diag2.push_back(*i);
        };
    };
    area = std::accumulate(fquad.area.begin(), fquad.area.end(), 0.0);
    coor centroid(0.0, 0.0, 0.0);
    for(size_t i = 0; i < sizeof(fquad.centroid); i++)
    {
        centroid += fquad.centroid[i] * fquad.area[i]; 
    };
    centroid = centroid / area;
    coor v1 =  nodes_[0] - nodes_[diag1.first];
    coor v2 =  nodes_[0] - nodes_[diag1.second];
    coor vnorm = v1.cross(v2);
    vnorm.norm();
    face_.fcentroid = centroid;
    face_.fnormal = vnorm;
    face_.farea = area;
};
void make_cell(minfo& mesh_, cinfo cell_, make<finfo>::map_int& faces_)
{
    coor centroid(0.0, 0.0, 0.0);
    double area = 0.0;
    for(std::pair<int, finfo> entry : faces_)
    {
        centroid += entry.second.fcentroid * entry.second.farea;
        area += entry.second.farea;
    };
    centroid = centroid/area;
    double vol = 0.0; coor parallel;
    for(std::pair<int, finfo> entry : faces_)
    {
        coor norm = entry.second.fnormal;
        coor vinner = entry.second.fcentroid - centroid;
        if(norm.dot(vinner) < 0)
        {
            norm = norm * -1;
        };
        vol += centroid.dot(norm) * entry.second.farea / 3;
    };
    cell_.ccentroid = centroid;
    cell_.cvolume = vol;
};
// scheme
scheme::scheme(std::string msh_file, user*& user_p)
{
    mshio::MshSpec spec = mshio::load_msh(msh_file);
    minfo mesh_(spec); pinfo prop_(mesh_, user_p); winfo wall_(mesh_);
    const char* temp_which = "fluid";
    std::string temp_which_str(temp_which);
    make<std::string>::vec temp_which_vec;
    temp_which_vec.push_back(temp_which_str);
    vinfo pressure_(temp_which_vec, mesh_, user_p->P_init);
    make<make<double>::map_int>::map_str face_value__;
    make<std::string>::vec face_value_keys = get_map_keys<std::string, make<make<double>::map_int>::map_str>(face_value__);
    for(std::pair<std::string, make<make<double>::map_int>::map_int> entry_face1 : user_p->face_source)
    {
        if(std::find(face_value_keys.begin(), face_value_keys.end(), entry_face1.first) == face_value_keys.end())
        {
            face_value__.insert({entry_face1.first, make<double>::map_int()});
        };
        for(std::pair<int, make<double>::map_int> entry_face2 : entry_face1.second)
        {
            face_value__[entry_face1.first].insert({entry_face2.first, entry_face2.second[0]});                  
        };
    };
    make<make<double>::map_str>::map_str cell_value__;
    make<std::string>::vec cell_value_keys = get_map_keys<std::string, make<make<double>::map_str>::map_str>(cell_value__);
    for(std::pair<std::string, make<make<double>::map_int>::map_str> entry_cell1 : user_p->cell_source)
    {
        if(std::find(cell_value_keys.begin(), cell_value_keys.end(), entry_cell1.first) == cell_value_keys.end())
        {
            cell_value__.insert({entry_cell1.first, make<double>::map_str()});
        };
        for(std::pair<std::string, make<double>::map_int> entry_cell2 : entry_cell1.second)
        {
            cell_value__[entry_cell1.first].insert({entry_cell2.first, entry_cell2.second[0]});                  
        };
    };
    this->mesh = new minfo(mesh_);
    this->prop = new pinfo(prop_);
    this->wall = new winfo(wall_);
    this->pressure = new vinfo(pressure_);
    this->source = new binfo(face_value__, cell_value__);
};
// minfo
minfo::minfo(mshio::MshSpec spec)
{
    this->make_minfo(spec);
};
void minfo::make_minfo(mshio::MshSpec spec)
{
    // make map int of physical group name
    make<std::string>::map_int phys_name_map;
    for(const auto& group: spec.physical_groups)
    {
        phys_name_map.insert({group.tag, group.name});
    };
    make<coor>::map_int nodes;
    make<finfo>::map_int faces;
    make<cinfo>::map_int cells;
    for(int i = 0; i < spec.nodes.num_entity_blocks; i++)
    {
        mshio::NodeBlock& block = spec.nodes.entity_blocks[i];
        for(int j = 0; j < block.num_nodes_in_block; j++)
        {
            coor node__(block.data[j*3+0], block.data[j*3+1], block.data[j*3+2]);
            nodes.insert({block.tags[j] - 1, node__});
        };
    };
    make<int>::map_int type_to_nface{{5,6}, {6,5}};
    int tag; int type; int n;
    for(int i = 0; i < spec.elements.num_entity_blocks; i++)
    {
        mshio::ElementBlock& block = spec.elements.entity_blocks[i];
        tag = block.entity_tag;
        type = block.element_type;
        n = mshio::nodes_per_element(type);
        switch(block.entity_dim)
        {
            case 2:
            for(int j = 0; j < block.num_elements_in_block; j++)
            {
                finfo face__;
                make<int>::vec fnode__{int(block.data[j * (n+1) - 1]), int(block.data[j * (n+1)]), int(block.data[j * (n+1) + 1])};
                face__.fnode = fnode__; face__.fboundary = phys_name_map[tag];
                faces.insert({block.data[j * (n+1)] - 1, face__});
            };
            case 3:
            for(int j = 0; j < block.num_elements_in_block; j++)
            {
                cinfo cell__;
                make<int>::vec cnode__; for(int k = 0; k < n; k++) {cnode__.push_back(block.data[j * (n+1) + k - 1]);};
                make<int>::vec cface__; 
                for(std::pair<int, finfo> entry : faces)
                {
                    int ctd = 0;
                    for(auto k = entry.second.fnode.begin(); k != entry.second.fnode.end(); k++)
                    {
                        if(std::find(cnode__.begin(), cnode__.end(), *k) != cnode__.end())
                        {
                            ctd += 1;
                        };
                    };
                    if(ctd == int(sizeof(entry.second.fnode)))
                    {
                        cface__.push_back(entry.first);
                    }
                    else if(int(sizeof(cface__)) == type_to_nface[type])
                    {
                        break;
                    };
                };
                cell__.cnode = cnode__; cell__.cface = cface__; cell__.cdomain = phys_name_map[tag];
                cells.insert({block.data[j * (n+1)] - 1, cell__});
            };
            default:
            return;
        };
    };
    this->nodes = nodes;
    this->faces = faces;
    this->cells = cells;
    this->make_id(spec, cells);
    this->make_template(nodes, faces, cells);
    this->make_size(nodes, faces, cells);
    this->make_geom(faces, cells);
    this->make_constants();
    make<std::string>::vec this_fid_str_keys = get_map_keys<std::string, make<make<make<int>::vec>::map_int>::map_str>(this->fid);
    if(std::find(this_fid_str_keys.begin(), this_fid_str_keys.end(), std::string("s2s")) != this_fid_str_keys.end())
    {
       this->make_clust(faces);
    };
};
void minfo::make_id(mshio::MshSpec& spec, make<cinfo>::map_int& cells_)
{
    // physical names to boundary and domain id
    make<make<make<int>::vec*>::map_int>::map_str fid__;
    make<make<make<int>::vec*>::map_str>::map_str cid__;
    make<make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>::map_str bid__; // store cell-face-<boundary name, unique-constant/iter value>
    make<make<make<int>::vec*>::vec>::map_int it_id_tag;
    make<int>::vec* temp_fid_p = NULL;
    make<int>::vec* temp_cid_p = NULL;
    for(const auto& group : spec.physical_groups)
    {
        make<std::string>::vec parsed__;
        std::stringstream ssin(group.name);
        std::string token;
        make<make<int>::vec*>::vec it_id_l;
        switch(group.dim)
        {
            case 2:
            while (ssin >> token)
            {
                if(is_unique(token))
                {
                    std::string first = token.substr(0, sizeof(token) - 2);
                    char num(token[sizeof(token) - 1]);
                    int second(num);
                    make<std::string>::vec fid_str_keys = get_map_keys<std::string, make<make<make<int>::vec*>::map_int>::map_str>(fid__);
                    if(std::find(fid_str_keys.begin(), fid_str_keys.end(), first) == fid_str_keys.end())
                    {
                        fid__.insert({first, make<make<int>::vec*>::map_int()});
                        fid__[first].insert({second, temp_fid_p});
                        fid__[first][second] = new make<int>::vec();
                    }
                    else
                    {
                        make<int>::vec fid_int_keys = get_map_keys<int, make<make<int>::vec*>::map_int>(fid__[first]);
                        if(std::find(fid_int_keys.begin(), fid_int_keys.end(), second) == fid_int_keys.end())
                        {
                            fid__[first].insert({second, temp_fid_p});
                            fid__[first][second] = new make<int>::vec();
                        };
                    };
                    it_id_l.insert(it_id_l.end(), fid__[first][second]);
                }
                else
                {  
                    make<std::string>::vec fid_str_keys = get_map_keys<std::string, make<make<make<int>::vec*>::map_int>::map_str>(fid__);
                    if(std::find(fid_str_keys.begin(), fid_str_keys.end(), token) == fid_str_keys.end())
                    {
                        fid__.insert({token, make<make<int>::vec*>::map_int()});
                        fid__[token].insert({0, temp_fid_p});
                        fid__[token][0] = new make<int>::vec;
                    };
                    it_id_l.insert(it_id_l.end(), fid__[token][0]);
                };
            };
            case 3:
            while (ssin >> token)
            {
                std::string first = token.substr(0, 5);
                std::string second = token.substr(6, sizeof(token) - 6);
                make<std::string>::vec cid_str1_keys = get_map_keys<std::string, make<make<make<int>::vec*>::map_str>::map_str>(cid__);
                if(std::find(cid_str1_keys.begin(), cid_str1_keys.end(), token) == cid_str1_keys.end())
                {
                    cid__.insert({first, make<make<int>::vec*>::map_str()});
                    cid__[first].insert({second, temp_cid_p});
                    cid__[first][second] = new make<int>::vec();
                }
                else
                {
                    make<std::string>::vec cid_str2_keys = get_map_keys<std::string, make<make<int>::vec*>::map_str>(cid__[first]);
                    if(std::find(cid_str2_keys.begin(), cid_str2_keys.end(), token) == cid_str2_keys.end())
                    {
                        cid__[first].insert({second, temp_cid_p});
                        cid__[first][second] = new make<int>::vec();
                    };
                };
                it_id_l.insert(it_id_l.end(), cid__[first][second]);
            };
        };
    };
    fid__.insert({"none", make<make<int>::vec*>::map_int()});
    fid__["none"].insert({0, temp_fid_p});
    fid__["none"][0] = new make<int>::vec();
    int tag; int type; int n;
    for(int i = 0; i < spec.elements.num_entity_blocks; i++)
    {
        mshio::ElementBlock& block = spec.elements.entity_blocks[i];
        tag = block.entity_tag;
        type = block.element_type;
        n = mshio::nodes_per_element(type);
        make<int>::vec it_id_tag_keys = get_map_keys<int, make<make<make<int>::vec*>::vec>::map_int>(it_id_tag);
        if(std::find(it_id_tag_keys.begin(), it_id_tag_keys.end(), tag) != it_id_tag_keys.end())
        {
            for(auto j = it_id_tag[tag].begin(); j != it_id_tag[tag].end(); j++)
            {
                make<int>::vec*& temp_ref = *j;
                for(int k = 0; k < block.num_elements_in_block; k++)
                {
                    temp_ref->push_back(block.data[k*(n+1)]);
                };
            };
        }
        else
        {
            for(int k = 0; k < block.num_elements_in_block; k++)
            {
                fid__["none"][0]->push_back(block.data[k*(n+1)]);
            };
        };
    };
    for(std::pair<std::string, make<make<int>::vec*>::map_int> entry_f1 : fid__)
    {
        for(auto it = entry_f1.second.begin(); it != entry_f1.second.end(); ++it)
        {
            make<int>::vec to_fid__;
            make<std::string>::vec this_fid_keys = get_map_keys<std::string, make<make<make<int>::vec>::map_int>::map_str>(this->fid);
            if(std::find(this_fid_keys.begin(), this_fid_keys.end(), std::string(entry_f1.first)) == this_fid_keys.end())
            {
                for(auto i = it->second->begin(); i != it->second->end(); i++)
                {
                    to_fid__.push_back(*i);
                };
                this->fid.insert({entry_f1.first, make<make<int>::vec>::map_int()});
                this->fid[entry_f1.first].insert({it->first, make<int>::vec(to_fid__)});
            }
            else
            {
                for(auto i = it->second->begin(); i != it->second->end(); i++)
                {
                    to_fid__.push_back(*i);
                };
                this->fid[entry_f1.first].insert({it->first, make<int>::vec(to_fid__)});
            };
        };
    };
    for(std::pair<std::string, make<make<int>::vec*>::map_str> entry_c1 : cid__)
    {
        for(auto it = entry_c1.second.begin(); it != entry_c1.second.end(); ++it)
        {
            make<int>::vec to_cid__;
            make<std::string>::vec this_cid_keys = get_map_keys<std::string, make<make<make<int>::vec>::map_str>::map_str>(this->cid);
            if(std::find(this_cid_keys.begin(), this_cid_keys.end(), std::string(entry_c1.first)) == this_cid_keys.end())
            {
                for(auto i = it->second->begin(); i != it->second->end(); i++)
                {
                    to_cid__.push_back(*i);
                };
                this->cid.insert({entry_c1.first, make<make<int>::vec>::map_str()});
                this->cid[entry_c1.first].insert({it->first, make<int>::vec(to_cid__)});
            }
            else
            {
                for(auto i = it->second->begin(); i != it->second->end(); i++)
                {
                    to_cid__.push_back(*i);
                };
                this->cid[entry_c1.first].insert({it->first, make<int>::vec(to_cid__)});
            };
        };
    };
    // boundary id
    // make<make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>::map_str bid__
    make<std::pair<std::string, int>>::vec temp_pair;
    make<std::string>::vec boundary_reserve{"noslip", "inlet", "outlet", "temp", "fflux", "hamb"};
    bid__.insert({"fluid", make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int()});
    for(std::pair<std::string, make<int>::vec> entry_cell : this->cid["fluid"])
    {
        for(auto i = entry_cell.second.begin(); i != entry_cell.second.end(); i++)
        {
            make<int>::vec& faces_ = cells_[*i].cface;
            for(auto j = faces_.begin(); j != faces_.end(); j++)
            {
                temp_pair = is_boundary(*j, this->fid, boundary_reserve);
                if(temp_pair != make<std::pair<std::string, int>>::vec())
                {
                    make<int>::vec bid_int_keys = get_map_keys<int, make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>(bid__["fluid"]);
                    if(std::find(bid_int_keys.begin(), bid_int_keys.end(), int(*i)) == bid_int_keys.end())
                    {
                        bid__["fluid"].insert({*i, make<make<std::pair<std::string, int>>::vec>::map_int()});
                    };
                    bid__["fluid"][*i].insert({*j, temp_pair});
                };
            };
        };
    };
    make<std::string>::vec this_cid_str_keys = get_map_keys<std::string, make<make<make<int>::vec>::map_str>::map_str>(this->cid);
    if(std::find(this_cid_str_keys.begin(), this_cid_str_keys.end(), std::string("solid")) != this_cid_str_keys.end())
    {
        bid__.insert({"solid", make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int()});
        make<std::string>::vec boundary_reserve{"temp", "sflux", "hamb"};
        for(std::pair<std::string, make<int>::vec> entry_cell : this->cid["solid"])
        {
            for(auto i = entry_cell.second.begin(); i != entry_cell.second.end(); i++)
            {
                make<int>::vec& faces_ = cells_[*i].cface;
                for(auto j = faces_.begin(); j != faces_.end(); j++)
                {
                    temp_pair = is_boundary(*j, this->fid, boundary_reserve);
                    if(temp_pair != make<std::pair<std::string, int>>::vec())
                    {
                        make<int>::vec bid_int_keys = get_map_keys<int, make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>(bid__["solid"]);
                        if(std::find(bid_int_keys.begin(), bid_int_keys.end(), int(*i)) == bid_int_keys.end())
                        {
                            bid__["solid"].insert({*i, make<make<std::pair<std::string, int>>::vec>::map_int()});
                        };
                        bid__["solid"][*i].insert({*j, temp_pair});
                    };
                };
            };
        };
    };
    this->bid = bid__;
};
void minfo::make_template(make<coor>::map_int& nodes, make<finfo>::map_int& faces, make<cinfo>::map_int& cells)
{
    make<make<sparse_input>::vec>::map_str cc_fc_input;
    make<make<sparse_input>::vec>::map_str fc_input;
    make<make<sparse_input>::vec>::map_str cc_input;
    // ids
    make<make<make<int>::vec>::map_str>::map_str& cid_ = this->cid;
    // fluid
    cc_fc_input.insert({"fluid", make<sparse_input>::vec()});
    fc_input.insert({"fluid", make<sparse_input>::vec()});
    cc_input.insert({"fluid", make<sparse_input>::vec()});
    int fneigh;
    for(std::pair<std::string, make<int>::vec> entry_fluid : cid_["fluid"])
    {
        for(int i = 0; i < sizeof(entry_fluid.second) - 1; i++)
        {
            for(auto j = cells[entry_fluid.second[i]].cface.begin(); j != cells[entry_fluid.second[i]].cface.end(); j++)
            {
                fc_input["fluid"].push_back(sparse_input(entry_fluid.second[i], *j, 0.0));
            };
            cc_input["fluid"].push_back(sparse_input(entry_fluid.second[i], entry_fluid.second[i], 0.0));
            for(int j = i; j < sizeof(entry_fluid.second); j++)
            {
                fneigh = is_neighbor(cells[entry_fluid.second[i]].cface, cells[entry_fluid.second[j]].cface);
                if(fneigh >= 0)
                {
                    cc_fc_input["fluid"].push_back(sparse_input(entry_fluid.second[i], entry_fluid.second[j], fneigh));
                    cc_fc_input["fluid"].push_back(sparse_input(entry_fluid.second[j], entry_fluid.second[i], fneigh));
                    cc_input["fluid"].push_back(sparse_input(entry_fluid.second[i], entry_fluid.second[j], 0.0));
                    cc_input["fluid"].push_back(sparse_input(entry_fluid.second[j], entry_fluid.second[i], 0.0));
                };
            };
        };
        int tail = *entry_fluid.second.end();
        for(auto i = cells[tail].cface.begin(); i != cells[tail].cface.end(); i++)
        {
            fc_input["fluid"].push_back(sparse_input(tail, *i, 0.0));
        };
        cc_input["fluid"].push_back(sparse_input(tail, tail, 0.0));
    };
    // solid + conj
    make<std::string>::vec cid_str_keys = get_map_keys<std::string, make<make<make<int>::vec>::map_str>::map_str>(cid_);
    if(std::find(cid_str_keys.begin(), cid_str_keys.end(), std::string("solid")) != cid_str_keys.end())
    {
        // solid
        cc_fc_input.insert({"solid", make<sparse_input>::vec()});
        fc_input.insert({"solid", make<sparse_input>::vec()});
        cc_input.insert({"solid", make<sparse_input>::vec()});
        // conj
        cc_fc_input.insert({"conj", make<sparse_input>::vec()});
        fc_input.insert({"conj", make<sparse_input>::vec()});
        cc_input.insert({"conj", make<sparse_input>::vec()});
        // wall
        cid_.insert({"misc", make<make<int>::vec>::map_str()});
        cid_["misc"].insert({"wall", make<int>::vec()});
        // iter
        make<std::string>::vec done;
        for(std::pair<std::string, make<int>::vec> entry_solid1 : cid_["solid"])
        {
            done.push_back(entry_solid1.first);
            // self connection
            for(std::size_t i = 0; i < sizeof(entry_solid1.second) - 1; i++)
            {
                for(auto k = cells[entry_solid1.second[i]].cface.begin(); k != cells[entry_solid1.second[i]].cface.end(); k++)
                {
                    fc_input["solid"].push_back(sparse_input(entry_solid1.second[i], *k, 0.0));
                };
                for(std::size_t j = i; j < sizeof(entry_solid1.second); j++)
                {
                    fneigh = is_neighbor(cells[entry_solid1.second[i]].cface, cells[entry_solid1.second[j]].cface);
                    if(fneigh >= 0)
                    {
                        cc_fc_input["solid"].push_back(sparse_input(entry_solid1.second[i], entry_solid1.second[j], fneigh));
                        cc_fc_input["solid"].push_back(sparse_input(entry_solid1.second[j], entry_solid1.second[i], fneigh));
                        cc_input["solid"].push_back(sparse_input(entry_solid1.second[i], entry_solid1.second[j], 0.0));
                        cc_input["solid"].push_back(sparse_input(entry_solid1.second[j], entry_solid1.second[i], 0.0));
                    };
                };
                int tail = *entry_solid1.second.end();
                for(auto k = cells[tail].cface.begin(); k != cells[tail].cface.end(); k++)
                {
                    fc_input["solid"].push_back(sparse_input(tail, *k, 0.0));
                };
                cc_input["solid"].push_back(sparse_input(0.0, tail, tail));
            };
            // solid-solid connection
            for(std::pair<std::string, make<int>::vec> entry_solid2 : cid_["solid"])
            {
                if(std::find(done.begin(), done.end(), std::string(entry_solid2.first)) == done.end())
                {
                    for(std::size_t i = 0; i < sizeof(entry_solid1.second); i++)
                    {
                        for(std::size_t j = 0; j < sizeof(entry_solid2.second); j++)
                        {
                            fneigh = is_neighbor(cells[entry_solid1.second[i]].cface, cells[entry_solid2.second[j]].cface);
                            if(fneigh >= 0)
                            {
                                cc_fc_input["solid"].push_back(sparse_input(entry_solid1.second[i], entry_solid2.second[j], fneigh));
                                cc_fc_input["solid"].push_back(sparse_input(entry_solid2.second[j], entry_solid1.second[i], fneigh));
                                fc_input["solid"].push_back(sparse_input(entry_solid1.second[i], fneigh, 0.0));
                                fc_input["solid"].push_back(sparse_input(entry_solid2.second[j], fneigh, 0.0));
                                cc_input["solid"].push_back(sparse_input(entry_solid1.second[i], entry_solid2.second[j], 0.0));
                                cc_input["solid"].push_back(sparse_input(entry_solid2.second[j], entry_solid1.second[i], 0.0));
                            };
                        };
                    };
                };
            };
            // conj connection
            for(std::pair<std::string, make<int>::vec> entry_fluid : cid_["fluid"])
            {
                for(std::size_t i = 0; i < sizeof(entry_solid1.second); i++)
                {
                    for(std::size_t j = 0; j < sizeof(entry_fluid.second); j++)
                    {
                        fneigh = is_neighbor(cells[entry_solid1.second[i]].cface, cells[entry_fluid.second[j]].cface);
                        if(fneigh >= 0)
                        {
                            cc_fc_input["conj"].push_back(sparse_input(entry_solid1.second[i], entry_fluid.second[j], fneigh));
                            cc_fc_input["conj"].push_back(sparse_input(entry_fluid.second[j], entry_solid1.second[i], fneigh));
                            fc_input["conj"].push_back(sparse_input(entry_solid1.second[i], fneigh, 0.0));
                            fc_input["conj"].push_back(sparse_input(entry_fluid.second[j], fneigh, 0.0));
                            cc_input["conj"].push_back(sparse_input(entry_solid1.second[i], entry_fluid.second[j], 0.0));
                            cc_input["conj"].push_back(sparse_input(entry_fluid.second[j], entry_solid1.second[i], 0.0));
                            cid_["misc"]["wall"].push_back(entry_fluid.second[j]);
                        };
                    };
                };
            };
        };
    };
    // template
    make<make<double>::sp_mat>::map_str fc__;
    for(std::pair<std::string, make<sparse_input>::vec> entry : fc_input)
    {
        fc__.insert({entry.first, make<double>::sp_mat()});
        fc__[entry.first].setFromTriplets(entry.second.begin(), entry.second.end());
    };
    this->fc = fc__;
    make<make<double>::sp_mat>::map_str cc__;
    for(std::pair<std::string, make<sparse_input>::vec> entry : cc_input)
    {
        cc__.insert({entry.first, make<double>::sp_mat()});
        cc__[entry.first].setFromTriplets(entry.second.begin(), entry.second.end());
    };
    this->cc = cc__;
    make<make<int>::sp_mat>::map_str cc_fc__;
    for(std::pair<std::string, make<sparse_input>::vec> entry : cc_fc_input)
    {
        cc_fc__.insert({entry.first, make<int>::sp_mat()});
        cc_fc__[entry.first].setFromTriplets(entry.second.begin(), entry.second.end());
    };
    this->cc_fc = cc_fc__;
};
void minfo::make_size(make<coor>::map_int& nodes, make<finfo>::map_int& faces, make<cinfo>::map_int& cells)
{
    for(std::pair<int, finfo> entry : faces)
    {
        finfo& face_ = entry.second;
        switch(sizeof(face_.fnode))
        {
            case 3:
            make_tri(face_, nodes);
            case 4:
            make_quad(face_, nodes);
            default:
            return;
        };
    };
    for(std::pair<int, cinfo> entry : cells)
    {
        minfo& mesh_ = *this;
        cinfo& cell_ = entry.second;
        make_cell(mesh_, cell_, faces);
    };
    make<make<double>::map_int>::map_str size__;
    size__.insert({"area", make<double>::map_int()}); size__.insert({"volume", make<double>::map_int()});
    for(std::pair<int, finfo> entry : faces)
    {
        size__["area"].insert({entry.first, entry.second.farea});
    };
    for(std::pair<int, cinfo> entry : cells)
    {
        size__["volume"].insert({entry.first, entry.second.cvolume});
    };
    this->size = size__;
};
void minfo::make_geom(make<finfo>::map_int& faces, make<cinfo>::map_int& cells)
{
    // refs
    make<make<double>::sp_mat>::map_str& cc_ = this->cc;
    make<make<double>::sp_mat>::map_str& fc_ = this->fc;
    make<make<int>::sp_mat>::map_str& cc_fc_ = this->cc_fc;
    // cc-centered
    make<axes>::map_str eCF_; coor eCF__;
    make<axes>::map_str dCF_; coor dCF__;
    for(std::pair<std::string, make<double>::sp_mat> entry : cc_)
    {
        eCF_.insert({entry.first, axes()});
        dCF_.insert({entry.first, axes()});
        for(int i = 0; i < cc_[entry.first].outerSize(); i++)
        {
            for(make<double>::sp_mat::InnerIterator it(cc_[entry.first], i); it; ++it)
            {
                cinfo& row_current = cells[it.row()];
                cinfo& col_current = cells[it.col()];
                dCF__ = row_current.ccentroid - col_current.ccentroid;
                dCF_[entry.first](it.row(), it.col(), dCF__);
                dCF__.norm();
                eCF__ = dCF__;
                eCF_[entry.first](it.row(), it.col(), eCF__);
            };
        };
    };
    // fc-centered
    make<axes>::map_str Sf_; coor Sf__;
    make<axes>::map_str eCf_; coor eCf__;
    make<axes>::map_str dCf_; coor dCf__;
    for(std::pair<std::string, make<double>::sp_mat> entry : fc_)
    {
        Sf_.insert({entry.first, axes()});
        eCf_.insert({entry.first, axes()});
        dCf_.insert({entry.first, axes()});
        for(int i = 0; i < fc_[entry.first].outerSize(); i++)
        {
            for(make<double>::sp_mat::InnerIterator it(fc_[entry.first], i); it; ++it)
            {
                cinfo& row_current = cells[it.row()];
                finfo& col_current = faces[it.col()];
                dCf__ = col_current.fcentroid - row_current.ccentroid;
                eCf__(dCf__.norm());
                Sf__ = col_current.fnormal * col_current.farea;
                if(Sf__.dot(eCf__) < 0)
                {
                    Sf__ = Sf__ * (-1);
                };
                Sf_[entry.first](it.row(), it.col(), Sf__);
                eCf_[entry.first](it.row(), it.col(), eCf__);
                dCf_[entry.first](it.row(), it.col(), dCf__);
            };
        };
    };
    // using cc_fc
    make<axes>::map_str Ef_; coor Ef__;
    make<axes>::map_str Tf_; coor Tf__;
    for(std::pair<std::string, make<int>::sp_mat> entry : cc_fc_)
    {
        Ef_.insert({entry.first, axes()});
        Tf_.insert({entry.first, axes()});
        for(int i = 0; i < cc_fc_[entry.first].outerSize(); i++)
        {
            for(make<int>::sp_mat::InnerIterator it(cc_fc_[entry.first], i); it; ++it)
            {
                Sf__ = Sf_[entry.first].axes_to_coor(it.row(), it.value());
                eCF__ = eCF_[entry.first].axes_to_coor(it.row(), it.col());
                Ef__ = eCF__ * (Sf__.dot(Sf__) / eCF__.dot(Sf__));
                Tf__ = Sf__ - Ef__;
                Ef_[entry.first](it.row(), it.value(), Ef__);
                Tf_[entry.first](it.row(), it.value(), Tf__);
            };
        };
    };
    make<make<axes>::map_str>::map_str geom__;
    geom__.insert({"Sf", make<axes>::map_str(Sf_)}); geom__.insert({"Ef", make<axes>::map_str(Ef_)});
    geom__.insert({"Tf", make<axes>::map_str(Tf_)}); geom__.insert({"eCf", make<axes>::map_str(eCf_)});
    geom__.insert({"eCF", make<axes>::map_str(eCF_)}); geom__.insert({"dCf", make<axes>::map_str(dCf_)});
    geom__.insert({"dCF", make<axes>::map_str(dCF_)});
    this->geom = geom__;
};
void minfo::make_constants()
{
    // refs
    make<make<int>::sp_mat>::map_str& cc_fc_ = this->cc_fc;
    make<make<double>::sp_mat>::map_str& fc_ = this->fc;
    make<make<double>::sp_mat>::map_str& cc_ = this->cc;
    make<axes>::map_str& Sf_ = this->geom["Sf"]; coor Sf__;
    make<axes>::map_str& eCF_ = this->geom["eCF"]; coor eCF__;
    make<axes>::map_str& dCf_ = this->geom["dCf"]; coor dCf__;
    make<axes>::map_str& dCF_ = this->geom["dCF"]; double dCF__;
    // fcg_conv_aC = ( 1 - ( eCF.dot(dCf) / (2*sqrt_sum(dCF)) ) ), fc
    // fcg_conv_aF = (eCF.dot(dCf) / (2*sqrt_sum(dCF))), cc
    // fcg_diff_aC = (eCF.dot(Sf) / sqrt_sum(dCF))
    // fcg_diff_aF = - (eCF.dot(Sf) / sqrt_sum(dCF))
    // gc = sqrt_sum(dFf) / (sqrt_sum(dFf) + sqrt_sum(dCf)) , cc
    make<make<double>::sp_mat>::map_str g_conv_aC_(fc_);
    make<make<double>::sp_mat>::map_str g_conv_aF_(cc_);
    make<make<double>::sp_mat>::map_str g_diff_aC_(fc_);
    make<make<double>::sp_mat>::map_str g_diff_aF_(cc_);
    make<make<double>::sp_mat>::map_str gc_(cc_);
    double g_conv__; double g_diff__; double gc__;
    for(std::pair<std::string, make<int>::sp_mat> entry : cc_fc_)
    {
        for(int i = 0; i < cc_fc_[entry.first].outerSize(); i++)
        {
            for(make<int>::sp_mat::InnerIterator it(cc_fc_[entry.first], i); it; ++it)
            {
                eCF__ = eCF_[entry.first].axes_to_coor(it.row(), it.col());
                dCf__ = dCf_[entry.first].axes_to_coor(it.row(), it.value());
                dCF__ = dCF_[entry.first].axes_to_val(it.row(), it.col());
                g_conv__ = eCF__.dot(dCf__) / (2*dCF__);
                g_diff__ = eCF__.dot(Sf__) / dCF__;
                gc__ = dCf_[entry.first].axes_to_val(it.col(), it.value()) /
                       (dCf_[entry.first].axes_to_val(it.col(), it.value()) +
                       dCf_[entry.first].axes_to_val(it.row(), it.value()));
                g_conv_aC_[entry.first].coeffRef(it.row(), it.value()) = 1 - g_conv__;
                g_conv_aF_[entry.first].coeffRef(it.row(), it.col()) = g_conv__;
                g_diff_aC_[entry.first].coeffRef(it.row(), it.value()) = g_diff__;
                g_diff_aF_[entry.first].coeffRef(it.row(), it.col()) = (-1) * g_diff__;
                gc_[entry.first].coeffRef(it.row(), it.col()) = gc__;
            };
        };
    };
    make<make<make<double>::sp_mat>::map_str>::map_str constants__;
    constants__.insert({"g_conv_aC", make<make<double>::sp_mat>::map_str(g_conv_aC_)});
    constants__.insert({"g_conv_aF", make<make<double>::sp_mat>::map_str(g_conv_aF_)});
    constants__.insert({"g_diff_aC", make<make<double>::sp_mat>::map_str(g_diff_aC_)});
    constants__.insert({"g_diff_aC", make<make<double>::sp_mat>::map_str(g_diff_aC_)});
    constants__.insert({"gc", make<make<double>::sp_mat>::map_str(gc_)});
    this->constants = constants__;
};
void minfo::make_clust(make<finfo>::map_int& faces)
{
    make<make<int>::vec>::map_int& fs2s_ = this->fid["s2s"];
    make<sparse_input>::vec fcs2s_input;
    make<sparse_input>::vec ccs2s_input;
    for(int i = 0; i < sizeof(fs2s_) - 1; i++)
    {
        for(int j = i; j < sizeof(fs2s_); j++)
        {
            ccs2s_input.push_back(sparse_input(i, j, 0.0));
        };
    };
    int tail = sizeof(fs2s_) - 1;
    ccs2s_input.push_back(sparse_input(0, 0, 0.0));
    ccs2s_input.push_back(sparse_input(tail, tail, 0.0));
    for(std::pair<int, make<int>::vec> entry : fs2s_)
    {
        for(auto i = entry.second.begin(); i != entry.second.end(); i++)
        {
            fcs2s_input.push_back(sparse_input(entry.first, *i, 0.0));
        };
    };
    make<double>::sp_mat fc__; fc__.setFromTriplets(fcs2s_input.begin(), fcs2s_input.end());
    make<double>::sp_mat cc__; cc__.setFromTriplets(ccs2s_input.begin(), ccs2s_input.end());
    // constants g_s2s (view)
    make<double>::sp_mat view_(cc__);
    coor norm1__; coor norm2__; coor v1__; coor v2__; double cos1__; double cos2__; double sq_sum__;
    for(int i = 0; i < view_.outerSize(); i++)
    {
        for(make<double>::sp_mat::InnerIterator it(view_, i); it; ++it)
        {
            double view__ = 0.0;
            if(it.row() != it.col())
            {
                double area_tot = 0.0; 
                for(auto k = fs2s_[it.row()].begin(); k != fs2s_[it.row()].end(); k++)
                {
                    area_tot += faces[*k].farea;
                    for(auto l = fs2s_[it.col()].begin(); l != fs2s_[it.col()].end(); l++)
                    {
                        norm1__ = faces[*k].fnormal; norm2__ = faces[*l].fnormal;
                        v1__ = faces[*l].fcentroid - faces[*k].fcentroid;
                        v2__ = faces[*k].fcentroid - faces[*l].fcentroid;
                        if(norm1__.dot(v1__) < 0)
                        {
                            norm1__ = norm1__ * (-1);
                        };
                        if(norm2__.dot(v2__) < 0)
                        {
                            norm2__ = norm2__ * (-1);  
                        };
                        cos1__ = norm1__.dot(v1__) / sqrt_sum(v1__);
                        cos2__ = norm2__.dot(v2__) / sqrt_sum(v2__);
                        sq_sum__ = pow(v1__(0), 2) + pow(v1__(1), 2) + pow(v1__(2), 2);
                        view__ += (cos1__ * cos2__) / (4 * std::atan(1.0) * sq_sum__);
                    };
                };
                view_.coeffRef(it.row(), it.col()) = view__ / area_tot;
            };
        };
    };
    for(int i = 0; i < sizeof(fs2s_); i++)
    {
        cc__.coeffRef(i, i) = 1.0;
    };
    this->fc.insert({"s2s", make<double>::sp_mat(fc__)});
    this->cc.insert({"s2s", make<double>::sp_mat(cc__)});
    this->constants.insert({"s2s", make<make<double>::sp_mat>::map_str()});
    this->constants["s2s"].insert({"view", make<double>::sp_mat(view_)});
};
// pinfo
pinfo::pinfo(minfo& mesh_, user*& user_p)
{
    this->make_pinfo(mesh_, user_p);
};
void pinfo::make_pinfo(minfo& mesh_, user*& user_p)
{
    make<make<make<int>::vec>::map_str>::map_str& cid_ = mesh_.cid;
    make<make<double>::sp_mat>::map_str& cc_ = mesh_.cc;
    for(std::pair<std::string, make<make<int>::vec>::map_str> entry_cid : cid_)
    {
        std::string check = entry_cid.first;
        if(check.compare(std::string("fluid")) == 0)
        {
            make<double>::map_int rho_({});
            make<double>::map_int miu_({});
            make<double>::map_int cp_({});
            for(std::pair<std::string, make<int>::vec> entry_fluid : entry_cid.second)
            {
                for(auto i = entry_fluid.second.begin(); i != entry_fluid.second.end(); i++)
                {
                    rho_.insert({*i, calc_fluid_prop("rho", user_p->P_init, user_p->T_init, user_p->W_init)});
                    miu_.insert({*i, calc_fluid_prop("miu", user_p->P_init, user_p->T_init, user_p->W_init)});
                    cp_.insert({*i, calc_fluid_prop("cp", user_p->P_init, user_p->T_init, user_p->W_init)});
                };
            };
            this->rho["fluid"] = rho_;
            this->miu["fluid"] = miu_;
            this->cp["fluid"] = cp_;
        }
        else if(check.compare(std::string("solid")) == 0)
        {
            make<double>::map_int k_;
            for(std::pair<std::string, make<int>::vec> entry_solid : entry_cid.second)
            {
                for(auto i = entry_solid.second.begin(); i != entry_solid.second.end(); i++)
                {
                    k_.insert({*i, user_p->solid_k[entry_solid.first]});
                };
            };
            this->k["solid"] = k_;
        };
    };
    make<std::string>::vec cc_keys = get_map_keys<std::string, make<make<double>::sp_mat>::map_str>(cc_);
    if(std::find(cc_keys.begin(), cc_keys.end(), std::string("s2s")) != cc_keys.end())
    {
        make<double>::map_int eps_;
        make<double>::map_int alpha_;
        for(int i = 0; i < cc_["s2s"].rows(); i++)
        {
            eps_.insert({i, user_p->s2s_eps[i]});
        };
        this->eps["s2s"] = eps_;
    };
};
// winfo
winfo::winfo(minfo& mesh_)
{
    this->make_winfo(mesh_);
};
void winfo::make_winfo(minfo& mesh_)
{
    make<make<make<int>::vec>::map_str>::map_str& cid_ = mesh_.cid;
    make<std::string>::vec cid_str1_keys = get_map_keys<std::string, make<make<make<int>::vec>::map_str>::map_str>(cid_);
    if(std::find(cid_str1_keys.begin(), cid_str1_keys.end(), std::string("misc")) != cid_str1_keys.end())
    {
        make<double>::map_int ts_;
        make<double>::map_int utau_;
        make<double>::map_int miut_;
        for(auto i = cid_["misc"]["wall"].begin(); i != cid_["misc"]["wall"].end(); i++)
        {
            ts_.insert({*i, 0.0}); utau_.insert({*i, 0.0}); miut_.insert({*i, 0.0});
        };
        this->ts = ts_;
        this->utau = utau_;
        this->miut = miut_;
    };
};
//binfo
binfo::binfo(make<make<double>::map_int>::map_str& face_value_in, make<make<double>::map_str>::map_str& cell_value_in):
             face_value(face_value_in), cell_value(cell_value_in) {};
// vinfo
vinfo::vinfo(make<std::string>::vec which, minfo& mesh_, double init_value)
{
    this->make_vinfo(which, mesh_, init_value);
};
vinfo vinfo::operator()(vinfo& other)
{
    this->cvalue = other.cvalue;
    this->fvalue = other.fvalue;
    this->prev_cvalue = other.prev_cvalue;
    this->prev_fvalue = other.prev_fvalue;
    this->cgrad = other.cgrad;
    this->fgrad = other.fgrad;
    this->prev_cgrad = other.prev_cgrad;
    this->prev_fgrad = other.prev_fgrad;
    return *this;
};
void vinfo::make_vinfo(make<std::string>::vec which, minfo& mesh_, double init_value)
{
    make<make<double>::map_int>::map_str cvalue_;
    make<make<double>::map_int>::map_str fvalue_;
    make<make<double>::map_int>::map_str prev_cvalue_;
    make<make<double>::map_int>::map_str prev_fvalue_;
    make<make<coor>::map_int>::map_str cgrad_;
    make<make<coor>::map_int>::map_str fgrad_;
    make<make<coor>::map_int>::map_str prev_cgrad_;
    make<make<coor>::map_int>::map_str prev_fgrad_;
    for(auto i = which.begin(); i != which.end(); i++)
    {
        make<double>::sp_mat& cc_ = mesh_.cc[*i];
        cvalue_.insert({*i, make<double>::map_int()});
        cgrad_.insert({*i, make<coor>::map_int()});
        for(int j = 0; j < cc_.rows(); j++)
        {
            cvalue_[*i].insert({j, init_value});
            cgrad_[*i].insert({j, coor(0.0, 0.0, 0.0)});
        };
        make<double>::sp_mat& fc_ = mesh_.fc[*i];
        fvalue_.insert({*i, make<double>::map_int()});
        fgrad_.insert({*i, make<coor>::map_int()});
        for(int j = 0; j < fc_.cols(); j++)
        {
            fvalue_[*i].insert({j, init_value});
            fgrad_[*i].insert({j, coor(0.0, 0.0, 0.0)});
        };
    };
    this->cvalue = cvalue_; this->prev_cvalue = cvalue_;
    this->fvalue = fvalue_; this->prev_fvalue = fvalue_;
    this->cgrad = cgrad_; this->prev_cgrad = cgrad_;
    this->fgrad = fgrad_; this->prev_fgrad = fgrad_;
};

// other export func
template <class V>
std::string vec_to_str(V vec_in)
{
    std::string get("'");
    for(auto i = vec_in.begin(); i != vec_in.end(); i++)
    {
        if(i != vec_in.end())
        {
            get += std::to_string(*i) + ",";
        }
        else
        {
            get += std::to_string(*i);
        };
    };
    get += "'";
    return get;
};
// export
exports::exports(std::string outputname, std::string meshname, scheme*& scheme_ref)
{
    this->output_name = outputname;
    this->mesh_name = meshname;
    int n_cells = 0;
    make<double>::comp_str empty_err_res;
    empty_err_res.insert({"u", make<double>::vec()});
    empty_err_res.insert({"v", make<double>::vec()});
    empty_err_res.insert({"w", make<double>::vec()});
    empty_err_res.insert({"k", make<double>::vec()});
    empty_err_res.insert({"e", make<double>::vec()});
    empty_err_res.insert({"T", make<double>::vec()});
    make<double>::comp_str empty_converged(empty_err_res);
    empty_converged.insert({"P", make<double>::vec()});
    minfo*& mesh_ = scheme_ref->mesh;    
    for(std::pair<std::string, make<make<int>::vec>::map_str> entry1 : mesh_->cid)
    {
        for(std::pair<std::string, make<int>::vec> entry2 : entry1.second)
        {
            for(auto i = entry2.second.begin(); i != entry2.second.end(); i++)
            {
                n_cells += 1;
                this->cell_export.insert({*i, empty_converged});
                this->face_export.insert({*i, empty_converged});
            };
        };
    };
    this->number_of_cells = n_cells;
    this->err_export = make<double>::comp_str(empty_err_res);
    this->res_export = make<double>::comp_str(empty_err_res);
    this->time_export = make<long long>::vec();
    this->iter_export = make<int>::vec();
};
void exports::update_export(scphghe*& scphghe_ref, scheme*& scheme_ref, make<double>::map_str err,
                           make<double>::map_str res, long long time_iter, int outer_iter)
{
    momentum*& u_ = scphghe_ref->solv_u->eq;
    momentum*& v_ = scphghe_ref->solv_v->eq;
    momentum*& w_ = scphghe_ref->solv_w->eq;
    turb_k*& k_ = scphghe_ref->solv_k->eq;
    turb_e*& e_ = scphghe_ref->solv_e->eq;
    energy*& energy_ = scphghe_ref->solv_energy->eq;
    make<make<double>::comp_str>::map_int& cell_export_ = this->cell_export;
    make<make<double>::comp_str>::map_int& face_export_ = this->face_export;
    // converged cell values
    for(std::pair<int, double> entry : scheme_ref->pressure->cvalue["fluid"])
    {
        cell_export_[entry.first]["P"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : u_->value->cvalue["fluid"])
    {
        cell_export_[entry.first]["u"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : v_->value->cvalue["fluid"])
    {
        cell_export_[entry.first]["v"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : w_->value->cvalue["fluid"])
    {
        cell_export_[entry.first]["w"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : k_->value->cvalue["fluid"])
    {
        cell_export_[entry.first]["k"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : e_->value->cvalue["fluid"])
    {
        cell_export_[entry.first]["e"].push_back(entry.second);
    };
    for(std::pair<std::string, make<double>::map_int> entry_str : energy_->value->cvalue)
    {
        for(std::pair<int, double> entry_int : entry_str.second)
        {
            cell_export_[entry_int.first]["T"].push_back(entry_int.second);
        };
    };
    // converged face values
    for(std::pair<int, double> entry : scheme_ref->pressure->fvalue["fluid"])
    {
        face_export_[entry.first]["P"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : u_->value->fvalue["fluid"])
    {
        face_export_[entry.first]["u"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : v_->value->fvalue["fluid"])
    {
        face_export_[entry.first]["v"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : w_->value->fvalue["fluid"])
    {
        face_export_[entry.first]["w"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : k_->value->fvalue["fluid"])
    {
        face_export_[entry.first]["k"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : e_->value->fvalue["fluid"])
    {
        face_export_[entry.first]["e"].push_back(entry.second);
    };
    for(std::pair<std::string, make<double>::map_int> entry_str : energy_->value->fvalue)
    {
        for(std::pair<int, double> entry_int : entry_str.second)
        {
            face_export_[entry_int.first]["T"].push_back(entry_int.second);
        };
    };
    // BiCGSTAB estimated numerical errors and residuals
    for(std::pair<std::string, double> entry : err)
    {
        this->err_export[entry.first].push_back(entry.second);
    };
    for(std::pair<std::string, double> entry : res)
    {
        this->res_export[entry.first].push_back(entry.second);
    };
    // time and total iter
    this->time_export.push_back(time_iter);
    this->iter_export.push_back(outer_iter);
};
void exports::export_to_sql(scheme*& scheme_ref)
{
    sqlite3* database_p;
    int exit = 0;
    char* message_err;
    // open
    exit = sqlite3_open(this->output_name.c_str(), &database_p);
    // sql table templates
    // time step-wise values stored as string with comma limiter
    std::string commons_template =  "CREATE TABLE COMMONS("
                                    "NAME       TEXT            NOT NULL,"
                                    "N_CELLS    INT            NOT NULL);";
    exit = sqlite3_exec(database_p, commons_template.c_str(), NULL, 0, &message_err);
    std::string node_template = "CREATE TABLE NODE("
                                "ID         INT     PRIMARY KEY NOT NULL,"
                                "X          REAL               NOT NULL,"
                                "Y          REAL               NOT NULL,"
                                "Z          REAL               NOT NULL);";
    exit = sqlite3_exec(database_p, node_template.c_str(), NULL, 0, &message_err);
    std::string face_template = "CREATE TABLE FACE("
                                "ID         INT     PRIMARY KEY NOT NULL,"
                                "BOUNDARY   TEXT                NOT NULL,"
                                "NODES      TEXT                NOT NULL,"
                                "AREA       REAL                NOT NULL,"
                                "P          TEXT                NOT NULL,"
                                "U          TEXT                NOT NULL,"
                                "V          TEXT                NOT NULL,"
                                "W          TEXT                NOT NULL,"
                                "K          TEXT                NOT NULL,"
                                "E          TEXT                NOT NULL,"
                                "T          TEXT                NOT NULL);";
    exit = sqlite3_exec(database_p, face_template.c_str(), NULL, 0, &message_err);
    std::string cell_template = "CREATE TABLE CELL("
                                "ID         INT     PRIMARY KEY NOT NULL,"
                                "DOMAIN     TEXT                NOT NULL,"
                                "NODES      TEXT                NOT NULL,"
                                "VOLUME     REAL                NOT NULL,"
                                "P          TEXT                NOT NULL,"
                                "U          TEXT                NOT NULL,"
                                "V          TEXT                NOT NULL,"
                                "W          TEXT                NOT NULL,"
                                "K          TEXT                NOT NULL,"
                                "E          TEXT                NOT NULL,"
                                "T          TEXT                NOT NULL);";
    exit = sqlite3_exec(database_p, cell_template.c_str(), NULL, 0, &message_err);
    std::string err_template =  "CREATE TABLE ERROR("
                                "U          TEXT                NOT NULL,"
                                "V          TEXT                NOT NULL,"
                                "W          TEXT                NOT NULL,"
                                "K          TEXT                NOT NULL,"
                                "E          TEXT                NOT NULL,"
                                "T          TEXT                NOT NULL);";
    exit = sqlite3_exec(database_p, err_template.c_str(), NULL, 0, &message_err);
    std::string res_template =  "CREATE TABLE RESIDUAL("
                                "U          TEXT                NOT NULL,"
                                "V          TEXT                NOT NULL,"
                                "W          TEXT                NOT NULL,"
                                "K          TEXT                NOT NULL,"
                                "E          TEXT                NOT NULL,"
                                "T          TEXT                NOT NULL);";
    exit = sqlite3_exec(database_p, res_template.c_str(), NULL, 0, &message_err);
    std::string time_iter_template = "CREATE TABLE TIME_ITER("
                                     "TIME       TEXT                NOT NULL,"
                                     "ITER       TEXT                NOT NULL);";
    exit = sqlite3_exec(database_p, time_iter_template.c_str(), NULL, 0, &message_err);
    // commons
    std::string commons_input("INSERT INTO COMMONS VALUES(");
    commons_input += "'" + this->mesh_name + "'" + ", ";
    commons_input += std::to_string(this->number_of_cells) + ");";
    exit = sqlite3_exec(database_p, commons_input.c_str(), NULL, 0, &message_err);
    // nodes
    for(std::pair<int, coor> entry : scheme_ref->mesh->nodes)
    {
        std::string nodes_input("INSERT INTO NODES VALUES(");
        nodes_input += std::to_string(entry.first) + ", ";
        nodes_input += std::to_string(entry.second(0)) + ", ";
        nodes_input += std::to_string(entry.second(1)) + ", ";
        nodes_input += std::to_string(entry.second(2)) + ");";
        exit = sqlite3_exec(database_p, nodes_input.c_str(), NULL, 0, &message_err);
    };
    // face
    for(std::pair<int, finfo> entry_face : scheme_ref->mesh->faces)
    {
        std::string faces_input("INSERT INTO FACES VALUES(");
        faces_input += std::to_string(entry_face.first) + ", ";
        faces_input += "'" + entry_face.second.fboundary + "'" + ", ";
        faces_input += vec_to_str<make<int>::vec>(entry_face.second.fnode) + ", ";
        faces_input += std::to_string(entry_face.second.farea) + ", ";
        faces_input += vec_to_str<make<double>::vec>(this->face_export[entry_face.first]["P"]) + ", ";
        faces_input += vec_to_str<make<double>::vec>(this->face_export[entry_face.first]["u"]) + ", ";
        faces_input += vec_to_str<make<double>::vec>(this->face_export[entry_face.first]["v"]) + ", ";
        faces_input += vec_to_str<make<double>::vec>(this->face_export[entry_face.first]["w"]) + ", ";
        faces_input += vec_to_str<make<double>::vec>(this->face_export[entry_face.first]["k"]) + ", ";
        faces_input += vec_to_str<make<double>::vec>(this->face_export[entry_face.first]["e"]) + ", ";
        faces_input += vec_to_str<make<double>::vec>(this->face_export[entry_face.first]["T"]) + ");";
        exit = sqlite3_exec(database_p, faces_input.c_str(), NULL, 0, &message_err);
    };
    // cell
    for(std::pair<int, cinfo> entry_cell : scheme_ref->mesh->cells)
    {
        std::string cells_input("INSERT INTO CELLS VALUES(");
        cells_input += std::to_string(entry_cell.first) + ", ";
        cells_input += "'" + entry_cell.second.cdomain + "'" + ", ";
        cells_input += vec_to_str<make<int>::vec>(entry_cell.second.cnode) + ", ";
        cells_input += std::to_string(entry_cell.second.cvolume) + ", ";
        cells_input += vec_to_str<make<double>::vec>(this->cell_export[entry_cell.first]["P"]) + ", ";
        cells_input += vec_to_str<make<double>::vec>(this->cell_export[entry_cell.first]["u"]) + ", ";
        cells_input += vec_to_str<make<double>::vec>(this->cell_export[entry_cell.first]["v"]) + ", ";
        cells_input += vec_to_str<make<double>::vec>(this->cell_export[entry_cell.first]["w"]) + ", ";
        cells_input += vec_to_str<make<double>::vec>(this->cell_export[entry_cell.first]["k"]) + ", ";
        cells_input += vec_to_str<make<double>::vec>(this->cell_export[entry_cell.first]["e"]) + ", ";
        cells_input += vec_to_str<make<double>::vec>(this->cell_export[entry_cell.first]["T"]) + ");";
        exit = sqlite3_exec(database_p, cells_input.c_str(), NULL, 0, &message_err);
    };
    // err
    std::string err_input("INSERT INTO ERROR VALUES(");
    err_input += vec_to_str<make<double>::vec>(this->err_export["u"]) + ", ";
    err_input += vec_to_str<make<double>::vec>(this->err_export["v"]) + ", ";
    err_input += vec_to_str<make<double>::vec>(this->err_export["w"]) + ", ";
    err_input += vec_to_str<make<double>::vec>(this->err_export["k"]) + ", ";
    err_input += vec_to_str<make<double>::vec>(this->err_export["e"]) + ", ";
    err_input += vec_to_str<make<double>::vec>(this->err_export["T"]) + ");";
    exit = sqlite3_exec(database_p, err_input.c_str(), NULL, 0, &message_err);
    // res
    std::string res_input("INSERT INTO RESIDUAL VALUES(");
    res_input += vec_to_str<make<double>::vec>(this->res_export["u"]) + ", ";
    res_input += vec_to_str<make<double>::vec>(this->res_export["v"]) + ", ";
    res_input += vec_to_str<make<double>::vec>(this->res_export["w"]) + ", ";
    res_input += vec_to_str<make<double>::vec>(this->res_export["k"]) + ", ";
    res_input += vec_to_str<make<double>::vec>(this->res_export["e"]) + ", ";
    res_input += vec_to_str<make<double>::vec>(this->res_export["T"]) + ");";
    exit = sqlite3_exec(database_p, res_input.c_str(), NULL, 0, &message_err);
    // time_iter
    std::string time_iter_input("INSERT INTO TIME_ITER VALUES(");
    time_iter_input += vec_to_str<make<long long>::vec>(this->time_export) + ", ";
    time_iter_input += vec_to_str<make<int>::vec>(this->iter_export) + ");";
    exit = sqlite3_exec(database_p, time_iter_input.c_str(), NULL, 0, &message_err);
    // close
    sqlite3_close(database_p);
};

// other linear func
template <class V>
void append_template(make<std::string>::vec which, V* eq, scheme*& scheme_ref, double init_value)
{
    // refs
    minfo& mesh_ = *scheme_ref->mesh;
    make<make<double>::sp_mat>::map_str& fc_ = mesh_.fc;
    make<make<double>::sp_mat>::map_str& cc_ = mesh_.cc;
    // iter
    vinfo value_(which, mesh_, init_value);
    make<make<make<double>::map_int>::map_str>::map_str gamma_; // cell
    gamma_.insert({"cell", make<make<double>::map_int>::map_str()});
    gamma_.insert({"face", make<make<double>::map_int>::map_str()});
    make<make<double>::sp_mat>::map_str lhs_cc_;
    make<make<double>::sp_mat>::map_str lhs_fc_;
    make<make<double>::sp_mat>::map_str rhs_cc_;
    make<make<double>::sp_mat>::map_str rhs_fc_;
    for(auto i = which.begin(); i != which.end(); i++)
    {
        lhs_cc_.insert({*i, make<double>::sp_mat(cc_[*i])});
        rhs_cc_.insert({*i, make<double>::sp_mat(cc_[*i].rows(), 1)});
        std::string check(*i);
        if(check.compare("s2s") != 0)
        {
            gamma_["cell"].insert({*i, make<double>::map_int(scheme_ref->prop->rho["cell"])});
            gamma_["face"].insert({*i, make<double>::map_int(scheme_ref->prop->rho["face"])});
            lhs_fc_.insert({*i, make<double>::sp_mat(fc_[*i])});
            rhs_fc_.insert({*i, make<double>::sp_mat(fc_[*i])});
        };
    };
    vinfo* value_p = new vinfo(value_);
    eq->value = value_p; 
    eq->lhs_cc = lhs_cc_; eq->rhs_cc = rhs_cc_;
    if(std::find(which.begin(), which.end(), "s2s") == which.end())
    {
        eq->gamma = gamma_; eq->lhs_fc = lhs_fc_; eq->rhs_fc = rhs_fc_;
    };
};
// linear
// virtual defines
void linear::update_linear() {};
void linear::calc_wall() {};
void linear::make_linear() {};
void linear::calc_gamma() {};
void linear::calc_lhs() {};
void linear::calc_rhs() {};
void linear::calc_bound_lhs() {};
void linear::calc_bound_rhs() {};
// constructors
pcorrect::pcorrect(scheme*& scheme_ref, user*& user_ref) {this->make_linear(scheme_ref, user_ref);};
momentum::momentum(scheme*& scheme_ref, int axis_in) {this->make_linear(scheme_ref, axis_in);};
turb_k::turb_k(scheme*& scheme_ref) {this->make_linear(scheme_ref);};
turb_e::turb_e(scheme*& scheme_ref) {this->make_linear(scheme_ref);};
energy::energy(scheme*& scheme_ref, user*& user_ref) {this->make_linear(scheme_ref, user_ref);};
s2s::s2s(scheme*& scheme_ref) {this->make_linear(scheme_ref);};
// make_linear
void pcorrect::make_linear(scheme*& scheme_ref, user*& user_ref)
{
    append_template(make<std::string>::vec{"fluid"}, this, scheme_ref, user_ref->P_init);
};
void momentum::make_linear(scheme*& scheme_ref, int axis_in)
{
    this->axis = axis_in;
    append_template(make<std::string>::vec{"fluid"}, this, scheme_ref, 0.0);
};
void turb_k::make_linear(scheme*& scheme_ref)
{
    append_template(make<std::string>::vec{"fluid"}, this, scheme_ref, 0.0);
};
void turb_e::make_linear(scheme*& scheme_ref)
{
    append_template(make<std::string>::vec{"fluid"}, this, scheme_ref, 0.0);
};
void energy::make_linear(scheme*& scheme_ref, user*& user_ref)
{
    make<std::string>::vec which;
    for(std::pair<std::string, make<double>::sp_mat> entry : scheme_ref->mesh->cc)
    {
        if(entry.first.compare("s2s") != 0)
        {
            which.push_back(entry.first);
        };
    };
    append_template(which, this, scheme_ref, user_ref->T_init);
};
void s2s::make_linear(scheme*& scheme_ref)
{
    append_template(make<std::string>::vec{"s2s"}, this, scheme_ref, 0.0);
};
// update_linear
void pcorrect::update_linear(scheme*& scheme_ref, momentum*& u_ref, momentum*& v_ref, momentum*& w_ref)
{
    this->calc_lhs(scheme_ref, u_ref, v_ref, w_ref);
    this->calc_rhs(scheme_ref, u_ref, v_ref, w_ref);
    //solve
};
void momentum::update_linear(scheme*& scheme_ref, turb_k*& k_ref, momentum*& v1_ref, momentum*& v2_ref)
{
    this->calc_gamma(scheme_ref);
    this->calc_lhs(scheme_ref, v1_ref, v2_ref);
    this->calc_rhs(scheme_ref, k_ref, v1_ref, v2_ref);
    //solve
};
void turb_k::update_linear(scheme*& scheme_ref, turb_e*& e_ref, momentum*& u_ref,
                           momentum*& v_ref, momentum*& w_ref)
{
    this->calc_gamma(scheme_ref);
    this->calc_lhs(scheme_ref, u_ref, v_ref, w_ref);
    this->calc_rhs(scheme_ref, e_ref, u_ref, v_ref, w_ref);
    //solve
};
void turb_e::update_linear(scheme*& scheme_ref, turb_k*& k_ref, momentum*& u_ref, momentum*& v_ref, momentum*& w_ref)
{
    this->calc_gamma(scheme_ref);
    this->calc_lhs(scheme_ref, u_ref, v_ref, w_ref);
    this->calc_rhs(scheme_ref, k_ref, u_ref, v_ref, w_ref);
    //solve
};
void energy::update_linear(scheme*& scheme_ref, turb_k*& k_ref, momentum*& u_ref, momentum*& v_ref, momentum*& w_ref, double AH)
{
    this->calc_gamma(scheme_ref, AH);
    this->calc_lhs(scheme_ref, k_ref, u_ref, v_ref, w_ref, AH);
    this->calc_rhs(scheme_ref, u_ref, v_ref, w_ref, AH);
    //solve
};
void s2s::update_linear(scheme*& scheme_ref, energy*& energy_ref)
{
    this->calc_lhs(scheme_ref);
    this->calc_rhs(scheme_ref, energy_ref);
    //solve
};
// update_correcttion
void pcorrect::update_correction(scheme*& scheme_ref, momentum*& u_ref, momentum*& v_ref, momentum*& w_ref)
{
    make<int>::sp_mat& cc_fc_ = scheme_ref->mesh->cc_fc["fluid"];
    make<double>::map_int& rho_f_ = scheme_ref->prop->rho["face"];
    make<double>::map_int& rho_c_ = scheme_ref->prop->rho["cell"];
    make<double>::map_int& vol_ = scheme_ref->mesh->size["volume"];
    axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
    axes& dCF_ = scheme_ref->mesh->geom["dCF"]["fluid"];
    make<double>::sp_mat& u_aC_ = u_ref->lhs_cc["fluid"];
    make<double>::sp_mat& v_aC_ = v_ref->lhs_cc["fluid"];
    make<double>::sp_mat& w_aC_ = w_ref->lhs_cc["fluid"];
    make<coor>::map_int& pcorrect_cgrad_ = this->value->cgrad["fluid"];
    make<coor>::map_int& pcorrect_fgrad_ = this->value->fgrad["fluid"];
    double g_pcorrect__; coor pcorrect_cgrad__; coor pcorrect_fgrad__;
    coor Sf__; coor Ef__; coor Tf__; coor dCF__;
    double u_Df__; double v_Df__; double w_Df__;
    double u_DC__; double v_DC__; double w_DC__;
    int row;
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
        {
            u_Df__ = (vol_[it.row()]/u_aC_.coeffRef(it.row(), it.row()) + vol_[it.col()]/u_aC_.coeffRef(it.col(), it.col())) / 2;
            v_Df__ = (vol_[it.row()]/v_aC_.coeffRef(it.row(), it.row()) + vol_[it.col()]/v_aC_.coeffRef(it.col(), it.col())) / 2;
            w_Df__ = (vol_[it.row()]/w_aC_.coeffRef(it.row(), it.row()) + vol_[it.col()]/w_aC_.coeffRef(it.col(), it.col())) / 2;
            Sf__ = Sf_.axes_to_coor(it.row(), it.value());
            dCF__ = dCF_.axes_to_coor(it.row(), it.col());
            pcorrect_fgrad__ = pcorrect_fgrad_[it.value()];
            g_pcorrect__ = (-1) * rho_f_[it.value()] * pcorrect_fgrad__.dot(Sf__);
            u_ref->value->fvalue["fluid"][it.value()] += u_Df__ * g_pcorrect__;
            v_ref->value->fvalue["fluid"][it.value()] += v_Df__ * g_pcorrect__;
            w_ref->value->fvalue["fluid"][it.value()] += w_Df__ * g_pcorrect__;
            row = it.row();
        };
        pcorrect_cgrad__ = pcorrect_cgrad_[row];
        u_DC__ = vol_[row]/u_aC_.coeffRef(row, row);
        v_DC__ = vol_[row]/v_aC_.coeffRef(row, row);
        w_DC__ = vol_[row]/v_aC_.coeffRef(row, row);
        u_ref->value->cvalue["fluid"][row] += (-1) * rho_c_[row] * u_DC__ * pcorrect_cgrad__(0);
        v_ref->value->cvalue["fluid"][row] += (-1) * rho_c_[row] * v_DC__ * pcorrect_cgrad__(1);
        w_ref->value->cvalue["fluid"][row] += (-1) * rho_c_[row] * w_DC__ * pcorrect_cgrad__(2);
        scheme_ref->pressure->cvalue["fluid"][row] += this->value->cvalue["fluid"][row];
    };
};
void s2s::update_source_s2s(scheme*& scheme_ref)
{
    for(std::pair<int, double> entry : this->value->cvalue["s2s"])
    {
        scheme_ref->source->face_value["s2s"][entry.first] = entry.second;
    };
};
//coef calc
//pcorrect
void pcorrect::calc_lhs(scheme*& scheme_ref, momentum*& u_ref, momentum*& v_ref,
                        momentum*& w_ref)
{
    make<double>::sp_mat& fc_ = this->lhs_fc["fluid"];
    make<double>::sp_mat& cc_ = this->lhs_cc["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref->mesh->cc_fc["fluid"];
    make<double>::map_int& rho_f_ = scheme_ref->prop->rho["face"];
    make<double>::map_int& vol_ = scheme_ref->mesh->size["volume"];
    axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
    axes& dCF_ = scheme_ref->mesh->geom["dCF"]["fluid"];
    make<double>::sp_mat& u_aC_ = u_ref->lhs_cc["fluid"];
    make<double>::sp_mat& v_aC_ = v_ref->lhs_cc["fluid"];
    make<double>::sp_mat& w_aC_ = w_ref->lhs_cc["fluid"];
    coor Sf__; coor Ef__; coor Tf__; coor dCF__;
    double u_Df__; double v_Df__; double w_Df__;
    double Dauf__; double DauC__ = 0.0;
    int row; double aC;
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        aC = 0.0;
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
        {
            u_Df__ = (vol_[it.row()]/u_aC_.coeffRef(it.row(), it.row()) + vol_[it.col()]/u_aC_.coeffRef(it.col(), it.col())) / 2;
            v_Df__ = (vol_[it.row()]/v_aC_.coeffRef(it.row(), it.row()) + vol_[it.col()]/v_aC_.coeffRef(it.col(), it.col())) / 2;
            w_Df__ = (vol_[it.row()]/w_aC_.coeffRef(it.row(), it.row()) + vol_[it.col()]/w_aC_.coeffRef(it.col(), it.col())) / 2;
            Sf__ = Sf_.axes_to_coor(it.row(), it.value());
            dCF__ = dCF_.axes_to_coor(it.row(), it.col());
            Dauf__ = (pow(u_Df__*Sf__(0) ,2) + pow(v_Df__*Sf__(1),2) + pow(w_Df__*Sf__(2),2)) /
                     (dCF__(0)*u_Df__*Sf__(0) + dCF__(1)*v_Df__*Sf__(1) + dCF__(2)*w_Df__*Sf__(2));
            DauC__ += (-1) * Dauf__;
            fc_.coeffRef(it.row(), it.value()) = rho_f_[it.value()]*Dauf__;
            cc_.coeffRef(it.row(), it.col()) = (-1)*rho_f_[it.value()]*Dauf__;
            row = it.row();
        };
        // boundary
        make<int>::vec bid_int_keys = get_map_keys<int, make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>(scheme_ref->mesh->bid["fluid"]);
        make<std::pair<std::string, int>>::vec fc_bound;
        if(std::find(bid_int_keys.begin(), bid_int_keys.end(), row) != bid_int_keys.end())
        {
            for(std::pair<int, make<std::pair<std::string, int>>::vec> entry : scheme_ref->mesh->bid["fluid"][row])
            {
                fc_bound = scheme_ref->mesh->bid["fluid"][row][entry.first];
                for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                {
                    this->calc_bound_lhs(row, entry.first, j->first, j->second, scheme_ref, DauC__);
                };
            };
        };
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; ++it)
        {
            row = it.row();
            aC += fc_.coeffRef(it.row(), it.col());
        };
        cc_.coeffRef(row, row) = aC;
    };
};
void pcorrect::calc_rhs(scheme*& scheme_ref, momentum*& u_ref, momentum*& v_ref,
                        momentum*& w_ref)
{
    make<double>::sp_mat& fc_ = this->rhs_fc["fluid"];
    make<double>::sp_mat& cc_ = this->rhs_cc["fluid"];
    make<double>::sp_mat& gc_ = scheme_ref->mesh->constants["gc"]["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref->mesh->cc_fc["fluid"];
    make<double>::map_int& rho_f_ = scheme_ref->prop->rho["face"];
    make<double>::map_int& vol_ = scheme_ref->mesh->size["volume"];
    axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
    axes& Ef_ = scheme_ref->mesh->geom["Ef"]["fluid"];
    axes& Tf_ = scheme_ref->mesh->geom["Tf"]["fluid"];
    axes& dCF_ = scheme_ref->mesh->geom["dCF"]["fluid"];
    make<double>::sp_mat& u_aC_ = u_ref->lhs_cc["fluid"];
    make<double>::sp_mat& v_aC_ = v_ref->lhs_cc["fluid"];
    make<double>::sp_mat& w_aC_ = w_ref->lhs_cc["fluid"];
    coor mom_f__; coor mom_C__; coor mom_F__;
    double gc__;
    coor prev_pgradf__; coor prev_pgradf_itr__; coor prev_pgradf_min__;
    coor Sf__; coor Ef__; coor Tf__; coor dCF__;
    double u_Df__; double v_Df__; double w_Df__;
    double Dauf__;
    int row; double bC;
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        bC = 0.0;
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
        {
            gc__ = gc_.coeffRef(it.row(), it.value());
            prev_pgradf__ = scheme_ref->pressure->prev_fgrad["fluid"][it.value()];
            prev_pgradf_itr__ = gc__ * scheme_ref->pressure->prev_cgrad["fluid"][it.row()] + (1 - gc__) * scheme_ref->pressure->prev_cgrad["fluid"][it.col()];
            prev_pgradf_min__ = prev_pgradf__ - prev_pgradf_itr__;
            u_Df__ = (vol_[it.row()]/u_aC_.coeffRef(it.row(), it.row()) + vol_[it.col()]/u_aC_.coeffRef(it.col(), it.col())) / 2;
            v_Df__ = (vol_[it.row()]/v_aC_.coeffRef(it.row(), it.row()) + vol_[it.col()]/v_aC_.coeffRef(it.col(), it.col())) / 2;
            w_Df__ = (vol_[it.row()]/w_aC_.coeffRef(it.row(), it.row()) + vol_[it.col()]/w_aC_.coeffRef(it.col(), it.col())) / 2;
            mom_C__(0) = u_ref->value->cvalue["fluid"][it.row()]; mom_C__(1) = v_ref->value->cvalue["fluid"][it.row()]; mom_C__(2) = w_ref->value->cvalue["fluid"][it.row()];
            mom_F__(0) = u_ref->value->cvalue["fluid"][it.col()]; mom_F__(1) = v_ref->value->cvalue["fluid"][it.col()]; mom_F__(2) = w_ref->value->cvalue["fluid"][it.col()];
            mom_f__ = gc__* (mom_C__) + (1 - gc__) * (mom_F__);
            Sf__ = Sf_.axes_to_coor(it.row(), it.value());
            Ef__ = Ef_.axes_to_coor(it.row(), it.value());
            Tf__ = Tf_.axes_to_coor(it.row(), it.value());
            dCF__ = dCF_.axes_to_coor(it.row(), it.col());
            Dauf__ = (pow(u_Df__*Sf__(0) ,2) + pow(v_Df__*Sf__(1),2) + pow(w_Df__*Sf__(2),2)) /
                     (dCF__(0)*u_Df__*Sf__(0) + dCF__(1)*v_Df__*Sf__(1) + dCF__(2)*w_Df__*Sf__(2));
            fc_.coeffRef(it.row(), it.value()) = (rho_f_[it.value()] * (mom_f__.dot(Ef__) + mom_f__.dot(Tf__)) -
                                                 Dauf__ * (prev_pgradf_min__.dot(Ef__) + prev_pgradf_min__.dot(Tf__)));
        };
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; ++it)
        {
            row = it.row();
            bC += fc_.coeffRef(it.row(), it.col());
        };
        cc_.coeffRef(row, 0) = bC;
    };
};
void pcorrect::calc_bound_lhs(int row, int col, std::string check, int unique, scheme*& scheme_ref, double DauC__)
{
    if(check.compare("ns") == 0)
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
        axes& Tf_ = scheme_ref->mesh->geom["Tf"]["fluid"];
        double pc__ = scheme_ref->pressure->cvalue["fluid"][row];
        coor prev_pgradc__ = scheme_ref->pressure->prev_cgrad["fluid"][row];
        coor prev_pgradf__ = scheme_ref->pressure->prev_fgrad["fluid"][col];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor Tf__ = Tf_.axes_to_coor(row, col);
        scheme_ref->pressure->fvalue["fluid"][col] = pc__ + (prev_pgradc__.dot(Sf__) - prev_pgradf__.dot(Tf__)) / DauC__;
    }
    else if(check.compare("in") == 0)
    {
        this->lhs_fc["fluid"].coeffRef(row, col) = scheme_ref->prop->rho["face"][col] * DauC__; 
    }
    else if(check.compare("out") == 0)
    {
        this->lhs_fc["fluid"].coeffRef(row, col) = scheme_ref->prop->rho["face"][col] * DauC__; 
    }
    else
    {
        return;
    };
};
//momentum
void momentum::calc_gamma(scheme*& scheme_ref)
{
    make<double>::map_int& gamma_c_ = this->gamma["fluid"]["cell"];
    make<double>::map_int& gamma_f_ = this->gamma["fluid"]["face"];
    make<double>::sp_mat& gc_ = scheme_ref->mesh->constants["gc"]["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref->mesh->cc_fc["fluid"];
    make<double>::map_int& miu_c_ = scheme_ref->prop->miu["cell"];
    make<double>::map_int& miut_ = scheme_ref->wall->miut;
    for(int i = 0; i < sizeof(gamma_c_); i++)
    {
        gamma_c_[i] = miu_c_[i] + miut_[i];
    };
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
        {
            gamma_f_[it.value()] = gamma_c_[it.row()] * gamma_c_[it.col()] / (gamma_c_[it.col()] *
                                   gc_.coeffRef(it.row(), it.col()) + gamma_c_[it.row()] * (1 -
                                   gc_.coeffRef(it.row(), it.col())));
        };
    };
};
void momentum::calc_wall(scheme*& scheme_ref, momentum*& v1_ref, momentum*& v2_ref)
{
    make<int>::vec& wall_cell_ = scheme_ref->mesh->cid["misc"]["wall"];
    make<double>::map_int& valuec_ = this->value->cvalue["fluid"];
    make<coor>::map_int& gradc_ = this->value->cgrad["fluid"];
    coor gradc__; double valuec__;
    coor v__; coor v_norm__; coor v_cross__;
    coor wall_perp__; coor wall_parallel__;
    double theta__;
    double d_plus__; double d_perp__; double u_tau__; double viu__;
    for(auto i = wall_cell_.begin(); i != wall_cell_.end(); i++)
    {
        v__(this->axis) = valuec_[*i];
        v__(v1_ref->axis) = v1_ref->value->cvalue["fluid"][*i];
        v__(v2_ref->axis) = v2_ref->value->cvalue["fluid"][*i];
        v_norm__ = v__ / v__.norm();
        wall_perp__ = scheme_ref->wall->wall_parallel[*i];
        v_cross__ = v_norm__.cross(wall_perp__);
        wall_parallel__ = wall_perp__.cross(v_cross__);
        wall_parallel__ = wall_parallel__ / wall_parallel__.norm();
        theta__ = std::acos(v__.dot(wall_parallel__) / v__.norm());
        if(theta__ > 90)
        {
            wall_parallel__ = (-1) * wall_parallel__;
        };
        gradc__ = gradc_[*i];
        valuec__ = valuec_[*i];
        u_tau__ = scheme_ref->wall->utau[*i];
        viu__ = scheme_ref->prop->miu["cell"][*i] / scheme_ref->prop->rho["cell"][*i];
        d_perp__ = pow(gradc__.squaredNorm() + 2*valuec__, 0.5) - gradc__.norm();
        d_plus__ = d_perp__ * u_tau__ / viu__;
        valuec_[*i] = ((std::log(d_plus__) / 0.41) + 5.25) * wall_parallel__(this->axis);
    };
};
void momentum::calc_lhs(scheme*& scheme_ref, momentum*& v1_ref, momentum*& v2_ref)
{
    make<double>::sp_mat& fc_ = this->lhs_fc["fluid"];
    make<double>::sp_mat& cc_ = this->lhs_cc["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref->mesh->cc_fc["fluid"];
    make<double>::sp_mat& fcg_conv_aC_ = scheme_ref->mesh->constants["g_conv_aC"]["fluid"]; 
    make<double>::sp_mat& fcg_conv_aF_ = scheme_ref->mesh->constants["g_conv_aF"]["fluid"];
    make<double>::sp_mat& fcg_diff_aC_ = scheme_ref->mesh->constants["g_diff_aC"]["fluid"];
    make<double>::sp_mat& fcg_diff_aF_ = scheme_ref->mesh->constants["g_diff_aF"]["fluid"];
    make<double>::map_int& rho_f_ = scheme_ref->prop->rho["face"];
    make<double>::sp_mat& rho_v_sf_ = scheme_ref->rho_v_sf;
    make<double>::map_int& gamma_f_ = this->gamma["fluid"]["face"];
    axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
    axes& Ef_ = scheme_ref->mesh->geom["Ef"]["fluid"];
    axes& Tf_ = scheme_ref->mesh->geom["Tf"]["fluid"];
    axes& eCF_ = scheme_ref->mesh->geom["eCF"]["fluid"];
    axes& eCf_ = scheme_ref->mesh->geom["eCf"]["fluid"];
    axes& dCF_ = scheme_ref->mesh->geom["dCF"]["fluid"];
    coor Sf__; coor Ef__; coor Tf__; coor eCF__; coor eCf__;
    double dCFval__;
    int row; double aC;
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        aC = 0.0;
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
        {
            fc_.coeffRef(it.row(), it.value()) = fcg_conv_aC_.coeffRef(it.row(), it.value()) *
                                                 rho_v_sf_.coeffRef(it.row(), it.value()) + 
                                                 fcg_diff_aC_.coeffRef(it.row(), it.value()) *
                                                 gamma_f_[it.value()];
            cc_.coeffRef(it.row(), it.col()) = fcg_conv_aF_.coeffRef(it.row(), it.value()) *
                                               rho_v_sf_.coeffRef(it.row(), it.value()) + 
                                               fcg_diff_aF_.coeffRef(it.row(), it.value()) *
                                               gamma_f_[it.value()];
            row = it.row();
        };
        // boundary
        make<int>::vec bid_int_keys = get_map_keys<int, make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>(scheme_ref->mesh->bid["fluid"]);
        make<std::pair<std::string, int>>::vec fc_bound;
        if(std::find(bid_int_keys.begin(), bid_int_keys.end(), row) != bid_int_keys.end())
        {
            for(std::pair<int, make<std::pair<std::string, int>>::vec> entry : scheme_ref->mesh->bid["fluid"][row])
            {
                fc_bound = scheme_ref->mesh->bid["fluid"][row][entry.first];
                for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                {
                    this->calc_bound_lhs(row, entry.first, j->first, j->second, scheme_ref, v1_ref, v2_ref);
                };
            };
        };
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; ++it)
        {
            row = it.row();
            aC += fc_.coeffRef(it.row(), it.col());
        };
        cc_.coeffRef(row, row) = aC;
    };
};
void momentum::calc_rhs(scheme*& scheme_ref, turb_k*& k_ref, momentum*& v1_ref,
                        momentum*& v2_ref)
{
    make<double>::sp_mat& fc_ = this->rhs_fc["fluid"];
    make<double>::sp_mat& cc_ = this->rhs_cc["fluid"];
    make<double>::sp_mat& gc_ = scheme_ref->mesh->constants["gc"]["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref->mesh->cc_fc["fluid"];
    make<double>::sp_mat& fcg_conv_aC_ = scheme_ref->mesh->constants["g_conv_aC"]["fluid"]; 
    make<double>::sp_mat& fcg_conv_aF_ = scheme_ref->mesh->constants["g_conv_aF"]["fluid"];
    make<double>::sp_mat& fcg_diff_aC_ = scheme_ref->mesh->constants["g_diff_aC"]["fluid"];
    make<double>::sp_mat& fcg_diff_aF_ = scheme_ref->mesh->constants["g_diff_aF"]["fluid"];
    make<double>::map_int& rho_f_ = scheme_ref->prop->rho["face"];
    make<double>::map_int& rho_c_ = scheme_ref->prop->rho["cell"];
    make<double>::sp_mat& rho_v_sf_ = scheme_ref->rho_v_sf;
    make<double>::map_int& gamma_f_ = this->gamma["face"]["fluid"];
    make<double>::map_int& vol_ = scheme_ref->mesh->size["volume"];
    double gc__;
    axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
    axes& Ef_ = scheme_ref->mesh->geom["Ef"]["fluid"];
    axes& Tf_ = scheme_ref->mesh->geom["Tf"]["fluid"];
    axes& eCF_ = scheme_ref->mesh->geom["eCF"]["fluid"];
    axes& eCf_ = scheme_ref->mesh->geom["eCf"]["fluid"];
    axes& dCf_ = scheme_ref->mesh->geom["dCf"]["fluid"];
    coor Sf__; coor Ef__; coor Tf__; coor eCF__; coor eCf__; coor dCf__;
    double Pf__; double kf__;
    coor gradf_itr__; coor gradc__;
    int row; double bC; coor bf_conv; coor bf_diff; double bf_body;
    make<std::pair<std::string, int>>::vec fc_bound;
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        bC = 0.0;
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
        {
            gc__ = gc_.coeffRef(it.row(), it.value());
            Pf__ = scheme_ref->pressure->fvalue["fluid"][it.value()];
            kf__ = k_ref->value->fvalue["fluid"][it.value()];
            gradc__ = this->value->cgrad["fluid"][it.row()];
            gradf_itr__ = gc__ * gradc__ + (1 - gc__) * this->value->cgrad["fluid"][it.col()];
            Ef__ = Ef_.axes_to_coor(it.row(), it.value());
            Tf__ = Tf_.axes_to_coor(it.row(), it.value());
            eCF__ = eCF_.axes_to_coor(it.row(), it.col());
            dCf__ = dCf_.axes_to_coor(it.row(), it.value());
            bf_conv = (gradf_itr__.dot(eCF__) * eCF__ - (gradc__ + gradf_itr__)) / 2;
            bf_diff = gradf_itr__ - (gradf_itr__.dot(eCF__) * eCF__);
            bf_body = Pf__ + (2 * rho_f_[it.value()] * kf__ / 3) * Sf__(this->axis);
            fc_.coeffRef(it.row(), it.value()) = (bf_conv.dot(dCf__) * rho_v_sf_.coeffRef(it.row(), it.value())) +
                                                 (bf_diff.dot(Ef__) * gamma_f_[it.value()]) + 
                                                 (bf_diff.dot(Tf__) * gamma_f_[it.value()]) + bf_body;
            row = it.row();
        };
        // boundary
        make<int>::vec bid_int_keys = get_map_keys<int, make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>(scheme_ref->mesh->bid["fluid"]);
        make<std::pair<std::string, int>>::vec fc_bound;
        if(std::find(bid_int_keys.begin(), bid_int_keys.end(), row) != bid_int_keys.end())
        {
            for(std::pair<int, make<std::pair<std::string, int>>::vec> entry : scheme_ref->mesh->bid["fluid"][row])
            {
                fc_bound = scheme_ref->mesh->bid["fluid"][row][entry.first];
                for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                {
                    this->calc_bound_lhs(row, entry.first, j->first, j->second, scheme_ref, v1_ref, v2_ref);
                };
            };
        };
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; ++it)
        {
            row = it.row();
            bC += fc_.coeffRef(it.row(), it.col());
            if(this->axis == 1)
            {
                bC += rho_c_[it.row()] * 9.81 * vol_[it.row()];
            };
        };
        cc_.coeffRef(row, row) = bC;
    };
};
void momentum::calc_bound_lhs(int row, int col, std::string check, int unique, scheme*& scheme_ref,
                              momentum*& v1_ref, momentum*& v2_ref)
{
    if(check.compare("ns") == 0)
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref->mesh->geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref->mesh->geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        double Sfval__ = scheme_ref->mesh->size["area"][col];
        double miu_f__ = scheme_ref->prop->miu["face"][col];
        double d_perp__ = dCf__.dot(Sf__) / Sfval__;
        this->lhs_fc["fluid"].coeffRef(row, col) = (miu_f__ * Sfval__ * (1 - pow(eCf__(this->axis), 2)))/ d_perp__;
    }
    else if(check.compare("out") == 0)
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref->mesh->geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref->mesh->geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        coor vf__;
        coor v0_gradc__; coor v1_gradc__; coor v2_gradc__;
        coor v0_gradf__; coor v1_gradf__; coor v2_gradf__;
        double v0_f__; double v1_f__; double v2_f__;
        v0_gradc__ = this->value->cgrad["fluid"][row];
        v1_gradc__ = v1_ref->value->cgrad["fluid"][row];
        v2_gradc__ = v2_ref->value->cgrad["fluid"][row];
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v1_gradf__ = v1_gradc__ - (v1_gradc__.dot(eCf__) * eCf__);
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v0_f__ = this->value->cvalue["fluid"][row] + (v0_gradf__.dot(dCf__));
        v1_f__ = v1_ref->value->cvalue["fluid"][row] + (v1_gradf__.dot(dCf__));
        v2_f__ = v2_ref->value->cvalue["fluid"][row] + (v2_gradf__.dot(dCf__));
        vf__(this->axis) = v0_f__; vf__(v1_ref->axis) = v1_f__; vf__(v2_ref->axis) = v2_f__; 
        this->lhs_fc["fluid"].coeffRef(row, col) = scheme_ref->prop->rho["face"][col] * (vf__.dot(Sf__)); 
    }
    else
    {
        return;
    };
};
void momentum::calc_bound_rhs(int row, int col, std::string check, int unique, scheme*& scheme_ref,
                              momentum*& v1_ref, momentum*& v2_ref)
{
    if(check.compare("ns") == 0)
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref->mesh->geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref->mesh->geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        double Sfval__ = scheme_ref->mesh->size["area"][col];
        double miu_f__ = scheme_ref->prop->miu["face"][col];
        double d_perp__ = dCf__.dot(Sf__) / Sfval__;
        this->rhs_fc["fluid"].coeffRef(row, col) = (miu_f__ * Sfval__ / d_perp__) * (
                                          this->value->fvalue["fluid"][col] * (1 - pow(eCf__(this->axis), 2)) +
                                          (v1_ref->value->cvalue["fluid"][row] - v1_ref->value->fvalue["fluid"][col]) * eCf__(v1_ref->axis) * eCf__(this->axis) +
                                          (v2_ref->value->cvalue["fluid"][row] - v2_ref->value->fvalue["fluid"][col]) * eCf__(v2_ref->axis) * eCf__(this->axis)
                                          ) - (scheme_ref->pressure->fvalue["fluid"][col] * Sf__(this->axis));
    }
    else if(check.compare("in") == 0)
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref->mesh->geom["eCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = (-1) * eCf_.axes_to_coor(row, col);
        coor rho_v_sf__ = scheme_ref->source->face_value[check][unique] * eCf__;
        this->rhs_fc["fluid"].coeffRef(row, col) = scheme_ref->prop->rho["face"][col] * rho_v_sf__.dot(Sf__); 
    }
    else if(check.compare("out") == 0)
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref->mesh->geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref->mesh->geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        coor vf__;
        coor v0_gradc__; coor v1_gradc__; coor v2_gradc__;
        coor v0_gradf__; coor v1_gradf__; coor v2_gradf__;
        double v0_f__; double v1_f__; double v2_f__;
        v0_gradc__ = this->value->cgrad["fluid"][row];
        v1_gradc__ = v1_ref->value->cgrad["fluid"][row];
        v2_gradc__ = v2_ref->value->cgrad["fluid"][row];
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v1_gradf__ = v1_gradc__ - (v1_gradc__.dot(eCf__) * eCf__);
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v0_f__ = this->value->cvalue["fluid"][row] + (v0_gradf__.dot(dCf__));
        v1_f__ = v1_ref->value->cvalue["fluid"][row] + (v1_gradf__.dot(dCf__));
        v2_f__ = v2_ref->value->cvalue["fluid"][row] + (v2_gradf__.dot(dCf__));
        vf__(this->axis) = v0_f__; vf__(v1_ref->axis) = v1_f__; vf__(v2_ref->axis) = v2_f__; 
        this->rhs_fc["fluid"].coeffRef(row, col) = (-1) * scheme_ref->prop->rho["face"][col] * (vf__.dot(Sf__)) * (v0_gradc__.dot(dCf__)) -
                                          (scheme_ref->pressure->fvalue["fluid"][col] * Sf__(this->axis));
    }
    else
    {
        return;
    };
};
//turb_k
void turb_k::calc_gamma(scheme*& scheme_ref)
{
    make<double>::map_int& gamma_c_ = this->gamma["cell"]["fluid"];
    make<double>::map_int& gamma_f_ = this->gamma["face"]["fluid"];
    make<double>::sp_mat& gc_ = scheme_ref->mesh->constants["gc"]["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref->mesh->cc_fc["fluid"];
    make<double>::map_int& miu_c_ = scheme_ref->prop->miu["cell"];
    make<double>::map_int& miut_ = scheme_ref->wall->miut;
    for(int i = 0; i < sizeof(gamma_c_); i++)
    {
        gamma_c_[i] = miu_c_[i] + miut_[i];
    };
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
        {
            gamma_f_[it.value()] = gamma_c_[it.row()] * gamma_c_[it.col()] / (gamma_c_[it.col()] *
                                   gc_.coeffRef(it.row(), it.col()) + gamma_c_[it.row()] * (1 -
                                   gc_.coeffRef(it.row(), it.col())));
        };
    };
};
void turb_k::calc_wall(scheme*& scheme_ref, turb_e*& e_ref)
{
    make<int>::vec& wall_cell_ = scheme_ref->mesh->cid["misc"]["wall"];
    make<double>::map_int& valuec_ = this->value->cvalue["fluid"];
    double Re_t__; double cmiu__;
    for(auto i = wall_cell_.begin(); i != wall_cell_.end(); i++)
    {
        Re_t__ = scheme_ref->prop->rho["cell"][*i] * pow(this->value->cvalue["fluid"][*i], 2) / (scheme_ref->prop->miu["cell"][*i] * e_ref->value->cvalue["fluid"][*i]);
        cmiu__ = 0.09 * std::exp(-3.4 / pow(1 + (Re_t__/50), 2));
        valuec_[*i] = 1 / pow(cmiu__, 0.5);
    };
};
void turb_k::calc_lhs(scheme*& scheme_ref, momentum*& u_ref, momentum*& v_ref, momentum*& w_ref)
{
    make<double>::sp_mat& fc_ = this->lhs_fc["fluid"];
    make<double>::sp_mat& cc_ = this->lhs_cc["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref->mesh->cc_fc["fluid"];
    make<double>::sp_mat& fcg_conv_aC_ = scheme_ref->mesh->constants["g_conv_aC"]["fluid"]; 
    make<double>::sp_mat& fcg_conv_aF_ = scheme_ref->mesh->constants["g_conv_aF"]["fluid"];
    make<double>::sp_mat& fcg_diff_aC_ = scheme_ref->mesh->constants["g_diff_aC"]["fluid"];
    make<double>::sp_mat& fcg_diff_aF_ = scheme_ref->mesh->constants["g_diff_aF"]["fluid"];
    make<double>::map_int& rho_f_ = scheme_ref->prop->rho["face"];
    make<double>::sp_mat& rho_v_sf_ = scheme_ref->rho_v_sf;
    make<double>::map_int& gamma_f_ = this->gamma["face"]["fluid"];
    axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
    axes& Ef_ = scheme_ref->mesh->geom["Ef"]["fluid"];
    axes& Tf_ = scheme_ref->mesh->geom["Tf"]["fluid"];
    axes& eCF_ = scheme_ref->mesh->geom["eCF"]["fluid"];
    axes& eCf_ = scheme_ref->mesh->geom["eCf"]["fluid"];
    axes& dCF_ = scheme_ref->mesh->geom["dCF"]["fluid"];
    coor Sf__; coor Ef__; coor Tf__; coor eCF__; coor eCf__;
    double dCFval__;
    int row; double aC;
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        aC = 0.0;
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
        {
            fc_.coeffRef(it.row(), it.value()) = fcg_conv_aC_.coeffRef(it.row(), it.value()) *
                                                 rho_v_sf_.coeffRef(it.row(), it.value()) + 
                                                 fcg_diff_aC_.coeffRef(it.row(), it.value()) *
                                                 gamma_f_[it.value()];
            cc_.coeffRef(it.row(), it.col()) = fcg_conv_aF_.coeffRef(it.row(), it.value()) *
                                               rho_v_sf_.coeffRef(it.row(), it.value()) + 
                                               fcg_diff_aF_.coeffRef(it.row(), it.value()) *
                                               gamma_f_[it.value()];
            row = it.row();
        };
        // boundary
        make<int>::vec bid_int_keys = get_map_keys<int, make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>(scheme_ref->mesh->bid["fluid"]);
        make<std::pair<std::string, int>>::vec fc_bound;
        if(std::find(bid_int_keys.begin(), bid_int_keys.end(), row) != bid_int_keys.end())
        {
            for(std::pair<int, make<std::pair<std::string, int>>::vec> entry : scheme_ref->mesh->bid["fluid"][row])
            {
                fc_bound = scheme_ref->mesh->bid["fluid"][row][entry.first];
                for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                {
                    this->calc_bound_lhs(row, entry.first, j->first, j->second, scheme_ref, u_ref, v_ref, w_ref);
                };
            };
        };
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; ++it)
        {
            row = it.row();
            aC += fc_.coeffRef(it.row(), it.col());
        };
        cc_.coeffRef(row, row) = aC;
    };
};
void turb_k::calc_rhs(scheme*& scheme_ref, turb_e*& e_ref, momentum*& u_ref,
                      momentum*& v_ref, momentum*& w_ref)
{
    make<double>::sp_mat& fc_ = this->lhs_fc["fluid"];
    make<double>::sp_mat& cc_ = this->lhs_cc["fluid"];
    make<double>::sp_mat& gc_ = scheme_ref->mesh->constants["gc"]["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref->mesh->cc_fc["fluid"];
    make<double>::sp_mat& fcg_conv_aC_ = scheme_ref->mesh->constants["g_conv_aC"]["fluid"]; 
    make<double>::sp_mat& fcg_conv_aF_ = scheme_ref->mesh->constants["g_conv_aF"]["fluid"];
    make<double>::sp_mat& fcg_diff_aC_ = scheme_ref->mesh->constants["g_diff_aC"]["fluid"];
    make<double>::sp_mat& fcg_diff_aF_ = scheme_ref->mesh->constants["g_diff_aF"]["fluid"];
    make<double>::map_int& rho_f_ = scheme_ref->prop->rho["face"];
    make<double>::map_int& rho_c_ = scheme_ref->prop->rho["cell"];
    make<double>::sp_mat& rho_v_sf_ = scheme_ref->rho_v_sf;
    make<double>::map_int& gamma_f_ = this->gamma["face"]["fluid"];
    make<double>::map_int& vol_ = scheme_ref->mesh->size["volume"];
    double gc__;
    axes& Sf_ = scheme_ref->mesh->geom["Sf"]["face"];
    axes& Ef_ = scheme_ref->mesh->geom["Ef"]["face"];
    axes& Tf_ = scheme_ref->mesh->geom["Tf"]["face"];
    axes& eCF_ = scheme_ref->mesh->geom["eCF"]["face"];
    axes& eCf_ = scheme_ref->mesh->geom["eCf"]["face"];
    axes& dCf_ = scheme_ref->mesh->geom["dCf"]["face"];
    coor Sf__; coor Ef__; coor Tf__; coor eCF__; coor eCf__; coor dCf__;
    double Pf__; double kf__;
    coor gradf_itr__; coor gradc__;
    int row; double bC; coor bf_conv; coor bf_diff;
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        bC = 0.0;
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
        {
            gc__ = gc_.coeffRef(it.row(), it.value());
            Pf__ = scheme_ref->pressure->fvalue["fluid"][it.value()];
            gradc__ = this->value->cgrad["fluid"][it.row()];
            gradf_itr__ = gc__ * gradc__ + (1 - gc__) * this->value->cgrad["fluid"][it.col()];
            Ef__ = Ef_.axes_to_coor(it.row(), it.value());
            Tf__ = Tf_.axes_to_coor(it.row(), it.value());
            eCF__ = eCF_.axes_to_coor(it.row(), it.col());
            dCf__ = dCf_.axes_to_coor(it.row(), it.value());
            bf_conv = (gradf_itr__.dot(eCF__) * eCF__ - (gradc__ + gradf_itr__)) / 2;
            bf_diff = gradf_itr__ - (gradf_itr__.dot(eCF__) * eCF__);
            fc_.coeffRef(it.row(), it.value()) = (bf_conv.dot(dCf__) * rho_v_sf_.coeffRef(it.row(), it.value())) +
                                                 (bf_diff.dot(Ef__) * gamma["face"]["fluid"][it.value()]) + 
                                                 (bf_diff.dot(Tf__) * gamma["face"]["fluid"][it.value()]);
        row = it.row();
        };
        // boundary
        make<int>::vec bid_int_keys = get_map_keys<int, make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>(scheme_ref->mesh->bid["fluid"]);
        make<std::pair<std::string, int>>::vec fc_bound;
        if(std::find(bid_int_keys.begin(), bid_int_keys.end(), row) != bid_int_keys.end())
        {
            for(std::pair<int, make<std::pair<std::string, int>>::vec> entry : scheme_ref->mesh->bid["fluid"][row])
            {
                fc_bound = scheme_ref->mesh->bid["fluid"][row][entry.first];
                for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                {
                    this->calc_bound_lhs(row, entry.first, j->first, j->second, scheme_ref, u_ref, v_ref, w_ref);
                };
            };
        };
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; ++it)
        {
            row = it.row();
            bC += fc_.coeffRef(it.row(), it.col());
        };
        double miut__ = scheme_ref->wall->miut[row];
        double phi_v__ = scheme_ref->phi_v[row];
        double rho_c__ = rho_c_[row];
        double e_c__ = e_ref->value->cvalue["fluid"][row];
        double vol__ = vol_[row];
        cc_.coeffRef(row, row) = bC + (((miut__ * phi_v__) - (rho_c__ * e_c__)) * vol__);
    };
};
void turb_k::calc_bound_lhs(int row, int col, std::string check, int unique, scheme*& scheme_ref,
                            momentum*& u_ref, momentum*& v_ref, momentum*& w_ref)
{
    if(check.compare("in") == 0)
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
        axes& dCf_ = scheme_ref->mesh->geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        this->lhs_fc["fluid"].coeffRef(row, col) = this->gamma["face"]["fluid"][col] * Sf__.norm() / dCf__.norm(); 
    }
    else if(check.compare("out") == 0)
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref->mesh->geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref->mesh->geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        // momentum
        coor vf__;
        coor v0_gradc__; coor v1_gradc__; coor v2_gradc__;
        coor v0_gradf__; coor v1_gradf__; coor v2_gradf__;
        double v0_f__; double v1_f__; double v2_f__;
        v0_gradc__ = u_ref->value->cgrad["fluid"][row];
        v1_gradc__ = v_ref->value->cgrad["fluid"][row];
        v2_gradc__ = w_ref->value->cgrad["fluid"][row];
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v1_gradf__ = v1_gradc__ - (v1_gradc__.dot(eCf__) * eCf__);
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v0_f__ = u_ref->value->cvalue["fluid"][row] + (v0_gradf__.dot(dCf__));
        v1_f__ = v_ref->value->cvalue["fluid"][row] + (v1_gradf__.dot(dCf__));
        v2_f__ = w_ref->value->cvalue["fluid"][row] + (v2_gradf__.dot(dCf__));
        vf__(u_ref->axis) = v0_f__; vf__(v_ref->axis) = v1_f__; vf__(w_ref->axis) = v2_f__; 
        this->lhs_fc["fluid"].coeffRef(row, col) = scheme_ref->prop->rho["face"][col] * (vf__.dot(Sf__)); 
    }
    else
    {
        return;
    };
}; 
void turb_k::calc_bound_rhs(int row, int col, std::string check, int unique, scheme*& scheme_ref,
                            momentum*& u_ref, momentum*& v_ref, momentum*& w_ref)
{
    if(check.compare("in") == 0)
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
        axes& dCf_ = scheme_ref->mesh->geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        coor vf__ = (-1) * scheme_ref->source->face_value[check][unique] * Sf__;
        double inlet_k = 1/2 * 0.01 * (vf__.dot(vf__));
        this->rhs_fc["fluid"].coeffRef(row, col) = (-1) * this->gamma["face"]["fluid"][col] * Sf__.norm() * inlet_k / dCf__.norm();
    }
    else if(check.compare("out") == 0)
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref->mesh->geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref->mesh->geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        // momentum
        coor vf__;
        coor v0_gradc__; coor v1_gradc__; coor v2_gradc__;
        coor v0_gradf__; coor v1_gradf__; coor v2_gradf__;
        double v0_f__; double v1_f__; double v2_f__;
        v0_gradc__ = u_ref->value->cgrad["fluid"][row];
        v1_gradc__ = v_ref->value->cgrad["fluid"][row];
        v2_gradc__ = w_ref->value->cgrad["fluid"][row];
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v1_gradf__ = v1_gradc__ - (v1_gradc__.dot(eCf__) * eCf__);
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v0_f__ = u_ref->value->cvalue["fluid"][row] + (v0_gradf__.dot(dCf__));
        v1_f__ = v_ref->value->cvalue["fluid"][row] + (v1_gradf__.dot(dCf__));
        v2_f__ = w_ref->value->cvalue["fluid"][row] + (v2_gradf__.dot(dCf__));
        vf__(u_ref->axis) = v0_f__; vf__(v_ref->axis) = v1_f__; vf__(w_ref->axis) = v2_f__;
        // turb_k
        coor k_gradc__; coor k_gradf__;
        k_gradc__ = this->value->cgrad["fluid"][row];
        k_gradf__ = k_gradc__ - (k_gradc__.dot(eCf__) * eCf__);
        this->rhs_fc["fluid"].coeffRef(row, col) = (-1) * scheme_ref->prop->rho["face"][col] * (vf__.dot(Sf__)) * (k_gradf__.dot(dCf__)); 
    }
    else
    {
        return;
    };
};
//turb_e
void turb_e::calc_gamma(scheme*& scheme_ref)
{
    make<double>::map_int& gamma_c_ = this->gamma["cell"]["fluid"];
    make<double>::map_int& gamma_f_ = this->gamma["face"]["fluid"];
    make<double>::sp_mat& gc_ = scheme_ref->mesh->constants["gc"]["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref->mesh->cc_fc["fluid"];
    make<double>::map_int& miu_c_ = scheme_ref->prop->miu["cell"];
    make<double>::map_int& miut_ = scheme_ref->wall->miut;
    for(int i = 0; i < sizeof(gamma_c_); i++)
    {
        gamma_c_[i] = miu_c_[i] + miut_[i] / 1.3;
    };
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
        {
            gamma_f_[it.value()] = gamma_c_[it.row()] * gamma_c_[it.col()] / (gamma_c_[it.col()] *
                                   gc_.coeffRef(it.row(), it.col()) + gamma_c_[it.row()] * (1 -
                                   gc_.coeffRef(it.row(), it.col())));
        };
    };
};
void turb_e::calc_wall(scheme*& scheme_ref)
{
    make<int>::vec& wall_cell_ = scheme_ref->mesh->cid["misc"]["wall"];
    make<double>::map_int& valuec_ = this->value->cvalue["fluid"];
    make<coor>::map_int& gradc_ = this->value->cgrad["fluid"];
    coor gradc__; double valuec__;
    double d_perp__; double u_tau__; double viu__;
    for(auto i = wall_cell_.begin(); i != wall_cell_.end(); i++)
    {
        gradc__ = gradc_[*i];
        valuec__ = valuec_[*i];
        u_tau__ = scheme_ref->wall->utau[*i];
        viu__ = scheme_ref->prop->miu["cell"][*i] / scheme_ref->prop->rho["cell"][*i];
        d_perp__ = pow(gradc__.squaredNorm() + 2*valuec__, 0.5) - gradc__.norm();
        valuec_[*i] = viu__ / (u_tau__ * 0.41 * d_perp__);
    };
};
void turb_e::calc_lhs(scheme*& scheme_ref, momentum*& u_ref, momentum*& v_ref, momentum*& w_ref)
{
    make<double>::sp_mat& fc_ = this->lhs_fc["fluid"];
    make<double>::sp_mat& cc_ = this->lhs_cc["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref->mesh->cc_fc["fluid"];
    make<double>::sp_mat& fcg_conv_aC_ = scheme_ref->mesh->constants["g_conv_aC"]["fluid"]; 
    make<double>::sp_mat& fcg_conv_aF_ = scheme_ref->mesh->constants["g_conv_aF"]["fluid"];
    make<double>::sp_mat& fcg_diff_aC_ = scheme_ref->mesh->constants["g_diff_aC"]["fluid"];
    make<double>::sp_mat& fcg_diff_aF_ = scheme_ref->mesh->constants["g_diff_aF"]["fluid"];
    make<double>::map_int& rho_f_ = scheme_ref->prop->rho["face"];
    make<double>::sp_mat& rho_v_sf_ = scheme_ref->rho_v_sf;
    make<double>::map_int& gamma_f_ = this->gamma["face"]["fluid"];
    axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
    axes& Ef_ = scheme_ref->mesh->geom["Ef"]["fluid"];
    axes& Tf_ = scheme_ref->mesh->geom["Tf"]["fluid"];
    axes& eCF_ = scheme_ref->mesh->geom["eCF"]["fluid"];
    axes& eCf_ = scheme_ref->mesh->geom["eCf"]["fluid"];
    axes& dCF_ = scheme_ref->mesh->geom["dCF"]["fluid"];
    coor Sf__; coor Ef__; coor Tf__; coor eCF__; coor eCf__;
    double dCFval__;
    int row; double aC;
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        aC = 0.0;
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
        {
            fc_.coeffRef(it.row(), it.value()) = fcg_conv_aC_.coeffRef(it.row(), it.value()) *
                                                 rho_v_sf_.coeffRef(it.row(), it.value()) + 
                                                 fcg_diff_aC_.coeffRef(it.row(), it.value()) *
                                                 gamma_f_[it.value()];
            cc_.coeffRef(it.row(), it.col()) = fcg_conv_aF_.coeffRef(it.row(), it.value()) *
                                               rho_v_sf_.coeffRef(it.row(), it.value()) + 
                                               fcg_diff_aF_.coeffRef(it.row(), it.value()) *
                                               gamma_f_[it.value()];
            row = it.row();
        };
        // boundary
        make<int>::vec bid_int_keys = get_map_keys<int, make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>(scheme_ref->mesh->bid["fluid"]);
        make<std::pair<std::string, int>>::vec fc_bound;
        if(std::find(bid_int_keys.begin(), bid_int_keys.end(), row) != bid_int_keys.end())
        {
            for(std::pair<int, make<std::pair<std::string, int>>::vec> entry : scheme_ref->mesh->bid["fluid"][row])
            {
                fc_bound = scheme_ref->mesh->bid["fluid"][row][entry.first];
                for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                {
                    this->calc_bound_lhs(row, entry.first, j->first, j->second, scheme_ref, u_ref, v_ref, w_ref);
                };
            };
        };
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; ++it)
        {
            row = it.row();
            aC += fc_.coeffRef(it.row(), it.col());
        };
        cc_.coeffRef(row, row) = aC;
    };
};
void turb_e::calc_rhs(scheme*& scheme_ref, turb_k*& k_ref, momentum*& u_ref, momentum*& v_ref, momentum*& w_ref)
{
    make<double>::sp_mat& fc_ = this->lhs_fc["fluid"];
    make<double>::sp_mat& cc_ = this->lhs_cc["fluid"];
    make<double>::sp_mat& gc_ = scheme_ref->mesh->constants["gc"]["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref->mesh->cc_fc["fluid"];
    make<double>::sp_mat& fcg_conv_aC_ = scheme_ref->mesh->constants["g_conv_aC"]["fluid"]; 
    make<double>::sp_mat& fcg_conv_aF_ = scheme_ref->mesh->constants["g_conv_aF"]["fluid"];
    make<double>::sp_mat& fcg_diff_aC_ = scheme_ref->mesh->constants["g_diff_aC"]["fluid"];
    make<double>::sp_mat& fcg_diff_aF_ = scheme_ref->mesh->constants["g_diff_aC"]["fluid"];
    make<double>::map_int& rho_f_ = scheme_ref->prop->rho["face"];
    make<double>::map_int& rho_c_ = scheme_ref->prop->rho["cell"];
    make<double>::sp_mat& rho_v_sf_ = scheme_ref->rho_v_sf;
    make<double>::map_int& gamma_f_ = this->gamma["face"]["fluid"];
    make<double>::map_int& vol_ = scheme_ref->mesh->size["volume"];
    double gc__;
    axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
    axes& Ef_ = scheme_ref->mesh->geom["Ef"]["fluid"];
    axes& Tf_ = scheme_ref->mesh->geom["Tf"]["fluid"];
    axes& eCF_ = scheme_ref->mesh->geom["eCF"]["fluid"];
    axes& eCf_ = scheme_ref->mesh->geom["eCf"]["fluid"];
    axes& dCf_ = scheme_ref->mesh->geom["dCf"]["fluid"];
    coor Sf__; coor Ef__; coor Tf__; coor eCF__; coor eCf__; coor dCf__;
    double Pf__; double kf__;
    coor gradf_itr__; coor gradc__;
    int row; double bC; coor bf_conv; coor bf_diff;
    double ce_1__ = 1.44; double ce_2__; double Re_t__;
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        bC = 0.0;
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
        {
            gc__ = gc_.coeffRef(it.row(), it.value());
            Pf__ = scheme_ref->pressure->fvalue["fluid"][it.value()];
            gradc__ = this->value->cgrad["fluid"][it.row()];
            gradf_itr__ = gc__ * gradc__ + (1 - gc__) * this->value->cgrad["fluid"][it.col()];
            Ef__ = Ef_.axes_to_coor(it.row(), it.value());
            Tf__ = Tf_.axes_to_coor(it.row(), it.value());
            eCF__ = eCF_.axes_to_coor(it.row(), it.col());
            dCf__ = dCf_.axes_to_coor(it.row(), it.value());
            bf_conv = (gradf_itr__.dot(eCF__) * eCF__ - (gradc__ + gradf_itr__)) / 2;
            bf_diff = gradf_itr__ - (gradf_itr__.dot(eCF__) * eCF__);
            fc_.coeffRef(it.row(), it.value()) = (bf_conv.dot(dCf__) * rho_v_sf_.coeffRef(it.row(), it.value())) +
                                                 (bf_diff.dot(Ef__) * gamma_f_[it.value()]) + 
                                                 (bf_diff.dot(Tf__) * gamma_f_[it.value()]);
            row = it.row();
        };
        // boundary
        make<int>::vec bid_int_keys = get_map_keys<int, make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>(scheme_ref->mesh->bid["fluid"]);
        make<std::pair<std::string, int>>::vec fc_bound;
        if(std::find(bid_int_keys.begin(), bid_int_keys.end(), row) != bid_int_keys.end())
        {
            for(std::pair<int, make<std::pair<std::string, int>>::vec> entry : scheme_ref->mesh->bid["fluid"][row])
            {
                fc_bound = scheme_ref->mesh->bid["fluid"][row][entry.first];
                for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                {
                    this->calc_bound_lhs(row, entry.first, j->first, j->second, scheme_ref, u_ref, v_ref, w_ref);
                };
            };
        };
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; ++it)
        {
            row = it.row();
            bC += fc_.coeffRef(it.row(), it.col());
        };
        double ts__ = scheme_ref->wall->ts[row];
        double miut__ = scheme_ref->wall->miut[row];
        double phi_v__ = scheme_ref->phi_v[row];
        double rho_c__ = rho_c_[row];
        double e_c__ = this->value->cvalue["fluid"][row];
        double vol__ = vol_[row];
        Re_t__ = scheme_ref->prop->rho["cell"][row] * pow(k_ref->value->cvalue["fluid"][row], 2) / (scheme_ref->prop->miu["cell"][row] * this->value->cvalue["fluid"][row]);
        ce_2__ = 1.92*(1 - 0.3*std::exp((-1) * pow(Re_t__, 2)));
        cc_.coeffRef(row, row) = bC + ((ce_1__ * miut__ * phi_v__ / ts__) - (ce_2__ * scheme_ref->prop->rho["cell"][row] *
                                 this->value->cvalue["fluid"][row] / ts__))*vol__;
    };
};
void turb_e::calc_bound_lhs(int row, int col, std::string check, int unique, scheme*& scheme_ref,
                            momentum*& u_ref, momentum*& v_ref, momentum*& w_ref)
{
    if(check.compare("in") == 0)
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
        axes& dCf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        this->lhs_fc["fluid"].coeffRef(row, col) = this->gamma["face"]["fluid"][col] * Sf__.norm() / dCf__.norm(); 
    }
    else if(check.compare("out") == 0)
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref->mesh->geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref->mesh->geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        // momentum
        coor vf__;
        coor v0_gradc__; coor v1_gradc__; coor v2_gradc__;
        coor v0_gradf__; coor v1_gradf__; coor v2_gradf__;
        double v0_f__; double v1_f__; double v2_f__;
        v0_gradc__ = u_ref->value->cgrad["fluid"][row];
        v1_gradc__ = v_ref->value->cgrad["fluid"][row];
        v2_gradc__ = w_ref->value->cgrad["fluid"][row];
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v1_gradf__ = v1_gradc__ - (v1_gradc__.dot(eCf__) * eCf__);
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v0_f__ = u_ref->value->cvalue["fluid"][row] + (v0_gradf__.dot(dCf__));
        v1_f__ = v_ref->value->cvalue["fluid"][row] + (v1_gradf__.dot(dCf__));
        v2_f__ = w_ref->value->cvalue["fluid"][row] + (v2_gradf__.dot(dCf__));
        vf__(u_ref->axis) = v0_f__; vf__(v_ref->axis) = v1_f__; vf__(w_ref->axis) = v2_f__; 
        this->lhs_fc["fluid"].coeffRef(row, col) = scheme_ref->prop->rho["face"][col] * (vf__.dot(Sf__)); 
    }
    else
    {
        return;
    };
}; 
void turb_e::calc_bound_rhs(int row, int col, std::string check, int unique, scheme*& scheme_ref,
                            momentum*& u_ref, momentum*& v_ref, momentum*& w_ref)
{
    if(check.compare("in") == 0)
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
        axes& dCf_ = scheme_ref->mesh->geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        coor vf__ = (-1) * scheme_ref->source->face_value[check][unique] * Sf__;
        make<double>::map_int& rho_f_ = scheme_ref->prop->rho["face"];
        make<double>::map_int& miu_f_ = scheme_ref->prop->miu["face"];
        make<double>::map_int& miut_ = scheme_ref->wall->miut;
        double rho_f__ = rho_f_[col]; double miut__ = miut_[row];
        double inlet_k = 1/2 * 0.01 * (vf__.dot(vf__));
        double inlet_e = 0.09 * rho_f__ * pow(inlet_k, 2) / miut__;
        this->rhs_fc["fluid"].coeffRef(row, col) = (-1) * this->gamma["face"]["fluid"][col] * Sf__.norm() * inlet_e / dCf__.norm(); 
    }
    else if(check.compare("out") == 0)
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref->mesh->geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref->mesh->geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        // momentum
        coor vf__;
        coor v0_gradc__; coor v1_gradc__; coor v2_gradc__;
        coor v0_gradf__; coor v1_gradf__; coor v2_gradf__;
        double v0_f__; double v1_f__; double v2_f__;
        v0_gradc__ = u_ref->value->cgrad["fluid"][row];
        v1_gradc__ = v_ref->value->cgrad["fluid"][row];
        v2_gradc__ = w_ref->value->cgrad["fluid"][row];
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v1_gradf__ = v1_gradc__ - (v1_gradc__.dot(eCf__) * eCf__);
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v0_f__ = u_ref->value->cvalue["fluid"][row] + (v0_gradf__.dot(dCf__));
        v1_f__ = v_ref->value->cvalue["fluid"][row] + (v1_gradf__.dot(dCf__));
        v2_f__ = w_ref->value->cvalue["fluid"][row] + (v2_gradf__.dot(dCf__));
        vf__(u_ref->axis) = v0_f__; vf__(v_ref->axis) = v1_f__; vf__(w_ref->axis) = v2_f__;
        // turb_e
        coor e_gradc__; coor e_gradf__;
        e_gradc__ = this->value->cgrad["fluid"][row];
        e_gradf__ = e_gradc__ - (e_gradc__.dot(eCf__) * eCf__);
        this->rhs_fc["fluid"].coeffRef(row, col) = (-1) * scheme_ref->prop->rho["face"][col] * (vf__.dot(Sf__)) * (e_gradf__.dot(dCf__)); 
    }
    else
    {
        return;
    };
};
//energy
void energy::calc_gamma(scheme*& scheme_ref, double AH)
{
    for(std::pair<std::string, make<double>::sp_mat> entry : this->lhs_cc)
    {
        make<double>::map_int& gamma_c_ = this->gamma["cell"][entry.first];
        make<double>::map_int& gamma_f_ = this->gamma["face"][entry.first];
        make<double>::sp_mat& gc_ = scheme_ref->mesh->constants["gc"][entry.first];
        make<int>::sp_mat& cc_fc_ = scheme_ref->mesh->cc_fc[entry.first];
        make<double>::map_int& k_ = scheme_ref->prop->k["face"];
        make<double>::map_int& cp_c_ = scheme_ref->prop->cp["cell"];
        make<double>::map_int& rho_c_ = scheme_ref->prop->rho["cell"];
        make<double>::map_int& miu_c_ = scheme_ref->prop->miu["cell"];
        make<double>::map_int& miut_ = scheme_ref->wall->miut;
        double Pr__;
        if(entry.first.compare("fluid") == 0)
        {
            for(int i = 0; i < cc_fc_.outerSize(); i++)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
                {
                    Pr__ = miu_c_[it.row()] / (rho_c_[it.row()] * calc_fluid_prop("alpha",
                           this->value->cvalue["fluid"][it.row()], scheme_ref->pressure->cvalue["fluid"][it.row()], AH));
                    gamma_c_[it.row()] = cp_c_[it.row()] * (miu_c_[it.row()] / Pr__ + miut_[it.row()] / Pr__ );
                    break;
                };
            };
            for(int i = 0; i < cc_fc_.outerSize(); i++)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
                {
                    gamma_f_[it.value()] = gamma_c_[it.row()] * gamma_c_[it.col()] / (gamma_c_[it.col()] *
                                        gc_.coeffRef(it.row(), it.col()) + gamma_c_[it.row()] * (1 -
                                        gc_.coeffRef(it.row(), it.col())));
                };
            };
        }
        else if(entry.first.compare("solid") == 0)
        {
            for(int i = 0; i < cc_fc_.outerSize(); i++)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
                {
                    gamma_c_[it.row()] = k_[it.row()];
                };
            };
            for(int i = 0; i < cc_fc_.outerSize(); i++)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
                {
                    gamma_f_[it.value()] = gamma_c_[it.row()] * gamma_c_[it.col()] / (gamma_c_[it.col()] *
                                        gc_.coeffRef(it.row(), it.col()) + gamma_c_[it.row()] * (1 -
                                        gc_.coeffRef(it.row(), it.col())));
                };
            };
        }
        else if(entry.first.compare("conj") == 0)
        {
            for(int i = 0; i < cc_fc_.outerSize(); i++)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
                {
                    if(std::find(scheme_ref->mesh->cid["misc"]["wall"].begin(), scheme_ref->mesh->cid["misc"]["wall"].end(), it.row())
                                 != scheme_ref->mesh->cid["misc"]["wall"].end())
                    {
                        // if fluid
                        gamma_f_[it.value()] = gamma_c_[it.row()] * gamma_c_[it.col()] / (gamma_c_[it.row()] *
                                                gc_.coeffRef(it.col(), it.row()) + gamma_c_[it.col()] * (1 -
                                                gc_.coeffRef(it.col(), it.row())));
                    }
                    else
                    {
                        // if solid
                        gamma_f_[it.value()] = gamma_c_[it.row()] * gamma_c_[it.col()] / (gamma_c_[it.col()] *
                                                gc_.coeffRef(it.row(), it.col()) + gamma_c_[it.row()] * (1 -
                                                gc_.coeffRef(it.row(), it.col())));
                    };

                };
            };
        };
    };
};
void energy::calc_wall(scheme*& scheme_ref, double AH)
{
    make<int>::vec& wall_cell_ = scheme_ref->mesh->cid["misc"]["wall"];
    make<double>::map_int& valuec_ = this->value->cvalue["fluid"];
    make<coor>::map_int& gradc_ = this->value->cgrad["fluid"];
    coor gradc__; double valuec__;
    double d_plus__; double d_perp__; double u_tau__; double viu__; double Pr__; double beta__;
    for(auto i = wall_cell_.begin(); i != wall_cell_.end(); i++)
    {
        gradc__ = gradc_[*i];
        valuec__ = valuec_[*i];
        u_tau__ = scheme_ref->wall->utau[*i];
        viu__ = scheme_ref->prop->miu["cell"][*i] / scheme_ref->prop->rho["cell"][*i];
        d_perp__ = pow(gradc__.squaredNorm() + 2*valuec__, 0.5) - gradc__.norm();
        d_plus__ = d_perp__ * u_tau__ / viu__;
        Pr__ = viu__ / calc_fluid_prop("alpha", this->value->cvalue["fluid"][*i], scheme_ref->pressure->cvalue["fluid"][*i], AH);
        beta__ = pow(3.85 * pow(Pr__, 1/3), 2)  + 2.12 * std::log(Pr__);
        valuec_[*i] = 2.12 * std::log(d_plus__) + beta__;
    };
};
void energy::calc_lhs(scheme*& scheme_ref, turb_k*& k_ref, momentum*& u_ref, momentum*& v_ref, momentum*& w_ref, double AH)
{
    for(std::pair<std::string, make<double>::sp_mat> entry : this->lhs_cc)
    {
        make<double>::sp_mat& fc_ = this->lhs_fc[entry.first];
        make<double>::sp_mat& cc_ = this->lhs_cc[entry.first];
        make<double>::sp_mat& gc_ = scheme_ref->mesh->constants["gc"][entry.first];
        make<int>::sp_mat& cc_fc_ = scheme_ref->mesh->cc_fc[entry.first];
        make<double>::sp_mat& fcg_conv_aC_ = scheme_ref->mesh->constants["g_conv_aC"][entry.first];
        make<double>::sp_mat& fcg_conv_aF_ = scheme_ref->mesh->constants["g_conv_aF"][entry.first];
        make<double>::sp_mat& fcg_diff_aC_ = scheme_ref->mesh->constants["g_diff_aC"][entry.first];
        make<double>::sp_mat& fcg_diff_aF_ = scheme_ref->mesh->constants["g_diff_aF"][entry.first];
        make<double>::map_int& rho_c_ = scheme_ref->prop->rho["cell"];
        make<double>::map_int& rho_f_ = scheme_ref->prop->rho["face"];
        make<double>::map_int& cp_c_ = scheme_ref->prop->cp["cell"];
        make<double>::sp_mat& rho_v_sf_ = scheme_ref->rho_v_sf;
        make<double>::map_int& gamma_f_ = this->gamma["face"][entry.first];
        axes& Sf_ = scheme_ref->mesh->geom["Sf"][entry.first];
        axes& Ef_ = scheme_ref->mesh->geom["Ef"][entry.first];
        axes& Tf_ = scheme_ref->mesh->geom["Tf"][entry.first];
        axes& eCF_ = scheme_ref->mesh->geom["eCF"][entry.first];
        axes& eCf_ = scheme_ref->mesh->geom["eCf"][entry.first];
        axes& dCF_ = scheme_ref->mesh->geom["dCF"][entry.first];
        coor Sf__; coor Ef__; coor Tf__; coor eCF__; coor eCf__; coor dCf__;
        double Pf__; double kf__;
        coor gradf_itr__; coor gradc__;
        double bC; coor bf_conv; coor bf_diff;
        double dCFval__; double cp_f__;
        int row; double aC;
        for(int i = 0; i < cc_.outerSize(); i++)
        {
            aC = 0.0;
            if(entry.first.compare("fluid") == 0)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
                {
                    cp_f__ = scheme_ref->prop->cp["face"][it.value()];
                    fc_.coeffRef(it.row(), it.value()) = cp_f__ * fcg_conv_aC_.coeffRef(it.row(), it.value()) *
                                                        rho_v_sf_.coeffRef(it.row(), it.value()) + 
                                                        fcg_diff_aC_.coeffRef(it.row(), it.value()) *
                                                        gamma_f_[it.value()];
                    cc_.coeffRef(it.row(), it.col()) = cp_f__ * fcg_conv_aF_.coeffRef(it.row(), it.value()) *
                                                    rho_v_sf_.coeffRef(it.row(), it.value()) + 
                                                    fcg_diff_aF_.coeffRef(it.row(), it.value()) *
                                                    gamma_f_[it.value()];
                    row = it.row();
                };
                // boundary
                make<std::string>::vec bid_str_keys = get_map_keys<std::string, make<make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>::map_str>(scheme_ref->mesh->bid);
                if(std::find(bid_str_keys.begin(), bid_str_keys.end(), "fluid") != bid_str_keys.end())
                {
                    make<int>::vec bid_int_keys = get_map_keys<int, make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>(scheme_ref->mesh->bid["fluid"]);
                    make<std::pair<std::string, int>>::vec fc_bound;
                    if(std::find(bid_int_keys.begin(), bid_int_keys.end(), row) != bid_int_keys.end())
                    {
                        for(std::pair<int, make<std::pair<std::string, int>>::vec> entry_bound : scheme_ref->mesh->bid["fluid"][row])
                        {
                            fc_bound = scheme_ref->mesh->bid["fluid"][row][entry_bound.first];
                            for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                            {
                                this->calc_bound_lhs(row, entry_bound.first, j->first, j->second, scheme_ref, "fluid", u_ref, v_ref, w_ref, AH);
                            };
                        };
                    };
                };
                for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; ++it)
                {
                    row = it.row();
                    aC += fc_.coeffRef(it.row(), it.col());
                };
                cc_.coeffRef(row, row) = aC;
            }
            else if(entry.first.compare("solid") == 0)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
                {
                    fc_.coeffRef(it.row(), it.value()) = fcg_diff_aC_.coeffRef(it.row(), it.value()) *
                                                        gamma_f_[it.value()];
                    cc_.coeffRef(it.row(), it.col()) = fcg_diff_aF_.coeffRef(it.row(), it.value()) *
                                                    gamma_f_[it.value()];
                    row = it.row();
                };
                // boundary
                make<std::string>::vec bid_str_keys = get_map_keys<std::string, make<make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>::map_str>(scheme_ref->mesh->bid);
                if(std::find(bid_str_keys.begin(), bid_str_keys.end(), "solid") != bid_str_keys.end())
                {
                    make<int>::vec bid_int_keys = get_map_keys<int, make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>(scheme_ref->mesh->bid["solid"]);
                    make<std::pair<std::string, int>>::vec fc_bound;
                    if(std::find(bid_int_keys.begin(), bid_int_keys.end(), row) != bid_int_keys.end())
                    {
                        for(std::pair<int, make<std::pair<std::string, int>>::vec> entry_bound : scheme_ref->mesh->bid["solid"][row])
                        {
                            fc_bound = scheme_ref->mesh->bid["solid"][row][entry_bound.first];
                            for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                            {
                                this->calc_bound_lhs(row, entry_bound.first, j->first, j->second, scheme_ref, "solid", u_ref, v_ref, w_ref, AH);
                            };
                        };
                    };
                };
                for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; ++it)
                {
                    row = it.row();
                    aC += fc_.coeffRef(it.row(), it.col());
                };
                cc_.coeffRef(row, row) = aC;
            }
            else if(entry.first.compare("conj") == 0)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
                {
                    if(std::find(scheme_ref->mesh->cid["misc"]["wall"].begin(), scheme_ref->mesh->cid["misc"]["wall"].end(), it.row())
                                 != scheme_ref->mesh->cid["misc"]["wall"].end())
                    {
                        double h_conj = rho_c_[it.row()] * cp_c_[it.row()] * pow(0.09, 0.25) * pow(k_ref->value->cvalue["fluid"][it.row()], 0.5) / this->value->cvalue["fluid"][it.row()];
                        double Sf_val__ = Sf_.axes_to_val(it.col(), it.value());
                        double hsf__ = h_conj * Sf_val__;
                        fc_.coeffRef(it.row(), it.value()) = (-1) * gamma_f_[it.value()] * hsf__ / (gamma_f_[it.value()] - hsf__);
                        cc_.coeffRef(it.row(), it.col()) = gamma_f_[it.value()] * hsf__ / (gamma_f_[it.value()] - hsf__);
                    }
                    else
                    {
                        double h_conj = rho_c_[it.col()] * cp_c_[it.col()] * pow(0.09, 0.25) * pow(k_ref->value->cvalue["fluid"][it.col()], 0.5) / this->value->cvalue["fluid"][it.col()];
                        double Sf_val__ = Sf_.axes_to_val(it.row(), it.value());
                        double hsf__ = h_conj * Sf_val__;
                        fc_.coeffRef(it.row(), it.value()) = gamma_f_[it.value()] * hsf__ / (gamma_f_[it.value()] - hsf__);
                        cc_.coeffRef(it.row(), it.col()) = (-1) * gamma_f_[it.value()] * hsf__ / (gamma_f_[it.value()] - hsf__);
                    };
                };
                for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; ++it)
                {
                    row = it.row();
                    aC += fc_.coeffRef(it.row(), it.col());
                };
                cc_.coeffRef(row, row) = aC;
            };
        };
    };
};
void energy::calc_rhs(scheme*& scheme_ref, momentum*& u_ref, momentum*& v_ref, momentum*& w_ref, double AH)
{
    for(std::pair<std::string, make<double>::sp_mat> entry : this->rhs_cc)
    {
        make<double>::sp_mat& fc_ = this->rhs_fc[entry.first];
        make<double>::sp_mat& cc_ = this->rhs_cc[entry.first];
        make<double>::sp_mat& gc_ = scheme_ref->mesh->constants["gc"][entry.first];
        make<int>::sp_mat& cc_fc_ = scheme_ref->mesh->cc_fc[entry.first];
        make<double>::sp_mat& fcg_conv_aC_ = scheme_ref->mesh->constants["g_conv_aC"][entry.first];
        make<double>::sp_mat& fcg_conv_aF_ = scheme_ref->mesh->constants["g_conv_aF"][entry.first];
        make<double>::sp_mat& fcg_diff_aC_ = scheme_ref->mesh->constants["g_diff_aC"][entry.first];
        make<double>::sp_mat& fcg_diff_aF_ = scheme_ref->mesh->constants["g_diff_aF"][entry.first];
        make<double>::map_int& rho_f_ = scheme_ref->prop->rho["face"];
        make<double>::sp_mat& rho_v_sf_ = scheme_ref->rho_v_sf;
        make<double>::map_int& gamma_f_ = this->gamma["face"][entry.first];
        make<double>::map_int& vol_ = scheme_ref->mesh->size["volume"];
        double gc__;
        axes& Sf_ = scheme_ref->mesh->geom["Sf"][entry.first];
        axes& Ef_ = scheme_ref->mesh->geom["Ef"][entry.first];
        axes& Tf_ = scheme_ref->mesh->geom["Tf"][entry.first];
        axes& eCF_ = scheme_ref->mesh->geom["eCF"][entry.first];
        axes& eCf_ = scheme_ref->mesh->geom["eCf"][entry.first];
        axes& dCf_ = scheme_ref->mesh->geom["dCf"][entry.first];
        coor Sf__; coor Ef__; coor Tf__; coor eCF__; coor eCf__; coor dCf__;
        double cp__;
        coor gradf_itr__; coor gradc__;
        int row; double bC; coor bf_conv; coor bf_diff; double bf_body;
        for(int i = 0; i < cc_.outerSize(); i++)
        {
            bC = 0.0;
            if(entry.first.compare("fluid") == 0)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
                {
                    gc__ = gc_.coeffRef(it.row(), it.value());
                    gradc__ = this->value->cgrad[entry.first][it.row()];
                    gradf_itr__ = gc__ * gradc__ + (1 - gc__) * this->value->cgrad[entry.first][it.col()];
                    Ef__ = Ef_.axes_to_coor(it.row(), it.value());
                    Tf__ = Tf_.axes_to_coor(it.row(), it.value());
                    eCF__ = eCF_.axes_to_coor(it.row(), it.col());
                    dCf__ = dCf_.axes_to_coor(it.row(), it.value());
                    bf_conv = (gradf_itr__.dot(eCF__) * eCF__ - (gradc__ + gradf_itr__)) / 2;
                    bf_diff = gradf_itr__ - (gradf_itr__.dot(eCF__) * eCF__);
                    fc_.coeffRef(it.row(), it.value()) = (bf_conv.dot(dCf__) * rho_v_sf_.coeffRef(it.row(), it.value())) +
                                                        (bf_diff.dot(Ef__) * gamma_f_[it.value()]) + 
                                                        (bf_diff.dot(Tf__) * gamma_f_[it.value()]);
                };
                // boundary
                make<std::string>::vec bid_str_keys = get_map_keys<std::string, make<make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>::map_str>(scheme_ref->mesh->bid);
                if(std::find(bid_str_keys.begin(), bid_str_keys.end(), "fluid") != bid_str_keys.end())
                {
                    make<int>::vec bid_int_keys = get_map_keys<int, make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>(scheme_ref->mesh->bid["fluid"]);
                    make<std::pair<std::string, int>>::vec fc_bound;
                    if(std::find(bid_int_keys.begin(), bid_int_keys.end(), row) != bid_int_keys.end())
                    {
                        for(std::pair<int, make<std::pair<std::string, int>>::vec> entry_bound : scheme_ref->mesh->bid["fluid"][row])
                        {
                            fc_bound = scheme_ref->mesh->bid["fluid"][row][entry_bound.first];
                            for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                            {
                                this->calc_bound_rhs(row, entry_bound.first, j->first, j->second, scheme_ref, "fluid", u_ref, v_ref, w_ref, AH);
                            };
                        };
                    };
                };
                // append total
                for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; ++it)
                {
                    row = it.row();
                    bC += fc_.coeffRef(it.row(), it.col());
                };
                // body term (fluid only)
                double miu_c__ = scheme_ref->prop->miu["cell"][row];
                double miut__ = scheme_ref->wall->miut[row];
                double phi_v__ = scheme_ref->phi_v[row];
                double vol__ = vol_[row];
                cc_.coeffRef(row, row) = bC + ((miu_c__ + miut__) * phi_v__ * vol__);
            }
            else if(entry.first.compare("solid") == 0)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; ++it)
                {
                    gc__ = gc_.coeffRef(it.row(), it.value());
                    gradc__ = this->value->cgrad[entry.first][it.row()];
                    gradf_itr__ = gc__ * gradc__ + (1 - gc__) * this->value->cgrad[entry.first][it.col()];
                    Ef__ = Ef_.axes_to_coor(it.row(), it.value());
                    Tf__ = Tf_.axes_to_coor(it.row(), it.value());
                    eCF__ = eCF_.axes_to_coor(it.row(), it.col());
                    bf_diff = gradf_itr__ - (gradf_itr__.dot(eCF__) * eCF__);
                    fc_.coeffRef(it.row(), it.value()) = (bf_diff.dot(Ef__) * gamma_f_[it.value()]) + 
                                                        (bf_diff.dot(Tf__) * gamma_f_[it.value()]);
                };
                // boundary
                make<std::string>::vec bid_str_keys = get_map_keys<std::string, make<make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>::map_str>(scheme_ref->mesh->bid);
                if(std::find(bid_str_keys.begin(), bid_str_keys.end(), "solid") != bid_str_keys.end())
                {
                    make<int>::vec bid_int_keys = get_map_keys<int, make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>(scheme_ref->mesh->bid["solid"]);
                    make<std::pair<std::string, int>>::vec fc_bound;
                    if(std::find(bid_int_keys.begin(), bid_int_keys.end(), row) != bid_int_keys.end())
                    {
                        for(std::pair<int, make<std::pair<std::string, int>>::vec> entry_bound : scheme_ref->mesh->bid["solid"][row])
                        {
                            fc_bound = scheme_ref->mesh->bid["solid"][row][entry_bound.first];
                            for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                            {
                                this->calc_bound_rhs(row, entry_bound.first, j->first, j->second, scheme_ref, "solid", u_ref, v_ref, w_ref, AH);
                            };
                        };
                    };
                };
                // append total
                for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; ++it)
                {
                    row = it.row();
                    bC += fc_.coeffRef(it.row(), it.col());
                };
                // soil source
                make<int>::vec soil_id = scheme_ref->mesh->cid["solid"]["soil"];
                if(std::find(soil_id.begin(), soil_id.end(), row) != soil_id.end())
                {
                    double vol__ = vol_[row];
                    bC += scheme_ref->source->cell_value["solid"]["soil"] * vol__;
                };
                cc_.coeffRef(row, row) = bC;
            };
        };
    };
};
void energy::calc_bound_lhs(int row, int col, std::string check, int unique, scheme*& scheme_ref, std::string domain__, momentum*& u_ref, momentum*& v_ref, momentum*& w_ref, double AH)
{
    if(check.compare("temp"))
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"][domain__];
        axes& dCf_ = scheme_ref->mesh->geom["dCf"][domain__];
        make<double>::map_int& gamma_f_ = this->gamma["face"][domain__];
        this->lhs_fc[domain__].coeffRef(row, col) = gamma_f_[col] * Sf_.axes_to_val(row, col) / dCf_.axes_to_val(row, col);
    }
    else if(check.compare("hamb"))
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"][domain__];
        axes& dCf_ = scheme_ref->mesh->geom["dCf"][domain__];
        make<double>::map_int& gamma_f_ = this->gamma["face"][domain__];
        // h_sky
        double g_hamb = gamma_f_[col] * Sf_.axes_to_val(row, col) / dCf_.axes_to_val(row, col);
        double T_f__ = this->value->fvalue[domain__][col];
        double T_amb__ = scheme_ref->source->face_value["hamb"][0]; double T_sky__ = 0.0552 * pow(T_amb__, 1.5);
        double eps_C__ = 0.8; // assumption for air
        double h_sky__ = 5.67 * pow(10, -8) * eps_C__ * (T_f__ + T_sky__) * (pow(T_f__, 2) + pow(T_sky__, 2)) * (T_f__ - T_sky__) / (T_f__ - T_amb__);
        // h_free_conv
        double T_film__ = (T_f__ + T_amb__) / 2;
        double miu_film__ = calc_fluid_prop("miu", T_film__, 101325, AH);
        double rho_film__ = calc_fluid_prop("rho", T_film__, 101325, AH);
        double viu_film__ = miu_film__ / rho_film__;
        double k_film__ = calc_fluid_prop("k", T_film__, 101325, AH);
        double alpha_film__ = calc_fluid_prop("alpha", T_film__, 101325, AH);
        double area__ = scheme_ref->mesh->size["area"][col];
        double charl__ = pow(area__, 0.5);
        double Ral__ = 9.81 * (T_f__ - T_amb__) * pow(charl__, 3) / (T_film__ * viu_film__ * alpha_film__);
        double NuN__ = 0.0;
        if(Ral__ < pow(10, 7))
        {
            NuN__ = 0.54 * pow(Ral__, 0.25);
        }
        else
        {
            NuN__ = 0.15 * pow(Ral__, 1/3);
        };
        double h_conv__ = NuN__ * k_film__ / charl__;
        double g_hamb__ = gamma_f_[col] * Sf_.axes_to_val(row, col) / dCf_.axes_to_val(row, col);
        this->lhs_fc[domain__].coeffRef(row, col) = g_hamb__ * (h_sky__ + h_conv__) * Sf_.axes_to_val(row, col) / (g_hamb__ + (h_sky__ + h_conv__) * Sf_.axes_to_val(row, col));
    }
    else if(check.compare("out") == 0)
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref->mesh->geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref->mesh->geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        // momentum
        coor vf__;
        coor v0_gradc__; coor v1_gradc__; coor v2_gradc__;
        coor v0_gradf__; coor v1_gradf__; coor v2_gradf__;
        double v0_f__; double v1_f__; double v2_f__;
        v0_gradc__ = u_ref->value->cgrad["fluid"][row];
        v1_gradc__ = v_ref->value->cgrad["fluid"][row];
        v2_gradc__ = w_ref->value->cgrad["fluid"][row];
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v1_gradf__ = v1_gradc__ - (v1_gradc__.dot(eCf__) * eCf__);
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v0_f__ = u_ref->value->cvalue["fluid"][row] + (v0_gradf__.dot(dCf__));
        v1_f__ = v_ref->value->cvalue["fluid"][row] + (v1_gradf__.dot(dCf__));
        v2_f__ = w_ref->value->cvalue["fluid"][row] + (v2_gradf__.dot(dCf__));
        vf__(u_ref->axis) = v0_f__; vf__(v_ref->axis) = v1_f__; vf__(w_ref->axis) = v2_f__; 
        this->lhs_fc["fluid"].coeffRef(row, col) = scheme_ref->prop->rho["face"][col] * (vf__.dot(Sf__)); 
    }
    else
    {
        return;
    };
};
void energy::calc_bound_rhs(int row, int col, std::string check, int unique, scheme*& scheme_ref, std::string domain__, momentum*& u_ref, momentum*& v_ref, momentum*& w_ref, double AH)
{
    if(check.compare("temp"))
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"][domain__];
        axes& dCf_ = scheme_ref->mesh->geom["dCf"][domain__];
        make<double>::map_int& gamma_f_ = this->gamma["face"][domain__];
        this->rhs_fc[domain__].coeffRef(row, col) = gamma_f_[col] * Sf_.axes_to_val(row, col) * scheme_ref->source->face_value[check][unique]/ dCf_.axes_to_val(row, col);
    }
    else if(check.compare("flux"))
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"][domain__];
        this->rhs_fc[domain__].coeffRef(row, col) = scheme_ref->source->face_value[check][unique] * Sf_.axes_to_val(row, col);
    }
    else if(check.compare("s2s"))
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"][domain__];
        this->rhs_fc[domain__].coeffRef(row, col) = scheme_ref->source->face_value[check][unique];
    }
    else if(check.compare("hamb"))
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"][domain__];
        axes& dCf_ = scheme_ref->mesh->geom["dCf"][domain__];
        make<double>::map_int& gamma_f_ = this->gamma["face"][domain__];
        // h_sky
        double g_hamb = gamma_f_[col] * Sf_.axes_to_val(row, col) / dCf_.axes_to_val(row, col);
        double T_f__ = this->value->fvalue[domain__][col];
        double T_amb__ = scheme_ref->source->face_value["hamb"][0]; double T_sky__ = 0.0552 * pow(T_amb__, 1.5);
        double eps_C__ = 0.8; // assumption for air
        double h_sky__ = 5.67 * pow(10, -8) * eps_C__ * (T_f__ + T_sky__) * (pow(T_f__, 2) + pow(T_sky__, 2)) * (T_f__ - T_sky__) / (T_f__ - T_amb__);
        // h_free_conv
        double T_film__ = (T_f__ + T_amb__) / 2;
        double miu_film__ = calc_fluid_prop("miu", T_film__, 101325, AH);
        double rho_film__ = calc_fluid_prop("rho", T_film__, 101325, AH);
        double viu_film__ = miu_film__ / rho_film__;
        double k_film__ = calc_fluid_prop("k", T_film__, 101325, AH);
        double alpha_film__ = calc_fluid_prop("alpha", T_film__, 101325, AH);
        double charl__ = pow(scheme_ref->mesh->size["area"][col], 0.5);
        double Ral__ = 9.81 * (T_f__ - T_amb__) * pow(charl__, 3) / (T_film__ * viu_film__ * alpha_film__);
        double NuN__ = 0.0;
        if(Ral__ < pow(10, 7))
        {
            NuN__ = 0.54 * pow(Ral__, 0.25);
        }
        else
        {
            NuN__ = 0.15 * pow(Ral__, 1/3);
        };
        double h_conv__ = NuN__ * k_film__ / charl__;
        double g_hamb__ = gamma_f_[col] * Sf_.axes_to_val(row, col) / dCf_.axes_to_val(row, col);
        this->rhs_fc[domain__].coeffRef(row, col) = (-1) * g_hamb__ * (h_sky__ + h_conv__) * Sf_.axes_to_val(row, col) * T_amb__ / (g_hamb__ + (h_sky__ + h_conv__) * Sf_.axes_to_val(row, col));
    }
    else if(check.compare("out") == 0)
    {
        axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref->mesh->geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref->mesh->geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        // momentum
        coor vf__;
        coor v0_gradc__; coor v1_gradc__; coor v2_gradc__;
        coor v0_gradf__; coor v1_gradf__; coor v2_gradf__;
        double v0_f__; double v1_f__; double v2_f__;
        v0_gradc__ = u_ref->value->cgrad["fluid"][row];
        v1_gradc__ = v_ref->value->cgrad["fluid"][row];
        v2_gradc__ = w_ref->value->cgrad["fluid"][row];
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v1_gradf__ = v1_gradc__ - (v1_gradc__.dot(eCf__) * eCf__);
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v0_f__ = u_ref->value->cvalue["fluid"][row] + (v0_gradf__.dot(dCf__));
        v1_f__ = v_ref->value->cvalue["fluid"][row] + (v1_gradf__.dot(dCf__));
        v2_f__ = w_ref->value->cvalue["fluid"][row] + (v2_gradf__.dot(dCf__));
        vf__(u_ref->axis) = v0_f__; vf__(v_ref->axis) = v1_f__; vf__(w_ref->axis) = v2_f__;
        // turb_e
        coor T_gradc__; coor T_gradf__;
        T_gradc__ = this->value->cgrad["fluid"][row];
        T_gradf__ = T_gradc__ - (T_gradc__.dot(eCf__) * eCf__);
        this->rhs_fc["fluid"].coeffRef(row, col) = (-1) * scheme_ref->prop->rho["face"][col] * (vf__.dot(Sf__)) * (T_gradf__.dot(dCf__)); 
    }
    else
    {
        return;
    };
};
// s2s
void s2s::calc_lhs(scheme*& scheme_ref)
{
    make<double>::sp_mat& cc_ = this->lhs_cc["s2s"];
    make<make<int>::vec>::map_int& s2s_f_l_ = scheme_ref->mesh->fid["s2s"];
    make<double>::sp_mat& view_ = scheme_ref->mesh->constants["view"]["s2s"];
    make<double>::map_int& rho_f_ = scheme_ref->prop->rho["face"];
    make<double>::map_int& area_ = scheme_ref->mesh->size["area"];
    double rho_f__; double area_clust__; make<int>::vec s2s_f_vec__;
    int row;
    for(int i = 0; i < cc_.outerSize(); i++)
    {
        rho_f__ = 0.0; area_clust__ = 0.0;
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(cc_, i); it; ++it)
        {
            cc_.coeffRef(it.row(), it.row()) = 1;
            s2s_f_vec__ = s2s_f_l_[it.row()];
            for(auto j = s2s_f_vec__.begin(); j != s2s_f_vec__.end(); j++)
            {
                rho_f__ += rho_f_[*j] * area_[*j];
                area_clust__ += area_[*j];
            };
            rho_f__ = rho_f__ / area_clust__;
            cc_.coeffRef(it.row(), it.col()) = rho_f__ * view_.coeffRef(it.row(), it.col());
        };
    };
};
void s2s::calc_rhs(scheme*& scheme_ref, energy*& energy_ref)
{
    make<double>::sp_mat& cc_ = this->rhs_cc["fluid"];
    make<make<int>::vec>::map_int& s2s_f_l_ = scheme_ref->mesh->fid["s2s"];
    make<double>::sp_mat& view_ = scheme_ref->mesh->constants["view"]["s2s"];
    make<double>::map_int& eps_s2s_ = scheme_ref->prop->eps["s2s"];
    make<double>::map_int& area_ = scheme_ref->mesh->size["area"];
    double eps_s2s__; double T_f__; double area_clust__; make<int>::vec s2s_f_vec__;
    int row;
    for(int i = 0; i < cc_.outerSize(); i++)
    {
        T_f__ = 0.0; area_clust__ = 0.0;
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(cc_, i); it; ++it)
        {
            cc_.coeffRef(it.row(), it.row()) = 1;
            for(auto j = s2s_f_vec__.begin(); j != s2s_f_vec__.end(); j++)
            {
                T_f__ += energy_ref->value->fvalue["fluid"][*j] * area_[*j];
                area_clust__ += area_[*j];
            };
            eps_s2s__ = eps_s2s_[it.row()];
            T_f__ = T_f__ / area_clust__;
            cc_.coeffRef(it.row(), it.col()) = eps_s2s__ * 5.67 * pow(10, -8) * pow(T_f__, 4);
        };
    };
};

// other solver func
template <class V>
double linear_face_itr_value(int C_id, int F_id, V* eq__, scheme*& scheme_ref, std::string domain__)
{
    double gc__ = scheme_ref->mesh->constants["gc"][domain__].coeffRef(C_id, F_id);
    double Cval__ = eq__->value->cvalue[domain__][C_id];
    double Fval__ = eq__->value->cvalue[domain__][F_id];
    return gc__ * Cval__ + (1 - gc__) * Fval__;
};
template <class V>
coor linear_face_itr_grad(int C_id, int F_id, V* eq__, scheme*& scheme_ref, std::string domain__)
{
    double gc__ = scheme_ref->mesh->constants["gc"][domain__].coeffRef(C_id, F_id);
    coor Cgrad__ = eq__->value->cgrad[domain__][C_id];
    coor Fgrad__ = eq__->value->cgrad[domain__][F_id];
    coor get = (Cgrad__ * gc__) + ((1 - gc__) * Fgrad__);
    return get;
};
template <class V>
double linear_face_itr_gamma(int C_id, int F_id, V* eq__, scheme*& scheme_ref, std::string domain__)
{
    double gc__ = scheme_ref->mesh->constants["gc"][domain__].coeffRef(C_id, F_id);
    double Cgamma__ = eq__->gamma["cell"][domain__][C_id];
    double Fgamma__ = eq__->gamma["cell"][domain__][F_id];
    return (gc__ * Fgamma__ + (1 - gc__) * Cgamma__) / (Cgamma__ * Fgamma__);
};
template <class V>
double quick_face_itr_value(int C_id, int F_id, V*& eq__, scheme*& scheme_ref, std::string domain__)
{
    int f_id; coor Cgrad__; coor Fgrad__; double Cval__;
    f_id = scheme_ref->mesh->cc_fc[domain__].coeffRef(C_id, F_id);
    coor cdCf__ = scheme_ref->mesh->geom["dCf"][domain__].axes_to_coor(C_id, f_id);
    Cgrad__ = eq__->value->cgrad[domain__][C_id];
    Fgrad__ = eq__->value->cgrad[domain__][F_id];
    coor fgrad__ = eq__->value->fgrad[domain__][f_id];
    Cval__ = eq__->value->cvalue[domain__][C_id];
    double add = fgrad__.dot(cdCf__);
    double get = Cval__ + add / 2;
    return get;
};
template <class V>
coor quick_face_itr_grad(int C_id, int F_id, V* eq__, scheme*& scheme_ref, std::string domain__)
{
    int f_id; coor Cgrad__; coor Fgrad__; double Cval__;
    f_id = scheme_ref->mesh->cc_fc[domain__].coeffRef(C_id, F_id);
    double gc__ = scheme_ref->mesh->constants["gc"][domain__].coeffRef(C_id, F_id);
    coor ceCF__ = scheme_ref->mesh->geom["eCF"][domain__].axes_to_coor(C_id, F_id);
    double vdCf__ = scheme_ref->mesh->geom["dCf"][domain__].axes_to_val(C_id,  f_id);
    Cgrad__ = eq__->value->cgrad[domain__][C_id];
    Fgrad__ = eq__->value->cgrad[domain__][F_id];
    Cval__ = eq__->value->cvalue[domain__][C_id];
    double Fval__ = eq__->value->cvalue[domain__][F_id];
    coor fgrad_itr__ = linear_face_itr_grad(C_id, F_id, eq__, scheme_ref, domain__);
    coor get = fgrad_itr__ + ((Fval__ - Cval__) / vdCf__) * ceCF__ - (Fgrad__.dot(ceCF__) * ceCF__);
    return get;
};
template <class V>
void least_square_itr(V* eq__, scheme*& scheme_ref, std::string domain__)
{
    // calculate cgrad with optimization [A][B] = [C]
    make<double>::sp_mat& lhs_cc_ = eq__->lhs_cc[domain__];
    axes& dCF_ = scheme_ref->mesh->geom["dCF"][domain__];
    double Cval__; double Fval__;
    coor cdCF__; double wk__;
    int F_id;
    for(int i = 0; i < lhs_cc_.outerSize(); i++)
    {
        Eigen::Matrix3d ls_lhs__; ls_lhs__.setZero();
        Eigen::Vector3d ls_rhs__; ls_rhs__.setZero();
        int row;
        for(make<double>::sp_mat::InnerIterator it(lhs_cc_, i); it; ++it)
        {
            Cval__ = eq__->value->cvalue[domain__][it.row()];
            Fval__ = eq__->value->cvalue[domain__][it.col()];
            cdCF__ = dCF_.axes_to_coor(it.row(), it.col());
            wk__ = 1 / dCF_.axes_to_val(it.row(), it.col());
            for(int j = 0; j < 2; j++)
            {
                for(int k = 0; k < 2; k++)
                {
                    ls_lhs__(i, j) += wk__ * cdCF__(i) * cdCF__(j);
                };
                ls_rhs__(i) += wk__ * cdCF__(i) * (Fval__ - Cval__);
            };
            row = it.row();
        };
        Eigen::Vector3d X = ls_lhs__.lu().solve(ls_rhs__);
        eq__->value->cgrad[domain__][row] = X.sparseView();
    };
};
template <class V>
void make_prev_lhs(solver<V>* solv__, scheme*& scheme_ref, int step_length__, bool is_init)
{
    make<double>::map_int& rho_c_ = scheme_ref->prop->rho["cell"];
    make<double>::map_int& vol_ = scheme_ref->mesh->size["volume"];
    if(is_init)
    {
        // rho * vol / step length
        for(int i = 0; i < solv__->prev_lhs.outerSize(); i++)
        {
            for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(solv__->prev_lhs, i); it; ++it)
            {
                if(it.row() == it.col())
                {
                    solv__->prev_lhs.coeffRef(it.row(), it.row()) = rho_c_[it.row()] * vol_[it.row()] / step_length__;
                };
            };
        };
    }
    else
    {
        // rho * vol / (2*step length)
        for(int i = 0; i < solv__->prev_lhs.outerSize(); i++)
        {
            for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(solv__->prev_lhs, i); it; ++it)
            {
                if(it.row() == it.col())
                {
                    solv__->prev_lhs.coeffRef(it.row(), it.row()) = rho_c_[it.row()] * vol_[it.row()] / (2 * step_length__);
                };
            };
        };
    };
};
template <class V>
void make_prev_rhs(solver<V>* solv__, scheme*& scheme_ref, int step_length__, bool is_init)
{
    make<double>::map_int& rho_c_ = scheme_ref->prop->rho["cell"];
    make<double>::map_int& vol_ = scheme_ref->mesh->size["volume"];
    make<double>::map_int& cvalue_ = solv__->eq->value->cvalue["fluid"];
    if(is_init)
    {
        // rho * vol * converged value / step length
        for(int i = 0; i < int(sizeof(solv__->prev_rhs)); i++)
        {
            make<double>::map_int::iterator it_cell = cvalue_.find(i);
            if(it_cell != cvalue_.end())
            {
                solv__->prev_rhs(i) = rho_c_[i] * vol_[i] * cvalue_[i] / step_length__;
            };
        };
    }
    else
    {
        // rho * vol * converged value / (2*step length)
        for(int i = 0; i < int(sizeof(solv__->prev_rhs)); i++)
        {
            make<double>::map_int::iterator it_cell = cvalue_.find(i);
            if(it_cell != cvalue_.end())
            {
                solv__->prev_rhs(i) = rho_c_[i] * vol_[i] * cvalue_[i] / (2*step_length__);
            };
        };
    };
};
template <class V>
void make_transient_lhs(solver<V>* solv__, double under_relax__)
{
    int ctd = 0;
    make<make<double>::sp_mat>::map_str& lhs_cc_ = solv__->eq->lhs_cc;
    make<double>::sp_mat& lhs_ = solv__->lhs;
    make<double>::sp_mat& prev_lhs_ = solv__->prev_lhs;
    for(std::pair<std::string, make<double>::sp_mat> entry_cc : lhs_cc_)
    {
        if(ctd == 0)
        {
            lhs_ = entry_cc.second + prev_lhs_;
        }
        else
        {
            lhs_ += entry_cc.second;
        };
        ctd += 1;
    };
    // under-relaxation
    for(int i = 0; i < lhs_.outerSize(); i++)
    {
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(lhs_, i); it; ++it)
        {
            if(it.row() == it.col())
            {
                lhs_.coeffRef(it.row(), it.row()) = it.value() / under_relax__;
            };
        };
    };
};
template <class V>
void make_transient_rhs(solver<V>* solv__, double under_relax__)
{
    int ctd = 0;
    make<make<double>::sp_mat>::map_str& rhs_cc_ = solv__->eq->rhs_cc;
    Eigen::VectorXd& rhs_ = solv__->rhs;
    Eigen::VectorXd& prev_rhs_ = solv__->prev_rhs;
    for(std::pair<std::string, make<double>::sp_mat> entry_cc : solv__->eq->rhs_cc)
    {
        if(ctd == 0)
        {
            for(int i = 0; i < int(sizeof(rhs_)); i++)
            {
                rhs_(i) = entry_cc.second.coeffRef(i, i) + prev_rhs_(i); 
            };
        }
        else
        {
            for(int i = 0; i < int(sizeof(rhs_)); i++)
            {
                rhs_(i) += entry_cc.second.coeffRef(i, i);
            };
        };
        ctd += 1;
    };
    // under-relaxation
    for(std::pair<std::string, make<double>::sp_mat> entry_rhs : rhs_cc_)
    {
        if(entry_rhs.first.compare("conj") != 0)
        {
            for(int i = 0; i < entry_rhs.second.outerSize(); i++)
            {
                for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(entry_rhs.second, i); it; ++it)
                {
                    if(it.row() == it.col())
                    {
                        rhs_(it.row()) += (1 - under_relax__) * rhs_.coeffRef(it.row(), it.row()) *
                                          solv__->eq->value->cvalue[entry_rhs.first][it.row()] / under_relax__;

                    };
                };
            };
        };
    };
}
template <class V>
Eigen::VectorXd get_new_values(solver<V>* solv__, std::string what, make<double>::map_str* err_p__)
{
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double, RowMajor>> solver;
    make<double>::sp_mat& lhs__ = solv__->lhs;
    Eigen::VectorXd& rhs__ = solv__->rhs;
    solver.compute(lhs__);
    Eigen::VectorXd x = solver.solve(rhs__);
    make<double>::map_str& err__ = *err_p__;
    err__[what] = solver.error();
    return x;
};
template <class V>
void update_values(solver<V>* solv__, scheme*& scheme_ref, std::string what, make<double>::map_str* err_p__)
{
    V*& eq__ = solv__->eq;
    Eigen::VectorXd new_cvalues__ = get_new_values(solv__, what, err_p__);
    // i.e. gradient computation
    // update cvalue -> cgrad (by iteration) -> fgrad -> fvalue
    for(std::pair<std::string, make<double>::sp_mat> entry : eq__->lhs_cc)
    {
        make<double>::map_int& cvalue_ = eq__->value->cvalue[entry.first];
        make<int>::sp_mat& cc_fc_ = scheme_ref->mesh->cc_fc[entry.first];
        for(std::pair<int, double> entry : eq__->value->cvalue[entry.first])
        {
            entry.second = new_cvalues__(entry.first);
        };
        least_square_itr(eq__, scheme_ref, entry.first);
        for(int i = 0; i < cc_fc_.outerSize(); i++)
        {
            for(make<int>::sp_mat::InnerIterator it(cc_fc_, i); it; ++it)
            {
                coor face_itr_grad__ = quick_face_itr_grad(it.row(), it.col(), eq__, scheme_ref, entry.first);
                eq__->value->fgrad[entry.first][it.value()](0) = face_itr_grad__(0);
                eq__->value->fgrad[entry.first][it.value()](1) = face_itr_grad__(1);
                eq__->value->fgrad[entry.first][it.value()](2) = face_itr_grad__(2);
                eq__->value->fvalue[entry.first][it.value()] = quick_face_itr_value(it.row(), it.col(), eq__, scheme_ref, entry.first);;
            };
        };
    };
};
void update_wall(scheme*& scheme_ref, momentum*& u_ref, momentum*& v_ref, momentum*& w_ref,
                 turb_k*& k_ref, turb_e*& e_ref, energy*& energy_ref)
{
    for(std::pair<int, double> entry : scheme_ref->wall->wall_dist)
    {
        // utau
        coor vc__(u_ref->value->cvalue["fluid"][entry.first], v_ref->value->cvalue["fluid"][entry.first], w_ref->value->cvalue["fluid"][entry.first]);
        double v__ = sqrt_sum(vc__);
        double dist__ = scheme_ref->wall->wall_dist[entry.first];
        double viu__ = scheme_ref->prop->miu["cell"][entry.first] / scheme_ref->prop->rho["cell"][entry.first];
        double u_star__ = pow(0.09, 0.25) * pow(k_ref->value->cvalue["fluid"][entry.first], 0.5);
        double utau_check__ = v__ / (5.25 + std::log(dist__) / 0.41);
        scheme_ref->wall->utau[entry.first] = std::max(utau_check__, 11.06);
        // ts
        double ts_vec1__ = k_ref->value->cvalue["fluid"][entry.first] / e_ref->value->cvalue["fluid"][entry.first];
        coor u_gradC__ = u_ref->value->cgrad["fluid"][entry.first]; coor v_gradC__ = v_ref->value->cgrad["fluid"][entry.first]; coor w_gradC__ = w_ref->value->cgrad["fluid"][entry.first];
        double St__ = pow(u_gradC__(0), 2) + pow(v_gradC__(1), 2) + pow(w_gradC__(2), 2) +
                    (pow(u_gradC__(1) + v_gradC__(0), 2) / 2) + (pow(u_gradC__(2) + w_gradC__(0), 2) / 2) + (pow(v_gradC__(2) + w_gradC__(1), 2) / 2);
        St__ = pow(St__, 0.5);
        double ts_vec2__ = 0.6 / (pow(6, 0.5) * 0.09 * St__);
        scheme_ref->wall->miut[entry.first] = std::min(ts_vec1__, ts_vec2__);
        // miut
        double rho_C__ = scheme_ref->prop->rho["cell"][entry.first];
        double miut__ = rho_C__ * 0.09 * k_ref->value->cvalue["fluid"][entry.first] * scheme_ref->wall->ts[entry.first];
        scheme_ref->wall->miut[entry.first] = miut__;
    };
};
void update_fluid_prop(scheme*& scheme_ref, energy*& energy_ref, double AH)
{
    double pval__; double Tval__;
    for(std::pair<std::string, make<double>::map_int> entry_str : scheme_ref->prop->rho)
    {
        for(std::pair<int, double> entry_int : entry_str.second)
        {
            if(entry_str.first.compare("cell") == 0)
            {
                pval__ = scheme_ref->pressure->cvalue["fluid"][entry_int.first];
                Tval__ = energy_ref->value->cvalue["fluid"][entry_int.first];
            }
            else
            {
                pval__ = scheme_ref->pressure->fvalue["fluid"][entry_int.first];
                Tval__ = energy_ref->value->fvalue["fluid"][entry_int.first];
            };
            scheme_ref->prop->rho[entry_str.first][entry_int.first] = calc_fluid_prop("rho", pval__, Tval__, AH);
            scheme_ref->prop->miu[entry_str.first][entry_int.first] = calc_fluid_prop("miu", pval__, Tval__, AH);
            scheme_ref->prop->cp[entry_str.first][entry_int.first] = calc_fluid_prop("cp", pval__, Tval__, AH);
        };
    };
};
void update_body_term(scheme*& scheme_ref, momentum*& u_ref, momentum*& v_ref, momentum*& w_ref)
{
    // rho_v_sf
    make<double>::sp_mat& rho_v_sf_ = scheme_ref->rho_v_sf;
    make<double>::map_int& rho_f_ = scheme_ref->prop->rho["face"];
    make<double>::map_int& u_fvalue_ = u_ref->value->fvalue["fluid"];
    make<double>::map_int& v_fvalue_ = v_ref->value->fvalue["fluid"];
    make<double>::map_int& w_fvalue_ = w_ref->value->fvalue["fluid"];
    axes& Sf_ = scheme_ref->mesh->geom["Sf"]["fluid"];
    coor Sf__;
    for(int i = 0; i < rho_v_sf_.outerSize(); i++)
    {
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(rho_v_sf_, i); it; ++it)
        {
            coor v__(u_fvalue_[it.col()], v_fvalue_[it.col()], w_fvalue_[it.col()]);
            Sf__ = Sf_.axes_to_coor(it.row(), it.col());
            rho_v_sf_.coeffRef(it.row(), it.col()) = rho_f_[it.col()] * (v__.dot(Sf__));
        };
    };
    // phi_v
    make<double>::map_int& phi_v_ = scheme_ref->phi_v;
    make<coor>::map_int& u_fgrad_ = u_ref->value->fgrad["fluid"];
    make<coor>::map_int& v_fgrad_ = v_ref->value->fgrad["fluid"];
    make<coor>::map_int& w_fgrad_ = w_ref->value->fgrad["fluid"];
    coor u_fgrad__; coor v_fgrad__; coor w_fgrad__;
    double phi_v__; 
    for(std::pair<int, double> entry : phi_v_)
    {
        u_fgrad__ = u_fgrad_[entry.first]; v_fgrad__ = v_fgrad_[entry.first]; w_fgrad__ = w_fgrad_[entry.first]; 
        phi_v__ = 2 * (pow(u_fgrad__[0], 2) + pow(v_fgrad__[1], 2) + pow(w_fgrad__[2], 2))
                  + pow(u_fgrad__[1] + v_fgrad__[0], 2) + pow(u_fgrad__[2] + w_fgrad__[0], 2)
                  + pow(v_fgrad__[2] + w_fgrad__[1], 2);
        phi_v_[entry.first] = phi_v__;
    };
};
// solver
template <class V>
solver<V>::solver(V* eq, scheme*& scheme_ref, int step_length)
{
    this->eq = eq;
    make<double>::sp_mat lhs__; make<double>::sp_mat rhs__;
    make<double>::sp_mat prev_lhs__; make<double>::sp_mat prev_rhs__;
    make<sparse_input>::vec prev_lhs_input__; make<sparse_input>::vec prev_rhs_input__;
    make<double>::map_int rho_c_ = scheme_ref->prop->rho["cell"];
    make<double>::map_int vol_ = scheme_ref->mesh->size["volume"];
    int ctd = 0;
    for(std::pair<std::string, make<double>::sp_mat> entry : this->eq->lhs_cc)
    {
        if(ctd == 0)
        {
            lhs__ = entry.second;
            rhs__ = this->eq->rhs_cc[entry.first];
            if(entry.first.compare("fluid") == 0)
            {
                for(int i = 0; i < lhs__.outerSize(); i++)
                {
                    for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(lhs__, i); it; ++it)
                    {
                        if(it.row() == it.col())
                        {
                            double rho_c__ = rho_c_[it.row()];
                            double vol__ = vol_[it.row()];
                            sparse_input temp_lhs__(it.row(), it.row(), rho_c__ * vol__ / step_length);
                            sparse_input temp_rhs__(it.row(), 0, rho_c__ * vol__ * eq->value->cvalue["fluid"][it.row()]/ step_length);
                            prev_lhs_input__.push_back(temp_lhs__); prev_rhs_input__.push_back(temp_rhs__);
                        };
                    };
                };
            };
        }
        else
        {
            lhs__ += entry.second;
            rhs__ += this->eq->rhs_cc[entry.first];
            if(entry.first.compare("fluid") == 0)
            {
                for(int i = 0; i < lhs__.outerSize(); i++)
                {
                    for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(lhs__, i); it; ++it)
                    {
                        if(it.row() == it.col())
                        {
                            double rho_c__ = rho_c_[it.row()];
                            double vol__ = vol_[it.row()];
                            sparse_input temp_lhs__(it.row(), it.row(), rho_c__ * vol__ / step_length);
                            sparse_input temp_rhs__(it.row(), 0, rho_c__ * vol__ * eq->value->cvalue["fluid"][it.row()]/ step_length);
                            prev_lhs_input__.push_back(temp_lhs__); prev_rhs_input__.push_back(temp_rhs__);
                        };
                    };
                };
            };
        };
        ctd += 1;
    };
    prev_lhs__.setFromTriplets(prev_lhs_input__.begin(), prev_lhs_input__.end());
    prev_rhs__.setFromTriplets(prev_rhs_input__.begin(), prev_rhs_input__.end());
    this->lhs = lhs__; this->rhs = rhs__;
    this->prev_lhs = prev_lhs__; this->prev_rhs = prev_rhs__;
};

// other scphghe func
template <class V>
double check_convergence(solver<V>* solv__)
{
    make<double>::sp_mat lhs__ = solv__->lhs;
    Eigen::VectorXd rhs__ = solv__->rhs;
    Eigen::VectorXd cvalue__ = Eigen::VectorXd(solv__->rhs.rows());
    for(std::pair<std::string, make<double>::map_int> entry1 : solv__->eq->value->cvalue)
    {
        for(std::pair<int, double> entry2 : entry1.second)
        {
            cvalue__(entry2.first) = entry2.second;
        };
    };
    Eigen::VectorXd x = lhs__ * cvalue__;
    x = rhs__ - x;
    double res = 0.0;
    for(auto i : x)
    {
        res += pow(i, 2);
    };
    res = pow(res / rhs__.size(), 0.5);
    return res;
};

// scphghe
scphghe::scphghe(scheme*& scheme_ref, user*& user_ref, double step_length_in, double under_relax_in,
                 double min_residual_in, int max_iter_in)
{
    this->make_scphghe(scheme_ref, user_ref, step_length_in, under_relax_in, min_residual_in, max_iter_in);
};
void scphghe::make_scphghe(scheme*& scheme_ref, user*& user_ref, double step_length_in, double under_relax_in,
                           double min_residual_in, int max_iter_in)
{
    this->step_length = step_length_in;
    this->under_relax = under_relax_in;
    this->max_iter = max_iter_in;
    this->min_residual = min_residual_in;
    this->current_time = 0;
    momentum* u_ = new momentum(scheme_ref, 0);
    this->solv_u = new solver(u_, scheme_ref, step_length_in);
    momentum* v_ = new momentum(scheme_ref, 1);
    this->solv_v = new solver(v_, scheme_ref, step_length_in);
    momentum* w_ = new momentum(scheme_ref, 2);
    this->solv_w = new solver(w_, scheme_ref, step_length_in);
    pcorrect* pcor_ = new pcorrect(scheme_ref, user_ref);
    this->solv_pcor = new solver(pcor_, scheme_ref, step_length_in);
    turb_k* k_turb_ = new turb_k(scheme_ref);
    this->solv_k = new solver(k_turb_, scheme_ref, step_length_in);
    turb_e* e_turb_= new turb_e(scheme_ref);
    this->solv_e = new solver(e_turb_, scheme_ref, step_length_in);
    energy* energy_ = new energy(scheme_ref, user_ref);
    this->solv_energy = new solver(energy_, scheme_ref, step_length_in);
    s2s* s2s_ = new s2s(scheme_ref);
    this->solv_s2s = new solver(s2s_, scheme_ref, step_length_in);
};
int scphghe::SIMPLE_loop(scheme*& scheme_ref, bool is_init, make<double>::map_str* err_p__, make<double>::map_str* res_p__)
{
    // momentum - pcorrect loop
    momentum*& u_ = this->solv_u->eq;
    momentum*& v_ = this->solv_v->eq;
    momentum*& w_ = this->solv_w->eq;
    turb_k*& k_ = this->solv_k->eq;
    pcorrect*& pcor_ = this->solv_pcor->eq;
    int ctrl = 0;
    int passes = 0;
    double res_u; double res_v; double res_w; double res_pcor;
    while(ctrl < this->max_iter)
    {
        passes += 1;
        // momentum
        u_->update_linear(scheme_ref, k_, v_, w_);
        make_transient_lhs<momentum>(this->solv_u, this->under_relax);
        make_transient_rhs<momentum>(this->solv_u, this->under_relax);
        v_->update_linear(scheme_ref, k_, u_, w_);
        make_transient_lhs<momentum>(this->solv_v, this->under_relax);
        make_transient_rhs<momentum>(this->solv_v, this->under_relax);
        w_->update_linear(scheme_ref, k_, u_, v_);
        make_transient_lhs<momentum>(this->solv_w, this->under_relax);
        make_transient_rhs<momentum>(this->solv_w, this->under_relax);
        update_values<momentum>(this->solv_u, scheme_ref, std::string("u"), err_p__);
        res_u = check_convergence<momentum>(this->solv_u);
        update_values<momentum>(this->solv_v, scheme_ref, std::string("v"), err_p__);
        res_v = check_convergence<momentum>(this->solv_v);
        update_values<momentum>(this->solv_w, scheme_ref, std::string("w"), err_p__);
        res_w = check_convergence<momentum>(this->solv_w);
        u_->calc_wall(scheme_ref, v_, w_);
        v_->calc_wall(scheme_ref, u_, w_);
        w_->calc_wall(scheme_ref, u_, v_);
        // pcor
        pcor_->update_linear(scheme_ref, u_, v_, w_);
        make_transient_lhs<pcorrect>(this->solv_pcor, this->under_relax);
        make_transient_rhs<pcorrect>(this->solv_pcor, this->under_relax);
        update_values<pcorrect>(this->solv_pcor, scheme_ref, std::string("pcor"), err_p__);
        res_pcor = check_convergence<pcorrect>(this->solv_pcor);
        pcor_->update_correction(scheme_ref, u_, v_, w_);
        // new cvalue, old lhs - rhs
        if(res_u <= this->min_residual && res_v <= this->min_residual &&
           res_w <= this->min_residual && res_pcor <= this->min_residual)
        {
            if(passes < 2)
            {
                make_prev_lhs<momentum>(this->solv_u, scheme_ref, this->step_length, is_init);
                make_prev_lhs<momentum>(this->solv_v, scheme_ref, this->step_length, is_init);
                make_prev_lhs<momentum>(this->solv_w, scheme_ref, this->step_length, is_init);
                make_prev_lhs<pcorrect>(this->solv_pcor, scheme_ref, this->step_length, is_init);
                make_prev_rhs<momentum>(this->solv_u, scheme_ref, this->step_length, is_init);
                make_prev_rhs<momentum>(this->solv_v, scheme_ref, this->step_length, is_init);
                make_prev_rhs<momentum>(this->solv_w, scheme_ref, this->step_length, is_init);
                make_prev_rhs<pcorrect>(this->solv_pcor, scheme_ref, this->step_length, is_init);
                make<double>::map_str& res__ = *res_p__;
                res__["u"] = res_u; res__["v"] = res_v; res__["w"] = res_w; res__["pcor"] = res_pcor;
                return ctrl;
            }
            else
            {
                passes = 0;
                make_prev_lhs<momentum>(this->solv_u, scheme_ref, this->step_length, is_init);
                make_prev_lhs<momentum>(this->solv_v, scheme_ref, this->step_length, is_init);
                make_prev_lhs<momentum>(this->solv_w, scheme_ref, this->step_length, is_init);
                make_prev_lhs<pcorrect>(this->solv_pcor, scheme_ref, this->step_length, is_init);
                ctrl += 1;
            };
        };
    };
    return 0;
};
int scphghe::turb_loop(scheme*& scheme_ref, bool is_init, make<double>::map_str* err_p__, make<double>::map_str* res_p__)
{
    momentum*& u_ = this->solv_u->eq;
    momentum*& v_ = this->solv_v->eq;
    momentum*& w_ = this->solv_w->eq;
    turb_k*& k_ = this->solv_k->eq;
    turb_e*& e_ = this->solv_e->eq;
    int ctrl = 0;
    int passes = 0;
    double res_k; double res_e;
    while (ctrl < this->max_iter)
    {
        passes += 1;
        // turb_k
        k_->update_linear(scheme_ref, e_, u_, v_, w_);
        make_transient_lhs<turb_k>(this->solv_k, this->under_relax);
        make_transient_rhs<turb_k>(this->solv_k, this->under_relax);
        update_values<turb_k>(this->solv_k, scheme_ref, "k", err_p__);
        res_k = check_convergence<turb_k>(this->solv_k);
        k_->calc_wall(scheme_ref, e_);
        // turb_e
        e_->update_linear(scheme_ref, k_, u_, v_, w_);
        make_transient_lhs<turb_e>(this->solv_e, this->under_relax);
        make_transient_rhs<turb_e>(this->solv_e, this->under_relax);
        update_values<turb_e>(this->solv_e, scheme_ref, "e", err_p__);
        res_e = check_convergence<turb_e>(this->solv_e);
        e_->calc_wall(scheme_ref);
        if(res_k <= this->min_residual && res_e <= this->min_residual)
        {
            if(passes < 2)
            {
                make_prev_lhs<turb_k>(this->solv_k, scheme_ref, this->step_length, is_init);
                make_prev_lhs<turb_e>(this->solv_e, scheme_ref, this->step_length, is_init);
                make_prev_rhs<turb_k>(this->solv_k, scheme_ref, this->step_length, is_init);
                make_prev_rhs<turb_e>(this->solv_e, scheme_ref, this->step_length, is_init);
                make<double>::map_str& res__ = *res_p__;
                res__["k"] = res_k; res__["e"] = res_e;
                return ctrl;
            }
            else
            {
                passes = 0;
                make_prev_lhs<turb_k>(this->solv_k, scheme_ref, this->step_length, is_init);
                make_prev_lhs<turb_e>(this->solv_e, scheme_ref, this->step_length, is_init);
                ctrl += 1;
            };
        };
    };
    return 0;
};
int scphghe::energy_loop(scheme*& scheme_ref, double AH, bool is_init, make<double>::map_str* err_p__, make<double>::map_str* res_p__)
{
    // s2s and energy
    momentum*& u_ = this->solv_u->eq;
    momentum*& v_ = this->solv_v->eq;
    momentum*& w_ = this->solv_w->eq;
    turb_k*& k_ = this->solv_k->eq;
    energy*& energy_ = this->solv_energy->eq;
    s2s*& s2s_ = this->solv_s2s->eq;
    int ctrl = 0;
    int passes = 0;
    double res_s2s; double res_energy;
    while(ctrl < this->max_iter)
    {
        passes += 1;
        // s2s
        s2s_->update_linear(scheme_ref, energy_);
        make_transient_lhs<s2s>(this->solv_s2s, this->under_relax);
        make_transient_rhs<s2s>(this->solv_s2s, this->under_relax);
        update_values<s2s>(this->solv_s2s, scheme_ref, "s2s", err_p__);
        res_s2s = check_convergence<s2s>(this->solv_s2s);
        s2s_->update_source_s2s(scheme_ref);
        // energy
        energy_->update_linear(scheme_ref, k_, u_, v_, w_, AH);
        make_transient_lhs<energy>(this->solv_energy, this->under_relax);
        make_transient_rhs<energy>(this->solv_energy, this->under_relax);
        update_values<energy>(this->solv_energy, scheme_ref, "T", err_p__);
        res_energy = check_convergence<energy>(this->solv_energy);
        energy_->calc_wall(scheme_ref, AH);
        if(res_s2s <= this->min_residual && res_energy <= this->min_residual)
        {
            if(passes < 2)
            {
                make_prev_lhs<energy>(this->solv_energy, scheme_ref, this->step_length, is_init);
                make_prev_rhs<energy>(this->solv_energy, scheme_ref, this->step_length, is_init);
                make<double>::map_str& res__ = *res_p__;
                res__["T"] = res_energy; res__["s2s"] = res_s2s;
                return ctrl;
            }
            else
            {
                passes = 0;
                make_prev_lhs<energy>(this->solv_energy, scheme_ref, this->step_length, is_init);
                ctrl += 1;
            };
        };
    };
    return 0;
};
void scphghe::iterate(scheme*& scheme_ref, exports*& export_ref, user*& user_ref)
{
    bool is_init = false;
    if(this->current_time == 0)
    {
        bool is_init = true;
    };
    make<double>::map_str empty__;
    empty__.insert({"u", 0.0}); empty__.insert({"v", 0.0}); empty__.insert({"w", 0.0});
    empty__.insert({"k", 0.0}); empty__.insert({"e", 0.0}); empty__.insert({"T", 0.0});
    make<double>::map_str* err_p = new make<double>::map_str(empty__);
    make<double>::map_str* res_p = new make<double>::map_str(empty__);
    momentum*& u_ = this->solv_u->eq;
    momentum*& v_ = this->solv_v->eq;
    momentum*& w_ = this->solv_w->eq;
    turb_k*& k_ = this->solv_k->eq;
    turb_e*& e_ = this->solv_e->eq;
    energy*& energy_ = this->solv_energy->eq;
    int ctrl = 0;
    int passes = 0;
    user_ref->update_source(this->current_time, scheme_ref);
    double AH = user_ref->W_init;
    std::cout << "Start time step iteration..." << this->current_time << std::endl; 
    while(ctrl < this->max_iter)
    {
        auto start_timer = std::chrono::high_resolution_clock::now();
        passes += SIMPLE_loop(scheme_ref, is_init, err_p, res_p);
        std::cout << "SIMPLE loop done. " << "[" << ctrl << "]" << std::endl;
        passes += turb_loop(scheme_ref, is_init, err_p, res_p);
        std::cout << "Turbulence loop done. " << "[" << ctrl << "]" << std::endl;
        passes += energy_loop(scheme_ref, AH, is_init, err_p, res_p);
        std::cout << "Energy loop done. " << "[" << ctrl << "]" << std::endl;
        if(passes < 5)
        {
            // export solution at time t
            std::cout << "Converged. " << "[" << ctrl << "]" << std::endl;
            update_fluid_prop(scheme_ref, energy_, AH);
            update_wall(scheme_ref, u_, v_, w_, k_, e_, energy_);
            update_body_term(scheme_ref, u_, v_, w_);
            this->current_time += 1;
            // update export
            auto stop_timer = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<sec>(stop_timer - start_timer);
            long long time__ = duration.count();
            std::cout << "Exporting time step... " << this->current_time << "exec. time: " << time__ << " s" << std::endl;
            // export_ref.update_export(*this, scheme_ref, *err_p, *res_p, time__, ctrl);
            return;
        }
        else
        {
            std::cout << "Updating scheme... " << "[" << ctrl << "]" << std::endl;
            update_fluid_prop(scheme_ref, energy_, AH);
            update_wall(scheme_ref, u_, v_, w_, k_, e_, energy_);
            update_body_term(scheme_ref, u_, v_, w_);
            passes = 0;
            ctrl += 1;
        };
    };
};
