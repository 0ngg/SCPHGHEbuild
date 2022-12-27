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
void user::update_source(int current_time, scheme* scheme_ref)
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
int is_neighbor(make<int>::vec v1, make<int>::vec v2)
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
make<std::pair<std::string, int>>::vec is_boundary(int fid, make<make<make<int>::vec>::map_int>::map_str fid_, make<std::string>::vec reserve)
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
void make_tri(finfo face_, make<coor>::map_int nodes_)
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
void make_quad(finfo face_, make<coor>::map_int nodes_)
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
void make_cell(minfo mesh_, cinfo cell_, make<finfo>::map_int faces_)
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
// scheme
scheme::scheme(std::string msh_file, user* user_p)
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
                cell__.cnode = cnode__; cell__.cface = cface__; 
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
pinfo::pinfo(minfo& mesh_, user* user_p)
{
    this->make_pinfo(mesh_, user_p);
};
void pinfo::make_pinfo(minfo& mesh_, user* user_p)
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
