// std lib
#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<string>
#include<utility>
#include<algorithm>
#include<memory>
#include<numeric>
#include<functional>
#include<chrono>
// third party
#include<mshio/mshio.h>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include"sqlite3.h"
#include"scphghe/temp_fwd_cfd.h"

using namespace Eigen;

namespace cfd
{
// UTIL CFD
typedef std::chrono::seconds sec;
typedef Eigen::Vector3d coor;
typedef Eigen::Triplet<double> sparse_input;
template<typename U>
struct make
{
    typedef std::vector<U> vec;
    typedef std::vector<std::vector<U>> std_mat;
    typedef std::map<int, U> map_int;
    typedef std::map<std::string, U> map_str;
    typedef Eigen::SparseMatrix<U, RowMajor> sp_mat;
    typedef std::pair<std::string, U> unique;
    typedef std::map<std::string, std::vector<U>> comp_str;
};
struct axes
{
    make<double>::sp_mat x;
    make<double>::sp_mat y;
    make<double>::sp_mat z;
    axes() {};
    axes(make<double>::sp_mat&);
    axes(const axes&);
    axes& operator()(int, int, coor&);
    axes& operator=(int);
    axes& operator+=(const axes&);
    coor axes_to_coor(int, int);
    double axes_to_val(int, int);
};

// STRUCT CFD
// user
class user
{
    public:
    double P_init;
    double T_init;
    double W_init;
    make<double>::map_str solid_k;
    make<double>::map_int s2s_eps;
    make<make<make<double>::map_int>::map_int>::map_str face_source; // .... - time step - value
    make<make<make<double>::map_int>::map_str>::map_str cell_source; // .... - time step - value
    user() {};
    user(double, double, double, std::string, std::string);
    void update_source(int, scheme*);
    private:
    void read_source_csv(std::string);
    void read_solid_prop_csv(std::string);
};

// scheme
class scheme
{
    public:
    minfo* mesh;
    pinfo* prop;
    winfo* wall;
    vinfo* pressure;
    binfo* source;
    make<double>::sp_mat rho_v_sf; // fluid, fc
    make<double>::map_int phi_v; // fluid, f
    scheme() {};
    scheme(std::string, user*);
};
struct finfo
{
    make<int>::vec fnode;
    coor fcentroid;
    coor fnormal;
    coor fparallel;
    double farea;
};
struct cinfo
{
    make<int>::vec cnode;
    make<int>::vec cface;
    coor ccentroid;
    double cvolume;
};
class minfo
{  
    public:
    make<coor>::map_int nodes;
    make<finfo>::map_int faces;
    make<cinfo>::map_int cells;
    make<make<axes>::map_str>::map_str geom; // Sf, Ef, Tf, eCf, eCF, dCf, dCF, parallel = -Sf.norm() from fluid perspective
    make<make<double>::map_int>::map_str size; // area, volume
    make<make<make<int>::vec>::map_int>::map_str fid; // store common (ns, inlet, outlet, conj) index 0, temp, flux, and s2s index unique id
    make<make<make<int>::vec>::map_str>::map_str cid; // store domain name unique id (derived from fluid/solid)
    make<make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>::map_str bid; // store cell-face-<boundary name, unique-constant/iter value>
    make<make<int>::sp_mat>::map_str cc_fc; // store neighbor cell to face alias id;
    make<make<double>::sp_mat>::map_str fc; // store face-cell sp_mat templates for fluid, solid, and conj (energy use)
    make<make<double>::sp_mat>::map_str cc; // store cell-cell sp_mat templates for fluid, solid, conj, and s2s (energy use)
    make<make<make<double>::sp_mat>::map_str>::map_str constants; // g_conv_aC, g_conv_aF, g_diff_aC, g_diff_aF, gc, view (cc, s2s)
    minfo() {};
    minfo(mshio::MshSpec);
    private:
    void make_minfo(mshio::MshSpec);
    void make_id(mshio::MshSpec&, make<cinfo>::map_int&);
    void make_template(make<coor>::map_int&, make<finfo>::map_int&, make<cinfo>::map_int&);
    void make_size(make<coor>::map_int&, make<finfo>::map_int&, make<cinfo>::map_int&);
    void make_geom(make<finfo>::map_int&, make<cinfo>::map_int&);
    void make_constants();
    void make_clust(make<finfo>::map_int&);
};
class pinfo
{
    public:
    make<make<double>::map_int>::map_str rho; // cell-face, fluid
    make<make<double>::map_int>::map_str miu; // cell-face, fluid 
    make<make<double>::map_int>::map_str cp; // cell-face, fluid
    make<make<double>::map_int>::map_str k; // cell-face, solid -> constant
    make<make<double>::map_int>::map_str eps; // face, solid(s2s -> glass/absorber, hamb -> glass/soil ids -> inline! (0.8)) -> constant
    pinfo() {};
    pinfo(minfo&, user*);
    private:
    void make_pinfo(minfo&, user*);
};
class winfo
{
    public:
    make<double>::map_int wall_dist; // dCf wall
    make<coor>::map_int wall_parallel; // cell parallel transformator
    make<double>::map_int ts; // cell, fluid
    make<double>::map_int utau; // cell, fluid
    make<double>::map_int miut; // cell, fluid
    winfo() {};
    winfo(minfo&);
    private:
    void make_winfo(minfo&);
};
struct vinfo
{
    make<make<double>::map_int>::map_str cvalue;
    make<make<double>::map_int>::map_str fvalue;
    make<make<double>::map_int>::map_str prev_cvalue;
    make<make<double>::map_int>::map_str prev_fvalue;
    make<make<coor>::map_int>::map_str cgrad;
    make<make<coor>::map_int>::map_str fgrad;
    make<make<coor>::map_int>::map_str prev_cgrad;
    make<make<coor>::map_int>::map_str prev_fgrad;
    vinfo() {};
    vinfo(make<std::string>::vec, minfo&, double);
    void make_vinfo(make<std::string>::vec, minfo&, double);
};
struct binfo
{
    /*momentum, pcorrect, turb_k, and turb_e*/ // -> fluid only
    // inlet value, environment velocity magnitude
    /*energy*/
    // isotemp(temp) value (unique), environment/ambient temperature
    // isoflux(flux) value (unique), environment irr q
    // s2s(flux) value (unique), s2s solver q // create inline program
    // hamb(hamb) value, ambient temperature
    make<make<double>::map_int>::map_str face_value; // boundary name-unique-constant/iter value
    make<make<double>::map_str>::map_str cell_value;
    binfo() {};
    binfo(make<make<double>::map_int>::map_str&, make<make<double>::map_str>::map_str&);
};
};