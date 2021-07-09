#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include <string>
class v3
{
    public:

    v3(){};
    v3(char* bin);
    v3(double x, double y, double z);
    double len();

    double m_x, m_y, m_z;


};

double dotProd(v3 vec1, v3 vec2);

v3 operator+(v3 left, v3 right);
v3 operator-(v3 left, v3 right);
v3 operator/(v3 left, double right);
v3 operator*(v3 left, double right);

class tri
{
    public:

    tri(){};
    tri(v3 p1, v3 p2, v3 p3);
    void draw();

    double distP(v3 &point, double &dot_out);
    double distP_naive(v3 point);
    double getSign(v3 point); //get side with respect to the normal


    v3 m_p[3];
    double m_vert_angle[3];

    v3 shP[3];

    v3 normal; //normal vector
    v3 center; //point of triangle centroid

    double **matrix; //transformation matrix (translation + rotation)
};

struct i3
{
public:
    int m_i[3];
};

#define NI_max 30
#define NJ_max 30
#define NK_max 30
class model
{
    public:

    void load(char *fname);

    std::vector<tri> m_tris;
    void draw();
    double distP(v3 point);
    double distP_naive(v3 point);
     double distP_fast(v3 point);

    void getMinMax();
    double m_xMin,m_xMax,m_yMin,m_yMax,m_zMin,m_zMax;

    static std::vector<int> m_grid[NI_max][NJ_max][NK_max]; //triangles in cells
    static i3 m_nearest_cell[NI_max][NJ_max][NK_max];// nearest non-empty cell to the current cell

    int m_ni,m_nj,m_nk;
    double w_x0,w_y0,w_z0,w_x1,w_y1,w_z1;
    double m_dx,m_dy,m_dz;

    void distributeTriangles();
    model(){};
    i3 getNearest(int i, int j, int k);

    int is_inside(int i,int j,int k);

    int cell_tri_overlap(int i, int j, int k, tri &trian);
};

#endif // MODEL_H
