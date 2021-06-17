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


    double m_x, m_y, m_z;

};

class tri
{
    public:

    tri(){};
    tri(v3 p1, v3 p2, v3 p3);
    void draw();


    v3 m_p[3];

    v3 normal;
};

class model
{
    public:

    void load(char *fname);

    std::vector<tri> m_tris;
    void draw();
    void getMinMax();
    double m_xMin,m_xMax,m_yMin,m_yMax,m_zMin,m_zMax;
    model(){};


};

#endif // MODEL_H