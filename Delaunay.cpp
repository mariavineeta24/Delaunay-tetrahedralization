//----------------------------------------------------------------------------------------------------------------------
/// @file Delaunay.cpp
/// @brief Class that has all functions necessary to compute the delaunay tetrahedralization
//----------------------------------------------------------------------------------------------------------------------

#include "Delaunay.h"
#include "ngl/Mat3.h"
#include "ngl/Mat4.h"
#include "Point4.h"

Delaunay::Delaunay(Tetrahedron* _tetrahedron)
{
    m_tetCount = 0;
    m_tetrahedron = _tetrahedron;
    m_tetrahedron->m_tetid = ++m_tetCount;
    m_tetrahedra.push_back(m_tetrahedron);
    m_third = new Tetrahedron();
    m_third = NULL;
    m_predicates = new Predicates();
    m_a = -1;
    m_b = -1;

}

Delaunay::~Delaunay()
{
    delete m_predicates;
    delete m_tetrahedron;
}

//----------------------------------------------------------------------------------------------------------------------
// This function computes the delaunay tetrahedralization by passing each point one at a time
//----------------------------------------------------------------------------------------------------------------------
std::vector<Tetrahedron*> Delaunay::compute(std::vector<ngl::Vec3> _points)
{

    ngl::Real tolerance = 0.000001;
    int flag = 0;
    Tetrahedron *oldT = m_tetrahedron;

    for(unsigned int i=0; i<_points.size(); ++i)
    {
        ngl::Vec3 point = _points[i];
        m_tetrahedron = walk(point,m_tetrahedron);

        // Degenerate case :Test if the point is exactly on the vertex
        for(int i=0; i<4; ++i)
        {
            double dist = distance(m_tetrahedron->getVertexData()[i],point);
            if(dist <= tolerance)
            {
                flag = 1;
                break;
            }
        }
        if(flag == 1)
        {
            flag = 0;
            continue;
        }
        m_tetrahedron = flip14(m_tetrahedron,point);
        Tetrahedron *tmp = new Tetrahedron();
        tmp = NULL;
        tmp = checkDelaunay();
        if(tmp!=NULL)
        {
            m_tetrahedron = tmp;
            tmp = NULL;
        }
    }

    // delete all tetrahedra that are modified
    for(unsigned int i=0; i<m_tetrahedra.size();++i)
    {
        if(m_tetrahedra[i]->m_modified == true)
        {
            m_tetrahedra.erase(m_tetrahedra.begin()+i);
            i--;
        }
        else
        {
            m_tetrahedra[i]->createVAO();
        }
    }

    // for deleting tetrahedra that contains the vertices of the big tetrahedron
    for(unsigned int i=0; i<m_tetrahedra.size(); ++i)
    {
        int flag = 0;
        for(int j=0; j<4; ++j)
        {
            for(int k=0; k<4; ++k)
            {
                if(m_tetrahedra[i]->getVertexData()[j] == oldT->getVertexData()[k])
                {
                    m_tetrahedra.erase(m_tetrahedra.begin()+i);
                    i--;
                    flag = 1;
                }
                if(flag == 1)
                {
                    break;
                }
            }
            if(flag == 1)
            {
                break;
            }
        }

    }

    return m_tetrahedra;
}

//----------------------------------------------------------------------------------------------------------------------
// This function computes the distance between two points
//----------------------------------------------------------------------------------------------------------------------
double Delaunay::distance(ngl::Vec3 _a, ngl::Vec3 _b)
{
    float deltaX = _b.m_x - _a.m_x;
    float deltaY = _b.m_y - _a.m_y;
    float deltaZ = _b.m_z - _a.m_z;

    double distance = std::sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);
    return distance;
}

//----------------------------------------------------------------------------------------------------------------------
// This function checks for delaunay criterion and performs the flips accordingly
//----------------------------------------------------------------------------------------------------------------------
Tetrahedron* Delaunay::checkDelaunay()
{
    int flipcase = -1;
    ngl::Vec3 p,d;
    Tetrahedron *t,*ta;
    Tetrahedron *next = new Tetrahedron();
    next = NULL;
    flipData f;
    int pid, did;

    while(!m_flipStack.empty())
    {
        // Find the next flip
        f = m_flipStack.top();
        m_flipStack.pop();

        if(f.m_t->m_modified != true)
        {
            // Set the tetrahedron and the point
            t = f.m_t;
            pid = f.m_ptPos;
            p = t->getVertexData()[pid];

            // Get the adjacent tetrahedron
            ta = t->m_neighbours[pid];

            // ta may be NULL for the first tetrahedron :
            if(ta != NULL)
            {
                // Get the apex of ta
                did = getApex(t,ta);
                d = ta->getVertexData()[did];

                // Check circumsphere check(inSphere) for d
                if(m_predicates->insphere3d(t,d) > 0)
                {
                    bool r = isCoplanar(t);
                    // Either apex(d) is non-existing,convex or concave or coplanar from p
                    flipcase = checkcase(t,ta,p,d);
                    if(flipcase == 3 && r == true)
                    {
                       //flipcase = 4;
                    }
                    switch(flipcase)
                    {
                        case 1 : // Both t and ta are convex
                                 next = flip23(t,pid,did);
                                 break;
                        case 2 : if(m_third != NULL)
                                 {
                                    next = flip32(t,ta,m_third,pid,did);
                                    m_third = NULL;
                                 }
                                 break;
                         case 3 : Tetrahedron* t3, *t4;
                                  if(config44(t,ta,p,m_cid,t3,t4) == true)
                                  {
                                    // The tetrahedra are in config44
                                    // next = flip44(t,ta,p,d,m_cid); - optional flip
                                  }
                                  break;
                         case 4 : next = flip23(t,pid,did);
                                  break;
                        default : std::cerr<<"Invalid flipcase!!!"<<std::endl;
                    }
                }
                else if(m_predicates->insphere3d(t,d)==0)
                {
                    // The point is on sphere.
                    // No special cases need to be done
                }
            else
            {
               // The point is outside sphere
            }
           }
        }
      }
    return next;
}

//----------------------------------------------------------------------------------------------------------------------
// This function performs Flip14
//----------------------------------------------------------------------------------------------------------------------
Tetrahedron* Delaunay::flip14(Tetrahedron* _tetra, ngl::Vec3 _point)
{
    std::vector<ngl::Vec3> verts = _tetra -> getVertexData();
    std::vector<ngl::Vec3> newverts1;
    newverts1.push_back(verts[0]);
    newverts1.push_back(verts[1]);
    newverts1.push_back(verts[2]);
    newverts1.push_back(_point);
    Tetrahedron *t1 = new Tetrahedron(newverts1);
    t1->m_tetid = ++m_tetCount;

    std::vector<ngl::Vec3> newverts2;
    newverts2.push_back(verts[0]);
    newverts2.push_back(verts[1]);
    newverts2.push_back(verts[3]);
    newverts2.push_back(_point);
    Tetrahedron *t2 = new Tetrahedron(newverts2);
    t2->m_tetid = ++m_tetCount;

    std::vector<ngl::Vec3> newverts3;
    newverts3.push_back(verts[0]);
    newverts3.push_back(verts[2]);
    newverts3.push_back(verts[3]);
    newverts3.push_back(_point);
    Tetrahedron *t3 = new Tetrahedron(newverts3);
    t3->m_tetid = ++m_tetCount;

    std::vector<ngl::Vec3> newverts4;
    newverts4.push_back(verts[1]);
    newverts4.push_back(verts[2]);
    newverts4.push_back(verts[3]);
    newverts4.push_back(_point);
    Tetrahedron *t4 = new Tetrahedron(newverts4);
    t4->m_tetid = ++m_tetCount;

    // setting up neighbours for t1
    t1->m_neighbours[0] = t4;
    t1->m_neighbours[1] = t3;
    t1->m_neighbours[2] = t2;
    t1->m_neighbours[3] = _tetra->m_neighbours[3];

    // setting up neighbours for t2
    t2->m_neighbours[0] = t4;
    t2->m_neighbours[1] = t3;
    t2->m_neighbours[2] = t1;
    t2->m_neighbours[3] = _tetra->m_neighbours[2];

    // setting up neighbours for t3
    t3->m_neighbours[0] = t4;
    t3->m_neighbours[1] = t2;
    t3->m_neighbours[2] = t1;
    t3->m_neighbours[3] = _tetra->m_neighbours[1];

    // setting up neighbours for t4
    t4->m_neighbours[0] = t3;
    t4->m_neighbours[1] = t2;
    t4->m_neighbours[2] = t1;
    t4->m_neighbours[3] = _tetra->m_neighbours[0];

    // updating nighbours of the new tetrahedra
    updateNeighbour(_tetra->m_neighbours[0],t4,_tetra);
    updateNeighbour(_tetra->m_neighbours[1],t3,_tetra);
    updateNeighbour(_tetra->m_neighbours[2],t2,_tetra);
    updateNeighbour(_tetra->m_neighbours[3],t1,_tetra);

    // Push new tetrahedra into the tetrahedra stack
    m_tetrahedra.push_back(t1);
    m_tetrahedra.push_back(t2);
    m_tetrahedra.push_back(t3);
    m_tetrahedra.push_back(t4);

    // Push new tetrahedra into the flip stack
    m_flipStack.push(createFlip(3,t1));
    m_flipStack.push(createFlip(3,t2));
    m_flipStack.push(createFlip(3,t3));
    m_flipStack.push(createFlip(3,t4));


    // Updating status of the old tetrahedra
    _tetra->m_modified = true;

    return t4;
}

//----------------------------------------------------------------------------------------------------------------------
// This function performs Flip44
//----------------------------------------------------------------------------------------------------------------------
Tetrahedron* Delaunay::flip44(Tetrahedron* _t1, Tetrahedron* _t2, ngl::Vec3 _p, ngl::Vec3 _d, int _cid)
{
    int pid = findPoint(_t1,_p);
    int did = findPoint(_t2,_d);
    ngl::Vec3 c = _t1->getVertexData()[_cid];
    int c2id = findPoint(_t2,c);
    Tetrahedron *t3 = _t1->m_neighbours[_cid];
    Tetrahedron *t4 = _t2->m_neighbours[c2id];
    findSharedEdge(_t1,_t2,t3);

    ngl::Vec3 v1 = _t1->getVertexData()[pid];
    ngl::Vec3 v2 = _t1->getVertexData()[(pid+1)%4];
    ngl::Vec3 v3 = _t1->getVertexData()[(pid+2)%4];
    ngl::Vec3 v4 = t3->getVertexData()[_cid];
    Tetrahedron *newt1 = new Tetrahedron(v1,v2,v3,v4);
    v1 = _t1->getVertexData()[pid];
    v2 = _t1->getVertexData()[(pid+1)%4];
    v3 = _t1->getVertexData()[(pid+3)%4];
    v4 = t3->getVertexData()[_cid];
    Tetrahedron *newt2 = new Tetrahedron(v1,v2,v3,v4);
    v1 = _t2->getVertexData()[did];
    v2 = _t2->getVertexData()[(did+1)%4];
    v3 = _t2->getVertexData()[(did+2)%4];
    v4 = t4->getVertexData()[c2id];
    Tetrahedron *newt3 = new Tetrahedron(v1,v2,v3,v4);
    v1 = _t2->getVertexData()[did];
    v2 = _t2->getVertexData()[(did+1)%4];
    v3 = _t2->getVertexData()[(did+3)%4];
    v4 = t4->getVertexData()[c2id];
    Tetrahedron *newt4 = new Tetrahedron(v1,v2,v3,v4);

    newt1->m_tetid = ++m_tetCount;
    newt2->m_tetid = ++m_tetCount;
    newt3->m_tetid = ++m_tetCount;
    newt4->m_tetid = ++m_tetCount;

    int a_t2 = findPoint(_t2,_t1->getVertexData()[m_a]);
    int b_t2 = findPoint(_t2,_t1->getVertexData()[m_b]);
    int a_t3 = findPoint(t3,_t1->getVertexData()[m_a]);
    int b_t3 = findPoint(t3,_t1->getVertexData()[m_b]);
    int a_t4 = findPoint(t4,_t1->getVertexData()[m_a]);
    int b_t4 = findPoint(t4,_t1->getVertexData()[m_b]);
    // neigbours for new t1
    newt1->m_neighbours[0] = newt3;
    newt1->m_neighbours[1] = t3->m_neighbours[b_t3];
    newt1->m_neighbours[2] = newt2;
    newt1->m_neighbours[3] = _t1->m_neighbours[m_b];

    // neigbours for new t2
    newt2->m_neighbours[0] = newt4;
    newt2->m_neighbours[1] = t3->m_neighbours[a_t3];
    newt2->m_neighbours[2] = newt1;
    newt2->m_neighbours[3] = _t1->m_neighbours[m_a];

    // neigbours for new t1
    newt3->m_neighbours[0] = newt1;
    newt3->m_neighbours[1] = t4->m_neighbours[b_t4];
    newt3->m_neighbours[2] = newt4;
    newt3->m_neighbours[3] = _t2->m_neighbours[b_t2];

    // neigbours for new t2
    newt4->m_neighbours[0] = newt2;
    newt4->m_neighbours[1] = t4->m_neighbours[a_t4];
    newt4->m_neighbours[2] = newt3;
    newt4->m_neighbours[3] = _t2->m_neighbours[a_t2];

    // update neighbours
    updateNeighbour(t3->m_neighbours[b_t3],newt1,t3);
    updateNeighbour(_t1->m_neighbours[m_b],newt1,_t1);
    updateNeighbour(t3->m_neighbours[a_t3],newt2,t3);
    updateNeighbour(_t1->m_neighbours[m_a],newt2,_t1);
    updateNeighbour(t4->m_neighbours[b_t4],newt3,t4);
    updateNeighbour(_t2->m_neighbours[b_t2],newt3,_t2);
    updateNeighbour(t4->m_neighbours[a_t4],newt4,t4);
    updateNeighbour(_t2->m_neighbours[a_t2],newt4,_t2);

    // Update flipstack
    m_flipStack.push(createFlip(0,newt1));
    m_flipStack.push(createFlip(0,newt2));
    m_flipStack.push(createFlip(0,newt3));
    m_flipStack.push(createFlip(0,newt4));

    // Push the new tetrahera created into tetrahedra stack
    m_tetrahedra.push_back(newt1);
    m_tetrahedra.push_back(newt2);
    m_tetrahedra.push_back(newt3);
    m_tetrahedra.push_back(newt4);

    _t1->m_modified = true;
    _t2->m_modified = true;
    t3->m_modified = true;
    t4->m_modified = true;

    return newt1;

}

//----------------------------------------------------------------------------------------------------------------------
// This function performs Flip23
//----------------------------------------------------------------------------------------------------------------------
Tetrahedron* Delaunay::flip23(Tetrahedron* _t, int _pid, int _did)
{
    ngl::Vec3 p = _t->getVertexData()[_pid];
    Tetrahedron *ta = _t->m_neighbours[_pid];
    ngl::Vec3 d = ta->getVertexData()[_did];

    int aid = (_pid+1)%4;
    int bid = (_pid+2)%4;
    int cid = (_pid+3)%4;

    ngl::Vec3 a = _t->getVertexData()[aid];
    ngl::Vec3 b = _t->getVertexData()[bid];
    ngl::Vec3 c = _t->getVertexData()[cid];

    int a2id = findPoint(ta,_t->getVertexData()[aid]);
    int b2id = findPoint(ta,_t->getVertexData()[bid]);
    int c2id = findPoint(ta,_t->getVertexData()[cid]);

    //create 3 new tetrahedra
    Tetrahedron *t1 = new Tetrahedron(p,a,b,d);
    Tetrahedron *t2 = new Tetrahedron(p,a,c,d);
    Tetrahedron *t3 = new Tetrahedron(p,b,c,d);

    t1->m_tetid = ++m_tetCount;
    t2->m_tetid = ++m_tetCount;
    t3->m_tetid = ++m_tetCount;

    //Setup neighbours for t1
    t1->m_neighbours[0] = ta->m_neighbours[c2id];
    t1->m_neighbours[1] = t3;
    t1->m_neighbours[2] = t2;
    t1->m_neighbours[3] = _t->m_neighbours[cid];


    //Setup neighbours for t2
    t2->m_neighbours[0] = ta->m_neighbours[b2id];
    t2->m_neighbours[1] = t3;
    t2->m_neighbours[2] = t1;
    t2->m_neighbours[3] = _t->m_neighbours[bid];

    //Setup neighbours for t3
    t3->m_neighbours[0] = ta->m_neighbours[a2id];
    t3->m_neighbours[1] = t2;
    t3->m_neighbours[2] = t1;
    t3->m_neighbours[3] = _t->m_neighbours[aid];

    //Setup neighbours of the neighbours
    updateNeighbour(_t->m_neighbours[aid],t3,_t);
    updateNeighbour(ta->m_neighbours[a2id],t3, ta);
    updateNeighbour(_t->m_neighbours[bid],t2,_t);
    updateNeighbour(ta->m_neighbours[b2id],t2,ta);
    updateNeighbour(_t->m_neighbours[cid],t1,_t);
    updateNeighbour(ta->m_neighbours[c2id],t1,ta);

    // Updating status of the two tetrahedra modified
     _t->m_modified = true;
     ta->m_modified = true;

     // Update flip stack
     m_flipStack.push(createFlip(0,t1));
     m_flipStack.push(createFlip(0,t2));
     m_flipStack.push(createFlip(0,t3));

     // Push the new tetrahera created into tetrahedra stack
     m_tetrahedra.push_back(t1);
     m_tetrahedra.push_back(t2);
     m_tetrahedra.push_back(t3);

     return t1;
}

//----------------------------------------------------------------------------------------------------------------------
// This function performs Flip32
//----------------------------------------------------------------------------------------------------------------------
Tetrahedron* Delaunay::flip32(Tetrahedron* _t1, Tetrahedron* _t2, Tetrahedron* _t3, int _pid, int _did)
{
    int c1 = -1;
    for(int i=0;i<4;++i)
    {
        if(i!=m_a && i!=_pid && i!=m_b)
        {
            c1 = i;
            break;
        }
    }
    if(c1 == -1)
    {
        std::cerr<<"Problem with c1!!!"<<std::endl;
    }

    // Create two new tetrahedra
    Tetrahedron* t1 = new Tetrahedron(_t1->getVertexData()[_pid], _t1->getVertexData()[m_a], _t1->getVertexData()[c1], _t2->getVertexData()[_did]);
    Tetrahedron* t2 = new Tetrahedron(_t1->getVertexData()[_pid], _t1->getVertexData()[m_b], _t1->getVertexData()[c1], _t2->getVertexData()[_did]);

    t1->m_tetid = ++m_tetCount;
    t2->m_tetid = ++m_tetCount;

    // Find indices for m_a in t2 and t3
    int a1 = findPoint(_t2,_t1->getVertexData()[m_a]);
    int a2 = findPoint(_t3,_t1->getVertexData()[m_a]);

    // Find indices for m_b in t2 and t3
    int b1 = findPoint(_t2,_t1->getVertexData()[m_b]);
    int b2 = findPoint(_t3,_t1->getVertexData()[m_b]);

    // Set neighbours for t1
    t1->m_neighbours[0] = _t2->m_neighbours[b1];
    t1->m_neighbours[1] = t2;
    t1->m_neighbours[2] = _t3->m_neighbours[b2];
    t1->m_neighbours[3] = _t1->m_neighbours[m_b];

    // Set neighbours for t2
    t2->m_neighbours[0] = _t2->m_neighbours[a1];
    t2->m_neighbours[1] = t1;
    t2->m_neighbours[2] = _t3->m_neighbours[a2];
    t2->m_neighbours[3] = _t1->m_neighbours[m_a];

    //Setup neighbours of the neighbours:
    updateNeighbour(_t1->m_neighbours[m_a],t2,_t1);
    updateNeighbour(_t1->m_neighbours[m_b],t1,_t1);
    updateNeighbour(_t2->m_neighbours[a1],t2,_t2);
    updateNeighbour(_t2->m_neighbours[b1],t1,_t2);
    updateNeighbour(_t3->m_neighbours[a2],t2,_t3);
    updateNeighbour(_t3->m_neighbours[b2],t1,_t3);

    _t1->m_modified = true;
    _t2->m_modified = true;
    _t3->m_modified = true;

    // Update flip stack
    m_flipStack.push(createFlip(0,t1));
    m_flipStack.push(createFlip(0,t2));

    // Push the new tetrahera created into tetrahedra stack
    m_tetrahedra.push_back(t1);
    m_tetrahedra.push_back(t2);

    return t1;

}

//----------------------------------------------------------------------------------------------------------------------
// This function checks if four tetrahedra are in config44
//----------------------------------------------------------------------------------------------------------------------
bool Delaunay::config44(Tetrahedron *_t1, Tetrahedron *_t2, ngl::Vec3 _p, int _cid, Tetrahedron* &_t3, Tetrahedron* &_t4)
{
    ngl::Vec3 c = _t1->getVertexData()[_cid];
    int c2id = findPoint(_t2,c);

    _t3 = _t1->m_neighbours[_cid];
    _t4 = _t2->m_neighbours[c2id];

    if(_t3!=NULL && _t4!=NULL)
    {
        int pid = findPoint(_t3,_p);
        if(_t3->m_neighbours[pid] == _t4)
        {
            return true;
        }
    }
    return false;
}

//----------------------------------------------------------------------------------------------------------------------
// This function retrieves the apex(vertex not shared) of the adjacent tetrahedra t
//----------------------------------------------------------------------------------------------------------------------
int Delaunay::getApex(Tetrahedron *_t1, Tetrahedron *_t2)
{
        //Checking for t1 between t2's neigbours
        for(int i=0;i<4;i++)
        {
                if(_t2->m_neighbours[i]==_t1)
                {
                        return i;
                }
        }
        return -1;
}

//----------------------------------------------------------------------------------------------------------------------
// This function checks if the tetrahedron vertices are coplanar
//----------------------------------------------------------------------------------------------------------------------
bool Delaunay::isCoplanar(Tetrahedron *_t)
{
    ngl::Vec3 a = _t->getVertexData()[0];
    ngl::Vec3 b = _t->getVertexData()[1];
    ngl::Vec3 c = _t->getVertexData()[2];
    ngl::Vec3 d = _t->getVertexData()[3];
    ngl::Mat4 m;
    m.m_00 = a.m_x;
    m.m_01 = a.m_y;
    m.m_02 = a.m_z;
    m.m_03 = 1;

    m.m_10 = b.m_x;
    m.m_11 = b.m_y;
    m.m_12 = b.m_z;
    m.m_13 = 1;

    m.m_20 = c.m_x;
    m.m_21 = c.m_y;
    m.m_22 = c.m_z;
    m.m_23 = 1;

    m.m_30 = d.m_x;
    m.m_31 = d.m_y;
    m.m_32 = d.m_z;
    m.m_33 = 1;

    double det = m.determinant();
    if(det == 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

//----------------------------------------------------------------------------------------------------------------------
// This function checks if a third tetrahedron exists, given t1 and t2
//----------------------------------------------------------------------------------------------------------------------
Tetrahedron* Delaunay::findthird(Tetrahedron *_t1, Tetrahedron *_t2, int _c1)
{
        int c2= findPoint(_t2, _t1->getVertexData()[_c1]);

        if(_t1->m_neighbours[_c1]==_t2->m_neighbours[c2])
        {
             return _t1->m_neighbours[_c1];
        }
        else
        {
            return NULL;
        }
}

Tetrahedron* Delaunay::findthird(Tetrahedron *_t1, Tetrahedron *_t2)
{
    for(int i=0;i<4;++i)
    {
        for(int j=0;j<4;++j)
        {
            if(_t1->m_neighbours[i]==_t2->m_neighbours[j])
            {
                return _t1->m_neighbours[i];
            }
        }
    }
    return NULL;
}

//----------------------------------------------------------------------------------------------------------------------
// This function finds the edge shared between two tetrahedron t1 and t2
//----------------------------------------------------------------------------------------------------------------------
void Delaunay::findSharedEdge(Tetrahedron *_t1, Tetrahedron *_t2, Tetrahedron* _t3)
{
    bool flag = true;
    int count = 0;
    for(int i=0;i<4;++i)
    {
       for(int j=0;j<4;++j)
        {
            if(_t1->getVertexData()[i] == _t2->getVertexData()[j])
            {
                for(int k=0;k<4;++k)
                {

                    if(_t1->getVertexData()[i] == _t3->getVertexData()[k])
                    {
                        if(count == 0)
                        {
                            m_a = i;
                            count++;
                            flag = true;
                            break;
                        }
                        else
                        {
                            m_b = i;
                            flag = true;
                            break;
                        }

                    }
                }
                if(flag == true)
                {
                    flag = false;
                    break;
                }
            }
        }
    }
}
//----------------------------------------------------------------------------------------------------------------------
// This function checks if the union of the two tetrahedron t1 and t2 are concave, convex or coplanar
//----------------------------------------------------------------------------------------------------------------------
int Delaunay::checkcase(Tetrahedron* _t1,Tetrahedron* _t2, ngl::Vec3 _point, ngl::Vec3 _d)
{
   bool concave = false;
   int case1, case2, case3;
   int pid = findPoint(_t1,_point);
   m_cid = -1;

   int aid = (pid+1)%4;
   int bid = (pid+2)%4;
   int cid = (pid+3)%4;
   ngl::Vec3 a =  _t1->getVertexData()[aid];
   ngl::Vec3 b = _t1->getVertexData()[bid];
   ngl::Vec3 c = _t1->getVertexData()[cid];

   case1 = oppSides(_point,a,b,c,_d);
   case2 = oppSides(_point,a,c,b,_d);
   case3 = oppSides(_point,b,c,a,_d);

   if(case1 == 2)
   {
       m_third = findthird(_t1,_t2,cid);
       m_a = aid;
       m_b = bid;
       concave = true;
       if(m_third!=NULL)
       {
           return 2;
       }
   }
   if(case2 == 2)
   {
       m_third = findthird(_t1,_t2,bid);
       m_a = aid;
       m_b = cid;
       concave = true;
       if(m_third!=NULL)
       {
           return 2;
       }
   }
   if(case3 == 2)
   {
       m_third = findthird(_t1,_t2,aid);
       m_a = bid;
       m_b = cid;
       concave = true;
       if(m_third!=NULL)
       {
           return 2;
       }
   }
   if(concave)
   {
       return 2;
   }
   if(case1==3)
   {
       m_cid = cid;
       return 3;
   }
   if(case2==3)
   {
       m_cid = bid;
       return 3;
   }
   if(case3==3)
   {
       m_cid = aid;
       return 3;
   }
   else
   {
       return 1;
   }

}

//----------------------------------------------------------------------------------------------------------------------
// This function finds the tetrahedron _findT in _oldT and updates it with _newT
//----------------------------------------------------------------------------------------------------------------------
void Delaunay::updateNeighbour(Tetrahedron *_oldt, Tetrahedron *_newt, Tetrahedron *_findt)
{

        int x=0;
        if(_oldt!=NULL){

                for(int i=0; i<4; i++){
                        if(_oldt->m_neighbours[i]==_findt){
                                _oldt->m_neighbours[i]=_newt;
                                x=1;
                        }
                }
                if(x==0){
                        std::cerr<<"Error : Neighbour not found!"<<std::endl;
                }
        }
}

//----------------------------------------------------------------------------------------------------------------------
// This function returns the position of the given point in the tetrahedron vertices
//----------------------------------------------------------------------------------------------------------------------
int Delaunay::findPoint(Tetrahedron* _t, ngl::Vec3 _point)
{
    for(int i=0; i<4; ++i)
    {
        if(_t->getVertexData()[i] == _point)
        {
            return i;
        }
    }
    return -1;
}

//----------------------------------------------------------------------------------------------------------------------
// This function updates the data required for flipping
//----------------------------------------------------------------------------------------------------------------------
flipData Delaunay::createFlip(int _pid, Tetrahedron* _t)
{
    flipData f;
    f.m_ptPos = _pid;
    f.m_t = _t;
    return f;
}

//----------------------------------------------------------------------------------------------------------------------
// This function finds if two points p1 and p2 are on either side of a plane defined by a,b,c
//----------------------------------------------------------------------------------------------------------------------
int Delaunay::oppSides(ngl::Vec3 _a, ngl::Vec3 _b, ngl::Vec3 _c, ngl::Vec3 _p1, ngl::Vec3 _p2)
{
    unsigned int result;
    double a,b;
    a = m_predicates->orient3d(_a,_b,_c,_p1);
    b = m_predicates->orient3d(_a,_b,_c,_p2);

    if((a>0 && b<0) || (a<0 && b>0))
    {
        result = 2;
    }
    if((a>0 && b>0) || (a<0 && b<0))
    {
        result = 1;
    }
    if(a==0 || b==0)
    {
        result = 3;
    }

    return result;
}

//----------------------------------------------------------------------------------------------------------------------
// This function finds the tetrahedron that contains p (WALK Algortihm)
//----------------------------------------------------------------------------------------------------------------------
Tetrahedron* Delaunay::walk(ngl::Vec3 _p, Tetrahedron *_tetra)
{
    float op1, op2;
    bool next=false;
    while(true)
    {
        for(unsigned int i=0; i<4; ++i)
        {
            op1 = m_predicates->orient3d(_tetra->getVertexData()[(i+1)%4],_tetra->getVertexData()[(i+2)%4],_tetra->getVertexData()[(i+3)%4],_p);
            op2 = m_predicates->orient3d(_tetra->getVertexData()[(i+1)%4],_tetra->getVertexData()[(i+2)%4],_tetra->getVertexData()[(i+3)%4],_tetra->getVertexData()[i]);

            if((op1>0 && op2<0) || (op1<0 && op2>0))
            {
                next = true;
                _tetra = _tetra->m_neighbours[i];
                break;
            }
        }
        if(next)
        {
            //Next Tetrahedron found
            next=false;
        }
        else
        {
            if(_tetra==NULL)
            {
               std::cerr<<"Tetrahedron not found!!!"<<std::endl;
            }
            return _tetra;
        }
    }
}





