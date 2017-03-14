//----------------------------------------------------------------------------------------------------------------------
/// @file Tetrahedron.cpp
/// @brief Class that defines the tetrahedron data-structure
//----------------------------------------------------------------------------------------------------------------------

#include "Tetrahedron.h"
#include "ngl/Random.h"
#include "ngl/Mat4.h"
#include "cstdlib"

Tetrahedron::Tetrahedron()
{
    m_verts.clear();
    m_modified = false;
    m_tetid = 0;
    m_vao = false;
}

Tetrahedron::Tetrahedron(ngl::Vec3 _a, ngl::Vec3 _b, ngl::Vec3 _c, ngl::Vec3 _d)
{
    m_verts.push_back(_a);
    m_verts.push_back(_b);
    m_verts.push_back(_c);
    m_verts.push_back(_d);
    m_modified = false;
    m_tetid = 0;
    m_vao = false;
    calculateFormula();
    findCenter();
    findRadius();
    createFaces();
    setColour();
}

Tetrahedron::Tetrahedron(std::vector<ngl::Vec3> _vertexArray)
{
    m_verts = _vertexArray;
    m_modified = false;
    m_tetid = 0;
    m_vao = false;
    calculateFormula();
    findCenter();
    findRadius();
    createFaces();
    setColour();
}

void Tetrahedron::setColour()
{
    ngl::Random *rng=ngl::Random::instance();
    m_tetraColour = rng->getRandomColour();
}

void Tetrahedron::calculateFormula()
{
    ngl::Vec3 a = m_verts[1] - m_verts[0];
    ngl::Vec3 b = m_verts[2] - m_verts[0];
    ngl::Vec3 c = m_verts[3] - m_verts[0];
    m_formulaN = (a.dot(a)*(b.cross(c)) + b.dot(b)*(c.cross(a)) + c.dot(c)*(a.cross(b)));
    m_formulaD = a.dot(b.cross(c));
}

void Tetrahedron::createFaces()
{
    Face f1,f2,f3,f4;
    f1.m_vert.push_back(0);
    f1.m_vert.push_back(1);
    f1.m_vert.push_back(2);
    m_face.push_back(f1);

    f2.m_vert.push_back(0);
    f2.m_vert.push_back(2);
    f2.m_vert.push_back(3);
    m_face.push_back(f2);

    f3.m_vert.push_back(0);
    f3.m_vert.push_back(3);
    f3.m_vert.push_back(1);
    m_face.push_back(f3);

    f4.m_vert.push_back(1);
    f4.m_vert.push_back(2);
    f4.m_vert.push_back(3);
    m_face.push_back(f4);
}

void Tetrahedron::findCenter()
{
    m_circumcenter = m_verts[0] + m_formulaN/(2*m_formulaD);
}

void Tetrahedron::findRadius()
{
    ngl::Vec3 a = m_verts[1] - m_verts[0];
    ngl::Vec3 b = m_verts[2] - m_verts[0];
    ngl::Vec3 c = m_verts[3] - m_verts[0];

    m_circumradius = m_formulaN.length()/(2*std::abs(m_formulaD));
}

Tetrahedron::~Tetrahedron()
{
    m_face.erase(m_face.begin(),m_face.end());
    glDeleteBuffers(1,&m_vboBuffers);
    if(m_vaoTetrahedron!=0)
    {
       delete m_vaoTetrahedron;
    }
}

void Tetrahedron::draw()
{
   if(m_vao == true)
   {
     m_vaoTetrahedron->bind();
     m_vaoTetrahedron->draw();
     m_vaoTetrahedron->unbind();
   }
}

void Tetrahedron::createVAO()
{
     // if we have already created a VBO just return.
     if(m_vao == true)
     {
         std::cout<<"VAO exist so returning\n";
         return;
     }

     // now we are going to process and pack the mesh into an ngl::VertexArrayObject
     std::vector <ngl::Vec3> vboMesh;
     ngl::Vec3 d;
     int loopFaceCount=3;

     // loop for each of the faces
     for(unsigned int i=0;i<4;++i)
     {
        // now for each triangle in the face (remember we ensured tri above)
        for(int j=0;j<loopFaceCount;++j)
        {
          // pack in the vertex data first
          d.m_x=m_verts[m_face[i].m_vert[j]].m_x;
          d.m_y=m_verts[m_face[i].m_vert[j]].m_y;
          d.m_z=m_verts[m_face[i].m_vert[j]].m_z;
          vboMesh.push_back(d);
        }
     }

     // first we grab an instance of our VOA
     m_vaoTetrahedron= ngl::VertexArrayObject::createVOA(GL_TRIANGLES);
     // next we bind it so it's active for setting data
     m_vaoTetrahedron->bind();
     m_tetrahedronSize=vboMesh.size();
     // now we have our data add it to the VAO, we need to tell the VAO the following
     // how much (in bytes) data we are copying
     // a pointer to the first element of data (in this case the address of the first element of the
     // std::vector
     m_vaoTetrahedron->setData(m_tetrahedronSize*sizeof(ngl::Vec3),vboMesh[0].m_x);
     // so we need to set the vertexAttributePointer so the correct size and type as follows
     // vertex is attribute 0 with x,y,z(3) parts of type GL_FLOAT, our complete packed data is
     // sizeof(vertData) and the offset into the data structure for the first x component is 0
     m_vaoTetrahedron->setVertexAttributePointer(0,3,GL_FLOAT,sizeof(ngl::Vec3),0);
     // now we have set the vertex attributes we tell the VAO class how many indices to draw when
     // glDrawArrays is called, in this case we use buffSize (but if we wished less of the sphere to be drawn we could
     // specify less (in steps of 3))
     m_vaoTetrahedron->setNumIndices(m_tetrahedronSize);
     // finally we have finished for now so time to unbind the VAO
     m_vaoTetrahedron->unbind();
     // indicate we have a vao now
     m_vao=true;
}



