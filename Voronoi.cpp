//----------------------------------------------------------------------------------------------------------------------
/// @file Voronoi.cpp
/// @brief Class that has all elements of a single voronoi cell
//----------------------------------------------------------------------------------------------------------------------

#include "include/Voronoi.h"
#include "ngl/Random.h"

Voronoi::Voronoi()
{
    m_vertices.clear();
    m_tetrahedra.clear();
    m_faces.clear();
}

Voronoi::Voronoi(std::vector<Tetrahedron*> _t)
{
    m_tetrahedra = _t;
    setCenter();
    setEdge();
    createVAOLines();
}

void Voronoi::setCenter()
{
    for(unsigned int i=0; i<m_tetrahedra.size();++i)
    {
        m_vertices.push_back(m_tetrahedra[i]->getCirCenter());
    }
}

void Voronoi::setVertex(ngl::Vec3 _v)
{
    m_vertices.push_back(_v);
}


void Voronoi::setEdge()
{
    for(unsigned int i=0; i<m_tetrahedra.size();++i)
    {
       m_tetrahedra[i]->m_modified = true;
       for(int j=0;j<4;++j)
       {
          Tetrahedron *temp = m_tetrahedra[i]->m_neighbours[j];
          if(temp!=NULL)
          {
            if(temp->m_modified!=true)
            {
                  m_edges.push_back(m_vertices[i]);
                  m_edges.push_back(temp->getCirCenter());
            }
          }
       }
    }
}


void Voronoi::createVAOLines()
{

        // first we grab an instance of our VOA
        m_vaoMesh= ngl::VertexArrayObject::createVOA(GL_LINES);
        // next we bind it so it's active for setting data
        m_vaoMesh->bind();
        // now we have our data add it to the VAO, we need to tell the VAO the following
        // how much (in bytes) data we are copying
        // a pointer to the first element of data (in this case the address of the first element of the
        // std::vector
        m_vaoMesh->setData(m_edges.size()*sizeof(ngl::Vec3),m_edges[0].m_x);
        // so we need to set the vertexAttributePointer so the correct size and type as follows
        // vertex is attribute 0 with x,y,z(3) parts of type GL_FLOAT, our complete packed data is
        // sizeof(vertData) and the offset into the data structure for the first x component is 0
        m_vaoMesh->setVertexAttributePointer(0,3,GL_FLOAT,sizeof(ngl::Vec3),0);
        // now we have set the vertex attributes we tell the VAO class how many indices to draw when
        // glDrawArrays is called, in this case we use buffSize (but if we wished less of the sphere to be drawn we could
        // specify less (in steps of 3))
        m_vaoMesh->setNumIndices(m_edges.size());
        // finally we have finished for now so time to unbind the VAO
        m_vaoMesh->unbind();
        // indicate we have a vao now
        m_vao=true;

}

void Voronoi::draw()
{
    if(m_vao == true)
    {
     m_vaoMesh->bind();
     m_vaoMesh->draw();
     m_vaoMesh->unbind();
    }
}

