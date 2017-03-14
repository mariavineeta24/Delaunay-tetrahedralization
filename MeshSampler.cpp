//----------------------------------------------------------------------------------------------------------------------
/// @file MeshSampler.cpp
/// @brief Class that has all the functions necessary for sampling the mesh
//----------------------------------------------------------------------------------------------------------------------

#include "MeshSampler.h"
#include "ngl/Random.h"
#include "ngl/ShaderLib.h"
#include "ngl/Util.h"
#include "ngl/VAOPrimitives.h"
#include "include/sdf/signed_distance_field_from_mesh.hpp"

MeshSampler::MeshSampler()
{
    m_surfacePointsBBox.clear();
    m_surfacePointsMesh.clear();
    m_rayStart.clear();
    m_rayEnd.clear();
    m_volumePoints.clear();
    m_vertTri.clear();
    m_tetrahedra.clear();
    m_objfilename = "models/prism.obj";
    m_density = 1;
    m_ptLocation = 0;
    m_method = 0;
}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::setObjFilename(
                                    const std::string _filename
                                )
{
    m_objfilename = _filename;
}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::setDensity(
                                int _density
                            )
{
    m_density = _density;
}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::setPointLocationType(
                                        int _type
                                      )
{
    m_ptLocation = _type;
}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::setMethod(
                             int _type
                           )
{
    m_method = _type;
}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::loadMesh()
{
    // load the mesh
    m_mesh = new ngl::Obj(m_objfilename);
    m_mesh->createVAO();
    // checking if the object is triangulated
    if(!m_mesh->isTriangular())
    {
        std::cout<<"Works only for triangulated meshes!!!"<<std::endl;
        exit(EXIT_FAILURE);
    }
    ngl::VAOPrimitives *prim = ngl::VAOPrimitives::instance();
    prim->createSphere("sphere",0.01,10);
}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::findPointsOnSurfaceBBox()
{
    // Obtaining list of vertices from Bbox
    ngl::Vec3 *vertex_array = m_mesh->getBBox().getVertexArray();

    std::vector<Point3> tri_bbox;
    Point3 p1(vertex_array[0],vertex_array[1],vertex_array[3]);
    tri_bbox.push_back(p1);
    Point3 p2(vertex_array[2],vertex_array[3],vertex_array[1]);
    tri_bbox.push_back(p2);
    Point3 p3(vertex_array[4],vertex_array[5],vertex_array[7]);
    tri_bbox.push_back(p3);
    Point3 p4(vertex_array[6],vertex_array[7],vertex_array[5]);
    tri_bbox.push_back(p4);

    Point3 p5(vertex_array[2],vertex_array[1],vertex_array[6]);
    tri_bbox.push_back(p5);
    Point3 p6(vertex_array[5],vertex_array[6],vertex_array[1]);
    tri_bbox.push_back(p6);
    Point3 p7(vertex_array[0],vertex_array[4],vertex_array[3]);
    tri_bbox.push_back(p7);
    Point3 p8(vertex_array[3],vertex_array[7],vertex_array[4]);
    tri_bbox.push_back(p8);

    Point3 p9(vertex_array[1],vertex_array[5],vertex_array[0]);
    tri_bbox.push_back(p9);
    Point3 p10(vertex_array[4],vertex_array[0],vertex_array[5]);
    tri_bbox.push_back(p10);
    Point3 p11(vertex_array[6],vertex_array[7],vertex_array[2]);
    tri_bbox.push_back(p11);
    Point3 p12(vertex_array[3],vertex_array[2],vertex_array[7]);
    tri_bbox.push_back(p12);

    for(unsigned int i=0;i<tri_bbox.size();++i)
    {
        for(int j=0;j<m_density;++j)
        {
            // Generating random points on the surface of the bbox
            ngl::Vec3 A,B,C;
            A = tri_bbox[i].m_a;
            B = tri_bbox[i].m_b;
            C = tri_bbox[i].m_c;

            ngl::Random *rng=ngl::Random::instance();
            float a,b,c;
            a = rng->randomPositiveNumber(1);
            b = rng->randomPositiveNumber(1);
            if(a+b > 1)
            {
               a = 1 - a;
               b = 1 - b;
            }
            c = 1 - a - b;
            ngl::Vec3 point;
            point.m_x = (a*A.m_x) + (b*B.m_x) + (c*C.m_x);
            point.m_y = (a*A.m_y) + (b*B.m_y) + (c*C.m_y);
            point.m_z = (a*A.m_z) + (b*B.m_z) + (c*C.m_z);
            m_surfacePointsBBox.push_back(point);
        }
    }
}

void MeshSampler::findPointsOnSurfaceMesh( Point3 _triangle )
{
    // Generating random points on the surface of the bbox
    ngl::Vec3 A,B,C;
    A = _triangle.m_a;
    B = _triangle.m_b;
    C = _triangle.m_c;
    for(int i=0;i<m_density;++i)
    {
        ngl::Random *rng=ngl::Random::instance();
        float a,b,c;
        a = rng->randomPositiveNumber(1);
        b = rng->randomPositiveNumber(1);
        if(a+b > 1)
        {
           a = 1 - a;
           b = 1 - b;
        }
        c = 1 - a - b;
        ngl::Vec3 point;
        point.m_x = (a*A.m_x) + (b*B.m_x) + (c*C.m_x);
        point.m_y = (a*A.m_y) + (b*B.m_y) + (c*C.m_y);
        point.m_z = (a*A.m_z) + (b*B.m_z) + (c*C.m_z);
        m_surfacePointsMesh.push_back(point);
    }
}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::SampleMesh_usingSDF()
{
    m_volumePoints.clear();
    int no_of_points = m_density * 20;
    // Signed distance field
    sdf::signed_distance_field_from_mesh s;
    s.load_from_file(m_objfilename);
    s.is_valid();
    const float* max_bound = s.maximum_bound();
    const float* min_bound = s.minimum_bound();

    // generating random points within the bounds of bounding box
    for(int i=0;i<no_of_points;++i)
    {
        ngl::Random *randomNumberGenerator = ngl::Random::instance();
        ngl::Vec3 pos(0.0,0.0,0.0);
        pos = randomNumberGenerator->getRandomVec3();
        float newX = min_bound[0] + ( (max_bound[0] - min_bound[0]) * (pos[0] + 1.f) / 2.f );
        float newY = min_bound[1] + ( (max_bound[1] - min_bound[1]) * (pos[1] + 1.f) / 2.f );
        float newZ = min_bound[2] + ( (max_bound[2] - min_bound[2]) * (pos[2] + 1.f) / 2.f );
        pos[0] = newX;
        pos[1] = newY;
        pos[2] = newZ;
        m_volumePoints.push_back(pos);
    }

    // rejecting points based on SDF
    for(unsigned int i=0; i<m_volumePoints.size(); ++i)
    {
        ngl::Vec3 point = m_volumePoints[i];
        if(s(point.m_x,point.m_y,point.m_z)>=0)
        {
           m_volumePoints.erase(m_volumePoints.begin()+i);
           i--;
         }
    }

}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::SampleMesh()
{
    m_surfacePointsBBox.clear();
    m_surfacePointsMesh.clear();
    m_rayStart.clear();
    m_rayEnd.clear();
    m_volumePoints.clear();
    m_vertTri.clear();
    m_hitPoints.clear();

    // Obtaining mesh data
    std::vector <ngl::Vec3> verts = m_mesh->getVertexList();
    std::vector <ngl::Face> faces = m_mesh->getFaceList();

    ngl::Vec3 d;
    unsigned int nFaces = faces.size();
    // looping through every face
    for(unsigned int i=0;i<nFaces;++i)
    {
        // looping trough each triangle in the face
        for(int j=0;j<3;++j)
        {
            d.m_x = verts[faces[i].m_vert[j]].m_x;
            d.m_y = verts[faces[i].m_vert[j]].m_y;
            d.m_z = verts[faces[i].m_vert[j]].m_z;
            m_vertTri.push_back(d);
        }

     }

    if(m_ptLocation != 0)
    {
        for(unsigned int i=0; i<m_vertTri.size(); i+=3)
        {
            findPointsOnSurfaceMesh(Point3(m_vertTri[i], m_vertTri[i+1], m_vertTri[i+2]));
        }
        m_points = m_surfacePointsMesh;
    }
    else
    {
        if(m_method == 0)
        {
            // using SDF
            SampleMesh_usingSDF();

        }
        else
        {
            // Ray intersection
            SampleMesh_rayIntersection();
        }
        m_points = m_volumePoints;
    }

}

void MeshSampler::SampleMesh_rayIntersection()
{
    findPointsOnSurfaceBBox();

        int size;
        size = 2 * m_density; // number of tris(4) * density
        for(int i=0; i<size; ++i)
        {
                    m_rayStart.push_back(m_surfacePointsBBox[i]);
                    m_rayEnd.push_back(m_surfacePointsBBox[i+size]);

        }
        for(int i=2*size; i<3*size;++i)
        {
                    m_rayStart.push_back(m_surfacePointsBBox[i]);
                    m_rayEnd.push_back(m_surfacePointsBBox[i+size]);

        }
        for(int i=4*size; i<5*size;++i)
        {
                    m_rayStart.push_back(m_surfacePointsBBox[i]);
                    m_rayEnd.push_back(m_surfacePointsBBox[i+size]);
        }

        rayTriangleIntersect();
        double length;
        int noOfPoints;
        for(unsigned int i=0;i<m_hitPoints.size();i+=2)
        {
            ngl::Vec3 diff = m_hitPoints[i] - m_hitPoints[i+1];
            length = sqrt(pow(diff.m_x,2) + pow(diff.m_y,2) + pow(diff.m_z,2));
            noOfPoints = ceil(length) * 2;
            for(int j=0;j<noOfPoints;++j)
            {
                ngl::Vec3 points;
                float t = j/(float)noOfPoints;
                points = ngl::lerp(m_hitPoints[i],m_hitPoints[i+1],t);
                m_volumePoints.push_back(points);
            }
        }

}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::rayTriangleIntersect()
{
    ngl::Vec3 tvec, pvec, qvec;
    float det, inv_det;

    // Calculating ray direction
    for(unsigned int i=0;i<m_rayStart.size();++i) // check for all rays
    {
        int count = 0;
        for(unsigned int j=0;j<m_vertTri.size();j+=3) // check a single ray with every triangle
        {
            ngl::Vec3 edge1, edge2;
            edge1 = m_vertTri[j+1] - m_vertTri[j];
            edge2 = m_vertTri[j+2] - m_vertTri[j];

            // ray direction vector
            ngl::Vec3 dir = m_rayStart[i] - m_rayEnd[i];
            pvec = dir.cross(edge2);
            det = edge1.dot(pvec);
            // if this is 0 no hit
            if (det > -0.00001f && det < 0.00001)
            {
                continue;
            }
            // get the inverse det
            inv_det = 1.0f / det;
            // calculate the 2nd vector
            tvec = m_rayStart[i] - m_vertTri[j];
            // get the dot product of this and inv det
            ngl::Real u,v,w;
            u = tvec.dot(pvec) * inv_det;
            // if out of range no hit
            if (u < -0.001f || u > 1.001f)
            {
                continue;
            }
            // check the 2nd vector edge
            qvec = tvec.cross(edge1);
            // get the dot product
            v = dir.dot(qvec) * inv_det;
            // if out of range no hit
            if (v < -0.001f || u + v > 1.001f)
            {
                continue;
            }
            // check the final value
            w = edge2.dot(qvec) * inv_det;
            // if greater than 0 no hit
            if (w >= 0)
            {
                continue;
            }

                // otherwise we are inside the triangle
                // so get the hit point
                // Reference : http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm#intersect_RayTriangle()
                // get intersect point of ray with triangle plane
                // calculate the normal
                ngl::Vec3 n = ngl::calcNormal(m_vertTri[j],m_vertTri[j+1],m_vertTri[j+2]);
                ngl::Vec3 hitPoint;
                float a = -n.dot(tvec);
                float b = n.dot(dir);
                float r=a/b;
                // intersect point of ray and plane
                hitPoint=m_rayStart[i] + r * dir;
                m_hitPoints.push_back(hitPoint);
                count++;
        }

        if(count % 2 != 0) // odd surfaces => ray doesnt intersect the mesh
        {
            for(int i=0;i<count;++i)
            {
                m_hitPoints.pop_back();
            }
        }
        else
        {
            if(count > 2)
            {
                int size = m_hitPoints.size() - count;
                for(int j=0;j<count;j+=2)
                {
                    for(int k=0;k<count;k+=2)
                    {
                        if(j!=k)
                        {
                            if((m_hitPoints[j+size] == m_hitPoints[k+size]) && (m_hitPoints[j+1+size]==m_hitPoints[k+1+size]))
                            {
                                m_hitPoints.erase(m_hitPoints.begin()+k+size);
                                m_hitPoints.erase(m_hitPoints.begin()+k+1+size);
                            }
                        }
                    }
                }
              float length[count];
              for(int k=0; k<count; ++k)
              {
                  ngl::Vec3 diff = m_rayStart[i] - m_hitPoints[size+k];
                  length[k] = sqrt(pow(diff.m_x,2) + pow(diff.m_y,2) + pow(diff.m_z,2));
               }

              for(int k=0;k<count;++k)
              {
                  for(int j=0;j<count;++j)
                  {
                      if(k!=j)
                      {
                        if(length[k] > length[j])
                        {
                            ngl::Vec3 temp = m_hitPoints[size+k];
                            m_hitPoints[size+k] = m_hitPoints[size+j];
                            m_hitPoints[size+j] = temp;
                            float t;
                            t = length[k];
                            length[k] = length[j];
                            length[j] = t;

                        }
                      }
                  }
               }

            }
        }


     }

}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::drawVolumePoints(ngl::TransformStack &_transformStack, ngl::Camera *_cam)
{
    ngl::VAOPrimitives *_prim=ngl::VAOPrimitives::instance();
    ngl::ShaderLib *shader=ngl::ShaderLib::instance();


    for(unsigned int i=0;i<m_volumePoints.size();++i)
    {
        _transformStack.pushTransform();
        {
            _transformStack.setPosition(m_volumePoints[i].m_x, m_volumePoints[i].m_y, m_volumePoints[i].m_z);

            loadMatricesToColourShader(_transformStack,_cam);
            shader->setShaderParam4f("Colour",1.0,0.0,0.0,1.0);
            _prim->draw("sphere");
        }
        _transformStack.popTransform();
    }
}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::drawSurfaceMeshPoints(ngl::TransformStack &_transformStack, ngl::Camera *_cam)
{
    ngl::VAOPrimitives *_prim=ngl::VAOPrimitives::instance();
    ngl::ShaderLib *shader=ngl::ShaderLib::instance();

    for(unsigned int i=0;i<m_surfacePointsMesh.size();++i)
    {
        _transformStack.pushTransform();
        {
            _transformStack.setPosition(m_surfacePointsMesh[i].m_x, m_surfacePointsMesh[i].m_y, m_surfacePointsMesh[i].m_z);

            loadMatricesToColourShader(_transformStack,_cam);
            shader->setShaderParam4f("Colour",1.0,0.0,0.0,1.0);
            _prim->draw("sphere");
        }
        _transformStack.popTransform();
    }
}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::drawBBox(ngl::TransformStack &_transformStack, ngl::Camera *_cam)
{
    ngl::ShaderLib *shader=ngl::ShaderLib::instance();
    // draw the mesh bounding box
    loadMatricesToColourShader(_transformStack,_cam);
    shader->setShaderParam4f("Colour",0,0,0,1.0);
    m_mesh->drawBBox();
}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::drawSurfacePointsBBox(ngl::TransformStack &_transformStack, ngl::Camera *_cam)
{
    ngl::VAOPrimitives *_prim=ngl::VAOPrimitives::instance();
    ngl::ShaderLib *shader=ngl::ShaderLib::instance();

    for(unsigned int i=0;i<m_surfacePointsBBox.size();++i)
    {
        _transformStack.pushTransform();
        {
            _transformStack.setPosition(m_surfacePointsBBox[i].m_x, m_surfacePointsBBox[i].m_y, m_surfacePointsBBox[i].m_z);
            loadMatricesToColourShader(_transformStack,_cam);
            shader->setShaderParam4f("Colour",0.0,1.0,0.0,1.0);
            _prim->draw("sphere");
        }
        _transformStack.popTransform();
    }
}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::drawHitPoints(ngl::TransformStack &_transformStack, ngl::Camera *_cam)
{
    ngl::VAOPrimitives *_prim=ngl::VAOPrimitives::instance();
    ngl::ShaderLib *shader=ngl::ShaderLib::instance();

    for(unsigned int i=0;i<m_hitPoints.size();++i)
    {
        _transformStack.pushTransform();
        {
            _transformStack.setPosition(m_hitPoints[i].m_x, m_hitPoints[i].m_y, m_hitPoints[i].m_z);
            loadMatricesToColourShader(_transformStack,_cam);
            shader->setShaderParam4f("Colour",0.0,0.0,1.0,0.5);
            _prim->draw("sphere");
        }
        _transformStack.popTransform();
    }
}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::drawVoronoiVertices(ngl::TransformStack &_transformStack, ngl::Camera *_cam)
{
    ngl::VAOPrimitives *_prim=ngl::VAOPrimitives::instance();
    ngl::ShaderLib *shader=ngl::ShaderLib::instance();

    for(unsigned int i=0;i<m_tetrahedra.size();++i)
    {
        _transformStack.pushTransform();
        {
            ngl::Vec3 pos(0.0,0.0,0.0);
            pos = m_tetrahedra[i]->getCirCenter();
            _transformStack.setPosition(pos.m_x, pos.m_y, pos.m_z);
            loadMatricesToColourShader(_transformStack,_cam);
            shader->setShaderParam4f("Colour",0.0,0.0,1.0,1.0);
            _prim->draw("sphere");
        }
        _transformStack.popTransform();
    }
}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::loadMatricesToColourShader(ngl::TransformStack &_transformStack, ngl::Camera *_cam)
{
  ngl::ShaderLib *shader=ngl::ShaderLib::instance();
  (*shader)["Colour"]->use();
  ngl::Mat4 MV;
  ngl::Mat4 MVP;

  MV=_transformStack.getCurrAndGlobal().getMatrix()*_cam->getViewMatrix() ;
  MVP=MV*_cam->getProjectionMatrix() ;
  shader->setShaderParamFromMat4("MVP",MVP);
}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::drawMesh(ngl::TransformStack &_transformStack, ngl::Camera *_cam)
{
    // grab an instance of the shader manager
    ngl::ShaderLib *shader=ngl::ShaderLib::instance();
    loadMatricesToColourShader(_transformStack,_cam);
    shader->setShaderParam4f("Colour",0.0,1.0,1.0,0.0);
    // Mesh is drawn only in Wireframe mode
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    // draw the mesh
    m_mesh->draw();
}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::delaunay()
{
   int rangeX = m_mesh->getBBox().width();
   int rangeY = m_mesh->getBBox().height();
   int rangeZ = m_mesh->getBBox().depth();
   int maxA = std::max(rangeX,rangeY);
   int range = std::max(maxA,rangeZ) * 10;
   std::vector<ngl::Vec3> vertex;
   vertex.push_back(ngl::Vec3(range, 0, -range));
   vertex.push_back(ngl::Vec3(0, range, -range));
   vertex.push_back(ngl::Vec3(-range, -range, -range));
   vertex.push_back(ngl::Vec3(0, 0, range));
   Tetrahedron *m_tetra;
   m_tetra = new Tetrahedron(vertex);

   Tetrahedron *_neigh = new Tetrahedron();
   _neigh = NULL;
   // Initialising empty neighbours for the first tetrahedron
   for(int i=0;i<4;++i)
   {
      m_tetra->m_neighbours[i] = _neigh;
   }

   m_tetra->createVAO();
   Delaunay *dt = new Delaunay(m_tetra);
   m_tetrahedra = dt->compute(m_points);

   m_voronoi = new Voronoi(m_tetrahedra);

}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::drawVoronoi(ngl::TransformStack &_transformStack, ngl::Camera *_cam)
{
    ngl::ShaderLib *shader=ngl::ShaderLib::instance();
    _transformStack.pushTransform();
    {
       shader->setShaderParam4f("Colour",0.0,0.0,0.0,0.0);
       loadMatricesToColourShader(_transformStack,_cam);
       glLineWidth(2);
       m_voronoi->draw();
       glLineWidth(1);
    }
   _transformStack.popTransform();
}

//----------------------------------------------------------------------------------------------------------------------
void MeshSampler::drawTetrahedron(ngl::TransformStack &_transformStack, ngl::Camera *_cam)
{
    // grab an instance of the shader manager
    ngl::ShaderLib *shader=ngl::ShaderLib::instance();

    for(unsigned int i=0; i<m_tetrahedra.size();++i)
    {
        _transformStack.pushTransform();
        {
            ngl::Colour c = m_tetrahedra[i]->getColour();
            shader->setShaderParam4f("Colour",c.m_r,c.m_g,c.m_b,c.m_a);
            loadMatricesToColourShader(_transformStack,_cam);
            glLineWidth(3);
            m_tetrahedra[i]->draw();
            glLineWidth(1);
        }
        _transformStack.popTransform();
    }

}
//----------------------------------------------------------------------------------------------------------------------


