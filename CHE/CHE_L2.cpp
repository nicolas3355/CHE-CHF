/**
* @file    CHE_L2.cpp
* @author  Marcos Lage         <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner      <thomas.lewiner@polytechnique.org>
* @author  Helio  Lopes        <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmidia
* @date    14/02/2006
*
* @brief  (Compact Half-Edge Structure - Level 2)
*/


#include <GL/glut.h>

#include "CHE_L2.hpp"



#include <ctime>

using namespace std;

//--------------------------------------------------//
vector<Vid> CHE_L2::R_00(const Vid v)
//--------------------------------------------------//
/** Computes the vertices in the star of a given vertex.*/
{
  vector<Vid> star;
  
  if( !v_valid(v) ) { 
    cout << "CHE_L2::Vertex Star ERROR: invalid vertex id" << endl;
    star.push_back(INV); 
    return star; 
  }

  //Gets the first incident half--edge
  HEid h = VH(v), h0 = h;
  
  //Walks through the star.  
  do 
    star.push_back( V(next(h)) );
  while( (O(h) != -1) && ( (h = next(O(h))) != h0) );

  return star;
}
//--------------------------------------------------//
vector<TRid> CHE_L2::R_02(const Vid v)
//--------------------------------------------------//
/** Computes the triangles in the star of a given vertex.*/
{
  vector<Vid> star;
  
  if( !v_valid(v) ) { 
    cout << "CHE_L2::Vertex Star ERROR: invalid vertex id" << endl;
    star.push_back(INV); 
    return star; 
  }

  //Gets the first incident half--edge
  HEid h = VH(v), h0 = h;
  
  //Walks through the star.  
  do 
    star.push_back( trig(h) );
  while( (O(h) != -1) && ( (h = next(O(h))) != h0) );

  return star;
}
//--------------------------------------------------//
void CHE_L2::compute_EH()
//--------------------------------------------------//
/** Computes the edges of the model.*/
{
  cout << "Pet_CHE::compute_EH...   " ;

  _EH.clear();

  for(HEid i=0; i<3*ntrig(); ++i)
	{
    if( !he_valid(i) ) continue;

		HEid he0 = i;
    HEid he1 = O(i);

    if( he1 == -1 ) _EH.insert(make_pair(he0, he1));
    else
    {
      if(he1 < he0) { HEid tmp = he0; he0 = he1; he1 = tmp; }
      if(_EH.find(he0) == _EH.end())
        _EH.insert(make_pair(he0, he1));
    }
	}
  cout << " " << (int)EH().size() << " edges found." << endl; 
}
//--------------------------------------------------//
void CHE_L2::compute_VH()
//--------------------------------------------------//
/** Computes the vertex Half-edge table.*/
{
  cout << "Pet_CHE::compute_VH...   " ;

   _VH.clear();
   _VH.resize( nvert(), -1 );

   for( HEid i=0; i<3*ntrig(); ++i)
   {
	   if(O(i) == -1) 
		   set_VH( V(i), i );
	   else
		   if( VH( V(i) ) == -1 ) set_VH( V(i), i );
   }
  cout << " done." << endl;
}
//--------------------------------------------------//
void CHE_L2::check()
//--------------------------------------------------//
/** Checks the mesh.*/
{
  CHE_L1::check();

  if(nvert() != static_cast<int>(  _VH.size())  ){ cout << "CHE_L2:: Erro nvert()!= VH.size" << endl; return;}

  for(int i=0; i<nvert(); ++i)
  {
    if( VH(i) >= 3*ntrig() )
    {
      cout << "CHE_L2:: Erro VH(" << i << ") >= 3*ntrig()." << endl;
      return;
    }
    if( VH(i) < 0 )
    {
      cout << "CHE_L2:: Erro VH(" << i << ") < 0." << endl;
      return;
    }
  }

  for(Ecit i=EH().begin(); i!=EH().end(); ++i)
  {
    if( (i->first >= 3*ntrig()) &&  (i->first < 0)   &&  (i->second != O(i->first)))
    {
      cout << "CHE_L2:: invalid edge." << endl;
      return;
    }
  }
}
//--------------------------------------------------//
void CHE_L2::draw_wire()
//--------------------------------------------------//
/** Draws the surface in wireframe with opengl.*/
{
  for(Ecit i= EH().begin(); i!= EH().end(); ++i)

	{

		const Vertex &v1 = G(V(   i->first     ));

		const Vertex &v2 = G(V( next(i->first) ));



		glColor3f(0.2,0.8,0.8);

		glBegin(GL_LINES);

			glNormal3d(v1.nx(),  v1.ny(),  v1.nz());

			glVertex3d(v1.x()  , v1.y(),   v1.z() );

			

			glNormal3d(v2.nx(),  v2.ny(),  v2.nz());

			glVertex3d(v2.x()  , v2.y(),   v2.z() );

		glEnd();

	}
}
//--------------------------------------------------//
void CHE_L2::read_ply( const char* file )
//--------------------------------------------------//
/** Gets a boundary compound. Sets _B*/
{
  // Stores Start && L2 time;
  clock_t start_time = static_cast<clock_t>(0.0),
                   L2_time = static_cast<clock_t>(0.0);

  start_time = clock();


  CHE_L1::read_ply( file );
	compute_EH();
	compute_VH();

  L2_time = clock();

  //cout << "L2 load time:" << static_cast<double>(L2_time-start_time)/static_cast<double>(CLOCKS_PER_SEC) << endl;

}
//------------------------------------------------//
