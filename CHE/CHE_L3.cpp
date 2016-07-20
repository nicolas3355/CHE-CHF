/**
* @file    CHE_L3.cpp
* @author  Marcos Lage         <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner      <thomas.lewiner@polytechnique.org>
* @author  Helio  Lopes        <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmidia
* @date    14/02/2006
*
* @brief  (Compact Half-Edge Structure - Level 3)
*/

#include <gl/glut.h>
#include <ctime>

#include "CHE_L3.hpp"

using namespace std;
//--------------------------------------------------//
void CHE_L3::compute_CH()
//--------------------------------------------------//
/** Computes the vertex Boundary Curves table.*/
{
  cout << "Pet_CHE::compute_CH...   " ;
   
  vector<bool> vst;
	vst.resize(3*ntrig(), false);
	
   _CH.clear();
   for(HEid he=0; he<3*ntrig(); ++he)
   {
	   if(O(he) == -1 && vst[he] == false)
	   {
		   HEid he0 = he;
		   _CH.push_back(he0); // add boundary curve representative
		   _ncurves++;              // increment number of curves
		   do                              // walk throght the boundary
		   {
			   // _O[he0] = -(ncurves() + 1); // set half--edge component
			   vst[he0] = true;                     // mark half--edge as visited

         while(O(next(he0)) >= 0) {
          //get the next half--edge in the boundary
			    he0 = O( next(he0) );
         }
         // real next boundary he.
         he0 = next(he0);
			}
			while(he0 != he);
      break;
     }
   }
  vst.clear();
  cout << " done." << endl;
}
//--------------------------------------------------//
void CHE_L3::check()
//--------------------------------------------------//
/** Checks the mesh.*/
{
  CHE_L2::check();

  if( ncurves() != static_cast<int>(_CH.size()) ) {
    cout << "ncurves != _CH.size()." << endl; 
    return;
  }

  for(HEid i=0; i<3*ntrig(); ++i){
    if(O(i) < -(ncurves()+1)) {
      cout << "O(" << i << ") < ncurves()." << endl;
      return;
    }
  }
}
//--------------------------------------------------//
void CHE_L3::draw_wire()
//--------------------------------------------------//
/** Draws the surface in wireframe with opengl.*/
{
  for(Ecit i= EH().begin(); i!= EH().end(); ++i)
	{
		const Vertex &v1 = G(V(   i->first     ));
		const Vertex &v2 = G(V( next(i->first) ));

    if(O(i->first) < 0) glColor3f(1.0,0.2,0.2);
    else glColor3f(0.6,0.6,0.6);

		glBegin(GL_LINES);
			glNormal3d(v1.nx(),  v1.ny(),  v1.nz());
			glVertex3d(v1.x()  , v1.y(),   v1.z() );
			
			glNormal3d(v2.nx(),  v2.ny(),  v2.nz());
			glVertex3d(v2.x()  , v2.y(),   v2.z() );
		glEnd();
	}
}
//--------------------------------------------------//
void CHE_L3::read_ply( const char* file )
//--------------------------------------------------//
/** Gets a boundary compound. Sets _CH*/
{
  // Stores Start && L3 time;
  clock_t start_time = static_cast<clock_t>(0.0),
                   L3_time = static_cast<clock_t>(0.0);

  start_time = clock();

  CHE_L2::read_ply( file );
	compute_CH();

  L3_time = clock();
  //cout << "L3 load time:" << static_cast<double>(L3_time-start_time)/static_cast<double>(CLOCKS_PER_SEC) << endl;
}
//-------------------------------------------------------------//
