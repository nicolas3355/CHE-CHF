/**
* @file    CHF_L2.cpp
* @author  Marcos Lage      <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner   <thomas.lewiner@polytechnique.org>
* @author  Hélio  Lopes     <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmídia
* @date    14/02/2006
*
* @brief   CHF: A Scalable topological data structure for tetrahedral meshes
* @brief  (Level 2)
*/
//--------------------------------------------------//
#include <set>
#include <map>
#include <stack>
#include <gl/glut.h>

#include "CHF_L2.hpp"    /**< Level 1 inheritance*/
#include "colorramp.h"	 /**< Gl color maps*/
												 
using namespace std;
//--------------------------------------------------//
vector<Vid> CHF_L2::R_00(const Vid v)
//--------------------------------------------------//
/** Computes the vertices in the star of a vertex.*/
{
  TEid t  = 0;
  HFid h = 0,  h1 = 0,  h2 = 0,  h3 = 0;

  stack<TEid> s;
	
  // vertex star.
  set<Vid> sstar;
  vector<Vid>  star;
 	
  // visited tetrahedra.
  set<TEid> tet;

  // Tests if the vertex v is valid.
  if( !v_valid(v) ) {	star.push_back(INV); return star; }

  // Finds the first tetrahedron that contains v.
  t = tetra( VH(v) );	s.push(t);
	
  // while the stack isn't empty:
  while( !s.empty() )
  {
    // picks the top tetrahedron.
    t = s.top();	s.pop();
			
	// If tetra was already visited, continue.
	pair<set<int>::iterator,bool> status = tet.insert(t);
	if( !status.second ) continue;
			
	for( int i = t<<2; i<(t+1)<<2; ++i )
	  if( V(i) ==  v )
	  {
	 	neighbors( i, h1, h2, h3 );
				
	  	// insert the vertex in the star
		sstar.insert( V(h1) );
		sstar.insert( V(h2) );
		sstar.insert( V(h3) );
		
		//includes the neighboors tetrahedra in the stack 
		h = O(h1) ; if( hf_valid(h) ) s.push( tetra(h) ) ;
		h = O(h2) ; if( hf_valid(h) ) s.push( tetra(h) ) ;
		h = O(h3) ; if( hf_valid(h) ) s.push( tetra(h) ) ;
		break;
	  }
  }
	
  // Copies in star
  for(set<Vid>::const_iterator iter = sstar.begin(); iter != sstar.end(); ++iter )
    star.push_back( *iter );
 	
  return star; 
}
//--------------------------------------------------//
vector<TEid> CHF_L2::R_03(const Vid v)
//--------------------------------------------------//
/** Computes the tetrahedrons in the star of a vertex.*/
{
  TEid  t = 0;
  HFid h = 0,  h1 = 0,  h2 = 0,  h3 = 0;

  stack<TEid> s;
  
  // Vertex star
  set<TEid> sstar;
  vector<TEid>  star;
 
  //Tests if the vertex v is valid.
  if( !v_valid(v) ) {	star.push_back(INV); return star;	}

  // Finds the first tetrahedron that contains v.
  t = tetra( VH(v) );	s.push(t);

  //while the stack isn't empty:
  while( !s.empty() )
  {
	//picks the top tetrahedron.
    t = s.top();	s.pop();
		
	// If tetra was already visited, continue.
	pair<set<int>::iterator,bool> status = sstar.insert(t);
	if( !status.second ) continue;
		
   	for( int i = t<<2; i<(t+1)<<2; ++i )
   	  if( V(i) ==  v )
	  {
	    neighbors( i, h1, h2, h3 );

 		//includes the neighboors tetrahedra in the stack 
 		h = O(h1) ; if( hf_valid(h) ) s.push( tetra(h) ) ;
		h = O(h2) ; if( hf_valid(h) ) s.push( tetra(h) ) ;
		h = O(h3) ; if( hf_valid(h) ) s.push( tetra(h) ) ;
		break;
	  }
  }
	  
  // Copies in star
  for(set<Vid>::const_iterator iter = sstar.begin(); iter != sstar.end(); ++iter )
	 star.push_back( *iter );

  return star;
}
//--------------------------------------------------//
vector<Vid> CHF_L2::R_10(const Vid a, const Vid b)
//--------------------------------------------------//
/** Computes the vertices in the star of a given edge.*/
{
  Eit pos;
  set<Vid> sstar;
  vector<Vid>   star;
	
  pair<HFid, HFid>  start, runcw;
	
  if( !v_valid(a) || !v_valid(b) || a ==b) { star.push_back(INV); return star; }
	
  //Gets the first edge
  Vid A = a, B =b;
  if(b<a) { Vid tmp= A; A=B; B=tmp; }
  pos = _EH.find(make_pair(A,B));
	
  if(pos == _EH.end()) {
    star.push_back(INV);
    return star;
  }
	
  HFid lives = pos->second, hf0 = 0, hf1 = 0, hf2 = 0;
  neighbors(lives, hf0, hf1, hf2);
	
  if      (V(hf0) == a) start  = runcw = make_pair(lives, hf0);
  else if (V(hf1) == a) start  = runcw = make_pair(lives, hf1);
  else if (V(hf2) == a)	start  = runcw = make_pair(lives, hf2);
	
  //Walk arround the edge --"Clockwise"
  do
  {
    sstar.insert( V(runcw.first) );
	runcw = matehe(runcw.first, runcw.second);
	runcw = radialhe(runcw.first, runcw.second);
  }
  while((runcw.first != -1) && (runcw.first != start.first));

  // Copies in star
  for(set<Vid>::const_iterator iter = sstar.begin(); iter != sstar.end(); ++iter )
 	 star.push_back( *iter );

  return star;
}
//--------------------------------------------------//
vector<TEid> CHF_L2::R_13(const Vid a, const Vid b)
//--------------------------------------------------//
/** Computes the vertices in the star of a given edge.*/
{
  Eit pos;

  set<TEid> sstar;
  vector<TEid>   star;
	
  pair<HFid, HFid>  start, runcw;
	
  if( !v_valid(a) || !v_valid(b) || a ==b) { star.push_back(INV); return star; }
	
  //Gets the first edge
  Vid A = a, B =b;
  if(b<a) { Vid tmp= A; A=B; B=tmp; }
    pos = _EH.find(make_pair(A,B));
	
  if(pos == _EH.end()) 
  {
	star.push_back(INV);
	return star;
  }
	
  HFid lives = pos->second, hf0 = 0, hf1 = 0, hf2 = 0;
  neighbors(lives, hf0, hf1, hf2);
	
  if      (V(hf0) == a) start  = runcw = make_pair(lives, hf0);
  else if (V(hf1) == a) start  = runcw = make_pair(lives, hf1);
  else if (V(hf2) == a)	start  = runcw = make_pair(lives, hf2);
	
  //Walk arround the edge --"Clockwise"
  do
  {
    sstar.insert( tetra(runcw.first) );
	runcw = matehe(runcw.first, runcw.second);
	runcw = radialhe(runcw.first, runcw.second);
  }
  while((runcw.first != -1) && (runcw.first != start.first));

  // Copies in star
  for(set<Vid>::const_iterator iter = sstar.begin(); iter != sstar.end(); ++iter )
    star.push_back( *iter );

  return star;
}
//--------------------------------------------------//
void CHF_L2::create_VH()
//--------------------------------------------------//
/** Creates the "half-face of a vertex" container */
{
  _VH.clear();
  _VH.resize( nvert(), -1 );

  for(HFid i=0; i<4*ntetra(); ++i)
  {
	HFid hf1 = 0,  hf2 = 0,  hf3 = 0;
	neighbors(i, hf1, hf2, hf3);
	if( O(i) == -1 )
	{
   	  set_VH( V(hf1), i );
	  set_VH( V(hf2), i );
	  set_VH( V(hf3), i );
	 }
	else
	{
      if( VH( V(hf1) ) == -1 ) set_VH( V(hf1), i );
      if( VH( V(hf2) ) == -1 ) set_VH( V(hf2), i );
      if( VH( V(hf3) ) == -1 ) set_VH( V(hf3), i );
	}
  }
  cout << "CHF_L3::create_VH: " << (unsigned)_VH.size() << " vertices found." << endl;
}
//--------------------------------------------------//
void CHF_L2::create_EH()
//--------------------------------------------------//
/** Creates the "half-face of a edge" map */
{
  Eit pos;
  Vid   a, b;
  HFid ha, hb;

  static const char comb1[6] = { 0,0,0,1,1,2 } ;
  static const char comb2[6] = { 1,2,3,2,3,3 } ;

  _EH.clear();

  for( TEid t=0; t< ntetra(); ++t)
  {
    if( !te_valid(t) ) continue;
		
	TEid te = t << 2; 
	for( Vid i=0; i< 6; ++i )
	{
	  ha = te | comb1[i] ;
	  hb = te | comb2[i] ;
			
	  a= V( ha ) ;
	  b= V( hb ) ;
			
	  if(b<a)
	  { 
	    Vid     tmp= a; a=b; b=tmp; 
	    HFid  htmp= ha; ha=hb; hb = htmp;
	  }

      pos = _EH.find( Eid(a,b) );
			
	  if( pos == _EH.end() )
	  {
      	HFid hf = te | comb1[5-i];
		if( O(hf) != -1 ) 
		{
		  HFid hftmp = te | comb2[5-i];
		  if(O(hftmp) == -1 && V(nexthe(hftmp, ha).second) == b)
		  hf =hftmp;
		}

		_EH.insert( make_pair( Eid(a,b), hf ) );
	  }
	  else if ( O(pos->second) != -1)
	  {
	    HFid hf = te | comb1[5-i] ;
		if( O( hf ) == -1 && V(nexthe(hf, ha).second) == b) pos->second = hf ;
		else
		{
		  hf = te | comb2[5-i] ;
		  if( O( hf ) == -1 && V(nexthe(hf, ha).second) == b) pos->second = hf;
		}
	  }
	}
  }
  cout << "CHF_L3::create_EH: " << (unsigned)_EH.size() << " edges found." << endl;
}
//--------------------------------------------------//
void CHF_L2::create_FH()
//--------------------------------------------------//
{
  _FH.clear();

  for(HFid i=0; i<4*ntetra(); ++i)
  {
	HFid hf0 = i;
    HFid hf1 = O(i);

    if( hf1 == -1 ) _FH.insert(make_pair(hf0, hf1));
    else
    {
      if(hf1 < hf0) { HFid tmp = hf0; hf0 = hf1; hf1 = tmp; }
      if(_FH.find(hf0) == _FH.end())
        _FH.insert(make_pair(hf0, hf1));
    }
  }
  cout << "CHF_L3::create_FH: " << (unsigned)_FH.size() << " faces found." << endl;
}
//--------------------------------------------------//
void CHF_L2::check()
//--------------------------------------------------//
/** Checks Level 2 structure */
{ 
  CHF_L1::check();

  for(Vid i=0; i< nvert(); ++i)
  {
    if(VH(i) == -1)
    {
      cout << "CHF_L2::Check ERRO : VH[" << i <<"] == -1 " << endl;
      return;
    }
  }

  for( Ecit i= EH().begin() ; i != EH().end(); ++i)
  {
    if(i->second == -1) 
    {
      cout << "CHF_L2::Check ERRO: EH[i]" << "->second == -1" << endl;
      return;
    }

	if(i->first.first >= nvert() || i->first.second>= nvert()) 
    {
   	  cout << "CHF_L2::Check ERRO :i->first.first  " << i->first.first << " >= nvert() " << endl;
      cout << "                    i->firsf.second " << i->first.second<< " >= nvert() " << endl;
      return;
    }
  }
  cout << "CHF_L2::Check Ok" << endl;
}
//--------------------------------------------------//
void CHF_L2::draw_wire(const int t)
//--------------------------------------------------//
/** Draws the mesh in wireframe */
{
  ColorRamp c;

  if( t == 5 ) //Draws the edges
  {
    for( Ecit i= EH().begin() ; i != EH().end(); ++i)
    {
	  const Vertex &v0 = G(i->first.first);
	  const Vertex &v1 = G(i->first.second);

	  c.set_GLcolor( 0.512, 1, COLOR_GRAYSCALE, 1 );

      glBegin( GL_LINES );
 	    glNormal3f(v0.nx(), v0.ny(), v0.nz());
 	    glVertex3f(v0.x() , v0.y() , v0.z());
 	    glNormal3f(v1.nx(), v1.ny(), v1.nz());
        glVertex3f(v1.x() , v1.y() , v1.z());
  	  glEnd();
	}
  }
  if( t == 6 ) //Draws the edges with boundary classification
  {
    for( Ecit i= EH().begin() ; i != EH().end(); ++i)
    {
	  const Vertex &v0 = G(i->first.first);
	  const Vertex &v1 = G(i->first.second);

	  if( O(i->second) == -1 )c.set_GLcolor( 0.512, 1, COLOR_GRAYSCALE, 1 );
	  else c.set_GLcolor( 0.012, 1, COLOR_AUTUMN, 1 );
			
      glBegin( GL_LINES );
 	    glNormal3f(v0.nx(), v0.ny(), v0.nz());
 	    glVertex3f(v0.x() , v0.y() , v0.z());
 	    glNormal3f(v1.nx(), v1.ny(), v1.nz());
        glVertex3f(v1.x() , v1.y() , v1.z());
  	  glEnd();
	}
  }

  if( t == 7 ) //Draws the boundary edges
  {
    for( Ecit i= EH().begin() ; i != EH().end(); ++i)
    {
	  const Vertex &v0 = G(i->first.first);
	  const Vertex &v1 = G(i->first.second);

	  if( O(i->second) == -1 )c.set_GLcolor( 0.512, 1, COLOR_GRAYSCALE, 1 );
	  else continue;
			
      glBegin( GL_LINES );
 	    glNormal3f(v0.nx(), v0.ny(), v0.nz());
 	    glVertex3f(v0.x() , v0.y() , v0.z());
 	    glNormal3f(v1.nx(), v1.ny(), v1.nz());
        glVertex3f(v1.x() , v1.y() , v1.z());
  	  glEnd();
    }
  }
  if( t != 5 && t != 6 && t != 7 )
    cout << "CHF_L2::draw_wire ERRO." << endl;
 }
//--------------------------------------------------//
void CHF_L2::draw_vert(const int t, const int p)
//--------------------------------------------------//
/** Draws the vertices of the mesh */
{
  ColorRamp c;

  if( t == 1 ) // Draw Verts
  {
    for(int i=0; i<nvert(); ++i)
    {
	  const Vertex &v1 = G(i);
      c.set_GLcolor( 0.512, 1, COLOR_RAINBOW, 1 );
		  
  	  if( !p )
      {
        glPushMatrix();
        glTranslatef( v1.x(), v1.y(), v1.z());
	    glutSolidSphere(0.008, 7, 7);
	    glPopMatrix();
      }
      else
      {
        glPointSize(2.7);
	    glBegin(GL_POINTS);
          glNormal3f( v1.nx(), v1.ny(), v1.nz() );
	      glVertex3f( v1.x() , v1.y() , v1.z() );
	    glEnd();
	  }
	}    
  }
  if( t == 2 ) // Draw verts with boundary classification
  {
    for(int i=0; i<nvert(); ++i)
	{
	  const Vertex &v1 = G(i);
      if( O( VH(i) ) == -1 ) c.set_GLcolor( 0.512, 1, COLOR_RAINBOW, 1 );
	  else c.set_GLcolor( 0.212, 1, COLOR_AUTUMN, 1 );
		  
 	  if( !p )
      {
   	    glPushMatrix();
        glTranslatef  ( v1.x(), v1.y(), v1.z());
    	glutSolidSphere(0.008, 6, 6);
  		glPopMatrix();
      }
      else
      {
        glPointSize(2.7);
		glBegin(GL_POINTS);
          glNormal3f( v1.nx(), v1.ny(), v1.nz() );
		  glVertex3f( v1.x() , v1.y(),  v1.z()  );
		glEnd();
      }
	}    
  }

  if( t == 3 ) // Draw boundary verts
  {
    for(int i=0; i<nvert(); ++i)
	{
	  const Vertex &v1 = G(i);
      if( O( VH(i) ) == -1 )c.set_GLcolor( 0.512, 1, COLOR_RAINBOW, 1 );
      else continue;
		  
      if( !p )
      {
 	    glPushMatrix();
        glTranslatef  ( v1.x(), v1.y(), v1.z());
	    glutSolidSphere(0.008, 6, 6);
	    glPopMatrix();
      }
      else
      {
        glPointSize(2.7);
		glBegin(GL_POINTS);
          glNormal3f( v1.nx(), v1.ny(), v1.nz() );
		  glVertex3f( v1.x() , v1.y(),  v1.z() );
		glEnd();
      }      
	}
  }

  if( t == 4 ) // Draw verts with scalar atributes
  {
    float maxf=0;
    for(int i=0; i<nvert(); ++i)
    {
 	  const Vertex &v1 = G(i);
      if(maxf < v1.f()) maxf = v1.f();
    }

	for(int i=0; i<nvert(); ++i)
	{
	  const Vertex &v1 = G(i);
      c.set_GLcolor( v1.f(), maxf, COLOR_RAINBOW, 1 ,true );
		  
      if( !p )
      {
 		glPushMatrix();
        glTranslatef  ( v1.x(), v1.y(), v1.z());
		glutSolidSphere(0.008, 6, 6);
        glPopMatrix();
      }
      else
      {
        glPointSize(2.7);
		  glBegin(GL_POINTS);
          glNormal3f( v1.nx(), v1.ny(), v1.nz() );
 		  glVertex3f(v1.x()  , v1.y(),   v1.z() );
		glEnd();
      }
	}    
  }
  if( t != 1 && t != 2 && t != 3 && t != 4)
    cout << "CHF_L2::draw_vert ERRO." << endl;
}
//--------------------------------------------------------------//