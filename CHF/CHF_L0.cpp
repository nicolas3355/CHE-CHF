/**
* @file    CHF_L0.cpp
* @author  Marcos Lage      <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner   <thomas.lewiner@polytechnique.org>
* @author  Hélio  Lopes     <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmídia
* @date    14/02/2006
*
* @brief   CHF: A Scalable topological data structure for tetrahedral meshes
* @brief  (Level 0)
*/
//--------------------------------------------------//
#include <set>
#include <ctime>
#include <fstream>
#include <algorithm>
#include <GL/glut.h>

#include "ply.h"        /**< PLY exportation*/
#include "fparser.h"    /**< Parses scalar Field*/  
#include "colorramp.h"  /**< Gl color maps*/
#include "CHF_L0.hpp"   /**< Level 0 inheritance*/

using namespace std;
//--------------------------------------------------//
vector<Vid> CHF_L0::R_00(const Vid v)
//--------------------------------------------------//
/** Computes the vertices in the star of a vertex.*/
{
      set<Vid> sstar;
   vector<Vid>  star;
  
  HFid  h1 = 0,  h2 = 0,  h3 = 0;
  if( !v_valid(v) ) { star.push_back(INV); return star; }

  for(HFid i=0; i<ntetra()<<2; ++i){   
	if( V(i) == v )
	{
      neighbors( i, h1, h2, h3 );
	
	  sstar.insert( V(h1) );   
	  sstar.insert( V(h2) );   
	  sstar.insert( V(h3) );
	}
  }

  for(set<Vid>::const_iterator iter = sstar.begin(); iter != sstar.end(); ++iter )
	star.push_back( *iter );
  
  return star;
}
//--------------------------------------------------//
vector<TEid> CHF_L0::R_03(const Vid v)
//--------------------------------------------------//
/** Computes the tetrahedrons in the star of a vertex.*/
{
  vector<Vid> star;
  if( !v_valid(v) ) { star.push_back(INV); return star; }

  for(HFid i=0; i<ntetra()<<2; ++i)
  { 
    if( V(i) == v ) 
	  star.push_back( tetra(i) );
  }
  return star;
}
//--------------------------------------------------//
vector<Vid> CHF_L0::R_10(const Vid a, const Vid b)
//--------------------------------------------------//
/** Computes the vertices in the star of a given edge.*/
{
     set<Vid> sstar;
  vector<Vid>   star;

  HFid  h1 = 0,  h2 = 0,  h3 = 0;
  if( !v_valid(a) || !v_valid(b) || a ==b) { star.push_back(INV); return star; }

  for(HFid i=0; i<ntetra()<<2; ++i){ 
    if( (V(i) == a) && (  (V(nexthf(i)) == b) ||  (V(midhf(i)) == b)  || (V(prevhf(i)) == b) )  ){
	  neighbors(i, h1, h2, h3);
  	  if(V(h1) != b) sstar.insert(V(h1));
  	  if(V(h2) != b) sstar.insert(V(h2));
  	  if(V(h3) != b) sstar.insert(V(h3));
    }
  }
  
  if(sstar.empty()) { star.push_back(INV); return star; }
  
  for(set<Vid>::const_iterator iter = sstar.begin(); iter != sstar.end(); ++iter )
	star.push_back( *iter );
 
  return star;
}
//--------------------------------------------------//
vector<TEid> CHF_L0::R_13(const Vid a, const Vid b)
//--------------------------------------------------//
/** Computes the vertices in the star of a given edge.*/
{
     set<Vid>  sstar;
  vector<Vid>   star;
 
  if( !v_valid(a) || !v_valid(b) || a ==b ) { star.push_back(INV); return star; }

  for(HFid i=0; i<ntetra()<<2; ++i)
  {
    if( (V(i) == a) && ( (V(nexthf(i)) == b) ||  (V(midhf(i)) == b)  || (V(prevhf(i)) == b) )  )
	{
      sstar.insert(tetra(i));
	}
  }
  
  if(sstar.empty()) { star.push_back(INV); return star; }

  for(set<Vid>::const_iterator iter = sstar.begin(); iter != sstar.end(); ++iter )
	star.push_back( *iter );

  return star;
}
//--------------------------------------------------//
vector<TEid> CHF_L0::R_33(const TEid t)
//--------------------------------------------------//
/** Computes the tetrahedrons incidents to a tetrahedron.*/
{
     set<Vid> sstar;
  vector<Vid>  star;
	
  Vid a = V(t<<2), b = V((t<<2) +1), c = V((t<<2)+2), d = V((t<<2)+3);
 
  if( !te_valid(t) ) { star.push_back(INV); return star; }
  
  for(HFid i=0; i<ntetra()<<2; ++i)
  {
    if( (V(i) == a) && ( (V(nexthf(i)) == b) ||  (V(midhf(i)) == b)  || (V(prevhf(i)) == b) ) && ( (V(nexthf(i)) == c) ||  (V(midhf(i)) == c)  || (V(prevhf(i)) == c) )  )
	{
   	  sstar.insert( tetra(i) );
	}
 	if( (V(i) == a) && ( (V(nexthf(i)) == c) ||  (V(midhf(i)) == c)  || (V(prevhf(i)) == c) ) && ( (V(nexthf(i)) == d) ||  (V(midhf(i)) == d)  || (V(prevhf(i)) == d) )  )
	{
	  sstar.insert( tetra(i) );
	}
 	if( (V(i) == a) && ( (V(nexthf(i)) == b) ||  (V(midhf(i)) == b)  || (V(prevhf(i)) == b) ) && ( (V(nexthf(i)) == d) ||  (V(midhf(i)) == d)  || (V(prevhf(i)) == d) )  )
	{
	  sstar.insert( tetra(i) );
	}
 	if( (V(i) == b) && ( (V(nexthf(i)) == c) ||  (V(midhf(i)) == c)  || (V(prevhf(i)) == c) ) && ( (V(nexthf(i)) == d) ||  (V(midhf(i)) == d)  || (V(prevhf(i)) == d) )  )
	{
	  sstar.insert( tetra(i) );
	}
  }

  for(set<Vid>::const_iterator iter = sstar.begin(); iter != sstar.end(); ++iter )
	star.push_back( *iter );

  return star;
}
//--------------------------------------------------//
void CHF_L0::bounding_box( float *min, float *max )
//--------------------------------------------------//
{
  float t_mx,t_Mx,t_my,t_My,t_mz,t_Mz;

  t_mx=t_Mx=_G[0].x();
  t_my=t_My=_G[0].y();
  t_mz=t_Mz=_G[0].z();

  for(Vid i=1; i<nvert(); ++i)
  {
    if(_G[i].x() < t_mx ) t_mx=_G[i].x();
	if(_G[i].x() > t_Mx ) t_Mx=_G[i].x();
																	
	if(_G[i].y() < t_my ) t_my=_G[i].y();
	if(_G[i].y() > t_My ) t_My=_G[i].y();
																	
	if(_G[i].z() < t_mz ) t_mz=_G[i].z();
	if(_G[i].z() > t_Mz ) t_Mz=_G[i].z();
  }

  min[0]=t_mx; min[1]=t_my; min[2]=t_mz;
  max[0]=t_Mx; max[1]=t_My; max[2]=t_Mz;
}
//--------------------------------------------------//
void CHF_L0::legalize_model  (float *min, float *max)
//--------------------------------------------------//
{
  float c[3];
  float l[3];
  float size=0;

  c[0]=(float)(max[0]+min[0])/2;
  c[1]=(float)(max[1]+min[1])/2;
  c[2]=(float)(max[2]+min[2])/2;

  for(Vid i=0; i<3; i++) 
  {
    l[i]=fabs(max[i]-min[i]);
    if(l[i]>size) size=l[i];
  }

	size /= 1.5;

  for(Vid i=0; i<nvert(); i++)
  {
    float tx = _G[i].x();
    float ty = _G[i].y();
    float tz = _G[i].z();

    tx -= c[0];
    if(size != 0) tx /= size;
    _G[i].set_x(tx);

    ty -= c[1];
    if(size != 0) ty /= size;
    _G[i].set_y(ty);

    tz -= c[2];
    if(size != 0) tz /= size;
    _G[i].set_z(tz);
  }

  for(Vid i=0; i<3; i++)
  {
    min[i] -= c[i];
    if(size != 0) min[i] /= size;
    max[i] -= c[i];
    if(size != 0) max[i] /= size;
  }
}
//--------------------------------------------------//
void CHF_L0::scalar_field(const char* eq)
//--------------------------------------------------//
/** Assigns a parsed formula on model vertices and creates a scalar field*/
{
  FunctionParser fparse;
  fparse.Parse(eq, "x,y,z" );
  
  float num[3], val;
  
  for(Vid i=0; i<nvert(); i++)
  {
    if( !v_valid( i ) ) continue ;
    num[0]= G(i).x();
    num[1]= G(i).y();
    num[2]= G(i).z();

    val= fparse.Eval(num);
    _G[i].set_f(val);
  }
}
//--------------------------------------------------//
void CHF_L0::check()
//--------------------------------------------------//
/** Checks the basic structure validation */
{
  if( static_cast<int>(_G.size()) != nvert()) 
  {
    cout << "CHF_L0::Check ERRO : _G.size (" << (unsigned) _G.size() << ") != nverts (" << nvert() << ")" << endl;
    return;
  }
  if( static_cast<int>(_V.size() )!= ntetra()<<2) 
  {
    cout << "CHF_L0::Check ERRO : _V.size	(" << (unsigned)_V.size() << ") != 4*ntetras (" << ntetra()<<2 << ")" << endl;
    return;
  }
  
  for(Vid i=0; i<ntetra()<<2; ++i)
  {
    if(V(i) >= nvert()) 
    {
      cout << "CHF_L0::Check ERRO : V(" << i << ") >= nvert (" << nvert() << ")" << endl;
      return;
	}

    if(V(i) == INV) 
    {
      cout << "CHF_L0::Check ERRO : V(" << i << ") Invalid." << endl;
      return;
    }
  }
  cout << "CHF_L0::Check Ok" << endl;
  return;
}
//--------------------------------------------------//
void CHF_L0::draw_wire(const int t)
//--------------------------------------------------//
/** Draws the mesh in wireframe */
{
   ColorRamp c;

  if( t==5 )
  {
    for(int i=0; i< ntetra(); ++i)
    {
      const Vertex &v0 = G(V(i<<2));
      const Vertex &v1 = G(V(i<<2 | 1));
      const Vertex &v2 = G(V(i<<2 | 2));
      const Vertex &v3 = G(V(i<<2 | 3));

      c.set_GLcolor( 0.512, 1, COLOR_GRAYSCALE, 1 );         
      
      glBegin( GL_LINES );
        glNormal3f(v0.nx(), v0.ny(), v0.nz());
        glVertex3f(v0.x() , v0.y() , v0.z());
        glNormal3f(v1.nx(), v1.ny(), v1.nz());
        glVertex3f(v1.x() , v1.y() , v1.z());

        glNormal3f(v0.nx(), v0.ny(), v0.nz());
        glVertex3f(v0.x() , v0.y() , v0.z());
        glNormal3f(v2.nx(), v2.ny(), v2.nz());
		glVertex3f(v2.x() , v2.y() , v2.z());

        glNormal3f(v0.nx(), v0.ny(), v0.nz());
        glVertex3f(v0.x() , v0.y() , v0.z());
        glNormal3f(v3.nx(), v3.ny(), v3.nz());
		glVertex3f(v3.x() , v3.y() , v3.z());

        glNormal3f(v1.nx(), v1.ny(), v1.nz());
        glVertex3f(v1.x() , v1.y() , v1.z());
        glNormal3f(v2.nx(), v2.ny(), v2.nz());
		glVertex3f(v2.x() , v2.y() , v2.z());

        glNormal3f(v2.nx(), v2.ny(), v2.nz());
        glVertex3f(v2.x() , v2.y() , v2.z());
        glNormal3f(v3.nx(), v3.ny(), v3.nz());
		glVertex3f(v3.x() , v3.y() , v3.z());

        glNormal3f(v3.nx(), v3.ny(), v3.nz());
        glVertex3f(v3.x() , v3.y() , v3.z());
        glNormal3f(v1.nx(), v1.ny(), v1.nz());
		glVertex3f(v1.x() , v1.y() , v1.z());
      glEnd();
    }
  }
  else 
    cout << "CHF_L0::draw_wire ERRO." << endl;
}
//--------------------------------------------------//
void CHF_L0::draw_vert(const int t, const int p)
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
  		
	  if( p == 0 ) {
	    glPushMatrix();
		glTranslatef( v1.x(), v1.y(), v1.z());
		glutSolidSphere(0.008, 6, 6);
		glPopMatrix();
	  }
	  else {     
	    glPointSize(2.7);
		glBegin(GL_POINTS);
		  glNormal3f( v1.nx(), v1.ny(), v1.nz() );
		  glVertex3f( v1.x() , v1.y(),  v1.z()  );
		glEnd();
	  }
	}    
  }
  if( t == 4 ) // Draw verts with scalar atributes
  {
    float maxf = G(0).f();
    for(int i=0; i<nvert(); ++i)
    {
 	  const Vertex &v1 = G(i);
      if(maxf < v1.f()) maxf = v1.f();
    }

	for(int i=0; i<nvert(); ++i)
	{
	  const Vertex &v1 = G(i);
      c.set_GLcolor( v1.f(), maxf, COLOR_RAINBOW, 1 , true);
 			
      if( p == 0 )
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
 	      glVertex3f( v1.x(),  v1.y(),  v1.z()  );
	    glEnd();
      }
	}   
  }
  if ( t != 1 && t != 4) cout << "CHF_L0::draw_vert ERRO." << endl;
}
//--------------------------------------------------//
/** PLY Data*/
//--------------------------------------------------//
/** PLY Point Data*/
typedef struct PlyPoint
{
	float x, y, z;
}PlyPoint;

/* list of property information for a PlyVertex */
PlyProperty plyvert_props[]  = { 
  {"x",  Float32, Float32, offsetof( PlyPoint, x  ), 0, 0, 0, 0},
  {"y",  Float32, Float32, offsetof( PlyPoint, y  ), 0, 0, 0, 0},
  {"z",  Float32, Float32, offsetof( PlyPoint, z  ), 0, 0, 0, 0},
};

typedef struct PlyFace
{
  unsigned char nverts;    /* number of Vertex indices in list */
  int *verts;              /* Vertex index list */
} PlyFace;

/* list of property information for a PlyFace */
PlyProperty plyface_props[]  = { 
{  "vertex_indices", Int32, Int32, offsetof( PlyFace,verts ),
   1, Uint8, Uint8, offsetof( PlyFace,nverts )}
};

//--------------------------------------------------//
void CHF_L0::read_ply( const char* fn )
//--------------------------------------------------//
/** Reads a tetrahedral mesh in the PLY file format.*/
{
	_G.clear();
	_V.clear();
  float xmin= 0, xmax= 0, ymin= 0, ymax=0, zmin=0, zmax= 0;

  // Stores Start && L1 && L2 time;
  clock_t start_time = static_cast<clock_t>(0.0),
                   L0_time = static_cast<clock_t>(0.0);

  start_time = clock();

  printf("CHF_L0::read_ply(%s)...", fn) ;

  /*** the Ply object ***/
  PlyFile  *in_ply;
  PlyFace     face;
  PlyPoint       v;

  /** loop indices */
  int   j, ntet, nv;
  int   i, elem_count;
  char  *elem_name;

  FILE *fp = fopen( fn, "r" );

  if( fp==NULL ) { printf(" file not found.\n" ); return; }

  in_ply   = ::read_ply (fp);

  /** gets the number of faces and vertices*/
  ntet = nv = 0 ;
  for ( i = 0; i < in_ply->num_elem_types; i++ )
  {
    elem_name = setup_element_read_ply ( in_ply, i, &elem_count );
    if ( equal_strings ( "vertex", elem_name ) )
      nv   = elem_count;
    if ( equal_strings ( "face",   elem_name ) )
      ntet = elem_count;
  }

  set_nvert (nv);
  set_ntetra(ntet);

  _G.resize(   nv   );
  _V.resize( ntet<<2 );

  /* examine each element type that is in the file (PlyVertex, PlyFace) */
  for ( i = 0; i < in_ply->num_elem_types; i++ )
  {
    /* prepare to read the i'th list of elements */
    elem_name = setup_element_read_ply ( in_ply, i, &elem_count );

    if ( equal_strings ( "vertex", elem_name ) )
    {
      /* set up for getting PlyVertex elements */
      setup_property_ply ( in_ply, &plyvert_props[0] );
      setup_property_ply ( in_ply, &plyvert_props[1] );
      setup_property_ply ( in_ply, &plyvert_props[2] );

      for ( j = 0; j < nv; ++j )
      { 
		    Vertex p;
        get_element_ply ( in_ply, ( void * ) &(v) );
		    p.set_x ( (float)v.x ) ;
		    p.set_y ( (float)v.y ) ;
		    p.set_z ( (float)v.z ) ;
		    
				p.set_nx( 0 ) ;
		    p.set_ny( 0 ) ;
		    p.set_nz( 0 ) ;

        if( xmin > p.x() ) xmin = p.x();
        if( xmax < p.x() ) xmax = p.x();
        
        if( ymin > p.y() ) ymin = p.y();
        if( ymax < p.y() ) ymax = p.y();
        
        if( zmin > p.z() ) zmin = p.z();
        if( zmax < p.z() ) zmax = p.z();

		    set_G(j, p);
      }
    }
    else if ( equal_strings ( "face", elem_name ) )
    {
      setup_property_ply ( in_ply, &plyface_props[0] );
      for ( j = 0; j < ntet; j++ )
      {
        get_element_ply ( in_ply, ( void * ) &face );
        if( face.nverts != 4 ) { printf("PLY importation: not a tetrahedral mesh\n"); }
        _V[j<<2]     = face.verts[0] ;
        _V[j<<2 | 1] = face.verts[1] ;
        _V[j<<2 | 2] = face.verts[2] ;
        _V[j<<2 | 3] = face.verts[3] ;
      }
    }
    else  /* all non-PlyVertex and non-PlyFace elements are grabbed here */
      get_other_element_ply ( in_ply );
  }
  close_ply ( in_ply );
  free_ply  ( in_ply );

	//Normalizes the model.
	float min[3], max[3];
	bounding_box( min, max );
	legalize_model( min, max);

  L0_time = clock();

  printf(" %d vertices and %d tetrahedrons found\n", nv, ntet ) ;

  //cout << "L0 load time:" << static_cast<double>(L0_time-start_time)/static_cast<double>(CLOCKS_PER_SEC) << endl;
}
//--------------------------------------------------//
void CHF_L0::write_ply( const char* file, bool bin )
//--------------------------------------------------//
/** Writes a 3D triangulated model in the PLY file format.*/
{
  printf("CHF_L0::write_ply(%s)...", file) ;

  PlyFile    *ply;
  FILE       *fp = fopen( file, "w" );

  PlyFace     face ;
  int         v[4] ;
  char       *elem_names[]  = { "vertex", "face" };
  ply = ::write_ply ( fp, 2, elem_names, bin? PLY_BINARY_LE : PLY_ASCII );

  /* describe what properties go into the PlyVertex elements */
  describe_element_ply  ( ply, "vertex", nvert() );
  describe_property_ply ( ply, &plyvert_props[0] );
  describe_property_ply ( ply, &plyvert_props[1] );
  describe_property_ply ( ply, &plyvert_props[2] );

  /* describe PlyFace properties (just list of PlyVertex indices) */
  describe_element_ply  ( ply, "face",  ntetra() );
  describe_property_ply ( ply, &plyface_props[0] );

  header_complete_ply ( ply );

  /* set up and write the PlyVertex elements */
  put_element_setup_ply ( ply, "vertex" );

  for( int i=0; i<nvert(); i++ )
  {
    Vertex &vt = _G[i] ;
	  PlyPoint p;

	  p.x =  (float)vt.x();
	  p.y =  (float)vt.y();
	  p.z =  (float)vt.z();

    put_element_ply ( ply, ( void * ) &(p) );
  }

  /* set up and write the PlyFace elements */
  put_element_setup_ply ( ply, "face" );
  face.nverts = 4 ;

  for( int i = 0 ; i < ntetra() ; ++i )
  {
    v[0] = _V[i<<2]     ;
    v[1] = _V[i<<2 | 1] ;
    v[2] = _V[i<<2 | 2] ;
    v[3] = _V[i<<2 | 3] ;
    face.verts = v ;
    put_element_ply ( ply, ( void * ) &face );
  }

  close_ply ( ply );
  free_ply ( ply );

  printf(" %d vertices and %d triangles written\n", nvert(), ntetra() ) ;
}
//--------------------------------------------------------------//
