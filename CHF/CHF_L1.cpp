/**
* @file    CHF_L1.cpp
* @author  Marcos Lage      <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner   <thomas.lewiner@polytechnique.org>
* @author  Hélio  Lopes     <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmídia
* @date    14/02/2006
*
* @brief   CHF: A Scalable topological data structure for tetrahedral meshes
* @brief  (Level 1)
*/
//--------------------------------------------------//
#include <set>
#include <map>
#include <stack>
#include <GL/glut.h>

#include "colorramp.h"  /**< Gl color maps*/
#include "CHF_L1.hpp"  /**< Level 1 inheritance*/
											
using namespace std;	
//--------------------------------------------------//
vector<Vid> CHF_L1::R_00(const Vid v)
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
  for(HFid i=0; i<ntetra()<<2; ++i)
    if( V(i) == v ) {
	  // includes the tetrahedron in the stack.
	  t = tetra(i);	s.push(t);
	  break;
	}
	
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
vector<TEid> CHF_L1::R_03(const Vid v)
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

  //Finds the first tetrahedron that contains v.
  for(HFid i=0; i<ntetra()<<2; ++i)
    if( V(i) == v ) {
	  //includes the tetrahedron in the stack.
	  t = tetra(i);	s.push(t);
	  break;
	}

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
vector<Vid> CHF_L1::R_10(const Vid a, const Vid b)
//--------------------------------------------------//
/** Computes the vertices in the star of a given edge.*/
{
  set<Vid>   sstar;
  vector<Vid> star;
	
  pair<HFid, HFid>  start = make_pair(INV, INV), runcw = make_pair(INV, INV), runcc = make_pair(INV, INV);
	
  if( !v_valid(a) || !v_valid(b) || a == b) { star.push_back(INV); return star; }
	
  //Search the first edge
  for(HFid i=0; i<4*ntetra(); ++i)
  {
    if((V(i) == a) && ( (V(nexthf(i)) == b) || (V(midhf(i)) == b) || (V(prevhf(i)) == b) ))
	{
	  HFid lives = 0;
			
	  if(V(nexthf(i)) != b &&  V(nexthe(nexthf(i), i).second) == b) lives = nexthf(i);
	  else if(V(midhf(i))  != b &&  V(nexthe(midhf(i),  i).second) == b) lives = midhf(i); 
	  else if(V(prevhf(i)) != b &&  V(nexthe(prevhf(i), i).second) == b) lives = prevhf(i); 
				
	  start  = runcw = make_pair(lives, i);
	  runcc = matehe( start.first, start.second );
	  break;
	}
  }

  if(start == make_pair(INV, INV)) { star.push_back(INV);  return star; }

  //Walk arround the edge --"Clockwise"
  do
  {
    sstar.insert( V(runcw.first) );
	runcw = matehe(runcw.first, runcw.second);
	runcw = radialhe(runcw.first, runcw.second);
  }
  while((runcw.first != -1) && (runcw.first != start.first));

  //if not 360º
  if(runcw.first == -1)
  {
    start  = runcc;

    //Walk arround the edge --"CounterClockwise"
    do
    {
      sstar.insert( V(runcc.first) );
	  runcc = matehe(runcc.first, runcc.second);
	  runcc = radialhe( runcc.first, runcc.second);
	}
	while( runcc.first != -1 );
  }
	
  // Copies in star
  for(set<Vid>::const_iterator iter = sstar.begin(); iter != sstar.end(); ++iter )
    star.push_back( *iter );

  return star;
}
//--------------------------------------------------//
vector<TEid> CHF_L1::R_13(const Vid a, const Vid b)
//--------------------------------------------------//
/** Computes the vertices in the star of a given edge.*/
{
  set<TEid> sstar;
  vector<TEid>   star;
	
  pair<HFid, HFid>  start, runcw, runcc;

  if( !v_valid(a) || !v_valid(b) || a ==b) { star.push_back(INV); return star; }
	
  //Search the first edge
  for(HFid i=0; i<4*ntetra(); ++i)
  {
    if((V(i) == a) && ( (V(nexthf(i)) == b) || (V(midhf(i)) == b) || (V(prevhf(i)) == b) ))
	{		
	  HFid lives = 0;
	
	  if(V(nexthf(i)) != b &&  V(nexthe(nexthf(i), i).second) == b) lives = nexthf(i);
	  else if(V(midhf(i))  != b &&  V(nexthe(midhf(i),  i).second) == b) lives = midhf(i); 
	  else if(V(prevhf(i)) != b &&  V(nexthe(prevhf(i), i).second) == b) lives = prevhf(i); 
					
	  start  = runcw = make_pair(lives, i);
	  runcc = matehe( start.first, start.second );
	  break;
	}
  }

  if(start == make_pair(INV, INV)) { star.push_back(INV); return star; }
	
  //Walk arround the edge --"Clockwise"
  do
  {
    sstar.insert( tetra(runcw.first) );
	runcw = matehe(runcw.first, runcw.second);
	runcw = radialhe(runcw.first, runcw.second);
  }
  while((runcw.first != -1) && (runcw.first != start.first));
	
  //if not 360º
  if(runcw.first == -1)
  {
    start  = runcc;

    //Walk arround the edge --"CounterClockwise"
	do
	{
	  sstar.insert( tetra(runcc.first) );
	  runcc = matehe(runcc.first, runcc.second);
	  runcc = radialhe( runcc.first, runcc.second);
	}
	while( runcc.first != -1 );		
  }
	
  // Copies in star
  for(set<Vid>::const_iterator iter = sstar.begin(); iter != sstar.end(); ++iter )
    star.push_back( *iter );

  return star;
}
//--------------------------------------------------//
vector<TEid> CHF_L1::R_33(const TEid t)
//--------------------------------------------------//
/** Computes the tetrahedrons incidents to a tetrahedron.*/
{
  vector<TEid> star;

  TEid t0 = 0,  t1 = 0,  t2 = 0,  t3 = 0;
  if( !te_valid(t) ) { star.push_back(INV); return star;}

  t0 = tetra( O(t<<2)      ); t1 = tetra( O((t<<2) + 1) ); 
  t2 = tetra( O((t<<2) + 2) ); t3 = tetra( O((t<<2) + 3) );	

  star.push_back(t0), star.push_back(t1);
  star.push_back(t2), star.push_back(t3);

  return star;
}
//--------------------------------------------------//
void CHF_L1::create_O()
//--------------------------------------------------//
/** Create adjacency relations between the tetrahedrons of the mesh*/
{
  TEid t;
  HFid n, a[3];

  map< pair<pair<Vid,Vid>,Vid>  ,HFid> adj;
  map< pair<pair<Vid,Vid>,Vid>  ,HFid>::iterator pos ;

  _O.clear();
  _O.resize(ntetra()<<2,-1);

  for ( HFid h= 0 ; h < ntetra()<<2 ; ++h )
  {
    if( !hf_valid( h ) ) continue ;

    neighbors( h, a[0], a[1], a[2] ) ;
    a[0] = V( a[0] ) ; a[1] = V( a[1] ) ; a[2] = V( a[2] ) ;

    for (HFid m = 0; m < 3; ++m)
    {
      HFid min = m;
      for (n = m+1; n < 3; n++)
        if (a[n] < a[min]) min = n;
      t = a[m]; a[m] = a[min]; a[min] = t; 
    }

    pos = adj.find( make_pair( make_pair(a[0],a[1]), a[2] ) ) ;

    if( pos != adj.end() )
    {
      _O[h] = pos->second ;
      _O[O(h)] = h ;
      adj.erase( pos ) ;
    }
    else
      adj[ make_pair( make_pair(a[0],a[1]), a[2] ) ]=h ;
  }
  adj.clear() ;
}
//--------------------------------------------------//
void CHF_L1::compute_normals()
//--------------------------------------------------//
{
  for(HFid i=0; i< 4*ntetra(); ++i)
  {
    if( !(hf_valid(i)) || O(i) != -1 ) continue;

	float norm[3] ;

	HFid h0 = 0,  h1 = 0, h2 = 0;

	neighbors( i, h0, h1, h2 );

	Vertex &v0 = _G[ V( h0 ) ];
	Vertex &v1 = _G[ V( h1 ) ];
	Vertex &v2 = _G[ V( h2 ) ];

	Vertex::normal( v0,v1,v2, norm );
	
	v0.set_nx( v0.nx()+norm[0] );
	v0.set_ny( v0.ny()+norm[1] );
	v0.set_nz( v0.nz()+norm[2] );
	
	v1.set_nx( v1.nx()+norm[0] );
	v1.set_ny( v1.ny()+norm[1] );
	v1.set_nz( v1.nz()+norm[2] );
	
	v2.set_nx( v2.nx()+norm[0] );
	v2.set_ny( v2.ny()+norm[1] );
	v2.set_nz( v2.nz()+norm[2] );
  }
}
//--------------------------------------------------//
void CHF_L1::change_orientation(TEid t)
//--------------------------------------------------//
/** Changes a tetrahedron orientation  */
{  
  /** Pick 2 half faces*/
  HFid h1 = t<<2 | 1 ;
  HFid h2 = t<<2 | 2 ;

  /** Exchange the geometrical vertices*/
  Vid v  =   V(h1)  ;
  set_V( h1, V(h2) );
  set_V( h2,   v   );

  /** Exchange mates*/
  HFid m =      O(h1)  ;
  set_O(   h1 , O(h2) );
  set_O(   h2 ,   m   );
  set_O( O(h1),   h1  );
  set_O( O(h2),   h2  );
}
//--------------------------------------------------//
const bool CHF_L1::orient_check(const HFid c,const HFid t)
//--------------------------------------------------//
/** Checks the orientation of a tetrahedron in relation to its opposite*/
{
  Vid v;
  HFid hc = 0, ht = 0;

  hc = nexthf(c);
  v  = V(hc);

  for(int i= 4*(tetra(t)); i<4*(tetra(t)+1) ; ++i)
    if( V(i) == v ) ht = i;

  return  V(nexthe(c, hc).second) == V(prevhe(t, ht).second)
    &&    V(prevhe(c, hc).second) == V(nexthe(t, ht).second);
}
//--------------------------------------------------//
void CHF_L1::orient()
//--------------------------------------------------//
/** Orients the mesh coherently*/
{
  if( ntetra() < 1 ) return ;

  /** Stack of half-faces to check*/
  stack<HFid> st ; 
  /** Vector of visited tetras*/
  vector<bool> visited ;  visited.resize( ntetra(), false ) ;

  /** orient tetra 0*/
  const Vertex &v0 = G( V(0) );  const Vertex &v1 = G( V(1) );
  const Vertex &v2 = G( V(2) );  const Vertex &v3 = G( V(3) );

  float test = Vertex::signed_tetra_volume(v0, v1, v2, v3);

  if( test < 0 ) change_orientation( 0 ) ;

  /** push half faces of 0*/
  st.push( 0 ) ;  st.push( 1 ) ;
  st.push( 2 ) ;  st.push( 3 ) ;
  visited[0] = true ;
  
  while( !st.empty() )
  {
    HFid hf = st.top() ;
    st.pop() ;

    /** avoid null face*/
    HFid mf = O(hf) ;
    if( mf == -1 )   continue ;

    /** avoid loops*/
    TEid t = tetra(mf) ;
    if( visited[t] ) continue ;

    /** repairs*/
    if( !orient_check( hf,mf ) ) change_orientation( t ) ;

    /** marks as visited*/
    visited[t] = true ;

    /** push half faces of t*/
    t <<= 2 ;
    st.push(   t   ) ;  st.push( t | 1 ) ;
    st.push( t | 2 ) ;  st.push( t | 3 ) ;
  }
}
//--------------------------------------------------//
void CHF_L1::check()
//--------------------------------------------------//
/** Checks adjacency structure */
{ 
  CHF_L0::check();
  
  for(HFid i=0; i<(ntetra()<<2); i++)
  {
    if(O(i)!= -1 && O(i) != INV)
    {
      if( i!=O(O(i)) )
      {
        cout << "CHF_L1::Check ERRO: "<< i << " != O( O(" << i << ") ) = " << O(O(i)) << endl;
        return;
      }
    }
  }

  for(HFid i=0; i<(ntetra()<<2); i++)
  {
    if(O(i)!= -1)
    {
      if( !orient_check( i, O(i) ) )
      {
        cout << "CHF_L1::Check ERRO: Tetrahedron "<< i/4 << " bad oriented " << endl;
        return;
      }
    }
  }
  cout << "CHF_L1::Check Ok" << endl;
  return;
}
//--------------------------------------------------//
void CHF_L1::draw_smooth(const int t)
//--------------------------------------------------//
/** Draws the boundary surface */
{
  if( t != 8 && t !=9 )
  {
    cout << "CHF_L1::draw_smooth ERRO." << endl;
	return;
  }

  ColorRamp c;
  c.set_GLcolor(0.6,1.0,COLOR_RAINBOW, 1.0);

  if ( t == 8 ) glShadeModel(GL_FLAT);
  if ( t == 9 ) glShadeModel(GL_SMOOTH);

  glBegin( GL_TRIANGLES );
  {
    for(int i=0; i< ntetra(); ++i)
	{
	  HFid  h0 = (i<<2),
	        h1 = (i<<2 | 1),
			h2 = (i<<2 | 2),
			h3 = (i<<2 | 3);
	    
	  const Vertex &v0 = G( V(h0) );
	  const Vertex &v1 = G( V(h1) );
	  const Vertex &v2 = G( V(h2) );
	  const Vertex &v3 = G( V(h3) );

	  if( O(h0) == -1 )
	  {
	    if (t == 8)
		{
	      float n[3];
		  Vertex::normal( v2, v3, v1, n );
		  glNormal3f( n[0], n[1], n[2]);
		}
		if( t == 9 ) glNormal3f( v2.nx(), v2.ny(), v2.nz());
		  glVertex3f( v2.x(),  v2.y(),  v2.z());
				
		if( t == 9 )glNormal3f( v3.nx(), v3.ny(), v3.nz());
		  glVertex3f( v3.x(),  v3.y(),  v3.z());
				
		if( t == 9 )glNormal3f( v1.nx(), v1.ny(), v1.nz());
		  glVertex3f( v1.x(),  v1.y(),  v1.z());
	  }

	  if( O(h1) == -1 )
	  {
        if (t == 8)
		{
		  float n[3];
		  Vertex::normal( v2, v0, v3, n );
		  glNormal3f( n[0], n[1], n[2]);
		}
		if( t == 9 )glNormal3f( v3.nx(), v3.ny(), v3.nz());
		  glVertex3f( v3.x(),  v3.y(),  v3.z());

		if( t == 9 )glNormal3f( v0.nx(), v0.ny(), v0.nz());
		  glVertex3f( v0.x(),  v0.y(),  v0.z());

		if( t == 9 )glNormal3f( v2.nx(), v2.ny(), v2.nz());
		  glVertex3f( v2.x(),  v2.y(),  v2.z());
	  }

	  if( O(h2) == -1 )
	  {
        if (t == 8)
		{
		  float n[3];
		  Vertex::normal( v1, v3, v0, n );
		  glNormal3f( n[0], n[1], n[2]);
		}
		if( t == 9 )glNormal3f( v1.nx(), v1.ny(), v1.nz());
		  glVertex3f( v1.x(),  v1.y(),  v1.z());
				
		if( t == 9 )glNormal3f( v3.nx(), v3.ny(), v3.nz());
		  glVertex3f( v3.x(),  v3.y(),  v3.z());
				
		if( t == 9 )glNormal3f( v0.nx(), v0.ny(), v0.nz());
		  glVertex3f( v0.x(),  v0.y(),  v0.z());
	  }

	  if( O(h3) == -1 )
	  {
        if (t == 8)
		{
		  float n[3];
		  Vertex::normal( v2, v1, v0, n );
		  glNormal3f( n[0], n[1], n[2]);
		}
		
		if( t == 9 )glNormal3f( v2.nx(), v2.ny(), v2.nz());
		  glVertex3f( v2.x(),  v2.y(),  v2.z());
				
		if( t == 9 )glNormal3f( v0.nx(), v0.ny(), v0.nz());
		  glVertex3f( v0.x(),  v0.y(),  v0.z());
				
		if( t == 9 )glNormal3f( v1.nx(), v1.ny(), v1.nz());
		  glVertex3f( v1.x(),  v1.y(),  v1.z());
	  }
	}
  }
  glEnd();
}
//--------------------------------------------------------------//
