/**
* @file    CHF_L3.cpp
* @author  Marcos Lage      <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner   <thomas.lewiner@polytechnique.org>
* @author  Hélio  Lopes     <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmídia
* @date    14/02/2006
*
* @brief   CHF: A Scalable topological data structure for tetrahedral meshes
* @brief  (Level 3)
*/
//--------------------------------------------------//
#include <map>
#include <gl/glut.h>

#include "CHF_L3.hpp"   /**< Level 2 inheritance*/
#include "colorramp.h"	/**< Gl color maps*/

using namespace std;
//--------------------------------------------------//
Bid CHF_L3::get_bS( Bid i )
//--------------------------------------------------//
/** Sets bS and returns the value.*/
{
  Bid b =  bS(i) ;
  if( b < 0 || b == i ) return b ;

  b = get_bS(b) ;
  set_bS( i,b );

  return b ;
}
//--------------------------------------------------//
void CHF_L3::create_bS()
//--------------------------------------------------//
/** Creates the bS container*/
{
  _bS.clear () ;
  _bS.resize( nvert() ) ;

  for(Vid v = 0 ; v < nvert() ; ++v )
  {
    if( !v_valid(v) ) continue;

    if( O( VH(v) ) < 0 )
      set_bS( v, v ) ;
    else
      set_bS( v, -1) ;
  }

  for( TRid j = 0 ; j < bntrig() ; ++j )
  {
    if( !tr_valid(j) ) continue;

    Bid b0 = get_bS( bV( 3*j ) ) ;
    Bid b1 = get_bS( bV(3*j+1) ) ;
    Bid b2 = get_bS( bV(3*j+2) ) ;
    if( b0 < 0 || b1 < 0 || b2 < 0 )
      continue ;
    if( b0 < b1 )
    {
      if( b0 < b2 )
      { // b0 min
        set_bS( b1, b0 ) ;
        set_bS( b2, b0 ) ;
      }
      else
      { // b2 min
        set_bS( b0, b2 ) ;
        set_bS( b1, b2 ) ;
      }
    }
    else
    {
      if( b1 < b2 )
      { // b1 min
        set_bS( b0, b1 ) ;
        set_bS( b2, b1 ) ;
      }
      else
      { // b2 min
        set_bS( b0, b2 ) ;
        set_bS( b1, b2 ) ;
      }
    }
  }

  _bnsurf = 0 ;
  map< Bid, Bid > corresp ;
  for( Vid v = 0 ; v < nvert() ; ++v )
  {
    if( !v_valid(v) ) continue;

    Bid b = get_bS( v ) ;
    if( b < 0 ) continue ;
    map< Bid,Bid >::iterator it = corresp.find( b ) ;
    if( it == corresp.end() ) corresp.insert( it, make_pair( b, _bnsurf++ ) ) ;
  }

  for( Vid v = 0 ; v < nvert() ; ++v )
  {
    if( !v_valid(v) ) continue;

    Bid b = bS(v) ;
    if( b < 0 ) continue ;
     set_bS( v, corresp[b]) ;
  }

	cout << "CHF_L3::create_bS: " << bnsurf() << " boundary surfaces found." << endl;
}
//--------------------------------------------------//
void CHF_L3::create_bV()
//--------------------------------------------------//
/** Creates the bV container*/
{
    _bV.clear();
	_bntrig = 0;

  for(HFid i=0; i< 4*ntetra(); ++i)
  {
    if( !hf_valid(i) ) continue;
    if(  O(i) != -1  ) continue ;

    HFid hf1 = 0, hf2 = 0, hf3 = 0;
    neighbors( i, hf1,hf2,hf3 ) ;

    _bV.push_back( V(hf1) );
    _bV.push_back( V(hf2) );
    _bV.push_back( V(hf3) );

		_bntrig ++;
  }
}
//--------------------------------------------------//
void CHF_L3::create_bO()
//--------------------------------------------------//
/** Creates the bO container*/
{
  Vid a,b;
  map<pair<int,int>,HFid> adj;
  map<pair<int,int>,HFid>::iterator pos ;

  _bO.clear();
  _bO.resize( 3*bntrig(), -1 );

  for ( HEid c = 0 ; c < 3*bntrig() ; ++c )
  {
    if( !he_valid(c) ) continue;

    a = bV( bnexthe(c)) ;  b = bV( bprevhe(c)) ;

    if( b < a) { Vid tmp = a ; a = b ; b = tmp ; }

    pos = adj.find( pair<int,int>(a,b) ) ;

    if( pos != adj.end() )
    {
      set_bO(c, pos->second) ;
      set_bO( bO(c), c) ;

      adj.erase( pos ) ;
    }
    else
    {
      adj[pair<int,int>(a,b)] = c ;
    }
  }
  adj.clear() ;
}
//--------------------------------------------------//
void CHF_L3::compute_normals()
//--------------------------------------------------//
{
  if( _bV.size() == 0 ) return;

  for( Vid i=0; i< nvert(); ++i )
  {
	_G[i].set_nx(0);
	_G[i].set_ny(0);
	_G[i].set_nz(0);
  }

  for(TRid i=0; i< bntrig(); ++i)
  {
    if( !(tr_valid(i)) ) continue;
		
	float norm[3] ;

	Vertex &v0 = _G[ bV( 3*i ) ];
	Vertex &v1 = _G[ bV( 3*i + 1 ) ];
	Vertex &v2 = _G[ bV( 3*i + 2 ) ];

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
const bool CHF_L3::borient_check(const HEid c, const HEid t)
//--------------------------------------------------//
{
  Vid v1c, v2c, v1t, v2t;
  
  v1c = bV( bprevhe(c) );
  v2c = bV( bnexthe(c) );
  
  v1t = bV( bprevhe(t) );
  v2t = bV( bnexthe(t) );

  if( v1c == v2t && v1t == v2c ) return true;
  return false;
}
//--------------------------------------------------//
void CHF_L3::check()
//--------------------------------------------------//
{
  CHF_L2::check();

  for(HEid i=0; i< 3*bntrig(); ++i)
  {
    if( bO(i) != -1 )
    {
      if( i!=bO(bO(i)) )
      {
        cout << "CHF_L3::Check ERRO: " << i <<" != bO(bO(" << i <<")). "<< endl;
        return;
      }
    }

    if( !borient_check( i, bO(i) ) )
    {
      cout << "CHF_L3::Check ERRO: Triangle " << (int)(i/3) << " bad oriented." << endl;
      return;
    }
    
    if( bV(i) >= nvert() )
    {
      cout << "CHF_L3::Check ERRO: bV(" << i << ") >= nvert(). " << endl;
      return;
    }
  }

  for(Vid i=0; i< nvert(); i++)
  {
    if( bS(i) >= bnsurf() )
    {
      cout << "CHF_L3::Check ERRO: bS(" << i << ") >= bnsurf(). " << endl;
      return;
    }
  }
  cout << "CHF_L3::Check Ok" << endl;
  return;
}
//--------------------------------------------------//
void CHF_L3::draw_smooth(const int t)
//--------------------------------------------------//
/** Draws the boundary surface */
{
  if( t != 8 && t != 9 && t != 10 ) 
  {
    cout << "CHF_L3::draw_smooth ERRO" << endl;
    return;
  }

  ColorRamp c;
	
  if ( t == 8 || t ==9 )  c.set_GLcolor(0.6,1.0,COLOR_RAINBOW, 1.0);
  if ( t == 8 )           glShadeModel(GL_FLAT);
  if ( t == 9 || t == 10 )glShadeModel(GL_SMOOTH);
	
  glBegin( GL_TRIANGLES );
  {
  	for(int i=0; i< bntrig(); ++i)
	{
	  HEid h0 = (3*i), h1 = (3*i + 1), h2 = (3*i + 2);

	  if ( t == 10 ) c.set_GLcolor( bS(bV(h0)), bnsurf() ,COLOR_RAINBOW, 1.0);

	  const Vertex &v0 = G( bV(h0) );
	  const Vertex &v1 = G( bV(h1) );
	  const Vertex &v2 = G( bV(h2) );

	  if( t == 8 )
	  {
		float n[3];
	 	Vertex::normal( v0, v1, v2, n );
	  	glNormal3f( n[0], n[1], n[2]);
	  }
    
	  if ( t == 9 || t == 10 )glNormal3f( v0.nx(), v0.ny(), v0.nz());
	  glVertex3f( v0.x(),  v0.y(),  v0.z() );
			
	  if ( t == 9 || t == 10 )glNormal3f( v1.nx(), v1.ny(), v1.nz());
	  glVertex3f( v1.x(),  v1.y(),  v1.z() );

	  if ( t == 9 || t == 10 )glNormal3f( v2.nx(), v2.ny(), v2.nz());
      glVertex3f( v2.x(),  v2.y(),  v2.z() );
	}
  }
  glEnd();
}
//--------------------------------------------------//
