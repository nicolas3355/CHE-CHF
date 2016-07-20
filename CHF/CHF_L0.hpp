/**
* @file    CHF_L0.hpp
* @author  Marcos Lage     <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner  <thomas.lewiner@polytechnique.org>
* @author  Hélio  Lopes    <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmídia
* @date    14/02/2006
*
* @brief   CHF: A Scalable topological data structure for tetrahedral meshes
* @brief  (Level1 0)
*/
//--------------------------------------------------//
#ifndef  _CHF_L0_HPP_
#define  _CHF_L0_HPP_

#include <cfloat>
#include <vector>
#include "Vertex.hpp"

/** \brief Invalid integer index */
#define INV -2

/** \brief Tetrahedron id type */
typedef int TEid;

/** \brief Half Face id type */
typedef int HFid;

/** \brief Invalid vertex*/
static Vertex V_INV;

/** \brief Next of a half-edge (i,j)*/
static const HFid _N[4][4] = { {-1,3,1,2},{2,-1,3,0},{3,0,-1,1},{1,2,0,-1} };

/** \brief standart namespace definiton*/
using namespace std;
//--------------------------------------------------//
/** CHF data-structure for tetrahedral meshes
  * \brief CHF: Level 0*/
class CHF_L0
//--------------------------------------------------//
{
//-- CHF_L0 protected data.--//        
protected:
  /** \brief Number of vertices in the model */
  Vid _nvert;
  
  /** \brief Number of tetrahedrons in the model */
  TEid _ntetra;
  
  /** \brief Vertex Table */
  vector<Vid>    _V;
  
  /** \brief Geometry Table */
  vector<Vertex> _G;

public:
  /** \brief Default constructor.*/
  CHF_L0(): _nvert(0), _ntetra(0) {};
  
  /** \brief First constructor.
    * \param nv   -  const Vid.
    * \param ntet -  const TEid.*/
  CHF_L0(const Vid nv, const TEid ntet): _nvert(nv), _ntetra(ntet){_V.resize(4*ntet, -1); _G.resize( nv );}

  /** \brief Copy constructor
    * \param h  -  const CHF_L0 object.*/
  CHF_L0(const CHF_L0& h): _nvert(h.nvert()), _ntetra(h.ntetra()) { _V=h._V; _G=h._G; }

  /** \brief Destructor.*/
  virtual ~CHF_L0(){ _V.clear(); _G.clear(); }

public:
  /** \brief Access to the number of vertices of the model */
  inline const  Vid  nvert () const{ return _nvert;  }
	
  /** \brief Access to the number of tetrahedrons of the model */
  inline const TEid  ntetra() const{ return _ntetra; }

  /** \brief Access to the vertex of a half-face 
    * \param h - const HFid */
  inline const    Vid  V( const HFid h ) const { if( !hf_valid(h) )  return    INV; return _V[h] ; }
	
  /** \brief Access to a specific vertex of the model 
    * \param v - const Vid */
  inline const Vertex  G( const Vid  v ) const { if(  !v_valid(v) )  return V_INV; return _G[v] ; }

  /** \brief Tests if a vertex is valid
    * \param v - const Vid */
  inline const bool  v_valid( const Vid  v ) const { return ( v >= 0 && v < nvert() && (Vertex)_G[v]!= V_INV ); }
	
  /** \brief Tests if a tetrahedron is valid
    * \param t - const TEid */
  inline const bool te_valid( const TEid t ) const { return ( t >= 0 && t < ntetra() && hf_valid(t<<2) && hf_valid(t<<2 | 1) && hf_valid(t<<2 | 2)  && hf_valid(t<<2 | 3) ) ; }
	
  /** \brief Tests if a half-face is valid
    * \param h - const HFid */
  inline const bool hf_valid( const HFid h ) const { return ( h >= 0 && h < 4*ntetra() && _V[h] != INV );}
 
protected:
  /** \brief Sets to the number of vertices in the model  
    * \param  nvert - const Vid*/
  inline const void set_nvert ( const    Vid  nvert   ){ _nvert  = nvert; }
	
  /** \brief Sets to the number of tetrahedrons in the model 
    * \param ntetra*/
  inline const void set_ntetra( const TEid ntetra ){ _ntetra = ntetra; }

  /** \brief Sets the vertex of a half-edge
    * \param v - const HEid  
    * \param v - const  Vid */
  inline const void set_V( const HFid h, const Vid v ) {  if( h>=0 && h<4*ntetra()) _V[h]=v ; }

  /** \breaf Sets the geometry of a vertex in the model 
    * \param v - const Vid  
    * \param p - const Vertex*/
  inline const void set_G( const  Vid v, Vertex p ) { if(v>=0 && v<nvert()) _G[v]=p ; }

  /** \brief Sets a vertex as invalid
    * \param v - const Vid */
  inline const void  v_invalid( const   Vid v ) { if (  v_valid(v) ) _G[v] = V_INV ; return; }
	
  /** \brief Sets a tetrahedron as invalid
    * \param t - const TEid */
  inline const void te_invalid( const  TEid t ) { if ( te_valid(t) ) { hf_invalid(t<<2);  hf_invalid((t<<2) | 1); hf_invalid(t<<2 | 2);  hf_invalid(t<<2 | 3); } return;  }
	
  /** \brief Sets a half-face as invalid
    * \param h - const HFid */
  inline const void hf_invalid( const  HFid h ) { if ( hf_valid(h) ) _V[h] = INV; return; }

public:
  /** \brief Accesses the tetrahedron of a half-face 
    * \param h - const HFid*/
  inline const TEid  tetra( const HFid h ) const { if( !hf_valid(h) ) return INV; return h>>2 ;}
  /** \brief Accesses the next of a half-face 
    * \param h - const HFid*/
  inline const HFid nexthf( const HFid h ) const { if( !hf_valid(h) ) return INV; return (h&(~3)) | ((h+1) & 3) ;}
  /** \brief Accesses the mid of a half-face 
    * \param h - const HFid*/
  inline const HFid  midhf( const HFid h ) const { if( !hf_valid(h) ) return INV; return (h&(~3)) | ((h+2) & 3) ;}
  /** \brief Accesses the prev of a half-face 
    * \param h - const HFid*/
  inline const HFid prevhf( const HFid h ) const { if( !hf_valid(h) ) return INV; return (h&(~3)) | ((h+3) & 3) ;}

  /** \brief Accesses the next of a half-edge he in a half-face h 
    * \param h  - const HFid
    * \param he - const HFid*/
  inline const pair<HFid,HFid> nexthe( const HFid hf, const HFid he ) const 
  { 
    if( !hf_valid(hf) || !hf_valid(he) || (hf>>2) != (he>>2) || he == hf ) { cout<< "CHF ERRO: nexthe" << endl; return make_pair(INV,INV); }
     return ( make_pair (hf, ( hf&(~3) ) | _N[ he & 3 ][ hf & 3 ] ) ) ;
  }
  /** \brief Accesses the prev of a half-edge he in a half-face h
    * \param h  - const HFid
    * \param he - const HFid*/
  inline const pair<HFid,HFid> prevhe( const HFid hf, const HFid he ) const 
  { 
    if( !hf_valid(hf) || !hf_valid(he) || (hf>>2) != (he>>2) || he == hf ) { cout<< "CHF ERRO: prevhe" << endl; return make_pair(INV,INV); } 
    return ( make_pair (hf, ( hf&(~3) ) | _N[ hf & 3 ][ he & 3 ] ) );
  }
  /** \brief Accesses the mate of a half-edge he in a half-face h 
    * \param h  - const HFid
    * \param he - const HFid*/
  inline const pair<HFid,HFid> matehe( const HFid hf, const HFid he ) const 
  { 
    if( !hf_valid(hf) || !hf_valid(he) || (hf>>2) != (he>>2) || he == hf ) { cout<< "CHF ERRO: matehe" << endl; exit(1);}//return make_pair(INV,INV); } 
    return ( make_pair(prevhe(hf,he).second, nexthe(hf,he).second) );
  }
  /** \brief Accesses the four half-faces of a tetrahedron 
    * \param h  - const HFid
    * \param h1 - HFid&
    * \param h2 - HFid&
    * \param h3 - HFid& */
   void neighbors (const HFid h, HFid &h1, HFid &h2, HFid &h3) const
   {
     if( !hf_valid(h) ) { h1 = h2 = h3 = INV;  return ; }
     TEid t = h & (~3) ;
     switch( h & 3 )
     {
		case 0 :  h1 = t | 1 ;  h2 = t | 2 ;  h3 = t | 3 ;  return ;
		case 1 :  h1 =   t   ;  h2 = t | 3 ;  h3 = t | 2 ;  return ;
		case 2 :  h1 = t | 3 ;  h2 =   t   ;  h3 = t | 1 ;  return ;
		case 3 :  h1 = t | 2 ;  h2 = t | 1 ;  h3 =   t   ;  return ;
	 }
   }

public:
  /** \brief Computes the vertices in the star of a vertex 
    * \param v - const Vid */
  virtual vector<Vid>  R_00( const Vid v  );
  /** \brief Computes the tetrahedrons in the star of a vertex 
    * \param v - const Vid */
  virtual vector<TEid> R_03( const Vid v  );
  /** \brief Computes the vertices in the star of an edge 
    * \param a - const Vid  
    * \param b - const Vid */
  virtual vector<Vid>  R_10( const Vid a, const Vid b  );
  /** \brief Computes the tetrahedrons in the star of an edge 
    * \param a - const Vid  
    * \param b - const Vid */
  virtual vector<TEid> R_13( const Vid a, const Vid b  );
  /** \brief Computes the tetrahedrons in the star of a tetrahedron 
    * \param r - const TEid */
  virtual vector<TEid> R_33( const TEid t );
 
  /** \brief Gets the model bouding_box
    * \param min - float*.
    * \param max - float*. */
  void bounding_box( float *min, float *max );
  /** \brief Legalizes the model
    * \param min - float*.
    * \param max - float*. */
  void legalize_model( float *min, float *max );
  /** \brief Assigns the scalar field
    * \param eq - const char*. */
  void scalar_field ( const char* eq );
  
  /** \brief Checks mesh validation */
  void check ();

  /** \brief Draws the mesh in wireframe 
    * \param in= true - const bool*/
  virtual void draw_wire ( const int t=4 );
  /** \brief Draws the vertices of the mesh 
    * \param t= 0 - const bool*/
  virtual void draw_vert ( const int t=1, const int p=false );
  
  /** \brief Reads the model from a file
    * \param fn - const char*. */
  void  read_ply ( const char* fn );

  /** \brief Writes the model in a file
    * \param fn - const char* 
    * \param bin = false - bool*/ 
  void  write_ply ( const char* fn, bool bin= false );
};
#endif
//--------------------------------------------------------------//
