/**
* @file    CHF_GLUI.cpp
* @author  Marcos Lage      <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner <thomas.lewiner@polytechnique.org>
* @author  Hélio  Lopes       <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmídia
* @version 0.1
* @date    16/05/2005
*
* @brief   CHF: A Scalable topological data structure for tetrahedral meshes
* @brief  (Test GLUI program)
*/

#include <GL/glui.h>
#include <GL/glut.h>

#include "gl2ps.h"
#include "CHF_L3.hpp"

#define         VERT  1
#define      IN_VERT  2
#define       B_VERT  3
#define        FIELD  4
#define         WIRE  5
#define      IN_WIRE  6
#define       B_WIRE  7
#define         FLAT  8
#define       SMOOTH  9
#define        BOUND  10

GLUI             *glui, *glui_r;

GLUI_Rotation    mouse_rt , *object_rt, *light_rt ;

GLUI_Translation mouse_xy ,  *object_xy;
GLUI_Translation mouse_zm, *object_zm;

GLUI_RadioGroup  *sgroup, *egroup, *vgroup;

GLUI_RadioButton *sb1, *sb2, *sb3,
							 *eb1, *eb2, *eb3,
							 *vb1, *vb2, *vb3, *vb4;

GLUI_Spinner *simplexid;

int fileid= 0;
int view_w=745,
	  view_h=650,
	  window_id,
	  ncalls,
      motion_type;

int  level    =    0;
int  sdraw  =    1,  edraw = 0,  vdraw = 0, dim =0;
int        srf = true,
           pnt = true,
          edg = true,
	   vtype = false,
	   estar = false;

float view_rotate[16]  = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
float         obj_pos[3]	= { 0.0, 0.0, 0.0 };
float                 scale  = 1.0;

char  file_ply[1024], file_eps[1024];
char *file_list[1024]= {         "test.ply",          "8tetra.ply",                "102.ply",                 "104.ply",
                                        "apparatus.ply",              "bird.ply",            "blunto.ply",                 "box.ply",
									  "bunny5188.ply",   "bunny9553.ply",  "bunny19047.ply",   "bunny48810.ply",
										     "deltao.ply",    "disconnect.ply",                "dog.ply",              "epcot.ply",
                                            "fandisk.ply",       "gargoyle.ply",               "gear.ply",     "hand28792.ply",
                                     "hand28793.ply",             "hand.ply",              "heart.ply",           "mesh1.ply",
										    "mesh2.ply",          "mesh3.ply",        "molecule.ply",             "pmdc.ply",
												"skull.ply",             "snail.ply",          "sphere0.ply",         "sphere1.ply",
											  "spine.ply",               "spx.ply",             "torus1.ply",           "torus2.ply",
												 "trex.ply",             "trice.ply"  };

CHF_L0 ch0;  CHF_L1 ch1;
CHF_L2 ch2;  CHF_L3 ch3;

GLfloat light0_rotation[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
GLfloat light1_rotation[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };

int eid = 0, nedges = 0;
vector<int> star;
//----------------------------------------------------------------//
vector<int> test_vstar(int eid,  int dim)
//----------------------------------------------------------------//
{
  vector<int> st;
  int VA[6] = {0,0,0,1,1,2};
  int VB[6] = {1,2,3,2,3,3};

  Vid va = 0, vb =0;

    if(level == 0) {
		   va = ch0.V(4*(eid/6) + VA[eid%6]);
		   vb = ch0.V(4*(eid/6) + VB[eid%6]);
		}
    if(level == 1) {
		   va = ch1.V(4*(eid/6) + VA[eid%6]);
		   vb = ch1.V(4*(eid/6) + VB[eid%6]);
		}
    if(level == 2) {
		   va = ch2.V(4*(eid/6) + VA[eid%6]);
		   vb = ch2.V(4*(eid/6) + VB[eid%6]);
		}
    if(level == 3) {
		   va = ch3.V(4*(eid/6) + VA[eid%6]);
		   vb = ch3.V(4*(eid/6) + VB[eid%6]);
		}

  if( dim == 0 ){
    if(level == 0) st = ch0.R_10(va, vb);
    if(level == 1) st = ch1.R_10(va, vb);
    if(level == 2) st = ch2.R_10(va, vb);
    if(level == 3) st = ch3.R_10(va, vb);
  }
  else{
    if(level == 0) st = ch0.R_13(va, vb);
    if(level == 1) st = ch1.R_13(va, vb);
    if(level == 2) st = ch2.R_13(va, vb);
    if(level == 3) st = ch3.R_13(va, vb);
  }

  return st;
}
//----------------------------------------------------------------//
void draw_vstar(int dim)
//----------------------------------------------------------------//
{
  //cout << "star size" << star.size() << endl;
  for(int i=0; i< static_cast<int>(star.size()); ++i)
 {
	//--- Vertex ---//
	if( dim == 0 )
	{
		  Vertex v1;

		  if(level == 0) v1 = ch0.G(star[i]);
		  if(level == 1) v1 = ch1.G(star[i]);
		  if(level == 2) v1 = ch2.G(star[i]);
		  if(level == 3) v1 = ch3.G(star[i]);

		  glColor3f(1.0,0.3,0.1);
		  glPointSize(2.0);
		  glBegin(GL_POINTS);
			  glVertex3d(v1.x() ,  v1.y() ,  v1.z() );
		  glEnd();
    }

    //--- tetras ---//
    else
	{
      Vertex v0, v1, v2, v3;

      if(level == 0)
	  {
        v0 = ch0.G( ch0.V(4*star[i])   );
        v1 = ch0.G( ch0.V(4*star[i]+1) );
        v2 = ch0.G( ch0.V(4*star[i]+2) );
        v3 = ch0.G( ch0.V(4*star[i]+3) );
      }
      if(level == 1)
	  {
        v0 = ch1.G( ch1.V(4*star[i])   );
        v1 = ch1.G( ch1.V(4*star[i]+1) );
        v2 = ch1.G( ch1.V(4*star[i]+2) );
        v3 = ch1.G( ch1.V(4*star[i]+3) );
     }
      if(level == 2)
	 {
        v0 = ch2.G( ch2.V(4*star[i])   );
        v1 = ch2.G( ch2.V(4*star[i]+1) );
        v2 = ch2.G( ch2.V(4*star[i]+2) );
        v3 = ch2.G( ch2.V(4*star[i]+3) );
      }
      if(level == 3)
	 {
        v0 = ch3.G( ch3.V(4*star[i])   );
        v1 = ch3.G( ch3.V(4*star[i]+1) );
        v2 = ch3.G( ch3.V(4*star[i]+2) );
        v3 = ch3.G( ch3.V(4*star[i]+3) );
     }

      glColor3f(1.0,0.0,0.0);
		  glBegin(GL_TRIANGLES);
			  // Triangle 0
	 		  glNormal3d(v1.nx(),  v1.ny(),  v1.nz());
			  glVertex3d(v1.x()  , v1.y(),   v1.z() );

			  glNormal3d(v2.nx(),  v2.ny(),  v2.nz());
			  glVertex3d(v2.x()  , v2.y(),   v2.z() );

			  glNormal3d(v3.nx(),  v3.ny(),  v3.nz());
			  glVertex3d(v3.x()  , v3.y(),   v3.z() );

			  // Triangle 1
	 		  glNormal3d(v0.nx(),  v0.ny(),  v0.nz());
			  glVertex3d(v0.x()  , v0.y(),   v0.z() );

			  glNormal3d(v2.nx(),  v2.ny(),  v2.nz());
			  glVertex3d(v2.x()  , v2.y(),   v2.z() );

			  glNormal3d(v1.nx(),  v1.ny(),  v1.nz());
			  glVertex3d(v1.x()  , v1.y(),   v1.z() );

			  // Triangle 2
	 		  glNormal3d(v0.nx(),  v0.ny(),  v0.nz());
			  glVertex3d(v0.x()  , v0.y(),   v0.z() );

			  glNormal3d(v1.nx(),  v1.ny(),  v1.nz());
			  glVertex3d(v1.x()  , v1.y(),   v1.z() );

			  glNormal3d(v3.nx(),  v3.ny(),  v3.nz());
			  glVertex3d(v3.x()  , v3.y(),   v3.z() );

			  // Triangle 3
			  glNormal3d(v0.nx(),  v0.ny(),  v0.nz());
			  glVertex3d(v0.x()  , v0.y(),   v0.z() );

			  glNormal3d(v3.nx(),  v3.ny(),  v3.nz());
			  glVertex3d(v3.x()  , v3.y(),   v3.z() );

			  glNormal3d(v2.nx(),  v2.ny(),  v2.nz());
			  glVertex3d(v2.x()  , v2.y(),   v2.z() );
		  glEnd();
    }
  }
}
//----------------------------------------------------------------//
void light_disable()
//----------------------------------------------------------------//
{
    glDisable(GL_POLYGON_OFFSET_FILL);
	glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_NORMALIZE);
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
}
//----------------------------------------------------------------//
void light_enable()
//----------------------------------------------------------------//
{
  GLfloat light0_ambient [] = {  0.2, 0.2, 0.2, 1.0 };
  GLfloat light0_diffuse    [] = {  0.2, 0.2, 0.2, 1.0 };
  GLfloat light0_specular [] = {  1.0, 1.0, 1.0, 1.0 };
  GLfloat light0_position  [] = { -1.0, 1.0, 0.0, 0.0 };

  GLfloat light1_ambient [] = { 0.2, 0.2, 0.2, 1.0 };
  GLfloat light1_diffuse    [] = { 0.2, 0.2, 0.2, 1.0 };
  GLfloat light1_specular [] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat light1_position  [] = { 1.0, 1.0, 0.0, 0.0 };

  GLfloat material_specular []   = { 0.2, 0.2, 0.2, 1.0 };
  GLfloat material_shininess[]  = { 46.8 } ;

  glLoadIdentity();
  glMultMatrixf( light0_rotation );
  glLightfv(GL_LIGHT0, GL_AMBIENT    , light0_ambient );
  glLightfv(GL_LIGHT0, GL_DIFFUSE     , light0_diffuse );
  glLightfv(GL_LIGHT0, GL_SPECULAR , light0_specular);
  glLightfv(GL_LIGHT0, GL_POSITION   , light0_position);

  glLoadIdentity();
  glMultMatrixf( light1_rotation );
  glLightfv(GL_LIGHT1, GL_AMBIENT    , light1_ambient );
  glLightfv(GL_LIGHT1, GL_DIFFUSE     , light1_diffuse );
  glLightfv(GL_LIGHT1, GL_SPECULAR , light1_specular);
  glLightfv(GL_LIGHT1, GL_POSITION   , light1_position);

  glMaterialfv(GL_FRONT, GL_SPECULAR , material_specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, material_shininess );

  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_NORMALIZE);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);

  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(0.5,0.5);
}
//----------------------------------------------------------------//
void start_config()
//----------------------------------------------------------------//
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1, 1, -1, 1, -1, 1);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef (obj_pos[0], obj_pos[1], obj_pos[2] );
	glMultMatrixf( view_rotate );
	glScalef(scale,scale,scale);

	glEnable (GL_DEPTH_TEST);
	glDepthFunc (GL_LESS);
	glClearDepth(10.0);
	glClear(GL_DEPTH_BUFFER_BIT);
}
//----------------------------------------------------------------//
void display()
//----------------------------------------------------------------//
{
	glClearColor(1.0,1.0,1.0,1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	light_enable();
	start_config();

	 if(estar){
	  	draw_vstar(dim);
	 }

	if(level == 0){
		if(pnt) {
			if(vdraw == 0) ch0.draw_vert(VERT,  vtype);
			if(vdraw == 3) ch0.draw_vert(FIELD, vtype);
		}
		if(edg && edraw == 0)
			ch0.draw_wire(WIRE);
	}
	if(level == 1){
		if(pnt) {
			if(vdraw == 0) ch1.draw_vert(VERT , vtype);
			if(vdraw == 3) ch1.draw_vert(FIELD, vtype);
		}
		if(edg && edraw == 0)
			ch1.draw_wire(WIRE);

		if(srf) {
			if(sdraw == 0) ch1.draw_smooth(FLAT);
			if(sdraw == 1) ch1.draw_smooth(SMOOTH);
		}
	}
	if(level == 2){
		if(pnt){
			if(vdraw == 0) ch2.draw_vert(VERT     , vtype);
			if(vdraw == 1) ch2.draw_vert(IN_VERT, vtype);
			if(vdraw == 2) ch2.draw_vert(B_VERT , vtype);
			if(vdraw == 3) ch2.draw_vert(FIELD    , vtype);
		}
		if(edg){
			if(edraw == 0) ch2.draw_wire(WIRE);
			if(edraw == 1) ch2.draw_wire(IN_WIRE);
			if(edraw == 2) ch2.draw_wire(B_WIRE);
		}
		if(srf){
			if(sdraw == 0) ch2.draw_smooth(FLAT);
			if(sdraw == 1) ch2.draw_smooth(SMOOTH);
		}
	}
	if(level == 3)
	{
		if(pnt){
			if(vdraw == 0) ch3.draw_vert(VERT     , vtype);
			if(vdraw == 1) ch3.draw_vert(IN_VERT, vtype);
			if(vdraw == 2) ch3.draw_vert(B_VERT , vtype);
			if(vdraw == 3) ch3.draw_vert(FIELD    , vtype);
		}
		if(edg){
			if(edraw == 0) ch3.draw_wire(WIRE);
			if(edraw == 1) ch3.draw_wire(IN_WIRE);
			if(edraw == 2) ch3.draw_wire(B_WIRE);
		}
		if(srf){
			if(sdraw == 0) ch3.draw_smooth(FLAT);
			if(sdraw == 1) ch3.draw_smooth(SMOOTH);
			if(sdraw == 2) ch3.draw_smooth(BOUND);
		}
	}

	light_disable();
	glutSwapBuffers();
}
//---------------------------------------------------------------//
void myGlutIdle( void )
//---------------------------------------------------------------//
{
  if ( glutGetWindow() != window_id)
  glutSetWindow(window_id);
  glutPostRedisplay();
}
//----------------------------------------------------------------//
void Reshape  (int w, int h)
//----------------------------------------------------------------//
{
	int tx, ty, tw, th;
	GLUI_Master.get_viewport_area( &tx, &ty, &tw, &th );
	glViewport( tx, ty, tw, th );
	mouse_rt.set_w( tw ) ;
	mouse_rt.set_h( th ) ;
	mouse_xy .set_w( tw ) ;
	mouse_xy .set_h( th ) ;
	mouse_zm .set_w( tw ) ;
	mouse_zm .set_h( th ) ;

	glutPostRedisplay();
}
//----------------------------------------------------------------//
void mouse(int button, int button_state, int x, int y )
//----------------------------------------------------------------//
{
	if         ( glutGetModifiers() & GLUT_ACTIVE_CTRL  ) motion_type = 1 ;
	else if( glutGetModifiers() & GLUT_ACTIVE_SHIFT ) motion_type = 2 ;
	else motion_type = 0 ;

	switch( motion_type ){
		case 0 :
			if ( button == GLUT_LEFT_BUTTON && button_state == GLUT_DOWN ){
				mouse_rt.init_live() ;
				mouse_rt.mouse_down_handler(x,y) ;
			}
			if ( button_state != GLUT_DOWN ){
				mouse_rt.mouse_up_handler(x,y,1) ;
			}
			object_rt->sync_live(0,1) ;
		break ;

		case 1 :
			if ( button == GLUT_LEFT_BUTTON && button_state == GLUT_DOWN ){
				mouse_zm.init_live() ;
				mouse_zm.glui = glui ;
				mouse_zm.mouse_down_handler(x,y) ;
				mouse_zm.glui = NULL ;
			}
			if ( button_state != GLUT_DOWN ){
				mouse_zm.glui = glui ;
				mouse_zm.mouse_up_handler(x,y,1) ;
				mouse_zm.glui = NULL ;
			}
			object_zm->sync_live(0,1) ;
		break ;

		case 2 :
			if ( button == GLUT_LEFT_BUTTON && button_state == GLUT_DOWN ){
				mouse_xy.init_live() ;
				mouse_xy.glui = glui ;
				mouse_xy.mouse_down_handler(x,y) ;
				mouse_xy.glui = NULL ;
			}
			if ( button_state != GLUT_DOWN ){
				mouse_xy.glui = glui ;
				mouse_xy.mouse_up_handler(x,y,1) ;
				mouse_xy.glui = NULL ;
			}
			object_xy->sync_live(0,1) ;
		break ;
	}
	ncalls = 0 ;
}
//----------------------------------------------------------------//
void motion(int x, int y )
//----------------------------------------------------------------//
{
	switch( motion_type ){
		case 0 :
			mouse_rt.glui = glui ;
			mouse_rt.iaction_mouse_held_down_handler(x,y,1);
			mouse_rt.glui = NULL ;
			if( ++ncalls > 10 ) { object_rt->sync_live(0,1) ;  ncalls = 0 ; }
		break ;

		case 1 :
			mouse_zm.glui = glui ;
			mouse_zm.iaction_mouse_held_down_handler(x,y,1);
			mouse_zm.glui = NULL ;
			if( ++ncalls > 10 ) { object_zm->sync_live(0,1) ;  ncalls = 0 ; }
		break ;

		case 2 :
			mouse_xy.glui = glui ;
			mouse_xy.iaction_mouse_held_down_handler(x,y,1);
			mouse_xy.glui = NULL ;
			if( ++ncalls > 10 ) { object_xy->sync_live(0,1) ;  ncalls = 0 ; }
		break ;
	}
	glutPostRedisplay() ;
}
//----------------------------------------------------------------//
void read_key(int key)
//----------------------------------------------------------------//
{
	int        size =0;
	int      state = GL2PS_OVERFLOW;
	GLint buffer = 0;

	switch(key){
		case 1:
			size = (int)strlen( file_ply );
			if( size != 0 ){
				if( file_ply[size-4] != '.' )
					strcat( file_ply, ".ply" ) ;
				if(level == 0)ch0.write_ply( file_ply );
				if(level == 1)ch1.write_ply( file_ply );
				if(level == 2)ch2.write_ply( file_ply );
				if(level == 3)ch3.write_ply( file_ply );
			}
		break;

		case 2:
			if( level == 0 ){
				ch0.read_ply(file_list[fileid]);
				ch0.scalar_field("x*x+y*y+z*z-0.1");
				//ch0.check();
				nedges = 6*ch0.ntetra();
			}
			if( level == 1 ){
				ch1.read_ply(file_list[fileid]);
				ch1.scalar_field("x*x+y*y+z*z-0.1");
				//ch1.check();
				nedges = 6*ch1.ntetra();
			}
			if( level == 2 ){
				ch2.read_ply(file_list[fileid]);
				ch2.scalar_field("x*x+y*y+z*z-0.1");
				//ch2.check();
				nedges = 6*ch2.ntetra();
			}
			if( level == 3 ){
				ch3.read_ply(file_list[fileid]);
				ch3.scalar_field("x*x+y*y+z*z-0.1");
				//ch3.check();
				nedges = 6*ch3.ntetra();
			}
		    simplexid->set_int_limits(0, nedges);
			read_key(6);
		break;

		case 3:
			size = (int)strlen( file_eps );
			if ( size != 0 ){
				if( file_eps[size-4] != '.' ) strcat( file_eps, ".eps" );
				FILE* fp = fopen( file_eps, "wb" );
				if( !fp ) break;

				GLint viewport[4];
				glutSetWindow(window_id);
				GLUI_Master.get_viewport_area( viewport + 0, viewport + 1, viewport + 2, viewport + 3 );
				buffer = viewport[2]*viewport[3];

				while( state == GL2PS_OVERFLOW) {
					buffer += buffer;
					gl2psBeginPage ( file_eps, "CHF: Example (by Marcos Lage)", viewport,
												GL2PS_EPS, GL2PS_BSP_SORT,
												GL2PS_SIMPLE_LINE_OFFSET | GL2PS_NO_PS3_SHADING | GL2PS_OCCLUSION_CULL,
												GL_RGBA, 0, NULL, 0, 0, 0, buffer,
												fp, NULL );
					display();
					state = gl2psEndPage();
				}
				fclose(fp);
			}
		break;

		case 4:
			if( srf )
				sgroup->enable();
			else
				sgroup->disable();

			if( edg )
				egroup->enable();
			else
				egroup->disable();

			if( pnt )
				vgroup->enable();
			else
				vgroup->disable();

			if (level == 0 ){
				sgroup->disable();
				eb2->disable();
				eb3->disable();
				vb2->disable();
				vb3->disable();
			}
			if (level == 1 ){
				sb3->disable();
				eb2->disable();
				eb3->disable();
				vb2->disable();
				vb3->disable();
			}
			if (level == 2 )
				sb3->disable();
		break;

		case 5:
			display();
		break;

		case 6:
		star.clear();
		star = test_vstar(eid, dim);
		break;

		case 9:
			object_rt->reset();
			object_xy->set_x(obj_pos[0]);
			object_xy->set_y(obj_pos[1]);
			object_zm->set_z(1.0);
		break;

		case 10:
			exit(1);
		break;
	}
}
//----------------------------------------------------------------//
int main(int argc, char **argv)
//----------------------------------------------------------------//
{
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(view_w,view_h);
	glutInitWindowPosition(250,150);
	window_id= glutCreateWindow("CHF_GLUI.cpp");

	glutDisplayFunc(display);
	glutMotionFunc ( motion );
	glutMouseFunc  ( mouse  );

	GLUI_Master.set_glutReshapeFunc ( Reshape );
	GLUI_Master.set_glutIdleFunc( myGlutIdle );

	glui  = GLUI_Master.create_glui_subwindow( window_id, GLUI_SUBWINDOW_BOTTOM );
	glui_r= GLUI_Master.create_glui_subwindow( window_id ,GLUI_SUBWINDOW_RIGHT );

	glui->set_main_gfx_window( window_id );
	glui_r->set_main_gfx_window( window_id );

	//----------------------------------------------------------//
	/*Trackball initialization */
	//----------------------------------------------------------//
	int tx,ty,tw,th ;
	GLUI_Master.get_viewport_area( &tx, &ty, &tw, &th );

	mouse_rt.set_w( tw ) ;
	mouse_rt.set_h( th ) ;
	mouse_rt.set_ptr_val( view_rotate );
	mouse_rt.init_live() ;
	mouse_rt.hidden = true;

	mouse_xy.set_speed(0.01f) ;
	mouse_xy.set_w( tw ) ;
	mouse_xy.set_h( th ) ;
	mouse_xy.set_ptr_val( obj_pos );
	mouse_xy.init_live() ;
	mouse_xy.trans_type = GLUI_TRANSLATION_XY;
	mouse_xy.float_array_size = 2 ;
	mouse_xy.hidden = true ;

	mouse_zm.set_speed(0.01f) ;
	mouse_zm.set_w( tw ) ;
	mouse_zm.set_h( th ) ;
	mouse_zm.set_ptr_val( &scale );
	mouse_zm.init_live() ;
	mouse_zm.trans_type = GLUI_TRANSLATION_Z ;
	mouse_zm.float_array_size = 1 ;
	mouse_zm.hidden = true ;

	//----------------------------------------------------------//
	/*Main Window	*/
	//----------------------------------------------------------//

	glui->add_checkbox("Surface", &srf, 4, read_key);
	sgroup = glui->add_radiogroup ( &sdraw );
	sb1  = glui->add_radiobutton_to_group( sgroup, "Flat" );
	sb2  = glui->add_radiobutton_to_group( sgroup, "Smooth" );
	sb3  = glui->add_radiobutton_to_group( sgroup, "Bound Surfaces" );

	glui->add_column(false);

	glui->add_checkbox("Edges", &edg, 4, read_key);
	egroup = glui->add_radiogroup ( &edraw );
	eb1  = glui->add_radiobutton_to_group( egroup, "Wireframe" );
	eb2  = glui->add_radiobutton_to_group( egroup, "Edge Classification" );
	eb3  = glui->add_radiobutton_to_group( egroup, "Bound Wireframe" );

	glui->add_column(false);

	glui->add_checkbox("Points", &pnt, 4, read_key);
	vgroup = 	glui->add_radiogroup ( &vdraw );
	vb1  = glui->add_radiobutton_to_group( vgroup, "Vertices" );
	vb2  = glui->add_radiobutton_to_group( vgroup, "Vertex Classification" );
	vb3  = glui->add_radiobutton_to_group( vgroup, "Bound Vertices" );
	vb4  = glui->add_radiobutton_to_group( vgroup, "Scalar Field" );

	glui->add_column(false);

	//glui->add_checkbox("glBlending", &vtype, 5, read_key);
	//glui->add_checkbox("glLines", &vtype, 5, read_key);
	glui->add_checkbox("glVertex", &vtype, 5, read_key);

	glui->add_column(false);

	object_rt = glui->add_rotation( "Rotation", view_rotate );

	glui->add_column(false);

	object_zm = glui->add_translation("Zoom", GLUI_TRANSLATION_Z, &scale );
	object_zm->set_speed( .01f );

	glui->add_column(false);

	object_xy = glui->add_translation("Translation", GLUI_TRANSLATION_XY, obj_pos );
	object_xy->set_speed( .01f );

	//----------------------------------------------------------//
	/*Second Window	*/
	//----------------------------------------------------------//

	GLUI_Panel *chflevel =
	glui_r->add_panel( "CHF Level", GLUI_PANEL_EMBOSSED );

	GLUI_RadioGroup* lgroup =
	glui_r->add_radiogroup_to_panel ( chflevel, &level, 4 , read_key );
	glui_r->add_radiobutton_to_group( lgroup, "CHF_L0" );
	glui_r->add_radiobutton_to_group( lgroup, "CHF_L1" );
	glui_r->add_radiobutton_to_group( lgroup, "CHF_L2" );
	glui_r->add_radiobutton_to_group( lgroup, "CHF_L3" );
	lgroup->set_alignment( GLUI_ALIGN_CENTER);

	GLUI_Listbox *listfile=
	glui_r->add_listbox("Select File" ,&fileid);

	for( int i=0; i<38; ++i )
		listfile->add_item( i, file_list[i] );

	listfile->set_alignment( GLUI_ALIGN_LEFT);
	listfile->set_w( 150 );

	glui_r->add_button( "Read PLY" , 2, read_key );

	GLUI_EditText *write_ply=
	glui_r->add_edittext( "File Name", GLUI_EDITTEXT_TEXT, file_ply );
	write_ply->set_alignment( GLUI_ALIGN_LEFT);
	write_ply->set_w( 150 );

	glui_r->add_button( "Write PLY" , 1, read_key );

	GLUI_EditText *write_eps=
	glui_r->add_edittext( "File Name", GLUI_EDITTEXT_TEXT, file_eps );
	write_eps->set_alignment( GLUI_ALIGN_LEFT);
	write_eps->set_w( 150 );

	glui_r->add_button( "Write EPS" , 3, read_key );

	//----------------------------------------------------------//
	/*Stars Panel*/
	//----------------------------------------------------------//
	GLUI_Panel *panel4 = glui_r->add_panel("Stars", GLUI_PANEL_EMBOSSED );
	panel4->set_alignment(GLUI_ALIGN_LEFT);

	glui_r->add_checkbox_to_panel(panel4, "eStars", &estar);

	GLUI_RadioGroup *group2 =
	glui_r->add_radiogroup_to_panel(panel4, &dim, 6, read_key);
		glui_r->add_radiobutton_to_group( group2, "R10" );
		glui_r->add_radiobutton_to_group( group2, "R13" );

	simplexid = glui_r->add_spinner_to_panel(panel4, "Edge id", GLUI_SPINNER_INT, &eid, 6, read_key );
	simplexid->set_int_limits(0, nedges, GLUI_LIMIT_WRAP);
	simplexid->set_speed(0.0001);
	//----------------------------------------------------------//

	light_rt = glui_r->add_rotation( "Light 0", light0_rotation  );
	light_rt = glui_r->add_rotation( "Light 1", light1_rotation  );

	glui_r->add_button( "Clear Position" , 9, read_key );
	glui_r->add_button( "Quit" , 10, read_key );

	read_key(4);
	glutMainLoop();
}
//--------------------------------------------------------------//
