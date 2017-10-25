#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

int main(int argc, char **argv)
{
  
  const int maxvert=10000; 
  const int maxedge=40000; 
  const int maxface=50000; 
  const int maxelem=20000;

  // 3D vertex coordinates
  float v_x_3d[maxvert];
  float v_y_3d[maxvert];
  float v_z_3d[maxvert];
  // for 3D nas meshes
  int renumber[maxvert];

  // 2D vertex coordinates
  float v_x[maxvert];
  float v_y[maxvert];

  // triangle vertices
  int t_v1[maxface];
  int t_v2[maxface];
  int t_v3[maxface];

  // triangle neighbours or BC flags
  int t_n1[maxface];
  int t_n2[maxface];
  int t_n3[maxface];

  // edge vertices
  int e_v1[maxedge];
  int e_v2[maxedge];

  // edge neighbours and BC flags
  int e_bc[maxedge];
  int e_n1[maxedge];
  int e_n2[maxedge];
  

  int lower_BC_flag=0;
  int upper_BC_flag=0;


  char nas_filename[256];
  char jk_filename[256]; 
  char keyword[256];
  
  sprintf(nas_filename,"%s.nas",argv[1]);
  sprintf(jk_filename,"%s.jk",argv[1]);
  
  FILE* nas_file = fopen(nas_filename,"r");

  if(nas_file == NULL){
    printf("\nCannot find input file %s \n", nas_filename);
    exit(-1); 
  }


  fscanf(nas_file,"%8s",keyword);

  while(strncmp(keyword,"GRID", 4) != 0) {
    fscanf(nas_file,"%8s\n",keyword);
    printf("%s\n", keyword);
  }

  printf("Found GRID keyword\n");

  int nr_vert_3d=0;

  do{
    
    int data1;
    fscanf(nas_file,"%8d",&data1);
    int data2;
    fscanf(nas_file,"%8d",&data2);
    float data3;
    fscanf(nas_file,"%8g",&data3);
    float data4;
    fscanf(nas_file,"%8g",&data4);
    float data5;
    fscanf(nas_file,"%8g",&data5);

    printf("%8s%8d%8d%8G%8G%8G\n", keyword, data1, data2, data3, data4, data5);


    nr_vert_3d++;
    v_x_3d[nr_vert_3d] = data3;
    v_y_3d[nr_vert_3d] = data4;
    v_z_3d[nr_vert_3d] = data5;
    
    fscanf(nas_file,"%8s\n",keyword);
    
  } while(strncmp(keyword,"CTRIA3", 6) != 0 && 
	  strncmp(keyword,"CQUAD4", 6) != 0);


  printf("Found CTRIA3 keyword\n");


  int nr_tria = 0;
  int nr_edge = 0;

  do{
    
    if(strncmp(keyword,"CTRIA3", 6) == 0){

      int data1;
      fscanf(nas_file,"%8d",&data1);
      int data2;
      fscanf(nas_file,"%8d",&data2);
      int data3;
      fscanf(nas_file,"%8d",&data3);
      int data4;
      fscanf(nas_file,"%8d",&data4);
      int data5;
      fscanf(nas_file,"%8d",&data5);

      printf("TRIA: %8s%8d%8d%8d%8d%8d\n", 
	     keyword, data1, data2, data3, data4, data5);

      // get one of vertices
      if(fabs(v_z_3d[data3]) < 1.e-9){

	// triangle in 2d plane
	t_v1[nr_tria] = data3;
	t_v2[nr_tria] = data4;
	t_v3[nr_tria] = data5;

	float x1_temp = v_x_3d[t_v2[nr_tria]] - v_x_3d[t_v1[nr_tria]];
	float y1_temp = v_y_3d[t_v2[nr_tria]] - v_y_3d[t_v1[nr_tria]];
	float x2_temp = v_x_3d[t_v3[nr_tria]] - v_x_3d[t_v1[nr_tria]];
	float y2_temp = v_y_3d[t_v3[nr_tria]] - v_y_3d[t_v1[nr_tria]];
	printf("vector product for triangle %d: %f\n", nr_tria,
	       x1_temp*y2_temp-y1_temp*x2_temp);
	if(x1_temp*y2_temp-y1_temp*x2_temp < 0.0){
	  printf("Clockwise orientation in JK file !!!\n");
	}


	if(lower_BC_flag==0){
	  lower_BC_flag = data2;
	}
	else{
	  if(data2!=lower_BC_flag){

	    printf("different BC flags for lower 2D plane!\n");
	    exit(-1);

	  }

	}

	nr_tria++;

      }
      else{

	if(upper_BC_flag==0){
	  upper_BC_flag = data2;
	}
	else{
	  if(data2!=upper_BC_flag){

	    printf("different BC flags for upper 2D plane!\n");
	    exit(-1);

	  }

	}

      }

    }
    else if(strncmp(keyword,"CQUAD4", 6) == 0){

      int data1;
      fscanf(nas_file,"%8d",&data1);
      int data2;
      fscanf(nas_file,"%8d",&data2);
      int data3;
      fscanf(nas_file,"%8d",&data3);
      int data4;
      fscanf(nas_file,"%8d",&data4);
      int data5;
      fscanf(nas_file,"%8d",&data5);
      int data6;
      fscanf(nas_file,"%8d",&data6);

      printf("QUAD: %8s%8d%8d%8d%8d%8d%8d\n", 
	     keyword, data1, data2, data3, data4, data5, data5);

      // quadrilaterals in 2D are represented by edges

      int e_vert = 0;
      if(fabs(v_z_3d[data3]) < 1.e-9){
	if(e_vert==0) e_v1[nr_edge] = data3;
	else  e_v2[nr_edge] = data3;
	e_vert++;
      }
      if(fabs(v_z_3d[data4]) < 1.e-9){
	if(e_vert==0)  e_v1[nr_edge] = data4;
	else e_v2[nr_edge] = data4;
	e_vert++;
      }
      if(fabs(v_z_3d[data5]) < 1.e-9){
	if(e_vert==0) e_v1[nr_edge] = data5;
	else  e_v2[nr_edge] = data5;
	e_vert++;
      }
      if(fabs(v_z_3d[data6]) < 1.e-9){
	if(e_vert==0) e_v1[nr_edge] = data6;
	else  e_v2[nr_edge] = data6;
	e_vert++;
      }

      e_bc[nr_edge] =   data2;
      e_n1[nr_edge] = -111173;
      e_n2[nr_edge] = -111173;

      if(e_vert==2){
	printf("edge %8d from quad %8d: v1 %8d, v2 %8d, bc %d\n",
	       nr_edge, data1, e_v1[nr_edge], e_v2[nr_edge], e_bc[nr_edge]);
      }
      else{
	printf("cannot create an edge from quad %d\n", data1);
      }

      nr_edge++;

    }

    fscanf(nas_file,"%8s\n",keyword);
    
  } while(strncmp(keyword,"CPENTA", 6) != 0);

  fclose(nas_file);

  // process vertices - squeeze 3D to 2D
  int nr_vert=0;
  int ivert_3d; int ivert; 
  for(ivert_3d=0; ivert_3d<nr_vert_3d; ivert_3d++){

    if(fabs(v_z_3d[ivert_3d]) < 1.e-9){

      v_x[nr_vert] = v_x_3d[ivert_3d];
      v_y[nr_vert] = v_y_3d[ivert_3d];

      renumber[ivert_3d] = nr_vert;


      printf("2D vertex: %4d    (old %4d), x = %12lf, y = %12lf\n",
	     nr_vert, ivert_3d, v_x[nr_vert], v_y[nr_vert] );

      nr_vert++;


    }

  }

  for(ivert_3d=0; ivert_3d<nr_vert_3d; ivert_3d++){

    if(fabs(v_z_3d[ivert_3d]) >= 1.e-9){

      int found = 0;
      for(ivert=0; ivert< nr_vert; ivert++){

	if(fabs(v_x[ivert] - v_x_3d[ivert_3d]) < 1.e-9 && 
	   fabs(v_y[ivert] - v_y_3d[ivert_3d]) < 1.e-9) {

	  renumber[ivert_3d] = ivert;
	  found = 1;
	  break;

	}
      }

      if(found==1){
	printf("3D vertex: %4d    (new %4d), x = %12lf, y = %12lf\n",
	       ivert_3d, ivert, v_x[ivert], v_y[ivert] );
      }
      else{
	printf("3D vertex: %4d,  x = %12lf, y = %12lf - not found in 2D\n",
	       ivert_3d, v_x_3d[ivert_3d], v_y_3d[ivert_3d] );
	exit(-1);

      }

    }
  
  }

  int itria; int iedge;

  // rewrite 3d data to 2d data
  for(itria=0; itria<nr_tria; itria++){

    printf("renumbering for triangle: old (%d, %d, %d) -> new (%d, %d, %d)\n",
	   t_v1[itria], t_v2[itria], t_v3[itria],
	   renumber[t_v1[itria]], renumber[t_v2[itria]], renumber[t_v3[itria]]);

    t_v1[itria] = renumber[t_v1[itria]];
    t_v2[itria] = renumber[t_v2[itria]];
    t_v3[itria] = renumber[t_v3[itria]];



  }

  for(iedge=0; iedge<nr_edge; iedge++){
    e_v1[iedge] = renumber[e_v1[iedge]];
    e_v2[iedge] = renumber[e_v2[iedge]];
  }

  // match edges and triangles to find neighbours
  for(itria=0; itria<nr_tria; itria++){

    t_n1[itria]=-111173;
    t_n2[itria]=-111173;
    t_n3[itria]=-111173;

  }


  for(itria=0; itria<nr_tria; itria++){

    int jtria;
    for(jtria=0; jtria<nr_tria; jtria++){

      if(itria!=jtria){

	if((t_v1[itria]==t_v1[jtria] && t_v2[itria]==t_v2[jtria]) ||
	   (t_v2[itria]==t_v1[jtria] && t_v1[itria]==t_v2[jtria])){
	  
	  if(t_n1[itria] == -111173){
	    t_n1[itria] = jtria;
	  }
	  else if(t_n1[itria] == jtria){
	    // OK - we visit the edge second time
	  }
	  else{
	    printf("too much neighbours for tria %d, edge 1\n", itria);
	    exit(-1);
	  }
	  
	  if(t_n1[jtria] == -111173){
	    t_n1[jtria] = itria;
	  }
	  else if(t_n1[jtria] == itria){
	    // OK - we visit the edge second time
	  }
	  else{
	    printf("too much neighbours for tria %d, edge 1\n", jtria);
	    exit(-1);
	  }
	  
	}

	if((t_v1[itria]==t_v2[jtria] && t_v2[itria]==t_v3[jtria]) ||
	   (t_v2[itria]==t_v2[jtria] && t_v1[itria]==t_v3[jtria])){
	  
	  if(t_n1[itria] == -111173){
	    t_n1[itria] = jtria;
	  }
	  else if(t_n1[itria] == jtria){
	    // OK - we visit the edge second time
	  }
	  else{
	    printf("too much neighbours for tria %d, edge 1\n", itria);
	    exit(-1);
	  }
	  
	  if(t_n2[jtria] == -111173){
	    t_n2[jtria] = itria;
	  }
	  else if(t_n2[jtria] == itria){
	    // OK - we visit the edge second time
	  }
	  else{
	    printf("too much neighbours for tria %d, edge 1\n", jtria);
	    exit(-1);
	  }
	  
	}

	if((t_v1[itria]==t_v1[jtria] && t_v2[itria]==t_v3[jtria]) ||
	   (t_v2[itria]==t_v1[jtria] && t_v1[itria]==t_v3[jtria])){
	  
	  if(t_n1[itria] == -111173){
	    t_n1[itria] = jtria;
	  }
	  else if(t_n1[itria] == jtria){
	    // OK - we visit the edge second time
	  }
	  else{
	    printf("too much neighbours for tria %d, edge 1\n", itria);
	    exit(-1);
	  }
	  
	  if(t_n3[jtria] == -111173){
	    t_n3[jtria] = itria;
	  }
	  else if(t_n3[jtria] == itria){
	    // OK - we visit the edge second time
	  }
	  else{
	    printf("too much neighbours for tria %d, edge 1\n", jtria);
	    exit(-1);
	  }
	  
	}

	if((t_v3[itria]==t_v1[jtria] && t_v2[itria]==t_v2[jtria]) ||
	   (t_v2[itria]==t_v1[jtria] && t_v3[itria]==t_v2[jtria])){
	  
	  if(t_n2[itria] == -111173){
	    t_n2[itria] = jtria;
	  }
	  else if(t_n2[itria] == jtria){
	    // OK - we visit the edge second time
	  }
	  else{
	    printf("too much neighbours for tria %d, edge 1\n", itria);
	    exit(-1);
	  }
	  
	  if(t_n1[jtria] == -111173){
	    t_n1[jtria] = itria;
	  }
	  else if(t_n1[jtria] == itria){
	    // OK - we visit the edge second time
	  }
	  else{
	    printf("too much neighbours for tria %d, edge 1\n", jtria);
	    exit(-1);
	  }
	  
	}

	if((t_v3[itria]==t_v2[jtria] && t_v2[itria]==t_v3[jtria]) ||
	   (t_v2[itria]==t_v2[jtria] && t_v3[itria]==t_v3[jtria])){
	  
	  if(t_n2[itria] == -111173){
	    t_n2[itria] = jtria;
	  }
	  else if(t_n2[itria] == jtria){
	    // OK - we visit the edge second time
	  }
	  else{
	    printf("too much neighbours for tria %d, edge 1\n", itria);
	    exit(-1);
	  }
	  
	  if(t_n2[jtria] == -111173){
	    t_n2[jtria] = itria;
	  }
	  else if(t_n2[jtria] == itria){
	    // OK - we visit the edge second time
	  }
	  else{
	    printf("too much neighbours for tria %d, edge 1\n", jtria);
	    exit(-1);
	  }
	  
	}

	if((t_v3[itria]==t_v1[jtria] && t_v2[itria]==t_v3[jtria]) ||
	   (t_v2[itria]==t_v1[jtria] && t_v3[itria]==t_v3[jtria])){
	  
	  if(t_n2[itria] == -111173){
	    t_n2[itria] = jtria;
	  }
	  else if(t_n2[itria] == jtria){
	    // OK - we visit the edge second time
	  }
	  else{
	    printf("too much neighbours for tria %d, edge 1\n", itria);
	    exit(-1);
	  }
	  
	  if(t_n3[jtria] == -111173){
	    t_n3[jtria] = itria;
	  }
	  else if(t_n3[jtria] == itria){
	    // OK - we visit the edge second time
	  }
	  else{
	    printf("too much neighbours for tria %d, edge 1\n", jtria);
	    exit(-1);
	  }
	  
	}

	if((t_v1[itria]==t_v1[jtria] && t_v3[itria]==t_v2[jtria]) ||
	   (t_v3[itria]==t_v1[jtria] && t_v1[itria]==t_v2[jtria])){
	  
	  if(t_n3[itria] == -111173){
	    t_n3[itria] = jtria;
	  }
	  else if(t_n3[itria] == jtria){
	    // OK - we visit the edge second time
	  }
	  else{
	    printf("too much neighbours for tria %d, edge 1\n", itria);
	    exit(-1);
	  }
	  
	  if(t_n1[jtria] == -111173){
	    t_n1[jtria] = itria;
	  }
	  else if(t_n1[jtria] == itria){
	    // OK - we visit the edge second time
	  }
	  else{
	    printf("too much neighbours for tria %d, edge 1\n", jtria);
	    exit(-1);
	  }
	  
	}

	if((t_v1[itria]==t_v2[jtria] && t_v3[itria]==t_v3[jtria]) ||
	   (t_v3[itria]==t_v2[jtria] && t_v1[itria]==t_v3[jtria])){
	  
	  if(t_n3[itria] == -111173){
	    t_n3[itria] = jtria;
	  }
	  else if(t_n3[itria] == jtria){
	    // OK - we visit the edge second time
	  }
	  else{
	    printf("too much neighbours for tria %d, edge 1\n", itria);
	    exit(-1);
	  }
	  
	  if(t_n2[jtria] == -111173){
	    t_n2[jtria] = itria;
	  }
	  else if(t_n2[jtria] == itria){
	    // OK - we visit the edge second time
	  }
	  else{
	    printf("too much neighbours for tria %d, edge 1\n", jtria);
	    exit(-1);
	  }
	  
	}

	if((t_v1[itria]==t_v1[jtria] && t_v3[itria]==t_v3[jtria]) ||
	   (t_v3[itria]==t_v1[jtria] && t_v1[itria]==t_v3[jtria])){
	  
	  if(t_n3[itria] == -111173){
	    t_n3[itria] = jtria;
	  }
	  else if(t_n3[itria] == jtria){
	    // OK - we visit the edge second time
	  }
	  else{
	    printf("too much neighbours for tria %d, edge 1\n", itria);
	    exit(-1);
	  }
	  
	  if(t_n3[jtria] == -111173){
	    t_n3[jtria] = itria;
	  }
	  else if(t_n3[jtria] == itria){
	    // OK - we visit the edge second time
	  }
	  else{
	    printf("too much neighbours for tria %d, edge 1\n", jtria);
	    exit(-1);
	  }
	  
	}

	
      }

    }

  }


  for(itria=0; itria<nr_tria; itria++){

    for(iedge=0; iedge<nr_edge; iedge++){

      if((t_v1[itria]==e_v1[iedge] && t_v2[itria]==e_v2[iedge]) ||
	 (t_v2[itria]==e_v1[iedge] && t_v1[itria]==e_v2[iedge])){
	if(t_n1[itria] == -111173){
	  t_n1[itria] = -e_bc[iedge];
	}
	else{
	  printf("too much neighbours for tria %d, edge 1\n", itria);
	  exit(-1);
	}
	if(e_n1[iedge] == -111173) e_n1[iedge]=itria;
	else if(e_n2[iedge] == -111173) e_n2[iedge]=itria;
	else {
	  printf("too much neighbours for edge %d\n", iedge);
	  exit(-1);
	}
      }

      if((t_v2[itria]==e_v1[iedge] && t_v3[itria]==e_v2[iedge]) || 
	 (t_v3[itria]==e_v1[iedge] && t_v2[itria]==e_v2[iedge]) ){
	if(t_n2[itria] == -111173){
	  t_n2[itria] = -e_bc[iedge];
	}
	else{
	  printf("too much neighbours for tria %d, edge 2\n", itria);
	  exit(-1);
	}
	if(e_n1[iedge] == -111173) e_n1[iedge]=itria;
	else if(e_n2[iedge] == -111173) e_n2[iedge]=itria;
	else {
	  printf("too much neighbours for edge %d\n", iedge);
	  exit(-1);
	}
      }

      if((t_v1[itria]==e_v1[iedge] && t_v3[itria]==e_v2[iedge]) || 
	 (t_v3[itria]==e_v1[iedge] && t_v1[itria]==e_v2[iedge]) ){
	if(t_n3[itria] == -111173){
	  t_n3[itria] = -e_bc[iedge];
	}
	else{
	  printf("too much neighbours for tria %d, edge 3\n", itria);
	  exit(-1);
	}
	if(e_n1[iedge] == -111173) e_n1[iedge]=itria;
	else if(e_n2[iedge] == -111173) e_n2[iedge]=itria;
	else {
	  printf("too much neighbours for edge %d\n", iedge);
	  exit(-1);
	}
      }

    }


  }
      
  for(itria=0; itria<nr_tria; itria++){

    printf("TRIA %4d: v1 %4d,    v2 %4d,    v3 %4d,     edges: %8d, %8d, %8d\n", 
	   itria, t_v1[itria], t_v2[itria], t_v3[itria],
	   t_n1[itria], t_n2[itria], t_n3[itria] );
    printf("x1 %f, y1 %f, x2 %f, y2 %f, x3 %f, y3 %f\n",
	   v_x[t_v1[itria]], v_y[t_v1[itria]],
	   v_x[t_v2[itria]], v_y[t_v2[itria]],
	   v_x[t_v3[itria]], v_y[t_v3[itria]]);


    float x1_temp = v_x[t_v2[itria]] - v_x[t_v1[itria]];
    float y1_temp = v_y[t_v2[itria]] - v_y[t_v1[itria]];
    float x2_temp = v_x[t_v3[itria]] - v_x[t_v1[itria]];
    float y2_temp = v_y[t_v3[itria]] - v_y[t_v1[itria]];
    printf("vector product for triangle %d: %f\n", itria,
	   x1_temp*y2_temp-y1_temp*x2_temp);
    if(x1_temp*y2_temp-y1_temp*x2_temp < 0.0){
      printf("Clockwise orientation in JK file !!!\n");
    }
    

  }


  for(iedge=0; iedge<nr_edge; iedge++){

	printf("edge %d: v1 %d, v2 %d, bc %d, neig: %d %d \n",
	       iedge, e_v1[iedge], e_v2[iedge], e_bc[iedge], 
	       e_n1[iedge], e_n2[iedge]);

  }


  printf("maxvert %d, maxedge %d, maxface %d, maxelem %d\n",
	 maxvert, maxedge, maxface, maxelem);

  float lower_z = 0.0;
  float upper_z = 0.1;
  int nr_layer = 1;
  printf("nr_vert %d, lower z %lf, upper z %lf, nr_layer %d, lower BC flag %d, upper BC flag %d\n",
	 nr_vert, lower_z, upper_z, nr_layer, lower_BC_flag, upper_BC_flag);

  FILE* jk_file = fopen(jk_filename,"w");

  fprintf(jk_file, "%d %d %d %d\n",
	  maxvert, maxedge, maxface, maxelem);

  fprintf(jk_file, "%d %g %g %d %d %d\n", 
	  nr_vert, lower_z, upper_z, nr_layer, lower_BC_flag, upper_BC_flag);

  for(ivert=0; ivert<nr_vert; ivert++){

    fprintf(jk_file, "%g  %g\n", v_x[ivert], v_y[ivert]);

  }

  fprintf(jk_file, "%d\n", nr_tria);

  for(itria=0; itria<nr_tria; itria++){

    // OFFSET 1 NUMBERING IN JK !!!!!!!!!!!!!!
    t_v1[itria]++;
    t_v2[itria]++;
    t_v3[itria]++;

    if(t_n1[itria]>=0) t_n1[itria]++;
    if(t_n2[itria]>=0) t_n2[itria]++;
    if(t_n3[itria]>=0) t_n3[itria]++;

  }

  for(itria=0; itria<nr_tria; itria++){

    // OFFSET 1 NUMBERING IN JK !!!!!!!!!!!!!!
    fprintf(jk_file, "1 %d %d %d\n", 
	   t_v1[itria], t_v2[itria], t_v3[itria]);
    fprintf(jk_file, "%d %d %d\n", 
	   t_n1[itria], t_n2[itria], t_n3[itria]);

  }

  fclose(jk_file);


  return(0);
}
