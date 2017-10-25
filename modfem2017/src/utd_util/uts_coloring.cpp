#ifndef _uts_coloring_
#define _uts_coloring_

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/shared_array_property_map.hpp>
#include <boost/graph/smallest_last_ordering.hpp>
#include <boost/graph/sequential_vertex_coloring.hpp>
#include "mmh_intf.h"
#include "uth_log.h"


int utr_generate_mesh_coloring( // return number of colors
        const int Mesh_id, // IN: mesh id to color
        int *Elem_colors)  // OUT: array of size number of active elems containing colors (numbers)
{

    mf_check_mem(Elem_colors);

    // Constructing and filling the graph.
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
    const int nElems = mmr_get_nr_elem(Mesh_id);
    Graph g(nElems);

    int nel=0;
    while( (nel=mmr_get_next_act_elem(Mesh_id,nel)) != 0) {
        int neigs[MMC_MAXELFAC+1]={0};
        mmr_el_eq_neig(Mesh_id,nel,neigs,NULL);
        for(int i=1; i < neigs[0]; ++i) {
            add_edge(nel,neigs[i],g);
        }
    }

    std::fill(Elem_colors,Elem_colors+nElems,0);
    boost::smallest_last_vertex_ordering(g,Elem_colors);

    std::vector<int> color_counts;
    int cur_color=-1;
    do {
        ++cur_color;
        color_counts.push_back( std::count(Elem_colors,Elem_colors+nElems,cur_color) );
    }
    while(color_counts[cur_color] > 0);

    mf_check_mem(Elem_colors);

    return cur_color;
}

#ifdef __cplusplus
extern "C"
{
#endif

struct element{
  int id;
  element* next;
};


int utr_first_fit_coloring_with_front( // return number of colors
        const int nElems, // IN: number of elements to color
        int ** elem2elem, // IN: array elements neighborhood
        int *Elem_colors)  // OUT: array of size number of active elems containing colors (numbers)
{
  int MAX_COLOR=500;
  int i,color_id,flag;
  int current_elem=0,nr_color=0,colored_elems=0,nr_colored;
  element* head=NULL;
  element* tail=NULL;
  bool* avail_elem = (bool *)malloc(nElems*sizeof(bool));
  
  for(i=0;i<nElems;++i){
    avail_elem[i]=true;
    Elem_colors[i]=-1;
  }
  
  avail_elem[current_elem]=false;
  Elem_colors[current_elem]=0;
  nr_colored=1;
  nr_color++;
  bool* avail_color = (bool *)malloc(MAX_COLOR*sizeof(bool));
  for(i=0;i<MAX_COLOR;++i){
    avail_color[i]=true;
  }
  
  while(nr_colored<nElems){
    for(i=1;i<=elem2elem[current_elem][0];++i){ //add neighbors to check list
      if(avail_elem[elem2elem[current_elem][i]]==true){
        element* nowy=new element;
        nowy->id=elem2elem[current_elem][i];
        nowy->next=NULL;
        avail_elem[elem2elem[current_elem][i]]=false;
        if(head==NULL){head=nowy;tail=nowy;}
        else {tail->next=nowy;tail=nowy;}
      }
    }
    if(head!=NULL){ //is neighbor
      current_elem=head->id;
      element* nowy=head;
      head=head->next;
      if(head==NULL)tail=NULL;
      delete nowy;
      avail_elem[current_elem]=false;
      for(i=1;i<=elem2elem[current_elem][0];++i){ //checking neighbors colors
        color_id=Elem_colors[elem2elem[current_elem][i]];
        if(color_id != -1){
          avail_color[color_id]=false;
        }
      }
      flag=0;
      for(color_id=0;color_id<nr_color;++color_id){ //set color with minimal id
        if(avail_color[color_id]==true && flag==0){
          Elem_colors[current_elem]=color_id;
          nr_colored++;
          flag=1;break;
        }
        avail_color[color_id]==true;
      }
      if(flag==0){  //set color with new id
        Elem_colors[current_elem]=nr_color;
        nr_color++;
        nr_colored++;
        if(nr_color>MAX_COLOR){
          avail_color = (bool*) realloc (avail_color, nr_color * sizeof(bool));
          for(color_id=0;color_id<nr_color;++color_id)avail_color[color_id]=true;
        }
      }else{
        for(color_id=0;color_id<nr_color;++color_id)avail_color[color_id]=true;
      };
    }
    else{ //isnt neighbor
      for(current_elem=0;current_elem<nElems && Elem_colors[current_elem] != -1;++current_elem);
      if(Elem_colors[current_elem] == -1){
        avail_elem[current_elem]=false;
        for(i=1;i<=elem2elem[current_elem][0];++i){
          color_id=Elem_colors[elem2elem[current_elem][i]];
          if(color_id != -1){
            avail_color[color_id]=false;
          }
        }
        flag=0;
        for(color_id=0;color_id<nr_color;++color_id){
          if(avail_color[color_id]==true && flag==0){
            Elem_colors[current_elem]=color_id;
            flag=1;
            nr_colored++;break;
          }
          avail_color[color_id]==true;
        }
        if(flag==0){
          Elem_colors[current_elem]=nr_color;
          nr_color++;
          nr_colored++;
          if(nr_color>MAX_COLOR){
            avail_color = (bool*) realloc (avail_color, nr_color * sizeof(bool));
            MAX_COLOR=nr_color;
            for(color_id=0;color_id<nr_color;++color_id)avail_color[color_id]=true;
          }
        }else{
          for(color_id=0;color_id<nr_color;++color_id)avail_color[color_id]=true;
        }
      }
    }
  }//while
  
  free(avail_elem);
  free(avail_color);
  
  return nr_color;
}

int utr_first_fit_coloring_without_front( // return number of colors
        const int nElems, // IN: number of elements to color
        int ** elem2elem, // IN: array elements neighborhood
        int *Elem_colors)  // OUT: array of size number of active elems containing colors (numbers)
{
  int MAX_COLOR=500;
  int i,color_id,flag;
  int current_elem=0,nr_color=0,nr_colored;
  bool* avail_color = (bool *)malloc(MAX_COLOR*sizeof(bool));
  
  
  for(current_elem=0;current_elem<nElems;++current_elem){
    Elem_colors[current_elem]=-1;
  }
  
  Elem_colors[0]=0;
  nr_colored=1;
  nr_color++;


  for(current_elem=0;current_elem<nElems;++current_elem){
    for(color_id=0;color_id<MAX_COLOR;++color_id){
      avail_color[color_id]=true;
    }
    for(i=1;i<=elem2elem[current_elem][0];++i){ //get neighbors color
        color_id=Elem_colors[elem2elem[current_elem][i]];
        if(color_id != -1){
          avail_color[color_id]=false;
        }
    }
    for(color_id=0,flag=0;color_id<nr_color;++color_id){ //set color with minimal id
        if(avail_color[color_id]==true && flag==0){
          Elem_colors[current_elem]=color_id;
          nr_colored++;
          flag=1;break;
        }
    }
    if(flag==0){  //set color with new id
      Elem_colors[current_elem]=nr_color;
      nr_color++;
      nr_colored++;
      if(nr_color>MAX_COLOR){
        avail_color = (bool*) realloc (avail_color, nr_color * sizeof(bool));
        MAX_COLOR=nr_color;
      }
    }
  }
  free(avail_color);
  
  return nr_color;
}


int utr_boost_sequential_coloring( // return number of colors
        const int nElems, // IN: number of elements to color
        int ** elem2elem, // IN: array elements neighborhood
        int *Elem_colors)  // OUT: array of size number of active elems containing colors (numbers)
{
  
  int i,j;
  int nr_color=0;
  mf_check_mem(Elem_colors);

  // Constructing and filling the graph.
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
  //const int nElems = mmr_get_nr_elem(Mesh_id);
  Graph g(nElems);

  
  for(i=0;i<nElems;i++){
    for(j=1;j<=elem2elem[i][0];j++)
      add_edge(i,elem2elem[i][j],g);
  }
  
  
  nr_color = boost::sequential_vertex_coloring(g, Elem_colors);
  
  return nr_color;
}


int utr_generate_int_ent_coloring( // return number of colors
        const int nElems, // IN: number of elements to color
        int ** elem2elem, // IN: array elements neighborhood
        int *Elem_colors)  // OUT: array of size number of active elems containing colors (numbers)
{
  //int nr_colors=utr_boost_sequential_coloring(nElems,elem2elem,Elem_colors);
  //int nr_colors=utr_first_fit_coloring_without_front(nElems,elem2elem,Elem_colors);
  int nr_colors=utr_first_fit_coloring_with_front(nElems,elem2elem,Elem_colors);

  
  /*kbw
  printf("nr_colors=%d\n",nr_colors);
  for(i=0;i<nElems;i++){
    printf("%d ",Elem_colors[i]);
  }
  /*kew*/
  
  
  return nr_colors;
}

#ifdef __cplusplus
}
#endif

#endif // uts_coloring_
