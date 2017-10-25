/*
 * fv_config.cpp
 *
 *  Created on: 10 lut 2015
 *      Author: pmaciol
 */

#include "fv_config.h"


size_t tilesize[2];
frameconfig_s * frames_p = NULL;
frameconfig_s * current_frame_p = NULL;
frameconfig_s last_frame_s;
cl_int lightcount;
cl_int2 tiles;
struct light_s * lightlist=NULL;
float density = 6.f;
char * outname=NULL;
char * scenefile=NULL;
int accelstruct;
int interactive=0;


const int nr_frames = 50;
const bool initgl = true;
float dg_density = 0.0f;
int img_width = DFLT_IMAGE_WIDTH;
int img_height = DFLT_IMAGE_HEIGHT;
int gridcount;
int cellcount;
// Warrnig: not normalized
float cut_planes[3][4] = {
		{-0.1f, 0.0f, 0.0f, 0.0f},
		{ 0.0f, 0.0f, 1.0f,-0.75f},
		{ 1.0f, 0.0f, 0.0f,-0.5f}
};

int current_plane_idx = 1;
int test_plane_cut = 0;
float fixed_light_pos[3] = {3.0f,3.0f,3.0f};

float bkg_color[3] = {0.0f, 0.0f, 0.0f};
float node_id_color[3] = {0.0f, 1.0f, 1.0f};
float elem_id_color[3] = {1.0f, 1.0f, 0.0f};


