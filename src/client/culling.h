/*
Minetest
Copyright (C) 2010-2020 celeron55, Perttu Ahola <celeron55@gmail.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/
#pragma once

#include <vector3d.h>
#include <aabbox3d.h>
#include <matrix4.h>
#include <IImage.h>
#include <ITexture.h>
#include <IVideoDriver.h>
#include <SColor.h>
#include <cmath>
#include <vector>
#include "irrlichttypes.h"
#include "irr_v3d.h"
#include "irr_aabb3d.h"

using irr::core::matrix4;
using irr::video::ITexture;
using irr::video::IImage;
using irr::video::IVideoDriver;
using irr::video::SColor;

struct CullBuffer {
	// size in 8x8 tiles
	s32 w;
	s32 h;

	// layer count
	s32 depth;
	
	// current sprite
	s32 spr_x;
	s32 spr_y;
	s32 spr_w;
	s32 spr_h;
	s32 spr_z;
	s32* span_l;		// span starts
	s32* span_r;		// span ends
	u64* spr_data; 		// tile buffer, must be WxH
	void clearSprite()
	{
		memset(spr_data + spr_y * w + spr_x, 0, (spr_h * w + spr_w) * sizeof(u64));
	}

	// layout: row, column, slice
	u64 *data = NULL;
	size_t data_length;

	// delete before commit
	ITexture* debug_tex;
	void draw(IVideoDriver* driver) {
		auto img = driver->createImage(video::ECF_R8G8B8,core::dimension2du(w,h));
		for( s32 y = 0; y < h; ++y ) {			
			for( s32 x = 0; x < w; ++x ) {
				float d = 0;
				img->setPixel(x,h-y-1, SColor(0xFF000000));
				for( s32 z = 0; z < depth; ++z ) {
					d = (float)(z+1) / (float)(depth);          
					u64 tile = data[(w*y+x)*depth+z];
					if (tile) {
						if ( tile == 2 ) {
							img->setPixel(x,h-y-1, SColor(0xFFFFFF00));							
						}else{
							img->setPixel(x,h-y-1, SColor(0xFF00FFFF).getInterpolated(SColor(0xFFA000A0),1-d));
						}
						break;
					}
				}
			}
		}
		auto size = driver->getCurrentRenderTargetSize();
		debug_tex = driver->addTexture("DEBUG_VIS",img);
		driver->draw2DImage(debug_tex, core::recti(0,0,size.Width,size.Height), core::recti(0,0,w,h));
		driver->removeTexture(debug_tex);
		img->drop();
	}
	
	void rasterizeSprite(int n_verts, v3f* polygon, aabb3f bounds);
	bool testSprite();
	void writeSprite();
	void clear()
	{
		memset(data, 0, data_length * sizeof(u64));
	}
	// Normalized depth.
	s32 depthToLayer(float d)
	{
		return (s32)(depth * sqrt(d));
	}
	CullBuffer()
	{
		data = NULL;
		debug_tex = NULL;
		span_l = NULL;
		span_r = NULL;
		spr_data = NULL;
	}
	CullBuffer(s32 sx, s32 sy, s32 d): w{sx}, h{sy}, depth{d}
	{
		data_length = w * h * d;
		debug_tex = NULL;
		data = new u64[data_length]();
		span_l = new s32[h];
		span_r = new s32[h];
		spr_data = new u64[w*h]();
	}
	~CullBuffer()
	{
		if (data)
			delete[] data;
		if (span_l)
			delete[] span_l;
		if (span_r)
			delete[] span_r;
		if (spr_data)
			delete[] spr_data;
	}
};

struct CullNode {
	matrix4 last_mvp;
	v3f last_box_pos;
	aabb3f bounds;
	std::vector<aabb3f> occluders;
	v3f ndc_verts[8];
	aabb3f ndc_bounds;
	bool touching_scr;
	
	void cacheClipSpaceData(matrix4 mvp);
	bool visible(matrix4 mvp, v3f box_pos, v3f eye_pos, CullBuffer& buffer);
};
