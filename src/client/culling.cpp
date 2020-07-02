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
#include "culling.h"
#include <cmath>
#include "log.h"
#include <cstdio>

void CullBuffer::writeSprite() {
	if (!spr_w || !spr_h)
		return;
	// OR every matching tile on this layer
	// with the sprite's corresponding tiles.
	for (s32 y=spr_y; y < spr_y+spr_h; ++y) {
		for (s32 x=spr_x; x < spr_x+spr_w; ++x) {
			data[(y * w + x) * depth + spr_z] |= spr_data[y * w + x];
		}
	}
}

bool CullBuffer::testSprite() {
	if (!spr_w || !spr_h)
		return false;
	if (spr_z == 0)
		return true;
	for (s32 y=spr_y; y < spr_y+spr_h; ++y) {
		for (s32 x=spr_x; x < spr_x+spr_w; ++x) {
			u64 tile = spr_data[y * w + x];
			if ( !tile ) { continue; } // next
			u64 layers = 0;
			for (s32 z = spr_z-1; z >= 0; --z)
				layers |= data[(y * w + x) * depth + z];
			if ( !layers ) { return true; } // pass
			tile &= ~layers;
			if ( tile ) { return true; } // pass
		}
	}
	return false; // fail
}

inline int isign(int x) {
	return (x > 0) - (x < 0);
}

inline float lerp( float x, float y, float t ) {
	return ( y - x ) * t + x;
}


// polygon and bounds are in clip space
void CullBuffer::rasterizeSprite(int n_verts, v3f* polygon, aabb3f bounds) {
	using std::abs;
	using std::min;
	using std::max;
	using std::floor;
	using std::printf;
	
	// NDC to raster
	bounds.MinEdge = bounds.MinEdge * 0.5 + v3f(0.5,0.5,0.5);
	bounds.MaxEdge = bounds.MaxEdge * 0.5 + v3f(0.5,0.5,0.5);
	
	spr_w = floor((bounds.MaxEdge.X - bounds.MinEdge.X) * w + 1);
	spr_h = floor((bounds.MaxEdge.Y - bounds.MinEdge.Y) * h + 1);
	spr_x = floor(bounds.MinEdge.X * w);
	spr_y = floor(bounds.MinEdge.Y * h);
	// Clip top and left
	if (spr_x < 0) {
		spr_w += spr_x;
		spr_x = 0;
	}
	if (spr_y < 0) {
		spr_h += spr_y;
		spr_y = 0;
	}
	// Clip bottom and right
	spr_w = max(0,min(spr_w, w - spr_x));
	spr_h = max(0,min(spr_h, h - spr_y));
	if (!spr_w || !spr_h)
		return; // empty sprite
	
	spr_z = (s32)((1 - min(1.0f,500 * max(0.0f,bounds.MaxEdge.Z*0.5f+bounds.MinEdge.Z*0.5f))) * (float)(depth-1));
	
	// Wipe sprite with FF:
	// Most of the surface area is inside the polygon,
	// so this takes care of that part.
	memset((void*)spr_data + spr_y * w * 8, 0, spr_h * w * 8);
	
	for (int i = spr_y; i < spr_y + spr_h; ++i) {
		// Set a negative span on each raster, so that min/max can be used
		// to shrinkwrap the real thing.
		span_l[i] = w;
		span_r[i] = -1;
		// span_l[i] = spr_x;
		// span_r[i] = spr_x + spr_w;
	}
	
	v3f pos = polygon[0] * 0.5 + v3f(0.5,0.5,0.5);
	for (int i = 1; i <= n_verts; ++i) {
		v3f next = polygon[i%n_verts] * 0.5 + v3f(0.5,0.5,0.5);
		// Clip the line against the screen edge; skip if it misses the rect
		v3f c_cur = pos;
		v3f c_nxt = next;
		bool next_out = next.X<0 || next.Y<0 || next.X>1 || next.Y>1;
		bool pos_out = pos.X<0 || pos.Y<0 || pos.X>1 || pos.Y>1;
		if (next_out || pos_out) {
		// if (0) {
			float d_nxt = 1e99;
			float d_cur = 1e99;
			// A part of the line is outside, clip it.
			
			// Turn it into mx+b form:
			float m = (next.Y - pos.Y) / (next.X - pos.X);
			float b = pos.Y - pos.X * m;
			
			// Solve for X=0, X=1, Y=0, Y=1
			float y0 = b;		// Left
			float y1 = m+b;		// Right
			float x0 = (-b)/m;	// Down
			float x1 = (1-b)/m;	// Up
			
			if ((x0 < 0 && x1 < 0) ||	// Miss below
				(y0 < 0 && y1 < 0) ||	// Miss to the left
				(x0 > 1 && x1 > 1) ||	// Miss above
				(y0 > 1 && y1 > 1) )	// Miss to the right
				{
					pos = next;
					continue; // No intersection
			}
			
			v3f clips[4] = {
				v3f(0, y0, 0),
				v3f(x0, 0, 0),
				v3f(1, y1, 0),
				v3f(x1, 1, 0)
			};
			
			#if __GNUC__ >= 8
				#pragma GCC unroll 4
			#endif
			for (int j=0; j<4; j++) {
				v3f p = clips[j];
				if (p.X < 0 || p.Y < 0 || p.X > 1 || p.Y > 1)
					continue;				
				if (next_out) {
					v3f d = c_nxt - p;
					float dist = d.dotProduct( d );
					if (dist < d_nxt) {
						c_nxt = p;
						d_nxt = dist;
					}
				}
				if (pos_out) {
					v3f d = c_cur - p;
					float dist = d.dotProduct( d );
					if (dist < d_cur) {
						c_cur = p;
						d_cur = dist;
					}
				}
			}
		}
		
		// Current tile
		int ct[2] = {c_cur.X * w, c_cur.Y * h};
		int st[2] = {ct[0], ct[1]};
		// Next tile
		int nt[2] = {c_nxt.X * w, c_nxt.Y * h};
		// Line bounds in tiles
		int lbx = abs(nt[0]-ct[0]);
		int lby = abs(nt[1]-ct[1]);
		// intercept form of the line (haha, again?)
		// but make it relative to a single (starting) tile
		// float m = ((c_nxt.Y - c_cur.Y) / (c_nxt.X - c_cur.X)) * ((float)w/(float)h);
		// float b = (c_cur.Y - c_cur.X * m) * h - ct[1];
		// ah, and now our old direction info is useless
		// let 0,1,2,3 represent down, left, up, right, then:
		// float t[4] = { 0, 0, 0, 0 }; // these are the intercepts for those edges
		if ( lbx+lby > w+h ) {
			warningstream << "Walking an impossibly long line: " <<
				c_cur.X<<","<<c_cur.Y<<" -> "<<
				c_nxt.X<<","<<c_nxt.Y<<"\n";
		}
		if (!lbx) { // straight vertical
			int dir = isign(nt[1] - ct[1]);
			while(ct[1] != nt[1]) {
				if (dir > 0){
					span_l[ct[1]] = ct[0];
				} else {
					span_r[ct[1]] = ct[0];
				}
				spr_data[ct[1] * w + ct[0]] = 2;
				ct[1] += dir;
			}
		} else if (!lby) { // straight horizontal
			int dir = isign(nt[0] - ct[0]);
			if (dir > 0){
				span_l[ct[1]] = ct[0];
				span_r[ct[1]] = nt[0];
			} else {
				span_l[ct[1]] = nt[0];
				span_r[ct[1]] = ct[0];
			}
			while(ct[0] != nt[0]) {
				spr_data[ct[1] * w + ct[0]] = 2;
				ct[0] += dir;
			}
		} else { // 2D modified bresenham
			int dist = lbx + lby;
			int j = 0;
			// Determine major dimension
			int smaj = lby > lbx;
			int smin = !smaj;
			// Determine major direction
			bool up = isign(nt[1] - st[1]) > 0;
			int dmaj = isign(nt[smaj] - ct[smaj]);
			int dmin = isign(nt[smin] - ct[smin]);
			float lmaj = abs(nt[smaj] - ct[smaj]);
			while (1) {
				// Rasterize the edge precisely by clipping each tile with
				// a binary bisection created from the current line.
				// This CSG-like construction handles all cases,
				// even the polygon being contained by the tile.
				// TODO: actually do it, and use a LUT for speed
				if (ct[1] < h && ct[0] < w && ct[1] >=0 && ct[0] >= 0) {
					spr_data[ct[1] * w + ct[0]] = 2;
				
					// expand span using current tile X
					if (up)
						span_l[ct[1]] = min(span_l[ct[1]],ct[0]);
					else
						span_r[ct[1]] = max(span_r[ct[1]],ct[0]);
				}
				// break;
				
				// reached dest?
				if(ct[0] == nt[0] && ct[1] == nt[1])
					break;
				
				// if ( ++j > dist ) { 
					// warningstream << "Too many iterations when walking: " << j << "/" << dist <<"\n";
					// break;
				// }
				
				// if at a crossing, update minor but not major,
				// that is, don't move diagonally,
				// so that our tiles all contain the polygon's edge
				float p = (ct[smaj] - st[smaj]) * dmaj / (float)lmaj;
				if (floor(lerp(st[smin],nt[smin],p)) != ct[smin])
					ct[smin] += dmin;
					// ct[smin] = floor(lerp(st[smin],nt[smin],p));
				else
					ct[smaj] += dmaj;
			}
		}
		
		pos = next;
	}
	
	// Clear the space outside of the polygon
	for (int i = spr_y; i < spr_y + spr_h; ++i) {
		// fill 0 to left
		// if (span_l[i] > spr_x)
			// memset((void*)(spr_data + w * i + spr_x), 0, (span_l[i] - spr_x) * 8);
		// fill right to w-1
		// if (span_r[i] < spr_x + spr_w - 1)
			// memset((void*)(spr_data + w * i + span_r[i] + 1), 0, (spr_x + spr_w - 1 - span_r[i]) * 8);
	}
}

// The planes of a box divide 3D space into 27 cells.
// Whichever cell the eye is in determines the set of edges
// that form the contour of the box, except for the middle cell,
// which is the box itself.
// The resulting polygons are wound clockwise.
int borderVerts[][6] = {
	// Refer to irrlicht's aabbox3d.h for the index layout.
	
	// Left Low
	{1,5,4,6,2,3},	// 0
	{1,0,4,6,2,3},	// 0-2
	{1,0,4,6,7,3},	// 2
	// Left Mid
	{4,0,2,3,1,5},	// 1-0
	{1,0,2,3,-1,-1},// -X
	{1,0,2,6,7,3},	// 3-2
	// Left High
	{7,5,4,0,2,3},	// 1
	{7,5,1,0,2,3},	// 1-3
	{7,5,1,0,2,6},	// 3
	
	// Center Low
	{0,1,5,4,6,2},	// 0-4
	{0,4,6,2,-1,-1},// -Y
	{0,4,6,7,3,2},	// 2-6
	// Center Mid
	{0,1,5,4,-1,-1},// -Z
	{8,8,8,8,8,8},	// Inside - instant pass, do not read
	{7,3,2,6,-1,-1},// +Z
	// Center High
	{1,3,7,5,4,0},	// 1-5
	{1,3,7,5,-1,-1},// +Y
	{1,3,2,6,7,5},	// 3-7
	
	// Right Low
	{1,5,7,6,2,0},	// 4
	{4,5,7,6,2,0},	// 4-6
	{4,5,7,3,2,0},	// 6
	// Right Mid
	{1,5,7,6,4,0},	// 4-5
	{5,7,6,4,-1,-1},	// +X
	{5,7,3,2,6,4},	// 6-7
	// Right High
	{1,3,7,6,4,0},	// 5
	{1,3,7,6,4,5},	// 5-7
	{1,3,2,6,4,5},	// 7
	// For some reason index 1 can't be last in a polygon. Odd.
};
// Pairs of closest/furthest vertices.
// Similar layout to the edge lists,
// but based on a 2x2x2 partitioning instead.
int closestVertex[][2] = {
	{0,7},
	{2,5},
	{1,6},
	{3,4},
	{4,3},
	{6,1},
	{5,2},
	{7,0},
};

inline v3f mul4x4p( matrix4 tx, v3f point ) {
	float r[4];
	tx.transformVect( r, point );
	float z = -r[2] / r[3];
	return v3f(
		// (r[2] >= -1) fixes a problem with XY
		// flipping past the near plane,
		// initially this was r[2] >= 0
		// change it back if things fail frustum when they shouldn't
		(z >= -1) * r[0] / r[3],
		(z >= -1) * r[1] / r[3],
		z
	);
}

void CullNode::cacheClipSpaceData( matrix4 mvp ) {
	bounds.getEdges( ndc_verts ); // lol irrlicht
	ndc_verts[0] = mul4x4p( mvp, ndc_verts[0] );
	ndc_bounds = aabb3f( ndc_verts[0] );
	touching_scr = ndc_verts[0].Z < -1;;
	for ( int i = 1; i < 8; i++ ) {
		ndc_verts[i] = mul4x4p( mvp, ndc_verts[i] );
		touching_scr |= ndc_verts[i].Z < -1;
		ndc_bounds.addInternalPoint( ndc_verts[i] );
	}
}

aabb3f ndc = aabb3f(v3f(-1,-1,-1),v3f(1,1,1));

bool CullNode::visible(matrix4 mvp, v3f box_pos, v3f eye_pos, CullBuffer& buffer) {
	// Frustum cull
	cacheClipSpaceData( mvp );
	if (
		// Sorted by distribution, roughly.
		ndc_bounds.MaxEdge.Y < -1 ||
		ndc_bounds.MinEdge.X > 1 ||
		ndc_bounds.MaxEdge.X < -1 ||
		ndc_bounds.MinEdge.Y > 1 ||
		// This should be after a frustum AABB broad step, so near/far planes
		// will not have much space behind them.
		ndc_bounds.MinEdge.Z > 1 ||
		ndc_bounds.MaxEdge.Z < -1
		) {
		return false; // Frustum fail.
	}
	
	// Frustum passed, do occlusion
	v3f box_min_w = box_pos + bounds.MinEdge;
	v3f box_max_w = box_pos + bounds.MaxEdge;
	if (ndc_bounds.MinEdge.Z < -1)
		return true; // No valid primitive for occlusion.
	int edge_set =
		9 * ((eye_pos.X > box_min_w.X) + (eye_pos.X > box_max_w.X)) +
		3 * ((eye_pos.Y > box_min_w.Y) + (eye_pos.Y > box_max_w.Y)) +
			((eye_pos.Z > box_min_w.Z) + (eye_pos.Z > box_max_w.Z));
	if ( edge_set == 13 ) // the camera is inside the box
		return true;
	if (touching_scr) {
		warningstream << "squeek!\n";
	}
	int n_verts = borderVerts[edge_set][4] > 0 ? 6 : 4;
	v3f verts[6];
	for ( int i = 0; i < n_verts; i++ ) {
		if (borderVerts[edge_set][i] < 0) { break; }
		verts[i] = ndc_verts[borderVerts[edge_set][i]];
	}
	buffer.rasterizeSprite(n_verts, verts, ndc_bounds);
	buffer.writeSprite();
	// ...
	// Occlusion passed, blit occluders if they exist
	return true;
}

