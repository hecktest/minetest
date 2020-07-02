/*
Minetest
Copyright (C) 2010-2013 celeron55, Perttu Ahola <celeron55@gmail.com>

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

#include "irrlichttypes_extrabloated.h"
#include "map.h"
#include "camera.h"
#include "map_batch.h"
#include <set>
#include <map>

// We need this to use a material hashmap.
template<> struct std::hash<irr::video::SMaterial>
{
	size_t operator()(const irr::video::SMaterial& self) const
	{
		return // Just the data relevant to Minetest.
			(size_t)self.getTexture(0) ^
			((size_t)self.getTexture(1) << 3) ^
			(size_t)self.ZWriteEnable << 11 ^
			(size_t)self.Wireframe << 9;
	}
};
// Make positions hashable.
template<> struct std::hash<v3s16>
{
	size_t operator()(const v3s16& self) const
	{
		return
			(((size_t)self.X + 0x8000) & 0xFFFF) ^
			rollst(((size_t)self.Z + 0x8000) & 0xFFFF, 16) ^
			rollst(((size_t)self.Y + 0x8000) & 0xFFFF, 32);
	}
};

struct MapDrawControl
{
	// Overrides limits by drawing everything
	bool range_all = false;
	// Wanted drawing range
	float wanted_range = 0.0f;
	// show a wire frame for debugging
	bool show_wireframe = false;
};

class Client;
class ITextureSource;

/*
	ClientMap

	This is the only map class that is able to render itself on screen.
*/

class ClientMap : public Map, public scene::ISceneNode
{
public:
	ClientMap(
			Client *client,
			MapDrawControl &control,
			s32 id
	);

	virtual ~ClientMap() = default;

	s32 mapType() const
	{
		return MAPTYPE_CLIENT;
	}

	void drop()
	{
		ISceneNode::drop();
	}

	void updateCamera(const v3f &pos, const v3f &dir, f32 fov, const v3s16 &offset, const aabb3f &box)
	{
		m_camera_position = pos;
		m_camera_box = box;
		m_camera_direction = dir;
		m_camera_fov = fov;
		m_camera_offset = offset;
	}

	/*
		Forcefully get a sector from somewhere
	*/
	MapSector * emergeSector(v2s16 p);

	//void deSerializeSector(v2s16 p2d, std::istream &is);

	/*
		ISceneNode methods
	*/

	virtual void OnRegisterSceneNode();

	virtual void render()
	{
		video::IVideoDriver* driver = SceneManager->getVideoDriver();
		driver->setTransform(video::ETS_WORLD, AbsoluteTransformation);
		renderMap(driver, SceneManager->getSceneNodeRenderPass());
	}

	virtual const aabb3f &getBoundingBox() const
	{
		return m_box;
	}

	void getBlocksInViewRange(v3s16 cam_pos_nodes,
		v3s16 *p_blocks_min, v3s16 *p_blocks_max);
	void updateDrawList();
	void renderMap(video::IVideoDriver* driver, s32 pass);

	int getBackgroundBrightness(float max_d, u32 daylight_factor,
			int oldvalue, bool *sunlight_seen_result);

	void renderPostFx(CameraMode cam_mode);

	// For debug printing
	// virtual void PrintInfo(std::ostream &out);

	const MapDrawControl & getControl() const { return m_control; }
	f32 getCameraFov() const { return m_camera_fov; }
private:
	Client *m_client;

	aabb3f m_box = aabb3f(-BS * 1000000, -BS * 1000000, -BS * 1000000,
		BS * 1000000, BS * 1000000, BS * 1000000);

	MapDrawControl &m_control;

	v3f m_camera_position = v3f(0,0,0);
	aabb3f m_camera_box;
	v3f m_camera_direction = v3f(0,0,1);
	f32 m_camera_fov = M_PI;
	v3s16 m_camera_offset;

	bool m_cache_trilinear_filter;
	bool m_cache_bilinear_filter;
	bool m_cache_anistropic_filter;
	
	std::unordered_map<video::SMaterial,std::unordered_map<v3s16, MapBatch*>*> m_batches;
	
	CullBuffer* m_vis_buffer;
	
	std::vector<MapBlock*> m_grabbed_blocks;
};
