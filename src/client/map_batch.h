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

#include "culling.h"
#include "irrlichttypes_bloated.h"
#include <vector>
#include <map>
#include <IMeshBuffer.h>
#include <ISceneManager.h>

using namespace irr;

// 4x4x4 blocks yields us up to 64 meshes being batched into one.
// This should result in pretty neat batch sizes.
static const constexpr int BATCH_STRIDE = 4;
static const constexpr int BATCH_COUNT = BATCH_STRIDE * BATCH_STRIDE * BATCH_STRIDE;

// Stores a mesh reference and tells us if it has been touched.
struct BatchMeshRef
{
	v3f	offset;
	// Last pointer, last changeID.
	scene::IMeshBuffer* mesh = NULL;
	u64 change_id = 0;
	// Set the above and tell us if they were different.
	bool setAndCheck(scene::IMeshBuffer* c_mesh);
	~BatchMeshRef() {
		if (mesh)
			mesh->drop();
	}
};

/*
	A MapBatch represents a bigger chunk of the world,
	only meant for VBO rendering. It's a batch of several
	mapblock meshes sharing a single material.
*/
class MapBatch
{
private:

	// Drop any mesh buffers and wipe the vector.
	void clearBuffers();
	
	// Remake mesh buffers from the sources.
	template<typename T> void rebuildDrawMeshesT();
	void rebuildDrawMeshes();
	
	// Add mesh if not null, and rebuild its bounds.
	void commitDrawMesh(scene::IMeshBuffer* mesh);
	
	video::E_VERTEX_TYPE vertex_type;
	
	video::IVideoDriver* driver;
	
public:
	u16		mesh_count = 0;
	
	CullNode occludee;
	
	// Upstream meshes.
	BatchMeshRef sources[BATCH_COUNT];
	
	// Dirty flag.
	bool dirty;
	
	
	// Buffers that are actually being displayed.
	// This can be empty, but you should just
	// kill the batch if it stays empty for too long.
	std::vector<scene::IMeshBuffer*> mesh_bufs;

	// Set a source mesh for index i.
	// If it has changed or been modified, set dirty.
	bool setSourceMesh(int i, scene::IMeshBuffer* mesh, v3f offset);

	// Check if the batch is dirty, and update if needed.
	bool update();
	
	MapBatch(v3s16 c_offset, scene::ISceneManager* scene) {
		driver = scene->getVideoDriver();
	}
	~MapBatch()
	{
		clearBuffers();
	}
};
