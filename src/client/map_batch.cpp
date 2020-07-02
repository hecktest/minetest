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

#include "map_batch.h"
#include "mesh.h"
#include "util/numeric.h"
#include <EHardwareBufferFlags.h>
#include <cstring>


//	BatchMeshRef

bool BatchMeshRef::setAndCheck(scene::IMeshBuffer* c_mesh) {
	u64 new_id = 0;
	if (c_mesh)
		new_id = (u64)c_mesh->getChangedID_Vertex() | ((u64)(c_mesh->getChangedID_Index()) << 32);
	
	bool meshChanged	= mesh != c_mesh;
	bool idChanged		= new_id != change_id;
	
	change_id = new_id;
	
	if(c_mesh)
		c_mesh->grab();
	if(mesh)
		mesh->drop();
	
	mesh = c_mesh;
	
	return meshChanged || idChanged;
}


//	MapBatch

bool MapBatch::setSourceMesh(int i, scene::IMeshBuffer* mesh, v3f offset)
{
	bool hadMesh = sources[i].mesh ? 1 : 0;
	
	if (mesh && mesh->getIndexCount() < 0)
		mesh = NULL; // We have no need for a mesh that's empty.
	
	if (mesh) {
		vertex_type = mesh->getVertexType();
		sources[i].offset = offset;
	}
	
	bool check = sources[i].setAndCheck(mesh);
	
	
	if (mesh && !hadMesh)
		mesh_count++;
	if (!mesh && hadMesh)
		mesh_count--;
	
	dirty = dirty || check;
	return check;
}

void MapBatch::commitDrawMesh(scene::IMeshBuffer* buf)
{
	if (buf) {
		buf->setHardwareMappingHint( scene::EHM_STATIC );
		buf->setDirty();
		mesh_bufs.push_back(buf);
	}
}

void MapBatch::clearBuffers()
{
	if (mesh_bufs.size() > 0) {
		for (auto buf : mesh_bufs) {
			driver->removeHardwareBuffer(buf);
			buf->drop();
		}
		mesh_bufs.clear();
	}
}

// where T: irr::video::S3DVertex
template<typename T> void MapBatch::rebuildDrawMeshesT()
{
	static const constexpr u32 max_verts = 1 << 16;
	clearBuffers();
	u32 done = 0;
	u32 vert_count = 0;
	u32 index_count = 0;
	std::vector<std::pair<v3f,scene::IMeshBuffer*>> current;
	bool batch_has_box = false;
	for (u32 i = 0; i < BATCH_COUNT; i++) {
		scene::IMeshBuffer* mesh = sources[i].mesh;
		if (mesh) {
			aabb3f box = mesh->getBoundingBox();
			box.MaxEdge += sources[i].offset;
			box.MinEdge += sources[i].offset;
			if (batch_has_box)
				occludee.bounds.addInternalBox(box);
			else {
				occludee.bounds = box;
				batch_has_box = true;
			}
			u32 next_vert_count = mesh->getVertexCount();
			u32 next_index_count = mesh->getIndexCount();
			done++;
			if (((vert_count + next_vert_count) >= max_verts ) || done == mesh_count) {
				if (done == mesh_count) {
					vert_count += next_vert_count;
					index_count += next_index_count;
					current.push_back(std::make_pair(sources[i].offset,sources[i].mesh));
				}
				
				scene::CMeshBuffer<T>* current_mesh = new scene::CMeshBuffer<T>;
				
				u32 vertex_offset = 0;
				
				current_mesh->Vertices.set_used(vert_count);
				current_mesh->Indices.set_used(index_count);
				
				register T* v = current_mesh->Vertices.pointer();
				register u16* n = current_mesh->Indices.pointer();
				
				bool has_box = false;
				for (auto c : current) {
					u32 c_vcount = c.second->getVertexCount();
					u32 c_icount = c.second->getIndexCount();
					u16* c_indices = c.second->getIndices();
					T* c_verts = (T*)c.second->getVertices();
					for (u32 j=0; j < c_vcount; j++) {
						*v = c_verts[j];
						v++->Pos += c.first;
					}
					for (u32 k=0; k < c_icount; k++)
						*(n++) = c_indices[k] + vertex_offset;
					
					vertex_offset += c_vcount;
					aabb3f box = c.second->getBoundingBox();
					box.MaxEdge += c.first;
					box.MinEdge += c.first;
					if (has_box)
						current_mesh->BoundingBox.addInternalBox(box);
					else {
						current_mesh->setBoundingBox(box);
						has_box = true;
					}
				}
				commitDrawMesh(current_mesh);
				vert_count = 0;
				index_count = 0;
				current.clear();
				if (done == mesh_count) {
					return;
				}
			}
			vert_count += next_vert_count;
			index_count += next_index_count;
			current.push_back(std::make_pair(sources[i].offset,sources[i].mesh));
		}
	}
}

void MapBatch::rebuildDrawMeshes() {
	switch (vertex_type) {
		case video::EVT_2TCOORDS:
			rebuildDrawMeshesT<video::S3DVertex2TCoords>();
			return;
		case video::EVT_TANGENTS:
			rebuildDrawMeshesT<video::S3DVertexTangents>();
			return;
		default:
			rebuildDrawMeshesT<video::S3DVertex>();
	}
}

bool MapBatch::update()
{
	bool r = dirty;
	if (dirty) {
		rebuildDrawMeshes();
		dirty = false;
	}
	return r;
}