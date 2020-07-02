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

#include "clientmap.h"
#include "client.h"
#include "mapblock_mesh.h"
#include <IMaterialRenderer.h>
#include <matrix4.h>
#include "mapsector.h"
#include "mapblock.h"
#include "profiler.h"
#include "settings.h"
#include "camera.h"               // CameraModes
#include "util/basic_macros.h"
#include <algorithm>
#include "client/renderingengine.h"
#include "client/mesh.h"

ClientMap::ClientMap(
		Client *client,
		MapDrawControl &control,
		s32 id
):
	Map(dout_client, client),
	scene::ISceneNode(RenderingEngine::get_scene_manager()->getRootSceneNode(),
		RenderingEngine::get_scene_manager(), id),
	m_client(client),
	m_control(control)
{
	/* TODO: Add a callback function so these can be updated when a setting
	 *       changes.  At this point in time it doesn't matter (e.g. /set
	 *       is documented to change server settings only)
	 *
	 * TODO: Local caching of settings is not optimal and should at some stage
	 *       be updated to use a global settings object for getting thse values
	 *       (as opposed to the this local caching). This can be addressed in
	 *       a later release.
	 */
	m_cache_trilinear_filter  = g_settings->getBool("trilinear_filter");
	m_cache_bilinear_filter   = g_settings->getBool("bilinear_filter");
	m_cache_anistropic_filter = g_settings->getBool("anisotropic_filter");
	m_vis_buffer = new CullBuffer(160,120,32);
}

MapSector * ClientMap::emergeSector(v2s16 p2d)
{
	// Check that it doesn't exist already
	MapSector *sector = getSectorNoGenerate(p2d);

	// Create it if it does not exist yet
	if (!sector) {
		sector = new MapSector(this, p2d, m_gamedef);
		m_sectors[p2d] = sector;
	}

	return sector;
}

void ClientMap::OnRegisterSceneNode()
{
	if(IsVisible)
	{
		SceneManager->registerNodeForRendering(this, scene::ESNRP_SOLID);
		SceneManager->registerNodeForRendering(this, scene::ESNRP_TRANSPARENT);
	}

	ISceneNode::OnRegisterSceneNode();
}

void ClientMap::getBlocksInViewRange(v3s16 cam_pos_nodes,
		v3s16 *p_blocks_min, v3s16 *p_blocks_max)
{
	v3s16 box_nodes_d = m_control.wanted_range * v3s16(1, 1, 1);
	// Define p_nodes_min/max as v3s32 because 'cam_pos_nodes -/+ box_nodes_d'
	// can exceed the range of v3s16 when a large view range is used near the
	// world edges.
	v3s32 p_nodes_min(
		cam_pos_nodes.X - box_nodes_d.X,
		cam_pos_nodes.Y - box_nodes_d.Y,
		cam_pos_nodes.Z - box_nodes_d.Z);
	v3s32 p_nodes_max(
		cam_pos_nodes.X + box_nodes_d.X,
		cam_pos_nodes.Y + box_nodes_d.Y,
		cam_pos_nodes.Z + box_nodes_d.Z);
	// Take a fair amount as we will be dropping more out later
	// Umm... these additions are a bit strange but they are needed.
	*p_blocks_min = v3s16(
			p_nodes_min.X / MAP_BLOCKSIZE - 3,
			p_nodes_min.Y / MAP_BLOCKSIZE - 3,
			p_nodes_min.Z / MAP_BLOCKSIZE - 3);
	*p_blocks_max = v3s16(
			p_nodes_max.X / MAP_BLOCKSIZE + 1,
			p_nodes_max.Y / MAP_BLOCKSIZE + 1,
			p_nodes_max.Z / MAP_BLOCKSIZE + 1);
}

void ClientMap::updateDrawList()
{
	// ScopeProfiler _(g_profiler, "CM::updateDrawList()", SPT_AVG);
	TimeTaker scope("updateDrawList scope");
	std::string prefix = "updateDrawList: ";
	
	for (int i=0; i < m_grabbed_blocks.size(); i++ )
		m_grabbed_blocks[i]->refDrop();
	m_grabbed_blocks.clear();
	
	v3f camera_position = m_camera_position;
	v3f camera_direction = m_camera_direction;
	f32 camera_fov = m_camera_fov;

	// Use a higher fov to accomodate faster camera movements.
	// Blocks are cropped better when they are drawn.
	// Or maybe they aren't? Well whatever.
	camera_fov *= 1.2;

	v3s16 cam_pos_nodes = floatToInt(camera_position, BS);
	v3s16 p_blocks_min;
	v3s16 p_blocks_max;
	getBlocksInViewRange(cam_pos_nodes, &p_blocks_min, &p_blocks_max);

	// Number of blocks with mesh in rendering range
	u32 blocks_in_range_with_mesh = 0;
	// Number of blocks occlusion culled
	u32 blocks_occlusion_culled = 0;

	// No occlusion culling when free_move is on and camera is
	// inside ground
	bool occlusion_culling_enabled = true;
	if (g_settings->getBool("free_move") && g_settings->getBool("noclip")) {
		MapNode n = getNode(cam_pos_nodes);
		if (n.getContent() == CONTENT_IGNORE ||
				m_nodedef->get(n).solidness == 2)
			occlusion_culling_enabled = false;
	}

	// Uncomment to debug occluded blocks in the wireframe mode
	// TODO: Include this as a flag for an extended debugging setting
	//if (occlusion_culling_enabled && m_control.show_wireframe)
	//    occlusion_culling_enabled = porting::getTimeS() & 1;
	/*
		Get animation parameters
	*/
	float animation_time = m_client->getAnimationTime();
	int crack = m_client->getCrackLevel();
	u32 daynight_ratio = m_client->getEnv().getDayNightRatio();
	
	u32 meshes_set = 0;
	TimeTaker gather("Block gathering");
	TimeTaker batch_notify("Batch updates");

	for (const auto &sector_it : m_sectors) {
		MapSector *sector = sector_it.second;
		v2s16 sp = sector->getPos();
		
		if (!m_control.range_all) {
			if (sp.X < p_blocks_min.X || sp.X > p_blocks_max.X ||
					sp.Y < p_blocks_min.Z || sp.Y > p_blocks_max.Z)
				continue;
		}

		MapBlockVect sectorblocks;
		sector->getBlocks(sectorblocks);

		/*
			Loop through blocks in sector
		*/

		u32 sector_blocks_drawn = 0;

		for (MapBlock *block : sectorblocks) {
			/*
				Compare block position to camera position, skip
				if not seen on display
			*/
						

			float range = 100000 * BS;
			if (!m_control.range_all)
				range = m_control.wanted_range * BS;

			// This block is in range. Reset usage timer.
			block->resetUsageTimer();
			// Add to set
			
			v3s16 block_pos = block->getPos();
			sector_blocks_drawn++;
			
			v3s16 batch_pos = v3s16(
				( block_pos.X + (block_pos.X < 0) ) / BATCH_STRIDE - (block_pos.X < 0),
				( block_pos.Y + (block_pos.Y < 0) ) / BATCH_STRIDE - (block_pos.Y < 0),
				( block_pos.Z + (block_pos.Z < 0) ) / BATCH_STRIDE - (block_pos.Z < 0)
			);
			// operator % breaks with negative numbers
			#define SMOD(l,r) (((l)%(r)+(r))%(r))
			v3s16 local_pos = v3s16(
				SMOD(block_pos.X, BATCH_STRIDE),
				SMOD(block_pos.Y, BATCH_STRIDE),
				SMOD(block_pos.Z, BATCH_STRIDE)
			);
			#undef SMOD
			
			// Swizzle the addressing so that even if a batch is forced
			// to split, the pieces will be spatially compact.
			// The next cell is always a neighbor of the last one.
			s16 batch_cell = local_pos.Y * BATCH_STRIDE * BATCH_STRIDE;
			
			if (local_pos.Y % 2)
				batch_cell += (BATCH_STRIDE - local_pos.Z - 1) * BATCH_STRIDE;
			else
				batch_cell += local_pos.Z * BATCH_STRIDE;
			
			if ((local_pos.Z + local_pos.Y) % 2)
				batch_cell += BATCH_STRIDE - local_pos.X - 1;
			else
				batch_cell += local_pos.X;
			
			
			MapBlockMesh *blockMesh = block->mesh;
			
			if (!blockMesh) {
				// for (auto layer = m_batches.begin(); layer != m_batches.end(); ++layer) {
					// for (auto batch = layer->second->begin(); batch != layer->second->end(); ++batch) {
						// if (batch->first == batch_pos) {
							// batch->second->setSourceMesh(batch_cell, NULL, v3f(0,0,0));
						// }
					// }
				// }
				continue;
			}
			
			block->refGrab();
			m_grabbed_blocks.push_back(block);
			// For budgeting animations, and for profiling
			/*
			u32 mesh_animate_count = 0;
			// Pretty random but this should work somewhat nicely
			bool faraway = d >= BS * 50;
			if (blockMesh->isAnimationForced() || !faraway ||
					mesh_animate_count < (m_control.range_all ? 200 : 50)) {

				bool animated = blockMesh->animate(faraway, animation_time,
					crack, daynight_ratio);
				if (animated)
					mesh_animate_count++;
			} else {
				blockMesh->decreaseAnimationForceTimer();
			}
			g_profiler->avg("renderMap(): animated meshes [#]", mesh_animate_count);
			*/
			
			for (int layer = 0; layer < MAX_TILE_LAYERS; layer++) {
				scene::IMesh *mesh = blockMesh->getMesh(layer);
				assert(mesh);

				u32 c = mesh->getMeshBufferCount();
				for (u32 i = 0; i < c; i++) {
					scene::IMeshBuffer *buf = mesh->getMeshBuffer(i);

					video::SMaterial& material = buf->getMaterial();
					// Validate material layer
					if (!m_batches[material]) {
						m_batches[material] = new std::unordered_map<v3s16, MapBatch*>;
					}

					// Validate block batch
					auto &map = *(m_batches[material]);
					if (!map[batch_pos]) {
						map[batch_pos] = new MapBatch(m_camera_offset, SceneManager);
					}
					
					// Set the source mesh for the current block batch at this cell
					if (map[batch_pos]->setSourceMesh(batch_cell, buf, intToFloat(local_pos * MAP_BLOCKSIZE * BS, 1)))
						meshes_set++;
				}
			}
		} // for blocks in sector
	} // foreach sectorblocks
	{
		g_profiler->avg(prefix + "mapblocks updated [#]", meshes_set);
		double batch_time = batch_notify.stop(true);
		g_profiler->avg(prefix + "batch notification [ms]", batch_time);
		g_profiler->avg(prefix + "block gathering [ms]", gather.stop(true) - batch_time);
		g_profiler->avg("CM::updateDrawList", scope.stop(true));
	}
}

struct MeshBufList
{
	video::SMaterial m;
	std::vector<scene::IMeshBuffer*> bufs;
};

struct MeshBufListList
{
	/*!
	 * Stores the mesh buffers of the world.
	 * The array index is the material's layer.
	 * The vector part groups vertices by material.
	 */
	std::vector<MeshBufList> lists[MAX_TILE_LAYERS];

	void clear()
	{
		for (auto &list : lists)
			list.clear();
	}

	void add(scene::IMeshBuffer *buf, u8 layer)
	{
		// Append to the correct layer
		std::vector<MeshBufList> &list = lists[layer];
		const video::SMaterial &m = buf->getMaterial();
		for (MeshBufList &l : list) {
			// comparing a full material is quite expensive so we don't do it if
			// not even first texture is equal
			if (l.m.TextureLayer[0].Texture != m.TextureLayer[0].Texture)
				continue;

			if (l.m == m) {
				l.bufs.push_back(buf);
				return;
			}
		}
		MeshBufList l;
		l.m = m;
		l.bufs.push_back(buf);
		list.push_back(l);
	}
};

void ClientMap::renderMap(video::IVideoDriver* driver, s32 pass)
{
	TimeTaker scope("renderMap scope");
	
	bool is_transparent_pass = pass == scene::ESNRP_TRANSPARENT;
	
	std::string prefix;
	if (pass == scene::ESNRP_SOLID)
		prefix = "renderMap(SOLID): ";
	else
		prefix = "renderMap(TRANSPARENT): ";
	
	/*
	if (!is_transparent_pass){
		// TODO move this out again
		u32 batches_remeshed = 0;
		TimeTaker preprocess("Batch processing");
		// Process batches:
		//  Delete empty
		//  Delete out of range
		//	Remesh
		//  Apply camera offset
		for (auto layer = m_batches.begin(); layer != m_batches.end(); ++layer) {
			for (auto batch = layer->second->begin(); batch != layer->second->end();) {

			}
		}
		g_profiler->avg(prefix + "batches remeshed [#]", batches_remeshed);
		g_profiler->avg(prefix + "batch processing [ms]", preprocess.stop(true));

	}
	*/
	
	
	u32 block_buffers = 0;
	u32 vertex_count = 0;
	u32 batches_drawn = 0;
	u32 culled = 0;
	double drawcall_time = 0;
	TimeTaker draw("Drawing mesh buffers");
	// Push model matrix
	core::matrix4 wmat_old = driver->getTransform(video::ETS_WORLD);
	core::matrix4 m; // Per-block model matrix
	core::matrix4 v = driver->getTransform(video::ETS_VIEW);
	core::matrix4 p = driver->getTransform(video::ETS_PROJECTION);
	core::matrix4 vp = p*v;
	v3f offset = intToFloat(m_camera_offset, BS);
	// Render all layers in order
	for (auto layer = m_batches.begin(); layer != m_batches.end(); ++layer) {
		video::SMaterial mat = layer->first;
		bool transparent = driver->getMaterialRenderer(mat.MaterialType)->isTransparent();
		if ( transparent == is_transparent_pass ) {
			mat.setFlag(video::EMF_TRILINEAR_FILTER, m_cache_trilinear_filter);
			mat.setFlag(video::EMF_BILINEAR_FILTER, m_cache_bilinear_filter);
			mat.setFlag(video::EMF_ANISOTROPIC_FILTER, m_cache_anistropic_filter);
			mat.setFlag(video::EMF_WIREFRAME, m_control.show_wireframe);
			driver->setMaterial(mat);
			for (auto batch = layer->second->begin(); batch != layer->second->end();) {
				v3f bpos = intToFloat(
					batch->first * MAP_BLOCKSIZE * BATCH_STRIDE +
					MAP_BLOCKSIZE * BATCH_STRIDE / 2, BS);
				float range = ( m_control.wanted_range + MAP_BLOCKSIZE * BATCH_STRIDE ) * BS;
				float sqrRange = range * range;
				float sqrDist = (bpos - m_camera_position).getLengthSQ();
				bool empty = batch->second->mesh_count == 0;
				v3f position = intToFloat(batch->first * BATCH_STRIDE * MAP_BLOCKSIZE, BS) - offset;
				if (empty|| (!m_control.range_all && sqrDist > sqrRange)) {
					delete batch->second;
					layer->second->erase(batch++);
				} else {
					m.setTranslation(position);
					// driver->draw3DBox(batch->second->occludee.bounds);
					batch->second->update();					
					if (batch->second->occludee.visible(vp*m,position,m_camera_position - offset,*m_vis_buffer)) {
						driver->setTransform(video::ETS_WORLD, m);
						for (auto mesh : batch->second->mesh_bufs) {
							auto vcount = mesh->getVertexCount();
							g_profiler->avg(prefix + "verts per batch [#]", vcount);
							vertex_count += vcount;
							batches_drawn++;
							block_buffers += batch->second->mesh_count;
							TimeTaker drawCall( "Drawcalls" );
							driver->drawMeshBuffer(mesh);
							drawcall_time += drawCall.stop(true);
						}
					} else {
						culled++;
					}
					++batch;
				}
			}
		}
	}
	// Pop model matrix
	driver->setTransform(video::ETS_WORLD, wmat_old);
	
	if (pass != scene::ESNRP_SOLID) {
		// Uncomment this to see my awful rasterizer at work
		// m_vis_buffer->draw( driver );
		m_vis_buffer->clear();
	}
	
	g_profiler->avg(prefix + "drawcalls [ms]", drawcall_time);
	g_profiler->avg(prefix + "frustum culled [#]", culled);
	g_profiler->avg(prefix + "draw setup [ms]", draw.stop(true) - drawcall_time);
	g_profiler->avg(prefix + "verts total [#]", vertex_count);
	g_profiler->avg(prefix + "batches drawn [#]", batches_drawn);
	g_profiler->avg(prefix + "batches drawn [#]", batches_drawn);
	g_profiler->avg("CM::renderMap", scope.stop(true) * 2);
}

static bool getVisibleBrightness(Map *map, const v3f &p0, v3f dir, float step,
	float step_multiplier, float start_distance, float end_distance,
	const NodeDefManager *ndef, u32 daylight_factor, float sunlight_min_d,
	int *result, bool *sunlight_seen)
{
	int brightness_sum = 0;
	int brightness_count = 0;
	float distance = start_distance;
	dir.normalize();
	v3f pf = p0;
	pf += dir * distance;
	int noncount = 0;
	bool nonlight_seen = false;
	bool allow_allowing_non_sunlight_propagates = false;
	bool allow_non_sunlight_propagates = false;
	// Check content nearly at camera position
	{
		v3s16 p = floatToInt(p0 /*+ dir * 3*BS*/, BS);
		MapNode n = map->getNode(p);
		if(ndef->get(n).param_type == CPT_LIGHT &&
				!ndef->get(n).sunlight_propagates)
			allow_allowing_non_sunlight_propagates = true;
	}
	// If would start at CONTENT_IGNORE, start closer
	{
		v3s16 p = floatToInt(pf, BS);
		MapNode n = map->getNode(p);
		if(n.getContent() == CONTENT_IGNORE){
			float newd = 2*BS;
			pf = p0 + dir * 2*newd;
			distance = newd;
			sunlight_min_d = 0;
		}
	}
	for (int i=0; distance < end_distance; i++) {
		pf += dir * step;
		distance += step;
		step *= step_multiplier;

		v3s16 p = floatToInt(pf, BS);
		MapNode n = map->getNode(p);
		if (allow_allowing_non_sunlight_propagates && i == 0 &&
				ndef->get(n).param_type == CPT_LIGHT &&
				!ndef->get(n).sunlight_propagates) {
			allow_non_sunlight_propagates = true;
		}

		if (ndef->get(n).param_type != CPT_LIGHT ||
				(!ndef->get(n).sunlight_propagates &&
					!allow_non_sunlight_propagates)){
			nonlight_seen = true;
			noncount++;
			if(noncount >= 4)
				break;
			continue;
		}

		if (distance >= sunlight_min_d && !*sunlight_seen && !nonlight_seen)
			if (n.getLight(LIGHTBANK_DAY, ndef) == LIGHT_SUN)
				*sunlight_seen = true;
		noncount = 0;
		brightness_sum += decode_light(n.getLightBlend(daylight_factor, ndef));
		brightness_count++;
	}
	*result = 0;
	if(brightness_count == 0)
		return false;
	*result = brightness_sum / brightness_count;
	/*std::cerr<<"Sampled "<<brightness_count<<" points; result="
			<<(*result)<<std::endl;*/
	return true;
}

int ClientMap::getBackgroundBrightness(float max_d, u32 daylight_factor,
		int oldvalue, bool *sunlight_seen_result)
{
	ScopeProfiler sp(g_profiler, "CM::getBackgroundBrightness", SPT_AVG);
	static v3f z_directions[50] = {
		v3f(-100, 0, 0)
	};
	static f32 z_offsets[sizeof(z_directions)/sizeof(*z_directions)] = {
		-1000,
	};

	if(z_directions[0].X < -99){
		for(u32 i=0; i<sizeof(z_directions)/sizeof(*z_directions); i++){
			// Assumes FOV of 72 and 16/9 aspect ratio
			z_directions[i] = v3f(
				0.02 * myrand_range(-100, 100),
				1.0,
				0.01 * myrand_range(-100, 100)
			).normalize();
			z_offsets[i] = 0.01 * myrand_range(0,100);
		}
	}

	int sunlight_seen_count = 0;
	float sunlight_min_d = max_d*0.8;
	if(sunlight_min_d > 35*BS)
		sunlight_min_d = 35*BS;
	std::vector<int> values;
	for(u32 i=0; i<sizeof(z_directions)/sizeof(*z_directions); i++){
		v3f z_dir = z_directions[i];
		core::CMatrix4<f32> a;
		a.buildRotateFromTo(v3f(0,1,0), z_dir);
		v3f dir = m_camera_direction;
		a.rotateVect(dir);
		int br = 0;
		float step = BS*1.5;
		if(max_d > 35*BS)
			step = max_d / 35 * 1.5;
		float off = step * z_offsets[i];
		bool sunlight_seen_now = false;
		bool ok = getVisibleBrightness(this, m_camera_position, dir,
				step, 1.0, max_d*0.6+off, max_d, m_nodedef, daylight_factor,
				sunlight_min_d,
				&br, &sunlight_seen_now);
		if(sunlight_seen_now)
			sunlight_seen_count++;
		if(!ok)
			continue;
		values.push_back(br);
		// Don't try too much if being in the sun is clear
		if(sunlight_seen_count >= 20)
			break;
	}
	int brightness_sum = 0;
	int brightness_count = 0;
	std::sort(values.begin(), values.end());
	u32 num_values_to_use = values.size();
	if(num_values_to_use >= 10)
		num_values_to_use -= num_values_to_use/2;
	else if(num_values_to_use >= 7)
		num_values_to_use -= num_values_to_use/3;
	u32 first_value_i = (values.size() - num_values_to_use) / 2;

	for (u32 i=first_value_i; i < first_value_i + num_values_to_use; i++) {
		brightness_sum += values[i];
		brightness_count++;
	}

	int ret = 0;
	if(brightness_count == 0){
		MapNode n = getNode(floatToInt(m_camera_position, BS));
		if(m_nodedef->get(n).param_type == CPT_LIGHT){
			ret = decode_light(n.getLightBlend(daylight_factor, m_nodedef));
		} else {
			ret = oldvalue;
		}
	} else {
		ret = brightness_sum / brightness_count;
	}

	*sunlight_seen_result = (sunlight_seen_count > 0);
	return ret;
}

void ClientMap::renderPostFx(CameraMode cam_mode)
{
	// Sadly ISceneManager has no "post effects" render pass, in that case we
	// could just register for that and handle it in renderMap().

	MapNode n = getNode(floatToInt(m_camera_position, BS));

	// - If the player is in a solid node, make everything black.
	// - If the player is in liquid, draw a semi-transparent overlay.
	// - Do not if player is in third person mode
	const ContentFeatures& features = m_nodedef->get(n);
	video::SColor post_effect_color = features.post_effect_color;
	if(features.solidness == 2 && !(g_settings->getBool("noclip") &&
			m_client->checkLocalPrivilege("noclip")) &&
			cam_mode == CAMERA_MODE_FIRST)
	{
		post_effect_color = video::SColor(255, 0, 0, 0);
	}
	if (post_effect_color.getAlpha() != 0)
	{
		// Draw a full-screen rectangle
		video::IVideoDriver* driver = SceneManager->getVideoDriver();
		v2u32 ss = driver->getScreenSize();
		core::rect<s32> rect(0,0, ss.X, ss.Y);
		driver->draw2DRectangle(post_effect_color, rect);
	}
}

// void ClientMap::PrintInfo(std::ostream &out)
// {
	// out<<"ClientMap: ";
// }


