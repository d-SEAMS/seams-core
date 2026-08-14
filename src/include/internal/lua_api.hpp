//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#ifndef SEAMS_LUA_API_H_
#define SEAMS_LUA_API_H_

#include <sol/sol.hpp>

/** @file lua_api.hpp
 *  @brief Lua registration of the d-SEAMS library for the yodaStruct front
 *   end. Each group registers one cohesive slice of the API; registerAll wires
 *   every group into the given state. New-style functions take and return
 *   plain Lua tables; the legacy names bound before this file existed keep
 *   their container-userdata semantics so older scripts keep running.
 */

namespace luaApi {

//! PointCloud usertype, trajectory readers and the file writers
void registerIO(sol::state &lua);

//! Neighbour list construction, by atom ID and by cloud index
void registerNeighbours(sol::state &lua);

//! Primitive ring enumeration and the RingUpdater usertype
void registerRings(sol::state &lua);

//! CHILL, CHILL+, Steinhardt and Voronoi-weighted order parameters
void registerOrder(sol::state &lua);

//! Template overlay, SOAP and Voronoi structure descriptors
void registerDescriptors(sol::state &lua);

//! Topological network criteria, clustering and selection analyses
void registerTopology(sol::state &lua);

//! Every registration group above
void registerAll(sol::state &lua);

} // namespace luaApi

#endif // SEAMS_LUA_API_H_
