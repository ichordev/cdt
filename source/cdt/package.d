/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

module cdt;

public import cdt.triangulation, cdt.utils;

import std.algorithm;

/+
pragma(inline,true) VerticesTriangles calculateTrianglesByVertex(
	const Triangle[] triangles, VertInd verticesSize,
){
	VerticesTriangles vertTris(verticesSize);
	for(TriInd iT(0); iT < triangles.size(); ++iT)
	{
		const VerticesArr3& vv = triangles[iT].vertices;
		for(VerticesArr3::const_iterator v = vv.begin(); v != vv.end(); ++v)
		{
			vertTris[*v].push_back(iT);
		}
	}
	return vertTris;
}
+/

struct DuplicatesInfo{
	size_t[] mapping; ///vertex index mapping
	size_t[] duplicates; ///duplicates' indices
}

DuplicatesInfo findDuplicates(V2)(scope const(V2)[] vertices) nothrow pure @safe{
	size_t[V2] uniqueVerts;
	auto di = DuplicatesInfo(new size_t[](vertices.length));
	size_t iOut = 0;
	foreach(size_t iIn, vert; vertices){
		if(auto ind = vert in uniqueVerts){
			//found a duplicate
			di.mapping[iIn] = *ind;
			di.duplicates ~= iIn;
		}else{
			uniqueVerts[vert] = iOut;
			di.mapping[iIn] = iOut++;
		}
	}
	return di;
}

void removeDuplicates(V2)(scope ref V2[] vertices, scope const size_t[] duplicates) nothrow pure @safe
in(duplicates.isStrictlyMonotonic){
	foreach(duplicate; duplicates)
		vertices = vertices.remove!(SwapStrategy.unstable)(duplicate);
}

void remapEdges(ref Edge[] edges, const size_t[] mapping) nothrow @nogc pure @safe{
	foreach(ref edge; edges)
		edge = Edge(mapping[edge.v1], mapping[edge.v2]);
}

DuplicatesInfo removeDuplicatesAndRemapEdges(V2)(scope ref V2[] vertices, scope ref Edge[] edges) nothrow pure @safe{
	DuplicatesInfo di = findDuplicates!V2(vertices);
	removeDuplicates(vertices, di.duplicates);
	remapEdges(edges, di.mapping);
	return di;
}
/+
template <typename T>
unordered_map<Edge, std::vector<VertInd> > EdgeToSplitVertices(
	const unordered_map<Edge, EdgeVec>& edgeToPieces,
	const std::vector<V2d<T> >& vertices)
{
	typedef std::pair<VertInd, T> VertCoordPair;
	struct ComparePred
	{
		bool operator()(const VertCoordPair& a, const VertCoordPair& b) const
		{
			return a.second < b.second;
		}
	} comparePred;

	unordered_map<Edge, std::vector<VertInd> > edgeToSplitVerts;
	typedef unordered_map<Edge, EdgeVec>::const_iterator It;
	for(It e2pIt = edgeToPieces.begin(); e2pIt != edgeToPieces.end(); ++e2pIt)
	{
		const Edge& e = e2pIt->first;
		const T dX = vertices[e.v2()].x - vertices[e.v1()].x;
		const T dY = vertices[e.v2()].y - vertices[e.v1()].y;
		const bool isX = std::abs(dX) >= std::abs(dY); // X-coord longer
		const bool isAscending =
			isX ? dX >= 0 : dY >= 0; // Longer coordinate ascends
		const EdgeVec& pieces = e2pIt->second;
		std::vector<VertCoordPair> splitVerts;
		// size is:  2[ends] + (pieces - 1)[split vertices] = pieces + 1
		splitVerts.reserve(pieces.size() + 1);
		typedef EdgeVec::const_iterator EIt;
		for(EIt pieceIt = pieces.begin(); pieceIt != pieces.end(); ++pieceIt)
		{
			const array<VertInd, 2> vv = {pieceIt->v1(), pieceIt->v2()};
			typedef array<VertInd, 2>::const_iterator VIt;
			for(VIt v = vv.begin(); v != vv.end(); ++v)
			{
				const T c = isX ? vertices[*v].x : vertices[*v].y;
				splitVerts.push_back(std::make_pair(*v, isAscending ? c : -c));
			}
		}
		// sort by longest coordinate
		std::sort(splitVerts.begin(), splitVerts.end(), comparePred);
		// remove duplicates
		splitVerts.erase(
			std::unique(splitVerts.begin(), splitVerts.end()),
			splitVerts.end());
		assert(splitVerts.size() > 2); // 2 end points with split vertices
		std::pair<Edge, std::vector<VertInd> > val =
			std::make_pair(e, std::vector<VertInd>());
		val.second.reserve(splitVerts.size());
		typedef typename std::vector<VertCoordPair>::const_iterator SEIt;
		for(SEIt it = splitVerts.begin() + 1; it != splitVerts.end() - 1; ++it)
		{
			val.second.push_back(it->first);
		}
		edgeToSplitVerts.insert(val);
	}
	return edgeToSplitVerts;
}
+/
