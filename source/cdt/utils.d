/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

///Utilities and helpers
module cdt.utils;

import cdt.predicates;

import std.algorithm, std.math, std.typecons;

///Index in triangle
alias Index = ubyte;
///Vertex index
alias VertInd = size_t;
///Triangle index
alias TriInd = size_t;

///Constant representing no valid value for index
enum invalidIndex = size_t.max;
///Constant representing no valid neighbour for a triangle
alias noNeighbour = invalidIndex;
///Constant representing no valid vertex for a triangle
alias noVertex = invalidIndex;

alias NeighboursArr3 = TriInd[3]; ///< array of three neighbours

///2D bounding box
struct Rect(T, V2=T[2]){
	V2 min = { V2 p; p[0] = p[1] =  T.max; return p; }(); ///Minimum rectangle coordinate
	V2 max = { V2 p; p[0] = p[1] = -T.max; return p; }(); ///Maximum rectangle coordinate
	
	///Envelop box around a point
	void envelopPoint(V2 p) nothrow @nogc pure @safe{
		this.min[0] = std.algorithm.comparison.min(p[0], this.min[0]);
		this.min[1] = std.algorithm.comparison.min(p[1], this.min[1]);
		this.max[0] = std.algorithm.comparison.max(p[0], this.max[0]);
		this.max[1] = std.algorithm.comparison.max(p[1], this.max[1]);
	}
}

///Bounding box of a collection of custom 2D points given coordinate getters
Rect!(T, V2) envelopBox(T, V2=T[2])(const V2[] list) nothrow @nogc pure @safe{
	Rect!(T, V2) box;
	foreach(item; list){
		box.envelopPoint(item);
	}
	return box;
}

alias Edge = Tuple!(size_t, "v1", size_t, "v2");
alias EdgeUSet = typeof(null)[Edge]; ///Hash table of edges
alias TriIndUSet = typeof(null)[TriInd]; ///Hash table of triangles
alias TriIndUMap = TriInd[TriInd]; ///Triangle hash map

/**
Triangulation triangle (anti-clockwise winding)
```
    v3
    /\
 n3/  \n2
  /____\
v1  n1  v2
```
*/
struct Triangle{
	VertInd[3] vertices; ///Triangle's three vertices
	NeighboursArr3 neighbours; ///Triangle's three neighbours
	
	/**
	Next triangle adjacent to a vertex (clockwise)
	Returns: pair of next triangle and the other vertex of a common edge
	*/
	auto next(VertInd i) const nothrow @nogc pure @safe
	in(vertices[0] == i || vertices[1] == i || vertices[2] == i){
		alias Ret = Tuple!(TriInd, VertInd);
		if(vertices[0] == i)
			return Ret(neighbours[0], vertices[1]);
		if(vertices[1] == i)
			return Ret(neighbours[1], vertices[2]);
		return Ret(neighbours[2], vertices[0]);
	}
	/**
	Previous triangle adjacent to a vertex (counter-clockwise)
	Returns: pair of previous triangle and the other vertex of a common edge
	*/
	auto prev(VertInd i) const nothrow @nogc pure @safe
	in(vertices[0] == i || vertices[1] == i || vertices[2] == i){
		alias Ret = Tuple!(TriInd, VertInd);
		if(vertices[0] == i)
			return Ret(neighbours[2], vertices[2]);
		if(vertices[1] == i)
			return Ret(neighbours[0], vertices[0]);
		return Ret(neighbours[1], vertices[1]);
	}
	
	bool containsVertex(VertInd i) const nothrow @nogc pure @safe =>
		vertices[].canFind(i);
}

///Advance vertex or neighbour index counter-clockwise
Index ccw(Index i) nothrow @nogc pure @safe =>
	Index((i + 1) % 3);
///Advance vertex or neighbour index clockwise
Index cw(Index i) nothrow @nogc pure @safe =>
	Index((i + 2) % 3);

///Location of point on a triangle
enum PtTriLocation{
	inside,
	outside,
	onEdge1,
	onEdge2,
	onEdge3,
	onVertex,
}

///Check if location is classified as on any of three edges
bool isOnEdge(PtTriLocation location) nothrow @nogc pure @safe =>
	location == PtTriLocation.onEdge1 ||
	location == PtTriLocation.onEdge2 ||
	location == PtTriLocation.onEdge3;

///Neighbour index from a on-edge location. Call only if located on the edge!
Index edgeNeighbour(PtTriLocation location) nothrow @nogc pure @safe
in(isOnEdge(location)) =>
	cast(Index)(location - PtTriLocation.onEdge1);

///Relative location of point to a line
enum PtLineLocation{
	left,
	right,
	onLine,
}

///Orient p against line v1-v2 2D: robust geometric predicate
T orient2D(T, V2=T[2])(V2 p, V2 v1, V2 v2) nothrow @nogc pure @safe =>
	orient2DAdaptive!T(v1[0], v1[1], v2[0], v2[1], p[0], p[1]);

///Check if point lies to the left of, to the right of, or on a line
PtLineLocation locatePointLine(T, V2=T[2])(V2 p, V2 v1, V2 v2, T orientationTolerance=0) nothrow @nogc pure @safe =>
	classifyOrientation!T(orient2D!(T, V2)(p, v1, v2), orientationTolerance);

///Classify value of orient2d predicate
PtLineLocation classifyOrientation(T)(T orientation, T orientationTolerance=T(0)) nothrow @nogc pure @safe{
	if(orientation < -orientationTolerance)
		return PtLineLocation.right;
	if(orientation > orientationTolerance)
		return PtLineLocation.left;
	return PtLineLocation.onLine;
}

///Check if point a lies inside of, outside of, or on an edge of a triangle
PtTriLocation locatePointTriangle(T, V2=T[2])(V2 p, V2 v1, V2 v2, V2 v3) nothrow @nogc pure @safe{
	PtTriLocation result = PtTriLocation.inside;
	PtLineLocation edgeCheck = locatePointLine!(T, V2)(p, v1, v2);
	if(edgeCheck == PtLineLocation.right)
		return PtTriLocation.outside;
	if(edgeCheck == PtLineLocation.onLine)
		result = PtTriLocation.onEdge1;
	edgeCheck = locatePointLine!(T, V2)(p, v2, v3);
	if(edgeCheck == PtLineLocation.right)
		return PtTriLocation.outside;
	if(edgeCheck == PtLineLocation.onLine){
		result = (result == PtTriLocation.inside) ?
			PtTriLocation.onEdge2 :
			PtTriLocation.onVertex;
	}
	edgeCheck = locatePointLine!(T, V2)(p, v3, v1);
	if(edgeCheck == PtLineLocation.right)
		return PtTriLocation.outside;
	if(edgeCheck == PtLineLocation.onLine){
		result = result == PtTriLocation.inside ?
			PtTriLocation.onEdge3 :
			PtTriLocation.onVertex;
	}
	return result;
}

pragma(inline,true) nothrow @nogc pure @safe{
	///Opposed neighbour index from vertex index
	Index opoNbr(Index vertIndex){
		switch(vertIndex){
			case 0: return 1;
			case 1: return 2;
			case 2: return 0;
			default: assert(0, "Invalid vertex index");
		}
	}
	
	///Opposed vertex index from neighbour index
	Index opoVrt(Index neighbourIndex){
		switch(neighbourIndex){
			case 0: return 2;
			case 1: return 0;
			case 2: return 1;
			default: assert(0, "Invalid neighbour index");
		}
	}
	
	///Index of triangle's neighbour opposed to a vertex
	Index opposedTriangleInd(ref const VertInd[3] vv, VertInd iVert)
	in(vv[0] == iVert || vv[1] == iVert || vv[2] == iVert){
		if(vv[0] == iVert)
			return 1;
		if(vv[1] == iVert)
			return 2;
		return 0;
	}
	
	///Index of triangle's neighbour opposed to an edge
	Index edgeNeighbourInd(ref const VertInd[3] vv, VertInd iVedge1, VertInd iVedge2)
	in(vv[0] == iVedge1 || vv[1] == iVedge1 || vv[2] == iVedge1)
	in(vv[0] == iVedge2 || vv[1] == iVedge2 || vv[2] == iVedge2)
	in(
		(vv[0] != iVedge1 && vv[0] != iVedge2) ||
		(vv[1] != iVedge1 && vv[1] != iVedge2) ||
		(vv[2] != iVedge1 && vv[2] != iVedge2)
	){
		/*
		*      vv[2]
		*       /\
		*  n[2]/  \n[1]
		*     /____\
		* vv[0] n[0] vv[1]
		*/
		if(vv[0] == iVedge1){
			if(vv[1] == iVedge2)
				return 0;
			return 2;
		}
		if(vv[0] == iVedge2){
			if(vv[1] == iVedge1)
				return 0;
			return 2;
		}
		return 1;
	}
	
	///Index of triangle's vertex opposed to a triangle
	Index opposedVertexInd(ref const NeighboursArr3 nn, TriInd iTopo)
	in(nn[0] == iTopo || nn[1] == iTopo || nn[2] == iTopo){
		if(nn[0] == iTopo)
			return 2;
		if(nn[1] == iTopo)
			return 0;
		return 1;
	}
	
	///If triangle has a given vertex return vertex-index
	Index vertexInd(ref const VertInd[3] vv, VertInd iV)
	in(vv[0] == iV || vv[1] == iV || vv[2] == iV){
		if(vv[0] == iV)
			return 0;
		if(vv[1] == iV)
			return 1;
		return 2;
	}
	
	///Given triangle and a vertex find opposed triangle
	TriInd opposedTriangle(ref const Triangle tri, VertInd iVert) =>
		tri.neighbours[opposedTriangleInd(tri.vertices, iVert)];
	
	///Given two triangles, return vertex of first triangle opposed to the second
	VertInd opposedVertex(ref const Triangle tri, TriInd iTopo) =>
		tri.vertices[opposedVertexInd(tri.neighbours, iTopo)];
	
	///Given triangle and an edge find neighbour sharing the edge
	TriInd edgeNeighbour(ref const Triangle tri, VertInd iVedge1, VertInd iVedge2) =>
		tri.neighbours[edgeNeighbourInd(tri.vertices, iVedge1, iVedge2)];
}

///Test if point lies in a circumscribed circle of a triangle
bool isInCircumCircle(T, V2=T[2])(V2 p, V2 v1, V2 v2, V2 v3) nothrow @nogc pure @safe =>
	inCircleAdaptive(v1[0], v1[1], v2[0], v2[1], v3[0], v3[1], p[0], p[1]) > T(0);

///Test if two vertices share at least one common triangle
bool verticesShareEdge(const TriInd[] aTris, const TriInd[] bTris) nothrow @nogc pure @safe{
	foreach(tri; aTris){
		if(bTris.canFind(tri))
			return true;
	}
	return false;
}

///Distance between two 2D points
T distance(T, V2=T[2])(V2 a, V2 b) nothrow @nogc pure @safe{
	const T dx = b[0] - a[0];
	const T dy = b[1] - a[1];
	return T(sqrt(cast(double)(dx * dx + dy * dy)));
}

///Check if any of triangle's vertices belongs to a super-triangle
bool touchesSuperTriangle(ref const Triangle t) nothrow @nogc pure @safe =>
	t.vertices[0] < 3 || t.vertices[1] < 3 || t.vertices[2] < 3;
