/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

module cdt.triangulation;

import cdt.kdtree, cdt.predicates, cdt.utils;

import std.algorithm, std.conv, std.math, std.meta, std.typecons;

/**
Enum of strategies specifying order in which a range of vertices is inserted

Note: VertexInsertionOrder.randomised will only randomise order of
inserting in triangulation, vertex indices will be preserved as
they were specified in the final triangulation
*/
enum VertexInsertionOrder{
	/**
	Automatic insertion order optimised for better performance
	breadth-first traversal of a Kd-tree for initial bulk-load,
	randomised for subsequent insertions
	*/
	auto_,
	///insert vertices in same order they are provided
	asProvided,
}

enum SuperGeometryType{
	superTriangle, ///conventional super-triangle
	custom, ///user-specified custom geometry (e.g., grid)
}

enum IntersectingConstraintEdges{
	notAllowed, ///constraint edge intersections are not allowed
	tryResolve, ///attempt to resolve constraint edge intersections
	/**
	No checks: slightly faster but less safe.
	User must provide a valid input without intersecting constraints.
	*/
	dontCheck,
}

/**
Type used for storing layer depths for triangles.

Note: LayerDepth should support 60K+ layers, which could be to much or
too little for some use cases. Use version `CDT_LayerDepth_UInt` or
`CDT_LayerDepth_UByte` to change its size.
*/
version(CDT_LayerDepth_UInt)
	alias LayerDepth = uint;
else version(CDT_LayerDepth_UByte)
	alias LayerDepth = ubyte;
else
	alias LayerDepth = ushort;
alias BoundaryOverlapCount = LayerDepth;

///Error thrown when triangulation modification is attempted after it was finalised
class FinalisedException: Exception{
	this(string file=__FILE__, size_t line=__LINE__) nothrow @nogc pure @safe{
		super(
			"Triangulation was finalised with `erase[...]` method. Further modification is not possible.",
			file, line,
		);
	}
}
///Error thrown when duplicate vertex is detected during vertex insertion
class DuplicateVertexException: Exception{
	this(VertInd v1, VertInd v2, string file=__FILE__, size_t line=__LINE__) nothrow pure @safe{
		super(
			text("Duplicate vertex detected: #",v1," is a duplicate of #",v2),
			file, line,
		);
	}
}
/**
Error thrown when intersecting constraint edges are detected, but
triangulation is not configured to attempt to resolve them
*/
class IntersectingConstraintsException: Exception{
	this(ref const Edge e1, ref const Edge e2) nothrow pure @safe{
		super(
			text(
				"Intersecting constraint edges detected: (",
				e1.v1,", ",e1.v2,") intersects (",
				e2.v1,", ",e2.v2,")"
			), file, line,
		);
	}
}

/**
Data structure representing a 2D constrained Delaunay triangulation

Params:
	T = type of vertex coordinates (e.g., `float`, `double`)
	V2 = vector for vertices (`T[2]` by default)
	TNearPointLocator = struct providing locating near point for efficiently
	inserting new points. Must provide methods `VertInd addInd(VertInd, const V2[])` and
	`VertInd nearestInd(V2, const V2[])`
*/
struct Triangulation(T, V2=T[2], TNearPointLocator=LocatorKDTree!(T, V2)){
	V2[] vertices; ///triangulation's vertices
	Triangle[] triangles; ///triangulation's triangles
	EdgeUSet fixedEdges; ///triangulation's constraints (fixed edges)
	/**
	Stores count of overlapping boundaries for a fixed edge. If no entry is
	present for an edge: no boundaries overlap.
	Note: map only has entries for fixed for edges that represent overlapping
	boundaries
	Note: needed for handling depth calculations and hole-removal in case of
	overlapping boundaries
	*/
	BoundaryOverlapCount[Edge] overlapCount;
	/**
	Stores list of original edges represented by a given fixed edge
	Note: map only has entries for edges where multiple original fixed edges
	overlap or where a fixed edge is a part of original edge created by
	conforming Delaunay triangulation vertex insertion
	*/
	Edge[][Edge] pieceToOriginals;
	
	private{
		TriInd[] m_dummyTris;
		TNearPointLocator m_nearPtLocator;
		size_t m_nTargetVerts = 0;
		SuperGeometryType m_superGeomType = SuperGeometryType.superTriangle;
		VertexInsertionOrder m_vertexInsertionOrder = VertexInsertionOrder.auto_;
		IntersectingConstraintEdges m_intersectingEdgesStrategy = IntersectingConstraintEdges.notAllowed;
		T m_minDistToConstraintEdge = T(0);
		TriInd[] m_vertTris; ///one triangle adjacent to each vertex
	}
	
	/**
	Params:
		vertexInsertionOrder = strategy used for ordering vertex insertions
	*/
	this(VertexInsertionOrder vertexInsertionOrder) nothrow @nogc pure @safe{
		m_vertexInsertionOrder = vertexInsertionOrder;
	}
	/**
	Params:
		vertexInsertionOrder = strategy used for ordering vertex insertions
		intersectingEdgesStrategy = strategy for treating intersecting constraint edges
		minDistToConstraintEdge = distance within which point is considered to be lying
			on a constraint edge. Used when adding constraints to the triangulation.
	*/
	this(
		VertexInsertionOrder vertexInsertionOrder,
		IntersectingConstraintEdges intersectingEdgesStrategy,
		T minDistToConstraintEdge,
	) nothrow @nogc pure @safe{
		m_vertexInsertionOrder = vertexInsertionOrder;
		m_intersectingEdgesStrategy = intersectingEdgesStrategy;
		m_minDistToConstraintEdge = minDistToConstraintEdge;
	}
	/**
	Params:
		vertexInsertionOrder = strategy used for ordering vertex insertions
		nearPtLocator = class providing locating near point for efficiently inserting
			new points
		intersectingEdgesStrategy = strategy for treating intersecting constraint edges
		minDistToConstraintEdge = distance within which point is considered to be lying
			on a constraint edge. Used when adding constraints to the triangulation.
	*/
	this(
		VertexInsertionOrder vertexInsertionOrder,
		ref TNearPointLocator nearPtLocator,
		IntersectingConstraintEdges intersectingEdgesStrategy,
		T minDistToConstraintEdge,
	) nothrow @nogc pure @safe{
		m_nearPtLocator = nearPtLocator;
		m_vertexInsertionOrder = vertexInsertionOrder;
		m_intersectingEdgesStrategy = intersectingEdgesStrategy;
		m_minDistToConstraintEdge = minDistToConstraintEdge;
	}
	/**
	Insert vertices into triangulation
	Params:
		vertices = array of vertices to insert
	*/
	void insertVertices(const V2[] vertices) @safe{
		if(isFinalised())
			throw new FinalisedException();
		
		const isFirstTime = this.vertices.length == 0;
		Rect!(T, V2) box;
		
		if(isFirstTime){ //called first time
			box = envelopBox!(T, V2)(vertices);
			addSuperTriangle(box);
		}
		tryIniNearestPointLocator();
		
		const VertInd nExistingVerts = this.vertices.length;
		const VertInd nVerts = nExistingVerts + vertices.length;
		//optimisation, try to pre-allocate tris
		triangles.reserve(triangles.length + 2 * nVerts);
		this.vertices.reserve(nVerts);
		m_vertTris.reserve(nVerts);
		foreach(vertex; vertices)
			addNewVertex(vertex, noNeighbour);
		
		switch(m_vertexInsertionOrder){
			case VertexInsertionOrder.asProvided:
				insertVerticesAsProvided(nExistingVerts);
				break;
			case VertexInsertionOrder.auto_:
				isFirstTime ?
					insertVerticesKDTreeBFS(nExistingVerts, box.min, box.max) :
					insertVerticesRandomised(nExistingVerts);
				break;
			default:
		}
	}
	/**
	Insert constraint edges into triangulation
	Note: Each fixed edge is inserted by deleting the triangles it crosses,
	followed by the triangulation of the polygons on each side of the edge.
	**No new vertices are inserted.**
	Note: If some edge appears more than once in the input this means that
	multiple boundaries overlap at the edge and impacts how hole detection
	algorithm of Triangulation.eraseOuterTrianglesAndHoles works.
	**Make sure there are no erroneous duplicates.**
	Params:
		edges = constraint edges
	 */
	void insertEdges(const Edge[] edges) pure @safe{
		if(isFinalised())
			throw new FinalisedException();
		
		TriangulatePseudoPolygonTask[] tppIterations;
		Edge[] remaining;
		foreach(e; edges){
			//+3 to account for super-triangle vertices
			const edge = Edge(e.v1 + m_nTargetVerts, e.v2 + m_nTargetVerts);
			insertEdge(edge, edge, remaining, tppIterations);
		}
		eraseDummies();
	}
	/**
	Ensure that triangulation conforms to constraints (fixed edges)
	Note: For each fixed edge that is not present in the triangulation its
	midpoint is recursively added until the original edge is represented by
	a sequence of its pieces. **New vertices are inserted.**
	Note: If some edge appears more than once the input this means that
	multiple boundaries overlap at the edge and impacts how hole detection
	algorithm of Triangulation.eraseOuterTrianglesAndHoles works.
	**Make sure there are no erroneous duplicates.**
	Params:
		edges edges to conform to
	*/
	void conformToEdges(const Edge[] edges) pure @safe{
		if(isFinalised())
			throw new FinalisedException();
		
		tryIniNearestPointLocator();
		//state shared between different runs for performance gains
		ConformToEdgeTask[] remaining;
		foreach(e; edges){
			//+3 to account for super-triangle vertices
			const edge = Edge(e.v1 + m_nTargetVerts, e.v2 + m_nTargetVerts);
			conformToEdge(edge, [edge], 0, remaining);
		}
		eraseDummies();
	}
	/**
	Erase triangles adjacent to super triangle
	Note: does nothing if custom geometry is used
	*/
	void eraseSuperTriangle() pure @safe{
		if(m_superGeomType != SuperGeometryType.superTriangle)
			return;
		//find triangles adjacent to super-triangle's vertices
		TriIndUSet toErase;
		foreach(iT, _; triangles){
			if(touchesSuperTriangle(triangles[iT]))
				toErase[iT] = null; //set insertion
		}
		finaliseTriangulation(toErase);
	}
	///Erase triangles outside of constrained boundary using growing
	void eraseOuterTriangles() pure @safe
	in(m_vertTris[0] != noNeighbour){
		//make dummy triangles adjacent to super-triangle's vertices
		const TriInd[] seed = [m_vertTris[0]];
		const TriIndUSet toErase = growToBoundary(seed);
		finaliseTriangulation(toErase);
	}
	/**
	Erase triangles outside of constrained boundary and auto-detected holes
	Note: detecting holes relies on layer peeling based on layer depth
	Note: supports overlapping or touching boundaries
	*/
	void eraseOuterTrianglesAndHoles() pure @safe{
		const LayerDepth[] triDepths = calculateTriangleDepths();
		TriIndUSet toErase;
		//toErase.reserve(triangles.length);
		foreach(iT, _; triangles){
			if(triDepths[iT] % 2 == 0)
				toErase[iT] = null; //insertion into set
		}
	}
	/**
	Call this method after directly setting custom super-geometry via
	vertices and triangles members
	*/
	void iniWithCustomSuperGeometry(){
		m_nearPtLocator = TNearPointLocator(vertices);
		m_nTargetVerts = vertices.length;
		m_superGeomType = SuperGeometryType.custom;
	}
	/**
	Check if the triangulation was finalised with `erase...` method and
	super-triangle was removed.
	Returns: `true` if triangulation is finalised, otherwise `false`
	*/
	bool isFinalised() const =>
		!m_vertTris.length && vertices.length;
	/**
	Calculate depth of each triangle in constraint triangulation.
	Supports overlapping boundaries.
	
	Perform depth peeling from super triangle to outermost boundary,
	then to next boundary and so on until all triangles are traversed.
	For example depth is:
		- 0 for triangles outside outermost boundary
		- 1 for triangles inside boundary but outside hole
		- 2 for triangles in hole
		- 3 for triangles in island and so on...
	Returns: vector where element at index i stores depth of i-th triangle
	*/
	LayerDepth[] calculateTriangleDepths() const pure @safe{
		auto triDepths = {
			auto ret = new LayerDepth[](triangles.length);
			foreach(ref item; ret)
				item = LayerDepth.max;
			return ret;
		}();
		TriInd[] seeds = [m_vertTris[0]]; //stack
		LayerDepth layerDepth = 0;
		LayerDepth deepestSeedDepth = 0;
		
		TriIndUSet[LayerDepth] seedsByDepth;
		do{
			const LayerDepth[TriInd] newSeeds = peelLayer(seeds, layerDepth, triDepths);
			
			seedsByDepth.remove(layerDepth);
			foreach(newSeed; newSeeds.byKeyValue){
				const triInd = newSeed.key;
				const depth = newSeed.value;
				deepestSeedDepth = max(deepestSeedDepth, depth);
				seedsByDepth.require(depth)[triInd] = null; //insertion into set
			}
			const TriIndUSet nextLayerSeeds = seedsByDepth.require(++layerDepth);
			seeds = nextLayerSeeds.keys;
		}while(seeds.length || deepestSeedDepth > layerDepth);
		
		return triDepths;
	}
	
	/**
	Advanced Advanced Triangulation Methods
	Advanced methods for manually modifying the triangulation from
	outside. Only to be used when you know what you are doing.
	*/
	
	/**
	Flip an edge between two triangles.
	Note: Advanced method for manually modifying the triangulation
	from outside. Please only use it if you know what you are doing.
	Params:
		iT = first triangle
		iTopo = second triangle
	```
	*            v4         | - old edge
	*           /|\         ~ - new edge
	*          / | \
	*      n3 /  T' \ n4
	*        /   |   \
	*       /    |    \
	* T -> v1~~~~~~~~~v3 <- Topo
	*       \    |    /
	*        \   |   /
	*      n1 \Topo'/ n2
	*          \ | /
	*           \|/
	*            v2
	```
	*/
	void flipEdge(TriInd iT, TriInd iTopo) nothrow pure @safe{
		Triangle* t = &triangles[iT];
		Triangle* tOpo = &triangles[iTopo];
		const TriInd[3] triNs = t.neighbours;
		const TriInd[3] triOpoNs = tOpo.neighbours;
		const VertInd[3] triVs = t.vertices;
		const VertInd[3] triOpoVs = tOpo.vertices;
		//find vertices and neighbours
		Index i = opposedVertexInd(t.neighbours, iTopo);
		const VertInd v1 = triVs[i];
		const VertInd v2 = triVs[ccw(i)];
		const TriInd n1 = triNs[i];
		const TriInd n3 = triNs[cw(i)];
		i = opposedVertexInd(tOpo.neighbours, iT);
		const VertInd v3 = triOpoVs[i];
		const VertInd v4 = triOpoVs[ccw(i)];
		const TriInd n4 = triOpoNs[i];
		const TriInd n2 = triOpoNs[cw(i)];
		//change vertices and neighbours
		*t    = Triangle([v4, v1, v3], [n3, iTopo, n4]);
		*tOpo = Triangle([v2, v3, v1], [n2, iT,    n1]);
		//adjust neighbouring triangles and vertices
		changeNeighbour(n1, iT, iTopo);
		changeNeighbour(n4, iTopo, iT);
		//only adjust adjacent triangles if triangulation is not finalised:
		//can happen when called from outside on an already finalised triangulation
		if(!isFinalised()){
			setAdjacentTriangle(v4, iT);
			setAdjacentTriangle(v2, iTopo);
		}
	}
	void flipEdge(
		TriInd iT, TriInd iTopo,
		VertInd v1, VertInd v2, VertInd v3, VertInd v4,
		TriInd n1, TriInd n2, TriInd n3, TriInd n4,
	) nothrow pure @safe{
		//change vertices and neighbours
		triangles[iT] = Triangle([v4, v1, v3], [n3, iTopo, n4]);
		triangles[iTopo] = Triangle([v2, v3, v1], [n2, iT, n1]);
		//adjust neighbouring triangles and vertices
		changeNeighbour(n1, iT, iTopo);
		changeNeighbour(n4, iTopo, iT);
		//only adjust adjacent triangles if triangulation is not finalised:
		//can happen when called from outside on an already finalised triangulation
		if(!isFinalised()){
			setAdjacentTriangle(v4, iT);
			setAdjacentTriangle(v2, iTopo);
		}
	}
	/**
	Remove triangles with specified indices.
	Adjust internal triangulation state accordingly.
	Params:
		removedTriangles = indices of triangles to remove
	*/
	void removeTriangles(const TriIndUSet removedTriangles) pure @safe{
		if(!removedTriangles.length)
			return;
		//remove triangles and calculate triangle index mapping
		TriIndUMap triIndMap;
		TriInd iTNew = 0;
		foreach(iT; 0..triangles.length){
			if(iT in removedTriangles)
				continue;
			triIndMap[iT] = iTNew;
			triangles[iTNew] = triangles[iT];
			iTNew++;
		}
		triangles = triangles[0..$-removedTriangles.length];
		//adjust triangles' neighbours
		foreach(iT; 0..triangles.length){
			Triangle* t = &triangles[iT];
			//update neighbours to account for removed triangles
			NeighboursArr3* nn = &t.neighbours;
			foreach(ref n; *nn){
				if(n in removedTriangles){
					n = noNeighbour;
				}else if(n != noNeighbour){
					n = triIndMap.require(n);
				}
			}
		}
	}
	///Access internal vertex adjacent triangles
	ref TriInd[] vertTrisInternal() nothrow @nogc pure @safe =>
		m_vertTris;
	
	private:
	
	void addSuperTriangle(Rect!(T, V2) box) nothrow pure @safe{
		m_nTargetVerts = 3;
		m_superGeomType = SuperGeometryType.superTriangle;
		
		const centre = {
			V2 p;
			p[0] = (box.min[0] + box.max[0]) / T(2);
			p[1] = (box.min[1] + box.max[1]) / T(2);
			return p;
		}();
		const T w = box.max[0] - box.min[0];
		const T h = box.max[1] - box.min[1];
		T ir = max(w, h); //incircle radius upper bound
		
		//Note: make sure radius is big enough. Constants chosen experimentally.
		//- for tiny bounding boxes: use 1.0 as the smallest radius
		//- multiply radius by 2.0 for extra safety margin
		ir = max(T(2) * ir, T(1));
		
		//Note: for very large floating point numbers rounding can lead to wrong
		//super-triangle coordinates. This is a very rare corner-case so the
		//handling is very primitive.
		{ //note: '<=' means '==' but avoids the warning
			while(centre.y <= centre.y - ir)
				ir *= T(2);
		}
		const T er = T(2) * ir; // excircle radius
		const T cos30Deg = T(0.8660254037844386); // note: (std::sqrt(3.0) / 2.0)
		const T shiftX = er * cos30Deg;
		const V2 posV1 = { V2 p = centre; p[0] -= shiftX; p[1] -= ir; return p; }();
		const V2 posV2 = { V2 p = centre; p[0] += shiftX; p[1] -= ir; return p; }();
		const V2 posV3 = { V2 p = centre;                 p[1] += er; return p; }();
		addNewVertex(posV1, 0);
		addNewVertex(posV2, 0);
		addNewVertex(posV3, 0);
		const superTri = Triangle(
			[0, 1, 2],
			[noNeighbour, noNeighbour, noNeighbour],
		);
		addTriangle(superTri);
		if(m_vertexInsertionOrder != VertexInsertionOrder.auto_){
			m_nearPtLocator = TNearPointLocator(vertices);
		}
	}
	void addNewVertex(V2 pos, TriInd iT) nothrow pure @safe{
		vertices ~= pos;
		m_vertTris ~= iT;
	}
	void insertVertex(VertInd iVert) pure @safe{
		const V2 v = vertices[iVert];
		const VertInd walkStart = m_nearPtLocator.nearestInd(v, vertices);
		insertVertex(iVert, walkStart);
		tryAddVertexToLocator(iVert);
	}
	void insertVertex(VertInd iVert, VertInd walkStart) pure @safe{
		const TriInd[2] trisAt = walkingSearchTrianglesAt(iVert, walkStart);
		TriInd[] triStack = //a stack
			trisAt[1] == noNeighbour ?
				insertVertexInsideTriangle(iVert, trisAt[0]) :
				insertVertexOnEdge(iVert, trisAt[0], trisAt[1]);
		ensureDelaunayByEdgeFlips(iVert, triStack);
	}
	void ensureDelaunayByEdgeFlips(VertInd iV1, ref TriInd[] triStack){
		TriInd iTopo, n1, n2, n3, n4;
		VertInd iV2, iV3, iV4;
		while(triStack.length){
			const TriInd iT = triStack[$-1]; //.top
			triStack = triStack[0..$-1]; //.pop()
			
			edgeFlipInfo(iT, iV1, iTopo, iV2, iV3, iV4, n1, n2, n3, n4);
			if(iTopo != noNeighbour && isFlipNeeded(iV1, iV2, iV3, iV4)){
				flipEdge(iT, iTopo, iV1, iV2, iV3, iV4, n1, n2, n3, n4);
				triStack ~= [iT, iTopo];
			}
		}
	}
	///Flip fixed edges and return a list of flipped fixed edges
	Edge[] insertVertexFlipFixedEdges(VertInd iV1) pure @safe{
		Edge[] flippedFixedEdges;
		
		const V2 v1 = vertices[iV1];
		const VertInd startVertex = m_nearPtLocator.nearestInd(v1, vertices);
		TriInd[2] trisAt = walkingSearchTrianglesAt(iV1, startVertex);
		TriInd[] triStack = trisAt[1] == noNeighbour ?
			insertVertexInsideTriangle(iV1, trisAt[0]) :
			insertVertexOnEdge(iV1, trisAt[0], trisAt[1]);
		
		TriInd iTopo, n1, n2, n3, n4;
		VertInd iV2, iV3, iV4;
		while(triStack.length){
			const TriInd iT = triStack[$-1]; //.top
			triStack = triStack[0..$-1]; //.pop()
			
			edgeFlipInfo(iT, iV1, iTopo, iV2, iV3, iV4, n1, n2, n3, n4);
			if(iTopo != noNeighbour && isFlipNeeded(iV1, iV2, iV3, iV4)){
				//if flipped edge is fixed, remember it
				const flippedEdge = Edge(iV2, iV4);
				if(fixedEdges.length && flippedEdge in fixedEdges){
					flippedFixedEdges ~= flippedEdge;
				}
				flipEdge(iT, iTopo, iV1, iV2, iV3, iV4, n1, n2, n3, n4);
				triStack ~= [iT, iTopo];
			}
		}
		tryAddVertexToLocator(iV1);
		return flippedFixedEdges;
	}
	///State for an iteration of triangulate pseudo-polygon
	alias TriangulatePseudoPolygonTask = Tuple!(size_t, size_t, TriInd, TriInd, Index);
	/**
	Insert an edge into constraint Delaunay triangulation
	Params:
		edge = edge to insert
		originalEdge = original edge inserted edge is part of
		remaining = parts of the edge that still need to be inserted
		tppIterations = stack to be used for storing iterations of
			triangulating pseudo-polygon
	Note: `remaining` and `tppIterations` are shared between different
	runs for performance gains (reducing memory allocations)
	*/
	void insertEdge(
		Edge edge, Edge originalEdge, ref Edge[] remaining,
		ref TriangulatePseudoPolygonTask[] tppIterations,
	) pure @safe{
		//use iteration over recursion to avoid stack overflows
		remaining = [edge];
		while(remaining.length){
			edge = remaining[$-1];
			remaining = remaining[0..$-1];
			insertEdgeIteration(edge, originalEdge, remaining, tppIterations);
		}
	}
	/**
	Insert an edge or its part into constraint Delaunay triangulation
	Params:
		edge = edge to insert
		originalEdge = original edge inserted edge is part of
		remainingStack = parts of the edge that still need to be inserted
		tppIterations = stack to be used for storing iterations of
			triangulating pseudo-polygon
	Note: `remainingStack` and `tppIterations` are shared between different
	runs for performance gains (reducing memory allocations)
	 */
	void insertEdgeIteration(
		Edge edge, Edge originalEdge, ref Edge[] remaining,
		ref TriangulatePseudoPolygonTask[] tppIterations,
	) pure @safe{
		const iA = edge.v1;
		auto  iB = edge.v2;
		if(iA == iB) //edge connects a vertex to itself
			return;
		if(hasEdge(iA, iB)){
			fixEdge(edge, originalEdge);
			return;
		}
		const V2 a = vertices[iA];
		const V2 b = vertices[iB];
		const T distanceTolerance = m_minDistToConstraintEdge == T(0) ?
			T(0) :
			m_minDistToConstraintEdge * distance!(T, V2)(a, b);
		
		TriInd iT;
		//Note: 'L' is left and 'R' is right of the inserted constraint edge
		VertInd iVL, iVR;
		AliasSeq!(iT, iVL, iVR) = intersectedTriangle(iA, a, b, distanceTolerance).expand;
		//if one of the triangle vertices is on the edge, move edge start
		if(iT == noNeighbour){
			const edgePart = Edge(iA, iVL);
			fixEdge(edgePart, originalEdge);
			remaining ~= Edge(iVL, iB);
			return;
		}
		Triangle t = triangles[iT];
		TriInd[] intersected = [iT];
		VertInd[] polyL = [iA, iVL], polyR = [iA, iVR];
		TriInd[Edge] outerTris;
		outerTris[Edge(iA, iVL)] = edgeNeighbour(t, iA, iVL);
		outerTris[Edge(iA, iVR)] = edgeNeighbour(t, iA, iVR);
		VertInd iV = iA;
		
		while(!t.containsVertex(iB)){
			const TriInd iTopo = opposedTriangle(t, iV);
			const Triangle tOpo = triangles[iTopo];
			const VertInd iVopo = opposedVertex(tOpo, iT);
			
			final switch(m_intersectingEdgesStrategy){
				case IntersectingConstraintEdges.notAllowed:
					if(Edge(iVL, iVR) in fixedEdges){
						//make sure to report original input edges in the exception
						Edge e1 = originalEdge;
						Edge e2 = Edge(iVL, iVR);
						if(auto edges = e2 in pieceToOriginals){
							e2 = (*edges)[0];
						}
						//don't count super-triangle vertices
						e1 = Edge(e1.v1 - m_nTargetVerts, e1.v2 - m_nTargetVerts);
						e2 = Edge(e2.v1 - m_nTargetVerts, e2.v2 - m_nTargetVerts);
						if(auto edges = e2 in pieceToOriginals){
							e2 = (*edges)[0];
						}
						throw new IntersectingConstraintsException(e1, e2);
					}
					break;
				case IntersectingConstraintEdges.tryResolve:
					if(Edge(iVL, iVR) !in fixedEdges)
						break;
					//split edge at the intersection of two constraint edges
					const V2 newV = intersectionPosition!(T, V2)(
						vertices[iA],  vertices[iB],
						vertices[iVL], vertices[iVR],
					);
					const VertInd iNewVert = splitFixedEdgeAt(Edge(iVL, iVR), newV, iT, iTopo);
					//TODO: is it's possible to re-use pseudo-polygons for inserting [iA, iNewVert] edge half?
					remaining ~= [Edge(iA, iNewVert), Edge(iNewVert, iB)];
					return;
				case IntersectingConstraintEdges.dontCheck:
					assert(Edge(iVL, iVR) !in fixedEdges);
					break;
			}
			
			switch(locatePointLine!(T, V2)(vertices[iVopo], a, b, distanceTolerance)){
				case PtLineLocation.left:
					const e = Edge(polyL[$-1], iVopo);
					outerTris.update(
						key: e,
						create: () => edgeNeighbour(tOpo, e.v1, e.v2),
						update: (TriInd _) => noNeighbour, //hanging edge detected
					);
					polyL ~= iVopo;
					iV = iVL;
					iVL = iVopo;
					break;
				case PtLineLocation.right:
					const e = Edge(polyR[$-1], iVopo);
					outerTris.update(
						key: e,
						create: () => edgeNeighbour(tOpo, e.v1, e.v2),
						update: (TriInd _) => noNeighbour, //hanging edge detected
					);
					polyR ~= iVopo;
					iV = iVR;
					iVR = iVopo;
					break;
				default: //encountered point on the edge
					iB = iVopo;
			}
			
			intersected ~= iTopo;
			iT = iTopo;
			t = triangles[iT];
		}
		outerTris[Edge(polyL[$-1], iB)] = edgeNeighbour(t, polyL[$-1], iB);
		outerTris[Edge(polyR[$-1], iB)] = edgeNeighbour(t, polyR[$-1], iB);
		polyL ~= iB;
		polyR ~= iB;
		
		assert(intersected.length);
		//make sure start/end vertices have a valid adjacent triangle
		//that is not intersected by an edge
		if(m_vertTris[iA] == intersected[0])
			pivotVertexTriangleCW(iA);
		if(m_vertTris[iB] == intersected[$-1])
			pivotVertexTriangleCW(iB);
		
		//remove intersected triangles
		foreach(inter; intersected)
			makeDummy(inter);
		
		{ //triangulate pseudo-polygons on both sides
			reverse(polyR);
			const TriInd iTL = addTriangle();
			const TriInd iTR = addTriangle();
			triangulatePseudoPolygon(polyL, outerTris, iTL, iTR, tppIterations);
			triangulatePseudoPolygon(polyR, outerTris, iTR, iTL, tppIterations);
		}
		
		if(iB != edge.v2){ //encountered point on the edge
			//fix edge part
			const edgePart = Edge(iA, iB);
			fixEdge(edgePart, originalEdge);
			remaining ~= Edge(iB, edge.v2);
			return;
		}else{
			fixEdge(edge, originalEdge);
		}
	}
	///State for iteration of conforming to edge
	alias ConformToEdgeTask = Tuple!(Edge, const Edge[], BoundaryOverlapCount);
	/**
	Conform Delaunay triangulation to a fixed edge by recursively inserting
	mid point of the edge and then conforming to its halves
	Params:
		edge = fixed edge to conform to
		originals = original edges that new edge is piece of
		overlaps = count of overlapping boundaries at the edge.
			Only used when re-introducing edge with overlaps > 0
		remaining = remaining edge parts to be conformed to
	Note: state of `remaining` is shared between different
	runs for performance gains (reducing memory allocations)
	 */
	void conformToEdge(
		Edge edge, const(Edge)[] originals,
		BoundaryOverlapCount overlaps,
		ref ConformToEdgeTask[] remaining,
	) pure @safe{
		//use iteration over recursion to avoid stack overflows
		remaining = [ConformToEdgeTask(edge, originals, overlaps)];
		while(remaining.length){
			AliasSeq!(edge, originals, overlaps) = remaining[$-1].expand;
			conformToEdgeIteration(edge, originals, overlaps, remaining);
		}
	}
	/**
	Iteration of conform to fixed edge.
	Params:
		edge = fixed edge to conform to
		originals = original edges that new edge is piece of
		overlaps = count of overlapping boundaries at the edge.
			Only used when re-introducing edge with overlaps > 0
		remaining = remaining edge parts
	Note: state of `remaining` is shared between different
	runs for performance gains (reducing memory allocations)
	*/
	void conformToEdgeIteration(
		Edge edge, const Edge[] originals,
		BoundaryOverlapCount overlaps,
		ref ConformToEdgeTask[] remaining,
	) pure @safe{
		const iA = edge.v1;
		auto  iB = edge.v2;
		if(iA == iB) //edge connects a vertex to itself
			return;
		
		if(hasEdge(iA, iB)){
			fixEdge(edge);
			if(overlaps > 0)
				overlapCount[edge] = overlaps;
			// avoid marking edge as a part of itself
			if(originals.length && edge != originals[0]){
				pieceToOriginals.require(edge).insertUnique(originals);
			}
			return;
		}
		
		const a = vertices[iA];
		const b = vertices[iB];
		const T distanceTolerance = m_minDistToConstraintEdge == T(0) ?
			T(0) :
			m_minDistToConstraintEdge * distance!(T, V2)(a, b);
		TriInd iT;
		VertInd iVLeft, iVRight;
		AliasSeq!(iT, iVLeft, iVRight) = intersectedTriangle(iA, a, b, distanceTolerance).expand;
		//if one of the triangle vertices is on the edge, move edge start
		if(iT == noNeighbour){
			const edgePart = Edge(iA, iVLeft);
			fixEdge(edgePart);
			if(overlaps > 0)
				overlapCount[edgePart] = overlaps;
			pieceToOriginals.require(edgePart).insertUnique(originals);
			remaining ~= ConformToEdgeTask(Edge(iVLeft, iB), originals, overlaps);
			return;
		}
		
		VertInd iV = iA;
		Triangle t = triangles[iT];
		while(!t.vertices[].canFind(iB)){
			const TriInd iTopo = opposedTriangle(t, iV);
			const Triangle tOpo = triangles[iTopo];
			const VertInd iVopo = opposedVertex(tOpo, iT);
			const V2 vOpo = vertices[iVopo];
			
			final switch(m_intersectingEdgesStrategy){
				case IntersectingConstraintEdges.notAllowed:
					if(Edge(iVLeft, iVRight) !in fixedEdges)
						break;
					//make sure to report original input edges in the exception
					Edge e1;
					if(auto edges = edge in pieceToOriginals){
						e1 = (*edges)[0];
					}else{
						e1 = edge;
					}
					auto e2 = Edge(iVLeft, iVRight);
					if(auto edges = e2 in pieceToOriginals){
						e2 = (*edges)[0];
					}
					//don't count super-triangle vertices
					e1 = Edge(e1.v1 - m_nTargetVerts, e1.v2 - m_nTargetVerts);
					e2 = Edge(e2.v1 - m_nTargetVerts, e2.v2 - m_nTargetVerts);
					throw new IntersectingConstraintsException(e1, e2);
				case IntersectingConstraintEdges.tryResolve:
					if(Edge(iVLeft, iVRight) !in fixedEdges)
						break;
					//split edge at the intersection of two constraint edges
					const V2 newV = intersectionPosition!(T, V2)(
						vertices[iA],     vertices[iB],
						vertices[iVLeft], vertices[iVRight],
					);
					const VertInd iNewVert = splitFixedEdgeAt(Edge(iVLeft, iVRight), newV, iT, iTopo);
					remaining ~= ConformToEdgeTask(Edge(iNewVert, iB), originals, overlaps);
					remaining ~= ConformToEdgeTask(Edge(iA, iNewVert), originals, overlaps);
					return;
				case IntersectingConstraintEdges.dontCheck:
					assert(Edge(iVLeft, iVRight) !in fixedEdges);
					break;
			}
			iT = iTopo;
			t = triangles[iT];
			
			switch(locatePointLine!(T, V2)(vOpo, a, b, distanceTolerance)){
				case PtLineLocation.left:
					iV = iVLeft;
					iVLeft = iVopo;
					break;
				case PtLineLocation.right:
					iV = iVRight;
					iVRight = iVopo;
					break;
				default: //encountered point on the edge
					iB = iVopo;
			}
		}
		
		//encountered one or more points on the edge: add remaining edge part
		if(iB != edge.v2){
			remaining ~= ConformToEdgeTask(Edge(iB, edge.v2), originals, overlaps);
		}
		
		//add mid-point to triangulation
		const VertInd iMid = vertices.length;
		const V2 start = vertices[iA];
		const V2 end = vertices[iB];
		{
			V2 newVert;
			newVert[0] = (start[0] + end[0]) / T(2);
			newVert[1] = (start[1] + end[1]) / T(2);
			addNewVertex(newVert, noNeighbour);
		}
		const Edge[] flippedFixedEdges = insertVertexFlipFixedEdges(iMid);
		
		remaining ~= ConformToEdgeTask(Edge(iMid, iB), originals, overlaps);
		remaining ~= ConformToEdgeTask(Edge(iA, iMid), originals, overlaps);
		
		//re-introduce fixed edges that were flipped and make sure overlap count is preserved
		foreach(flippedFixedEdge; flippedFixedEdges){
			fixedEdges.remove(flippedFixedEdge);
			
			BoundaryOverlapCount prevOverlaps = 0;
			if(auto olap = flippedFixedEdge in overlapCount){
				prevOverlaps = *olap;
				overlapCount.remove(flippedFixedEdge);
			}
			//override overlapping boundaries count when re-inserting an edge
			Edge[] prevOriginals;
			if(auto orig = flippedFixedEdge in pieceToOriginals){
				prevOriginals = *orig;
			}else{
				prevOriginals = [flippedFixedEdge];
			}
			remaining ~= ConformToEdgeTask(flippedFixedEdge, prevOriginals, prevOverlaps);
		}
	}
	/**
	Returns:
		- intersected triangle index
		- index of point on the left of the line
		- index of point on the right of the line
	If left point is right on the line: no triangle is intersected:
		- triangle index is no-neighbour (invalid)
		- index of point on the line
		- index of point on the right of the line
	*/
	auto intersectedTriangle(
		VertInd iA, V2 a, V2 b, T orientationTolerance=T(0),
	) const pure @safe{
		alias Ret = Tuple!(TriInd, VertInd, VertInd);
		const TriInd startTri = m_vertTris[iA];
		TriInd iT = startTri;
		do{
			const Triangle t = triangles[iT];
			const Index i = vertexInd(t.vertices, iA);
			const VertInd iP2 = t.vertices[ccw(i)];
			const T orientP2 = orient2D!(T, V2)(vertices[iP2], a, b);
			const PtLineLocation locP2 = classifyOrientation(orientP2);
			if(locP2 == PtLineLocation.right){
				const VertInd iP1 = t.vertices[cw(i)];
				const T orientP1 = orient2D!(T, V2)(vertices[iP1], a, b);
				const PtLineLocation locP1 = classifyOrientation(orientP1);
				if(locP1 == PtLineLocation.onLine){
					return Ret(noNeighbour, iP1, iP1);
				}
				if(locP1 == PtLineLocation.left){
					if(orientationTolerance){
						T closestOrient;
						VertInd iClosestP;
						if(abs(orientP1) <= abs(orientP2)){
							closestOrient = orientP1;
							iClosestP = iP1;
						}else{
							closestOrient = orientP2;
							iClosestP = iP2;
						}
						if(classifyOrientation(closestOrient, orientationTolerance) == PtLineLocation.onLine){
							return Ret(noNeighbour, iClosestP, iClosestP);
						}
					}
					return Ret(iT, iP1, iP2);
				}
			}
			iT = t.next(iA)[0];
		}while(iT != startTri);
		
		throw new Exception("Could not find vertex triangle intersected by an edge.");
	}
	/**
	Insert point into triangle: split into 3 triangles:
		- create 2 new triangles
		- re-use old triangle for the 3rd
	Returns: indices of three resulting triangles
	```
	*             v3
	*           / | \
	*          /  |  \ <-- original triangle (t)
	*         /   |   \
	*     n3 /    |    \ n2
	*       /newT2|newT1\
	*      /      v      \
	*     /    __/ \__    \
	*    /  __/       \__  \
	*   / _/      t'     \_ \
	* v1 ___________________ v2
	*            n1
	```
	*/
	TriInd[] insertVertexInsideTriangle(VertInd v, TriInd iT) nothrow pure @safe{
		const TriInd iNewT1 = addTriangle();
		const TriInd iNewT2 = addTriangle();
		
		Triangle* t = &triangles[iT];
		const VertInd[3] vv = t.vertices;
		const TriInd[3] nn = t.neighbours;
		const VertInd v1 = vv[0], v2 = vv[1], v3 = vv[2];
		const TriInd n1 = nn[0], n2 = nn[1], n3 = nn[2];
		//make two new triangles and convert current triangle to 3rd new triangle
		triangles[iNewT1] = Triangle([v2, v3, v], [n2, iNewT2, iT]);
		triangles[iNewT2] = Triangle([v3, v1, v], [n3, iT, iNewT1]);
		*t = Triangle([v1, v2, v], [n1, iNewT1, iNewT2]);
		//adjust adjacent triangles
		setAdjacentTriangle(v, iT);
		setAdjacentTriangle(v3, iNewT1);
		//change triangle neighbour's neighbours to new triangles
		changeNeighbour(n2, iT, iNewT1);
		changeNeighbour(n3, iT, iNewT2);
		//return newly added triangles
		return [iT, iNewT1, iNewT2];
	}
	/**
	Inserting a point on the edge between two triangles
	```
	*  T1 (top)        v1
	*                 /|\
	*            n1 /  |  \ n4
	*             /    |    \
	*           /  T1' | Tnew1\
	*         v2-------v-------v4
	*           \  T2' | Tnew2/
	*             \    |    /
	*            n2 \  |  / n3
	*                 \|/
	* T2 (bottom)      v3
	```
	Returns: indices of four resulting triangles
	*/
	TriInd[] insertVertexOnEdge(VertInd v, TriInd iT1, TriInd iT2) nothrow pure @safe{
		const TriInd iTNew1 = addTriangle();
		const TriInd iTNew2 = addTriangle();
		
		Triangle* t1 = &triangles[iT1];
		Triangle* t2 = &triangles[iT2];
		Index i = opposedVertexInd(t1.neighbours, iT2);
		const VertInd v1 = t1.vertices[i];
		const VertInd v2 = t1.vertices[ccw(i)];
		const TriInd n1 = t1.neighbours[i];
		const TriInd n4 = t1.neighbours[cw(i)];
		i = opposedVertexInd(t2.neighbours, iT1);
		const VertInd v3 = t2.vertices[i];
		const VertInd v4 = t2.vertices[ccw(i)];
		const TriInd n3 = t2.neighbours[i];
		const TriInd n2 = t2.neighbours[cw(i)];
		//add new triangles and change existing ones
		*t1 = Triangle([v, v1, v2], [iTNew1, n1, iT2]);
		*t2 = Triangle([v, v2, v3], [iT1, n2, iTNew2]);
		triangles[iTNew1] = Triangle([v, v4, v1], [iTNew2, n4, iT1]);
		triangles[iTNew2] = Triangle([v, v3, v4], [iT2, n3, iTNew1]);
		//adjust adjacent triangles
		setAdjacentTriangle(v, iT1);
		setAdjacentTriangle(v4, iTNew1);
		//adjust neighbouring triangles and vertices
		changeNeighbour(n4, iT1, iTNew1);
		changeNeighbour(n3, iT2, iTNew2);
		//return newly added triangles
		return [iT1, iTNew2, iT2, iTNew1];
	}
	TriInd[2] trianglesAt(V2 pos) const pure @safe{
		TriInd[2] output = [noNeighbour, noNeighbour];
		foreach(i; 0..triangles.length){
			const Triangle t = triangles[i];
			const PtTriLocation loc = locatePointTriangle!(T, V2)(
				pos,
				vertices[t.vertices[0]],
				vertices[t.vertices[1]],
				vertices[t.vertices[2]],
			);
			if(loc == PtTriLocation.outside)
				continue;
			output[0] = i;
			if(isOnEdge(loc))
				output[1] = t.neighbours[edgeNeighbour(loc)];
			return output;
		}
		throw new Exception("No triangle was found at position");
	}
	TriInd[2] walkingSearchTrianglesAt(VertInd iV, VertInd startVertex) const pure @safe{
		const V2 v = vertices[iV];
		TriInd[2] output = [noNeighbour, noNeighbour];
		const TriInd iT = walkTriangles(startVertex, v);
		
		//Finished walk, locate point in current triangle
		const Triangle t = triangles[iT];
		const V2 v1 = vertices[t.vertices[0]];
		const V2 v2 = vertices[t.vertices[1]];
		const V2 v3 = vertices[t.vertices[2]];
		const PtTriLocation loc = locatePointTriangle!(T, V2)(v, v1, v2, v3);
		switch(loc){
			case PtTriLocation.outside:
				throw new Exception("No triangle was found at position");
			case PtTriLocation.onVertex:
				const VertInd iDupe = v1 == v ?
					t.vertices[0] :   v2 == v ?
						t.vertices[1] :
							t.vertices[2];
				throw new DuplicateVertexException(iV - m_nTargetVerts, iDupe - m_nTargetVerts);
			default:
		}
		
		output[0] = iT;
		if(isOnEdge(loc))
			output[1] = t.neighbours[edgeNeighbour(loc)];
		return output;
	}
	TriInd walkTriangles(VertInd startVertex, V2 pos) const nothrow pure @safe{
		//begin walk in search of triangle at pos
		TriInd currTri = m_vertTris[startVertex];
		bool found = false;
		SplitMix64RandGen rng;
		while(!found){
			found = true;
			const Triangle t = triangles[currTri];
			// stochastic offset to randomize which edge we check first
			const offset = cast(Index)(rng() % 3);
			for(Index ii=0; ii < 3; ii++){
				const i = cast(Index)((ii + offset) % 3);
				const V2 vStart = vertices[t.vertices[i]];
				const V2 vEnd = vertices[t.vertices[ccw(i)]];
				const TriInd iN = t.neighbours[i];
				if(
					locatePointLine!(T, V2)(pos, vStart, vEnd) == PtLineLocation.right &&
					iN != noNeighbour
				){
					found = false;
					currTri = iN;
					break;
				}
			}
		}
		return currTri;
	}
	/**
	Given triangle and its vertex find opposite triangle and the
	other three vertices and surrounding neighbours
	```
	*                       v4         original edge: (v1, v3)
	*                      /|\   flip-candidate edge: (v,  v2)
	*                    /  |  \
	*              n3  /    |    \  n4
	*                /      |      \
	* new vertex--> v1    T | Topo  v3
	*                \      |      /
	*              n1  \    |    /  n2
	*                    \  |  /
	*                      \|/
	*                       v2
	```
	*/
	void edgeFlipInfo(
		TriInd iT, VertInd iV1,
		ref TriInd iTopo, ref VertInd iV2,
		ref VertInd iV3, ref VertInd iV4,
		ref TriInd n1, ref TriInd n2,
		ref TriInd n3, ref TriInd n4,
	) const nothrow @nogc pure @safe{
		/*
		*       v[2]
		*        / \
		*   n[2]/   \n[1]
		*      /_____\
		* v[0]  n[0]  v[1]
		*/
		const Triangle t = triangles[iT];
		if(t.vertices[0] == iV1){
			iV2 = t.vertices[1];
			iV4 = t.vertices[2];
			n1 = t.neighbours[0];
			n3 = t.neighbours[2];
			iTopo = t.neighbours[1];
		}else if(t.vertices[1] == iV1){
			iV2 = t.vertices[2];
			iV4 = t.vertices[0];
			n1 = t.neighbours[1];
			n3 = t.neighbours[0];
			iTopo = t.neighbours[2];
		}else{
			iV2 = t.vertices[0];
			iV4 = t.vertices[1];
			n1 = t.neighbours[2];
			n3 = t.neighbours[1];
			iTopo = t.neighbours[0];
		}
		if(iTopo == noNeighbour)
			return;
		const Triangle tOpo = triangles[iTopo];
		if(tOpo.neighbours[0] == iT){
			iV3 = tOpo.vertices[2];
			n2 = tOpo.neighbours[1];
			n4 = tOpo.neighbours[2];
		}else if(tOpo.neighbours[1] == iT){
			iV3 = tOpo.vertices[0];
			n2 = tOpo.neighbours[2];
			n4 = tOpo.neighbours[0];
		}else{
			iV3 = tOpo.vertices[1];
			n2 = tOpo.neighbours[0];
			n4 = tOpo.neighbours[1];
		}
	}
	/**
	Handles super-triangle vertices.
	Super-tri points are not infinitely far and influence the input points
	Three cases are possible:
		1. If one of the opposed vertices is super-tri: no flip needed
		2. One of the shared vertices is super-tri:
			check if on point is same side of line formed by non-super-tri
			vertices as the non-super-tri shared vertex
		3. None of the vertices are super-tri: normal circumcircle test
	```
	*                       v4         original edge: (v2, v4)
	*                      /|\   flip-candidate edge: (v1, v3)
	*                    /  |  \
	*                  /    |    \
	*                /      |      \
	* new vertex--> v1      |       v3
	*                \      |      /
	*                  \    |    /
	*                    \  |  /
	*                      \|/
	*                       v2
	```
	*/
	bool isFlipNeeded(VertInd iV1, VertInd iV2, VertInd iV3, VertInd iV4) const nothrow @nogc pure @safe{
		if(Edge(iV2, iV4) in fixedEdges)
			return false; //flip not needed if the original edge is fixed
		const V2 v1 = vertices[iV1];
		const V2 v2 = vertices[iV2];
		const V2 v3 = vertices[iV3];
		const V2 v4 = vertices[iV4];
		if(m_superGeomType == SuperGeometryType.superTriangle){
			//If flip-candidate edge touches super-triangle in-circumference
			//test has to be replaced with orient2D test against the line formed
			//by two non-artificial vertices (that don't belong to super-triangle)
			alias locatePointLineTV2 = locatePointLine!(T, V2);
			if(iV1 < 3){ //flip-candidate edge touches super-triangle
				//does original edge also touch super-triangle?
				if(iV2 < 3)
					return locatePointLineTV2(v2, v3, v4) == locatePointLineTV2(v1, v3, v4);
				if(iV4 < 3)
					return locatePointLineTV2(v4, v2, v3) == locatePointLineTV2(v1, v2, v3);
				return false; //original edge does not touch super-triangle
			}
			if(iV3 < 3){ //flip-candidate edge touches super-triangle
				//does original edge also touch super-triangle?
				if(iV2 < 3){
					return locatePointLineTV2(v2, v1, v4) == locatePointLineTV2(v3, v1, v4);
				}
				if(iV4 < 3){
					return locatePointLineTV2(v4, v2, v1) == locatePointLineTV2(v3, v2, v1);
				}
				return false; // original edge does not touch super-triangle
			}
			//flip-candidate edge does not touch super-triangle
			if(iV2 < 3)
				return locatePointLineTV2(v2, v3, v4) == locatePointLineTV2(v1, v3, v4);
			if(iV4 < 3)
				return locatePointLineTV2(v4, v2, v3) == locatePointLineTV2(v1, v2, v3);
		}
		return isInCircumCircle!(T, V2)(v1, v2, v3, v4);
	}
	void changeNeighbour(TriInd iT, TriInd oldNeighbour, TriInd newNeighbour) nothrow @nogc pure @safe{
		if(iT == noNeighbour)
			return;
		NeighboursArr3* nn = &triangles[iT].neighbours;
		assert((*nn)[0] == oldNeighbour || (*nn)[1] == oldNeighbour || (*nn)[2] == oldNeighbour);
		if((*nn)[0] == oldNeighbour)
			(*nn)[0] = newNeighbour;
		else if((*nn)[1] == oldNeighbour)
			(*nn)[1] = newNeighbour;
		else
			(*nn)[2] = newNeighbour;
	}
	void changeNeighbour(TriInd iT, VertInd iVEdge1, VertInd iVEdge2, TriInd newNeighbour) nothrow @nogc pure @safe
	in(iT != noNeighbour){
		Triangle* t = &triangles[iT];
		t.neighbours[edgeNeighbourInd(t.vertices, iVEdge1, iVEdge2)] = newNeighbour;
	}
	void triangulatePseudoPolygon(
		const VertInd[] poly, TriInd[Edge] outerTris, TriInd iT, TriInd iN,
		ref TriangulatePseudoPolygonTask[] iterations,
	) pure @safe
	in(poly.length > 2){
		//`iterations` is a stack
		//Note: uses iteration instead of recursion to avoid stack overflows
		iterations = [TriangulatePseudoPolygonTask(0, poly.length - 1, iT, iN, 0)];
		while(iterations.length){
			size_t iA, iB;
			TriInd iParent;
			Index iInParent;
			AliasSeq!(iA, iB, iT, iParent, iInParent) = iterations[$-1].expand; //.top
			iterations = iterations[0..$-1]; //.pop()
			assert(iB - iA > 1);
			assert(iT != noNeighbour);
			assert(iParent != noNeighbour);
			Triangle* t = &triangles[iT];
			//find Delaunay point
			const size_t iC = findDelaunayPoint(poly, iA, iB);
			const VertInd a = poly[iA];
			const VertInd b = poly[iB];
			const VertInd c = poly[iC];
			//split pseudo-polygon in two parts and triangulate them
			//Note: second part needs to be pushed on stack first to be processed first
			
			//second part: points after the Delaunay point
			if(iB - iC > 1){
				const TriInd iNext = addTriangle();
				iterations ~= TriangulatePseudoPolygonTask(iC, iB, iNext, iT, 1);
			}else{ //pseudo-poly is reduced to a single outer edge
				const outerEdge = Edge(b, c);
				const TriInd outerTri = outerTris.require(outerEdge);
				if(outerTri != noNeighbour){
					assert(outerTri != iT);
					t.neighbours[1] = outerTri;
					changeNeighbour(outerTri, c, b, iT);
				}else
					outerTris[outerEdge] = iT;
			}
			//first part: points before the Delaunay point
			if(iC - iA > 1){
				//add next triangle and add another iteration
				const TriInd iNext = addTriangle();
				iterations ~= TriangulatePseudoPolygonTask(iA, iC, iNext, iT, 2);
			}else{
				//pseudo-poly is reduced to a single outer edge
				const outerEdge = Edge(c, a);
				const TriInd outerTri = outerTris.require(outerEdge);
				if(outerTri != noNeighbour){
					assert(outerTri != iT);
					t.neighbours[2] = outerTri;
					changeNeighbour(outerTri, c, a, iT);
				}else
					outerTris[outerEdge] = iT;
			}
			/*
			Finalise the triangle
			Note: only when triangle is finalised to we add it as a neighbour
			to parent to maintain triangulation topology consistency
			*/
			triangles[iParent].neighbours[iInParent] = iT;
			t.neighbours[0] = iParent;
			t.vertices = [a, b, c];
			setAdjacentTriangle(c, iT);
		}
	}
	size_t findDelaunayPoint(
		const VertInd[] poly,
		size_t iA, size_t iB,
	) const nothrow @nogc pure @safe
	in(iB - iA > 1)
	out(o; o > iA && o < iB, "point is not between ends"){
		const V2 a = vertices[poly[iA]];
		const V2 b = vertices[poly[iB]];
		size_t output = iA + 1;
		const(V2)* c = &vertices[poly[output]]; //caching for better performance
		for(size_t i = iA + 1; i < iB; ++i){
			const(V2)* v = &vertices[poly[i]];
			if(isInCircumCircle!(T, V2)(*v, a, b, *c)){
				output = i;
				c = v; //change pointer
			}
		}
		return output;
	}
	///Note: invalidates iterators!
	TriInd addTriangle(Triangle t) nothrow pure @safe{
		if(!m_dummyTris.length){
			triangles ~= t;
			return triangles.length-1;
		}
		const TriInd nxtDummy = m_dummyTris[$-1];
		m_dummyTris = m_dummyTris[0..$-1];
		triangles[nxtDummy] = t;
		return nxtDummy;
	}
	///Note: invalidates triangle iterators!
	TriInd addTriangle() nothrow pure @safe{
		if(!m_dummyTris.length){
			const dummy = Triangle(
				[noVertex,    noVertex,    noVertex],
				[noNeighbour, noNeighbour, noNeighbour],
			);
			triangles ~= dummy;
			return triangles.length-1;
		}
		const TriInd nxtDummy = m_dummyTris[$-1];
		m_dummyTris = m_dummyTris[0..$-1];
		return nxtDummy;
	}
	/**
	Remove super-triangle (if used) and triangles with specified indices.
	Adjust internal triangulation state accordingly.
	Params:
		removedTriangles = indices of triangles to remove
	*/
	void finaliseTriangulation(const TriIndUSet removedTriangles) pure @safe{
		eraseDummies();
		m_vertTris = [];
		//remove super-triangle
		if(m_superGeomType == SuperGeometryType.superTriangle){
			vertices = vertices[3..$];
			//Edge re-mapping
			{ //fixed edges
				EdgeUSet updatedFixedEdges;
				foreach(edge; fixedEdges.keys){
					updatedFixedEdges[remapNoSuperTriangle(edge)] = null; //set insertion
				}
				fixedEdges = updatedFixedEdges;
			}
			{ //overlap count
				BoundaryOverlapCount[Edge] updatedOverlapCount;
				foreach(item; overlapCount.byKeyValue){
					updatedOverlapCount[remapNoSuperTriangle(item.key)] = item.value;
				}
				overlapCount = updatedOverlapCount;
			}
			{ //split edges mapping
				Edge[][Edge] updatedPieceToOriginals;
				foreach(ref ee; pieceToOriginals.byKeyValue){
					auto eeVal = ee.value;
					foreach(ref eeItem; eeVal){
						eeItem = remapNoSuperTriangle(eeItem);
					}
					updatedPieceToOriginals[remapNoSuperTriangle(ee.key)] = eeVal;
				}
				pieceToOriginals = updatedPieceToOriginals;
			}
		}
		//remove other triangles
		removeTriangles(removedTriangles);
		//adjust triangle vertices: account for removed super-triangle
		if(m_superGeomType == SuperGeometryType.superTriangle){
			foreach(ref t; triangles){
				foreach(ref v; t.vertices){
					v -= 3;
				}
			}
		}
	}
	TriIndUSet growToBoundary(const(TriInd)[] seeds) const pure @safe{
		TriIndUSet traversed;
		while(seeds.length){
			//`seeds` is a stack
			const TriInd iT = seeds[$-1]; //.top
			seeds = seeds[0..$-1]; //.pop()
			
			traversed[iT] = null; //this inserts `iT` into the `traversed` set
			const Triangle t = triangles[iT];
			for(Index i=0; i < 3; ++i){
				const opEdge = Edge(t.vertices[ccw(i)], t.vertices[cw(i)]);
				if(opEdge in fixedEdges)
					continue;
				const TriInd iN = t.neighbours[opoNbr(i)];
				if(iN != noNeighbour && iN !in traversed){
					seeds ~= iN;
				}
			}
		}
		return traversed;
	}
	void fixEdge(Edge edge) pure @safe{
		if(edge in fixedEdges){
			//if it's already in the set, increase overlap count:
			overlapCount.require(edge)++;
		}else{
			fixedEdges[edge] = null; //add `edge` to set
		}
	}
	void fixEdge(Edge edge, Edge originalEdge) pure @safe{
		fixEdge(edge);
		if(edge != originalEdge)
			pieceToOriginals.require(edge).insertUnique(originalEdge);
	}
	/**
	Split existing constraint (fixed) edge
	Params:
		edge = fixed edge to split
		iSplitVert = index of the vertex to be used as a split vertex
	*/
	void splitFixedEdge(Edge edge, VertInd iSplitVert) pure @safe{
		//split constraint (fixed) edge that already exists in triangulation
		const half1 = Edge(edge.v1, iSplitVert);
		const half2 = Edge(iSplitVert, edge.v2);
		//remove the edge that and add its halves
		fixedEdges.remove(edge);
		fixEdge(half1);
		fixEdge(half2);
		//maintain overlaps
		if(auto overlaps = edge in overlapCount){
			overlapCount.require(half1) += *overlaps;
			overlapCount.require(half2) += *overlaps;
			overlapCount.remove(edge);
		}
		//maintain piece-to-original mapping
		Edge[] newOriginals;
		if(auto originalEdge = edge in pieceToOriginals){
			//edge being split was split before: pass-through originals
			newOriginals = *originalEdge;
			pieceToOriginals.remove(edge);
		}else{
			newOriginals = [edge];
		}
		pieceToOriginals.require(half1).insertUnique(newOriginals);
		pieceToOriginals.require(half2).insertUnique(newOriginals);
	}
	/**
	Add a vertex that splits an edge into the triangulation
	Params:
		splitVert = position of split vertex
		iT = index of a first triangle adjacent to the split edge
		iTopo = index of a second triangle adjacent to the split edge
			(opposed to the first triangle)
	Returns: index of a newly added split vertex
	*/
	VertInd addSplitEdgeVertex(V2 splitVert, TriInd iT, TriInd iTopo) nothrow pure @safe{
		//add a new point on the edge that splits an edge in two
		const VertInd iSplitVert = vertices.length;
		addNewVertex(splitVert, noNeighbour);
		TriInd[] triStack = insertVertexOnEdge(iSplitVert, iT, iTopo);
		tryAddVertexToLocator(iSplitVert);
		ensureDelaunayByEdgeFlips(iSplitVert, triStack);
		return iSplitVert;
	}
	/**
	Split fixed edge and add a split vertex into the triangulation
	Params:
		edge = fixed edge to split
		splitVert = position of split vertex
		iT = index of a first triangle adjacent to the split edge
		iTopo = index of a second triangle adjacent to the split edge
			(opposed to the first triangle)
	Returns: index of a newly added split vertex
	*/
	VertInd splitFixedEdgeAt(Edge edge, V2 splitVert, TriInd iT, TriInd iTopo) pure @safe{
		const VertInd iSplitVert = addSplitEdgeVertex(splitVert, iT, iTopo);
		splitFixedEdge(edge, iSplitVert);
		return iSplitVert;
	}
	/**
	Flag triangle as dummy
	Note: Advanced method for manually modifying the triangulation from
	outside. Please call it when you know what you are doing.
	Params:
		iT = index of a triangle to flag
	*/
	void makeDummy(TriInd iT) nothrow pure @safe{
		m_dummyTris ~= iT;
	}
	/**
	Erase all dummy triangles
	Note: Advanced method for manually modifying the triangulation from
	outside. Please call it when you know what you are doing.
	*/
	void eraseDummies() pure @safe{
		if(!m_dummyTris.length)
			return;
		const dummySet = {
			TriIndUSet ret;
			foreach(dummyTri; m_dummyTris)
				ret[dummyTri] = null; //set insertion
			return ret;
		}();
		TriIndUMap triIndMap = [noNeighbour: noNeighbour];
		TriInd iTNew = 0;
		foreach(TriInd iT; 0..triangles.length){
			if(iT in dummySet)
				continue;
			triIndMap[iT] = iTNew;
			triangles[iTNew] = triangles[iT];
			iTNew++;
		}
		triangles = triangles[$-dummySet.length..$];
		//remap adjacent triangle indices for vertices
		foreach(ref iT; m_vertTris){
			iT = triIndMap.require(iT);
		}
		//remap neighbour indices for triangles
		foreach(ref t; triangles)
			foreach(ref iN; t.neighbours)
				iN = triIndMap.require(iN);
		//clear dummy triangles
		m_dummyTris = [];
	}
	/**
	Depth-peel a layer in triangulation, used when calculating triangle
	depths
	
	It takes starting seed triangles, traverses neighbouring triangles,
	and assigns given layer depth to the traversed triangles. Traversal
	is blocked by constraint edges. Triangles behind constraint edges
	are recorded as seeds of next layer and returned from the function.
	Params:
		seeds = indices of seed triangles
		layerDepth = current layer's depth to mark triangles with
		triDepths = depths of triangles
	Returns: triangles of the deeper layers that are adjacent to the
	peeled layer. To be used as seeds when peeling deeper layers.
	 */
	LayerDepth[TriInd] peelLayer(
		TriInd[] seeds, LayerDepth layerDepth, ref LayerDepth[] triDepths,
	) const nothrow pure @safe{
		LayerDepth[TriInd] behindBoundary;
		while(seeds.length){ //seeds is a stack
			const TriInd iT = seeds[$-1]; //.top
			seeds = seeds[0..$-1]; //.pop()
			
			triDepths[iT] = min(triDepths[iT], layerDepth);
			behindBoundary.remove(iT);
			const Triangle t = triangles[iT];
			for(Index i=0; i < 3; i++){
				const opEdge = Edge(t.vertices[ccw(i)], t.vertices[cw(i)]);
				const TriInd iN = t.neighbours[opoNbr(i)];
				if(iN == noNeighbour || triDepths[iN] <= layerDepth)
					continue;
				if(opEdge in fixedEdges){
					if(auto overlaps = opEdge in overlapCount){
						behindBoundary[iN] = cast(LayerDepth)(layerDepth + *overlaps + 1);
					}else{
						behindBoundary[iN] = cast(LayerDepth)(layerDepth + 1);
					}
					continue;
				}
				seeds ~= iN;
			}
		}
		return behindBoundary;
	}
	void insertVerticesAsProvided(VertInd superGeomVertCount) pure @safe{
		foreach(VertInd iV; superGeomVertCount..vertices.length)
			insertVertex(iV);
	}
	void insertVerticesRandomised(VertInd superGeomVertCount) pure @safe{
		size_t vertexCount = vertices.length - superGeomVertCount;
		auto ii = new VertInd[](vertexCount);
		foreach(ref item; ii)
			item = superGeomVertCount++;
		
		ii.randomShuffle();
		foreach(item; ii){
			insertVertex(item);
		}
	}
	void insertVerticesKDTreeBFS(VertInd superGeomVertCount, V2 boxMin, V2 boxMax) @safe{
		//calculate original indices
		const VertInd vertexCount = vertices.length - superGeomVertCount;
		if(vertexCount <= 0)
			return;
		auto ii = new VertInd[](vertexCount);
		foreach(ref item; ii)
			item = superGeomVertCount++;
		
		alias QueueItem = Tuple!(VertInd[], V2, V2, VertInd);
		QueueItem[] queue;
		queue.reserve(maxQueueLengthBFSKDTree(vertexCount));
		queue ~= QueueItem(ii, boxMin, boxMax, 0);
		
		VertInd[] range;
		V2 newBoxMin, newBoxMax;
		VertInd parent, mid;
		
		while(queue.length){
			AliasSeq!(range, boxMin, boxMax, parent) = queue[0].expand;
			queue = queue[1..$];
			assert(range.length);
			
			if(range.length == 1){
				insertVertex(range[0], parent);
				continue;
			}
			const size_t midInd = range.length / 2;
			if(boxMax[0] - boxMin[0] >= boxMax[1] - boxMin[1]){
				topN!((a, b) => vertices[a][0] < vertices[b][0])(range, midInd);
				mid = range[midInd];
				const T split = vertices[mid].x;
				newBoxMin.x = split;
				newBoxMin.y = boxMin.y;
				newBoxMax.x = split;
				newBoxMax.y = boxMax.y;
			}else{
				topN!((a, b) => vertices[a][1] < vertices[b][1])(range, midInd);
				mid = range[midInd];
				const T split = vertices[mid].y;
				newBoxMin.x = boxMin.x;
				newBoxMin.y = split;
				newBoxMax.x = boxMax.x;
				newBoxMax.y = split;
			}
			insertVertex(mid, parent);
			if(midInd > 0)
				queue ~= QueueItem(range[0..midInd], boxMin, newBoxMax, mid);
			if(midInd+1 < range.length)
				queue ~= QueueItem(range[midInd+1..$], newBoxMin, boxMax, mid);
		}
	}
	TriInd[2] edgeTriangles(VertInd a, VertInd b) const nothrow pure @safe{
		const TriInd triStart = m_vertTris[a];
		assert(triStart != noNeighbour);
		TriInd iT = triStart, iTNext = triStart;
		VertInd iV = noVertex;
		do{
			const Triangle t = triangles[iT];
			AliasSeq!(iTNext, iV) = t.next(a);
			assert(iTNext != noNeighbour);
			if(iV == b){
				return [iT, iTNext];
			}
			iT = iTNext;
		}while(iT != triStart);
		return [invalidIndex, invalidIndex];
	}
	bool hasEdge(VertInd a, VertInd b) const nothrow pure @safe =>
		edgeTriangles(a, b)[0] != invalidIndex;
	
	void setAdjacentTriangle(VertInd v, TriInd t) nothrow pure @safe
	in(t != noNeighbour)
	out(; triangles[t].vertices[0] == v || triangles[t].vertices[1] == v || triangles[t].vertices[2] == v){
		m_vertTris[v] = t;
	}
	void pivotVertexTriangleCW(VertInd v) nothrow pure @safe
	in(m_vertTris[v] != noNeighbour)
	out(; m_vertTris[v] != noNeighbour)
	out(; triangles[m_vertTris[v]].vertices[0] == v || triangles[m_vertTris[v]].vertices[1] == v || triangles[m_vertTris[v]].vertices[2] == v){
		m_vertTris[v] = triangles[m_vertTris[v]].next(v)[0];
	}
	///Add vertex to nearest-point locator if locator is initialised
	void tryAddVertexToLocator(VertInd v){
		if(m_nearPtLocator.length) //only if locator is initialised already
			m_nearPtLocator.addInd(v, vertices);
	}
	///Perform lazy initialisation of nearest-point locator after the Kd-tree
	///BFS bulk load if necessary
	void tryIniNearestPointLocator(){
		if(vertices.length && !m_nearPtLocator.length)
			m_nearPtLocator = TNearPointLocator(vertices);
	}
}

///pseudo-random number generator
struct SplitMix64RandGen{
	ulong m_state = 0;
	ulong opCall() nothrow @nogc pure @safe{
		ulong z = (m_state += 0x9E3779B97F4A7C15UL);
		z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL;
		z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL;
		return z ^ (z >> 31);
	}
}

void randomShuffle(T)(T[] range) nothrow @nogc pure @safe{
	SplitMix64RandGen rng;
	foreach_reverse(i, ref item; range[1..$]){
		swap(item, range[rng() % (i + 1)]);
	}
}

///Remap removing super-triangle: subtract 3 from vertices
pragma(inline,true) Edge remapNoSuperTriangle(Edge e) nothrow @nogc pure @safe =>
	Edge(VertInd(e.v1 - 3), VertInd(e.v2 - 3));

///Add element to 'to' if not already in 'to'
void insertUnique(T)(ref T[] to, T elem) nothrow pure @safe{
	if(!to.canFind(elem)){
		to ~= elem;
	}
}
///Add elements of 'from' that are not present in 'to' to 'to'
void insertUnique(T)(ref T[] to, const T[] from) nothrow pure @safe{
	to.reserve(to.length + from.length);
	foreach(item; from){
		insertUnique!T(to, item);
	}
}

private T lerp(T)(T a, T b, T t) nothrow @nogc pure @safe =>
	T((T(1) - t) * a + t * b);

///Precondition: ab and cd intersect normally
V2 intersectionPosition(T, V2=T[2])(V2 a, V2 b, V2 c, V2 d) nothrow @nogc pure @safe{
	//note: for better accuracy we interpolate x and y separately
	//on a segment with the shortest x/y-projection correspondingly
	const T aCD = orient2DAdaptive!T(c[0], c[1], d[0], d[1], a[0], a[1]);
	const T bCD = orient2DAdaptive!T(c[0], c[1], d[0], d[1], b[0], b[1]);
	const T tAB = aCD / (aCD - aCD);
	
	const T cAB = orient2DAdaptive!T(a[0], a[1], b[0], b[1], c[0], c[1]);
	const T dAB = orient2DAdaptive!T(a[0], a[1], b[0], b[1], d[0], d[1]);
	const T tCD = cAB / (cAB - cAB);
	
	V2 ret;
	ret[0] = abs(a[0] - b[0]) < abs(c[0] - d[0]) ?
		lerp(a.x, b.x, tAB) :
		lerp(c.x, d.x, tCD);
	ret[1] = abs(a[1] - b[1]) < abs(c[1] - d[1]) ?
		lerp(a.y, b.y, tAB) :
		lerp(c.y, d.y, tCD);
	return ret;
}

/**
Since KD-tree bulk load builds a balanced tree the maximum length of a
queue can be pre-calculated: it is calculated as size of a completely
filled tree layer plus the number of the nodes on a completely filled
layer that have two children.
*/
pragma(inline,true) size_t maxQueueLengthBFSKDTree(size_t vertexCount) nothrow @nogc pure @safe{
	const int filledLayerPow2 = cast(int)(floor(log2(cast(double)vertexCount)) - 1);
	const size_t nodesInFilledTree      = cast(size_t)((2 ^^ filledLayerPow2 + 1) - 1);
	const size_t nodesInLastFilledLayer = cast(size_t)((2 ^^ filledLayerPow2));
	const size_t nodesInLastLayer = vertexCount - nodesInFilledTree;
	return nodesInLastLayer >= nodesInLastFilledLayer ?
		nodesInLastFilledLayer + nodesInLastLayer - nodesInLastFilledLayer :
		nodesInLastFilledLayer;
}
