/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 * Contribution of original implementation:
 * Andre Fecteau <andre.fecteau1@gmail.com> */

module cdt.kdtree;

import cdt.utils;

import std.algorithm, std.typecons;

enum NodeSplitDirection{
	x, y,
}

/**
Simple tree structure with alternating half splitting nodes

- Tree to incrementally add points to the structure.
- Get the nearest point to a given input.
- Does not check for duplicates, expect unique points.

Params:
	TCoordType = type used for storing point coordinates
	TPointType = type used for storing points
	numVerticesInLeaf = the number of points per leaf
	initialStackDepth = initial size of stack depth for nearest query. Should be at least 1
	stackDepthIncrement = increment of stack depth for nearest query when stack depth is reached
**/
struct KDTree(
	TCoordType, TPointType=TCoordType[2],
	size_t numVerticesInLeaf,
	size_t initialStackDepth,
	size_t stackDepthIncrement,
){
	alias CoordType = TCoordType;
	alias PointType = TPointType;
	alias PointIndex = VertInd;
	alias ValueType = Tuple!(PointType, PointIndex);
	alias NodeIndex = VertInd;
	alias ChildrenType = NodeIndex[2];
	
	private{
		NodeIndex m_root = 0;
		Node[] m_nodes = [Node(0)];
		NodeSplitDirection m_rootDir = NodeSplitDirection.x;
		PointType m_min = { PointType p; p[0] = p[1] = -CoordType.max; return p; }();
		PointType m_max = { PointType p; p[0] = p[1] =  CoordType.max; return p; }();
		VertInd m_size = 0;
		
		bool m_isRootBoxIni = false;
		
		//used for nearest query
		struct NearestTask{
			NodeIndex node;
			PointType min, max;
			NodeSplitDirection dir;
			CoordType distSq;
		}
		//allocated at construction (not in the 'nearest' method) for better performance
		NearestTask[] m_tasksStack = new NearestTask[](initialStackDepth);
	}
	
	///Stores kd-tree node data
	struct Node{
		ChildrenType children; ///two children if not leaf; {0,0} if leaf
		PointIndex[] data; ///points' data if leaf
		
		this(int _){
			data.reserve(numVerticesInLeaf);
		}
		
		///Children setter for convenience
		void setChildren(NodeIndex c1, NodeIndex c2) nothrow @nogc pure @safe{
			children[0] = c1;
			children[1] = c2;
		}
		///Check if node is a leaf (has no valid children)
		bool isLeaf() const nothrow @nogc pure @safe =>
			children[0] == children[1];
	}
	
	///Constructor with bounding box known in advance
	this(PointType min, PointType max) nothrow @nogc pure @safe{
		m_min = min;
		m_max = max;
		m_isRootBoxIni = true;
	}
	
	@property VertInd length() const nothrow @nogc pure @safe =>
		m_size;
	
	/**
	Insert a point into kd-tree
	Note: external point-buffer is used to reduce kd-tree's memory footprint
	Params:
		iPoint = index of point in external point-buffer
		points = external point-buffer
	*/
	void insert(PointIndex iPoint, const PointType[] points) nothrow pure @safe{
		m_size++;
		
		//if point is outside root, extend tree by adding new roots
		const pos = points[iPoint];
		while(!isInsideBox(pos, m_min, m_max)){
			extendTree(pos);
		}
		//now insert the point into the tree
		NodeIndex node = m_root;
		PointType minVal = m_min;
		PointType maxVal = m_max;
		NodeSplitDirection dir = m_rootDir;
		
		//below: initialised only to suppress warnings
		auto newDir = NodeSplitDirection.x;
		CoordType mid = CoordType(0);
		PointType newMin, newMax;
		while(true){
			if(m_nodes[node].isLeaf()){
				//add point if capacity is not reached
				PointIndex[]* pd = &m_nodes[node].data;
				if(pd.length < numVerticesInLeaf){
					*pd ~= iPoint;
					return;
				}
				//initialise bbox first time the root capacity is reached
				if(!m_isRootBoxIni){
					iniRootBox(points);
					minVal = m_min;
					maxVal = m_max;
				}
				//split a full leaf node
				calcSplitInfo(minVal, maxVal, dir, mid, newDir, newMin, newMax);
				const NodeIndex c1 = addNewNode(), c2 = addNewNode();
				Node* n = &m_nodes[node];
				n.setChildren(c1, c2);
				PointIndex[]* c1Data = &m_nodes[c1].data;
				PointIndex[]* c2Data = &m_nodes[c2].data;
				
				//move node's points to children
				foreach(inds; n.data){
					if(whichChild(points[inds], mid, dir) == 0){
						*c1Data ~= inds;
					}else{
						*c2Data ~= inds;
					}
				}
				n.data = [];
			}else{
				calcSplitInfo(minVal, maxVal, dir, mid, newDir, newMin, newMax);
			}
			//add the point to a child
			const size_t iChild = whichChild(points[iPoint], mid, dir);
			iChild == 0 ? (maxVal = newMax) : (minVal = newMin);
			node = m_nodes[node].children[iChild];
			dir = newDir;
		}
	}
	
	/**
	Query kd-tree for a nearest neighbour point
	Note: external point-buffer is used to reduce kd-tree's memory footprint
	Params:
		point query point position
		points external point-buffer
	*/
	ValueType nearest(PointType point, const PointType[] points) nothrow pure @safe{
		ValueType output;
		
		int iTask = -1;
		CoordType minDistSq = CoordType.max;
		m_tasksStack[++iTask] = NearestTask(m_root, m_min, m_max, m_rootDir, minDistSq);
		
		while(iTask != -1){
			const NearestTask t = m_tasksStack[iTask--];
			if(t.distSq > minDistSq)
				continue;
			const Node n = m_nodes[t.node];
			if(n.isLeaf()){
				foreach(inds; n.data){
					const PointType p = points[inds];
					const CoordType distSq = (a, b){
						const CoordType[2] d = [b[0] - a[0], b[1] - a[1]];
						return d[0]*d[0] + d[1]*d[1];
					}(point, p);
					if(distSq < minDistSq){
						minDistSq = distSq;
						output[0] = p;
						output[1] = inds;
					}
				}
			}else{
				CoordType mid = 0;
				NodeSplitDirection newDir;
				PointType newMin, newMax;
				calcSplitInfo(t.min, t.max, t.dir, mid, newDir, newMin, newMax);
				
				const CoordType distToMid =
					t.dir == NodeSplitDirection.x ?
						point.x - mid :
						point.y - mid;
				const CoordType toMidSq = distToMid * distToMid;
				
				const size_t iChild = whichChild(point, mid, t.dir);
				if(iTask + 2 >= cast(int)m_tasksStack.length){
					m_tasksStack.length += stackDepthIncrement;
				}
				//node containing point should end up on top of the stack
				if(iChild == 0){
					m_tasksStack[++iTask] = NearestTask(n.children[1], newMin, t.max,  newDir, toMidSq);
					m_tasksStack[++iTask] = NearestTask(n.children[0], t.min,  newMax, newDir, toMidSq);
				}else{
					m_tasksStack[++iTask] = NearestTask(n.children[0], t.min,  newMax, newDir, toMidSq);
					m_tasksStack[++iTask] = NearestTask(n.children[1], newMin, t.max,  newDir, toMidSq);
				}
			}
		}
		return output;
	}
	
	private:
	///Add a new node and return its index in nodes buffer
	NodeIndex addNewNode(){
		const NodeIndex newNodeIndex = m_nodes.length;
		m_nodes ~= Node(0);
		return newNodeIndex;
	}
	
	/**
	Test which child point belongs to after the split
	Returns: 0 if first child, 1 if second child
	*/
	size_t whichChild(PointType point, CoordType split, NodeSplitDirection dir) const nothrow @nogc pure @safe =>
		dir == NodeSplitDirection.x ? point.x > split : point.y > split;
	
	///Calculate split location, direction, and children boxes
	static void calcSplitInfo(
		PointType min, PointType max,
		NodeSplitDirection dir,
		ref CoordType midOut,
		ref NodeSplitDirection newDirOut,
		ref PointType newMinOut,
		ref PointType newMaxOut,
	) nothrow @nogc pure @safe{
		newMaxOut = max;
		newMinOut = min;
		switch(dir){
			case NodeSplitDirection.x:
				midOut = (min.x + max.x) / 2;
				newDirOut = NodeSplitDirection.y;
				newMinOut.x = midOut;
				newMaxOut.x = midOut;
				return;
			case NodeSplitDirection.y:
				midOut = (min.y + max.y) / 2;
				newDirOut = NodeSplitDirection.x;
				newMinOut.y = midOut;
				newMaxOut.y = midOut;
				return;
			default:
		}
	}
	
	///Test if point is inside a box
	static bool isInsideBox(PointType p, PointType min, PointType max) nothrow @nogc pure @safe =>
		p.x >= min.x && p.x <= max.x && p.y >= min.y && p.y <= max.y;
	
	///Extend a tree by creating new root with old root and a new node as children
	void extendTree(PointType point) nothrow pure @safe{
		const NodeIndex newRoot = addNewNode();
		const NodeIndex newLeaf = addNewNode();
		switch(m_rootDir){
			case NodeSplitDirection.x:
				m_rootDir = NodeSplitDirection.y;
				point.y < m_min.y ?
					m_nodes[newRoot].setChildren(newLeaf, m_root) :
					m_nodes[newRoot].setChildren(m_root,  newLeaf);
				if(point.y < m_min.y)
					m_min.y -= m_max.y - m_min.y;
				else if(point.y > m_max.y)
					m_max.y += m_max.y - m_min.y;
				break;
			case NodeSplitDirection.y:
				m_rootDir = NodeSplitDirection.x;
				point.x < m_min.x ?
					m_nodes[newRoot].setChildren(newLeaf, m_root) :
					m_nodes[newRoot].setChildren(m_root,  newLeaf);
				if(point.x < m_min.x)
					m_min.x -= m_max.x - m_min.x;
				else if(point.x > m_max.x)
					m_max.x += m_max.x - m_min.x;
				break;
			default:
		}
		m_root = newRoot;
	}
	
	///Calculate root's box enclosing all root points
	void iniRootBox(const PointType[] points) nothrow @nogc pure @safe{
		const PointIndex[] data = m_nodes[m_root].data;
		m_min = points[data[0]];
		m_max = m_min;
		foreach(inds; data){
			const PointType p = points[inds];
			m_min[0] = min(m_min[0], p[0]);
			m_min[1] = min(m_min[1], p[1]);
			m_max[0] = max(m_max[0], p[0]);
			m_max[1] = max(m_max[1], p[1]);
		}
		//make sure bounding box does not have a zero size by adding padding:
		//zero-size bounding box cannot be extended properly
		const CoordType padding = 1;
		if(m_min[0] == m_max[0]){
			m_min[0] -= padding;
			m_max[0] += padding;
		}
		if(m_min[1] == m_max[1]){
			m_min[1] -= padding;
			m_max[1] += padding;
		}
		m_isRootBoxIni = true;
	}
}

///KD-tree holding points
struct LocatorKDTree(
	TCoordType, TPointType=TCoordType[2],
	size_t numVerticesInLeaf=32,
	size_t initialStackDepth=32,
	size_t stackDepthIncrement=32,
){
	private{
		alias Tree = KDTree!(TCoordType, TPointType, numVerticesInLeaf, initialStackDepth, stackDepthIncrement);
		Tree m_kdTree;
	}
	///Initialise KD-tree with points
	this(const TPointType[] points) nothrow pure @safe{
		TPointType minVal = points[0];
		TPointType maxVal = minVal;
		foreach(point; points){
			minVal[0] = min(minVal[0], point[0]);
			minVal[1] = min(minVal[1], point[1]);
			maxVal[0] = max(maxVal[0], point[0]);
			maxVal[1] = max(maxVal[1], point[1]);
		}
		m_kdTree = Tree(minVal, maxVal);
		for(VertInd i=0; i < points.length; ++i){
			m_kdTree.insert(i, points);
		}
	}
	///Add point to KD-tree
	void addInd(const VertInd i, const TPointType[] points) nothrow pure @safe{
		m_kdTree.insert(i, points);
	}
	///Find nearest point using R-tree
	VertInd nearestInd(
		ref const TPointType pos,
		const TPointType[] points,
	) nothrow pure @safe =>
		m_kdTree.nearest(pos, points)[1];
	
	VertInd length() const nothrow @nogc pure @safe =>
		m_kdTree.length;
}
