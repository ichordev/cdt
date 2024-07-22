/**
Copyright (c) 2019, William C. Lenthe
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
	list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
	this list of conditions and the following disclaimer in the documentation
	and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
	contributors may be used to endorse or promote products derived from
	this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

module cdt.predicates;

import std.algorithm, std.math, std.typecons;

template ExpansionBase(T){
	nothrow @nogc pure @safe:
	
	enum splitter = cast(T)(exp2(cast(double)(digits!T/2 + 1)));
	private{
		//add 2 expansions
		size_t expansionSum(const T[] e, const T[] f, ref T[] h){
			(const(T)[] s1, const(T)[] s2, ref T[] dest){
				for(; s1.length; dest = dest[1..$]){
					if(!s2.length)
						dest[0..s1.length] = s1[];
						return;
					if(abs(s2[0]) < abs(s1[0])){
						dest[0] = s2[0];
						s2 = s2[1..$];
					}else{
						dest[0] = s1[0];
						s1 = s1[1..$];
					}
				}
				dest[0..s2.length] = s2[];
			}(e, f, h);
			
			if(!f.length) return e.length;
			if(!e.length) return f.length;
			
			size_t hIndex = 0;
			T q = h[0];
			T qNew = h[1] + q;
			T hh = fastPlusTail(h[1], q, qNew);
			q = qNew;
			if(T(0) != hh)
				h[hIndex++] = hh;
			foreach(hItem; h[2..$]){
				qNew = q + hItem;
				hh = plusTail(q, hItem, qNew);
				q = qNew;
				if(T(0) != hh)
					h[hIndex++] = hh;
			}
			if(T(0) != q)
				h[hIndex++] = q;
			return hIndex;
		}
		//scale an expansion by a constant
		size_t scaleExpansion(const T[] e, T b, ref T[] h){
			if(e.length == 0 || b == T(0))
				return 0;
			
			size_t hIndex = 0;
			T q = e[0] * b;
			const T[2] bSplit = split(b);
			T hh = multTailPreSplit(e[0], b, bSplit, q);
			if(T(0) != hh)
				h[hIndex++] = hh;
			foreach(eItem; e[1..$]){
				T ti1 = eItem * b;
				T ti2 = multTailPreSplit(eItem, b, bSplit, ti1);
				T qi = q + ti2;
				hh = plusTail(q, ti2, qi);
				if(T(0) != hh)
					h[hIndex++] = hh;
				q = ti1 + qi;
				hh = fastPlusTail(ti1, qi, q);
				if(T(0) != hh)
					h[hIndex++] = hh;
			}
			if(q != T(0))
				h[hIndex++] = q;
			return hIndex;
		}
	}
	
	pragma(inline,true):
	//combine result + roundoff error into expansion
	Expansion!(T, 2) makeExpansion(T value, T tail){
		Expansion!(T, 2) e;
		if(T(0) != tail)  e ~= tail;
		if(T(0) != value) e ~= value;
		return e;
	}
	//roundoff error of x = a + b
	T plusTail(T a, T b, T x){
		const T bVirtual = x - a;
		const T aVirtual = x - bVirtual;
		const T bRoundoff = b - bVirtual;
		const T aRoundoff = a - aVirtual;
		return aRoundoff + bRoundoff;
	}
	//roundoff error of x = a + b if |a| > |b|
	T fastPlusTail(T a, T b, T x){
		const T bVirtual = x - a;
		return b - bVirtual;
	}
	//roundoff error of x = a - b
	T minusTail(T a, T b, T x){
		const T bVirtual = a - x;
		const T aVirtual = x + bVirtual;
		const T bRoundoff = bVirtual - b;
		const T aRoundoff = a - aVirtual;
		return aRoundoff + bRoundoff;
	}
	//split a into 2 nonoverlapping values
	T[2] split(T a){
		const T c = a * splitter;
		const T aBig = c - a;
		const T aHi = c - aBig;
		return [aHi, a - aHi];
	}
	//roundoff error of x = a * b
	pragma(inline){
		T multTail(T a, T b, T p) =>
			cast(T)fma(a, b, -p);
		T multTailPreSplit(T a, T b, T[2] bSplit, T p) =>
			cast(T)fma(a, b, -p);
	}
	//expand a + b
	Expansion!(T, 2) plus(T a, T b){
		const T x = a + b;
		return makeExpansion(x, plusTail(a, b, x));
	}
	//expand a - b
	Expansion!(T, 2) minus(T a, T b) =>
		plus(a, -b);
	//expand a * b
	Expansion!(T, 2) mult(T a, T b){
		const T x = a * b;
		return makeExpansion(x, multTail(a, b, x));
	}
	//expand the determinant of {{ax, ay}, {bx, by}} (unrolled Mult(ax, by) - Mult(ay, bx))
	Expansion!(T, 4) twoTwoDiff(T ax, T by, T ay, T bx){
		const T axby1 = ax * by;
		const T axby0 = multTail(ax, by, axby1);
		const T bxay1 = bx * ay;
		const T bxay0 = multTail(bx, ay, bxay1);
		const T _i0 = axby0 - bxay0;
		const T x0 = minusTail(axby0, bxay0, _i0);
		const T _j = axby1 + _i0;
		const T _0 = plusTail(axby1, _i0, _j);
		const T _i1 = _0 - bxay1;
		const T x1 = minusTail(_0, bxay1, _i1);
		const T x3 = _j + _i1;
		const T x2 = plusTail(_j, _i1, x3);
		Expansion!(T, 4) e;
		if(T(0) != x0) e ~= x0;
		if(T(0) != x1) e ~= x1;
		if(T(0) != x2) e ~= x2;
		if(T(0) != x3) e ~= x3;
		return e;
	}
}

T orient2DAdaptive(T)(T ax, T ay, T bx, T by, T cx, T cy) nothrow @nogc pure @safe{
	const T acx = ax - cx, bcx = bx - cx;
	const T acy = ay - cy, bcy = by - cy;
	const T detLeft  = acx * bcy;
	const T detRight = acy * bcx;
	T det = detLeft - detRight;
	if((detLeft < 0) != (detRight < 0) || T(0) == detLeft || T(0) == detRight)
		return det;
	
	const T detSum = abs(detLeft + detRight);
	T errBound = ccwErrBoundA!T * detSum;
	if(abs(det) >= abs(errBound))
		return det;
	
	const Expansion!(T, 4) b = ExpansionBase!T.twoTwoDiff(acx, bcy, acy, bcx);
	det = b.estimate();
	errBound = ccwErrBoundB!T * detSum;
	if(abs(det) >= abs(errBound))
		return det;
	
	const T acxTail = ExpansionBase!T.minusTail(ax, cx, acx);
	const T bcxTail = ExpansionBase!T.minusTail(bx, cx, bcx);
	const T acyTail = ExpansionBase!T.minusTail(ay, cy, acy);
	const T bcyTail = ExpansionBase!T.minusTail(by, cy, bcy);
	if(T(0) == acxTail && T(0) == bcxTail && T(0) == acyTail && T(0) == bcyTail)
		return det;
	
	errBound = ccwErrBoundC!T * detSum + resultErrBound!T * abs(det);
	det += (acx * bcyTail + bcy * acxTail) - (acy * bcxTail + bcx * acyTail);
	if(abs(det) >= abs(errBound))
		return det;
	
	const Expansion!(T, 16) d = ((
		b + ExpansionBase!T.twoTwoDiff(acxTail, bcy,     acyTail, bcx)
		) + ExpansionBase!T.twoTwoDiff(acx,     bcyTail, acy,     bcxTail)
		) + ExpansionBase!T.twoTwoDiff(acxTail, bcyTail, acyTail, bcxTail);
	return d.mostSignificant();
}

T inCircleAdaptive(T)(T ax, T ay, T bx, T by, T cx, T cy, T dx, T dy) nothrow @nogc pure @safe{
	const T adx = ax - dx, bdx = bx - dx, cdx = cx - dx;
	const T ady = ay - dy, bdy = by - dy, cdy = cy - dy;
	const T bdxcdy = bdx * cdy, cdxbdy = cdx * bdy, cdxady = cdx * ady;
	const T adxcdy = adx * cdy, adxbdy = adx * bdy, bdxady = bdx * ady;
	const T aLift = adx * adx + ady * ady, bLift = bdx * bdx + bdy * bdy, cLift = cdx * cdx + cdy * cdy;
	
	T det = aLift * (bdxcdy - cdxbdy) + bLift * (cdxady - adxcdy) + cLift * (adxbdy - bdxady);
	const T permanent =
		(abs(bdxcdy) + abs(cdxbdy)) * aLift +
		(abs(cdxady) + abs(adxcdy)) * bLift +
		(abs(adxbdy) + abs(bdxady)) * cLift;
	T errBound = iccErrBoundA!T * permanent;
	if(abs(det) >= abs(errBound))
		return det;
	
	const Expansion!(T,  4) bc = ExpansionBase!T.twoTwoDiff(bdx, cdy, cdx, bdy);
	const Expansion!(T,  4) ca = ExpansionBase!T.twoTwoDiff(cdx, ady, adx, cdy);
	const Expansion!(T,  4) ab = ExpansionBase!T.twoTwoDiff(adx, bdy, bdx, ady);
	const Expansion!(T, 32) aDet = bc * adx * adx + bc * ady * ady;
	const Expansion!(T, 32) bDet = ca * bdx * bdx + ca * bdy * bdy;
	const Expansion!(T, 32) cDet = ab * cdx * cdx + ab * cdy * cdy;
	const Expansion!(T, 96) fin1 = aDet + bDet + cDet;
	det = fin1.estimate();
	errBound = iccErrBoundB!T * permanent;
	if(abs(det) >= abs(errBound))
		return det;
	
	const T adxTail = ExpansionBase!T.minusTail(ax, dx, adx);
	const T adyTail = ExpansionBase!T.minusTail(ay, dy, ady);
	const T bdxTail = ExpansionBase!T.minusTail(bx, dx, bdx);
	const T bdyTail = ExpansionBase!T.minusTail(by, dy, bdy);
	const T cdxTail = ExpansionBase!T.minusTail(cx, dx, cdx);
	const T cdyTail = ExpansionBase!T.minusTail(cy, dy, cdy);
	if(
		T(0) == adxTail && T(0) == bdxTail && T(0) == cdxTail &&
		T(0) == adyTail && T(0) == bdyTail && T(0) == cdyTail
	)
		return det;
	
	errBound = iccErrBoundC!T * permanent + resultErrBound!T * abs(det);
	det += ((adx * adx + ady * ady) * ((bdx * cdyTail + cdy * bdxTail) - (bdy * cdxTail + cdx * bdyTail))
		+   (bdx * cdy - bdy * cdx) *  (adx * adxTail + ady * adyTail) * T(2))
		+  ((bdx * bdx + bdy * bdy) * ((cdx * adyTail + ady * cdxTail) - (cdy * adxTail + adx * cdyTail))
		+   (cdx * ady - cdy * adx) *  (bdx * bdxTail + bdy * bdyTail) * T(2))
		+  ((cdx * cdx + cdy * cdy) * ((adx * bdyTail + bdy * adxTail) - (ady * bdxTail + bdx * adyTail))
		+   (adx * bdy - ady * bdx) *  (cdx * cdxTail + cdy * cdyTail) * T(2));
	if(abs(det) >= abs(errBound))
		return det;
	return inCircleExact(ax, ay, bx, by, cx, cy, dx, dy);
}

T inCircleExact(T)(T ax, T ay, T bx, T by, T cx, T cy, T dx, T dy) nothrow @nogc pure @safe{
	const Expansion!(T, 4) ab = ExpansionBase!T.twoTwoDiff(ax, by, bx, ay);
	const Expansion!(T, 4) bc = ExpansionBase!T.twoTwoDiff(bx, cy, cx, by);
	const Expansion!(T, 4) cd = ExpansionBase!T.twoTwoDiff(cx, dy, dx, cy);
	const Expansion!(T, 4) da = ExpansionBase!T.twoTwoDiff(dx, ay, ax, dy);
	const Expansion!(T, 4) ac = ExpansionBase!T.twoTwoDiff(ax, cy, cx, ay);
	const Expansion!(T, 4) bd = ExpansionBase!T.twoTwoDiff(bx, dy, dx, by);
	
	const Expansion!(T, 12) abc = ab + bc - ac;
	const Expansion!(T, 12) bcd = bc + cd - bd;
	const Expansion!(T, 12) cda = cd + da + ac;
	const Expansion!(T, 12) dab = da + ab + bd;
	
	const Expansion!(T, 96) aDet = bcd * ax *  ax + bcd * ay *  ay;
	const Expansion!(T, 96) bDet = cda * bx * -bx + cda * by * -by;
	const Expansion!(T, 96) cDet = dab * cx *  cx + dab * cy *  cy;
	const Expansion!(T, 96) dDet = abc * dx * -dx + abc * dy * -dy;
	
	const Expansion!(T, 384) deter = (aDet + bDet) + (cDet + dDet);
	return deter.mostSignificant();
}

struct Expansion(T, size_t n){
	nothrow @nogc pure @safe:
	private{
		T[n] m_array;
		size_t m_size;
	}
	
	@property length() const =>
		m_size;
	@property empty() const =>
		m_size == 0;
	
	void opOpAssign(string op: "~")(const T v){
		m_array[m_size++] = v;
	}
	
	//estimates of expansion value
	@property mostSignificant() const =>
		empty ? T(0) : m_array[m_size-1];
	@property estimate() const =>
		sum(m_array[0..m_size]);
	
	Expansion!(T, n+m) opBinary(string op: "+", size_t m)(const Expansion!(T, m) f) const{
		Expansion!(T, n+m) h;
		auto hRef = h.m_array[];
		h.m_size = ExpansionBase!T.expansionSum(m_array[0..m_size], f.m_array[0..f.length], hRef);
		return h;
	}
	auto negate(){
		foreach(ref item; m_array){
			item = -item;
		}
		return this;
	}
	Expansion!(T, n+m) opBinary(string op: "-", size_t m)(Expansion!(T, m) f) const =>
		this.opBinary!"+"(f.negate());
	
	Expansion!(T, 2*n) opBinary(string op: "*")(T b) const{
		Expansion!(T, 2*n) h;
		auto hRef = h.m_array[];
		h.m_size = ExpansionBase!T.scaleExpansion(m_array[0..length], b, hRef);
		return h;
	}
}

template digits(T){
	static if(__traits(isFloating, T))
		enum digits = T.mant_dig;
	else
		enum digits = T.sizeof * 8;
}
enum epsilon(T) = cast(T)exp2(cast(double)-digits!T);
enum resultErrBound(T) = (T( 3) + T(  8) * epsilon!T) * epsilon!T;
enum ccwErrBoundA(T)   = (T( 3) + T( 16) * epsilon!T) * epsilon!T;
enum ccwErrBoundB(T)   = (T( 2) + T( 12) * epsilon!T) * epsilon!T;
enum ccwErrBoundC(T)   = (T( 9) + T( 64) * epsilon!T) * epsilon!T * epsilon!T;
enum iccErrBoundA(T)   = (T(10) + T( 96) * epsilon!T) * epsilon!T;
enum iccErrBoundB(T)   = (T( 4) + T( 48) * epsilon!T) * epsilon!T;
enum iccErrBoundC(T)   = (T(44) + T(576) * epsilon!T) * epsilon!T * epsilon!T;
