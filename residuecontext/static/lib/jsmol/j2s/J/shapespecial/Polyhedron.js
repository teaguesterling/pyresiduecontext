Clazz.declarePackage ("J.shapespecial");
Clazz.load (null, "J.shapespecial.Polyhedron", ["java.lang.Boolean", "$.Float", "java.util.Hashtable", "JU.AU", "$.BS", "$.Measure", "$.P3", "$.V3", "JM.Atom", "JS.SV", "JU.C", "$.Elements", "$.Escape", "$.Logger", "$.Node", "$.Normix", "$.Point3fi"], function () {
c$ = Clazz.decorateAsClass (function () {
this.info = null;
this.id = null;
this.center = null;
this.centralAtom = null;
this.vertices = null;
this.triangles = null;
this.faces = null;
this.nVertices = 0;
this.collapsed = false;
this.bsFlat = null;
this.distanceRef = 0;
this.normals = null;
this.normixes = null;
this.smiles = null;
this.smarts = null;
this.polySmiles = null;
this.pointGroup = null;
this.pointGroupFamily = null;
this.volume = null;
this.visible = true;
this.isFullyLit = false;
this.isValid = true;
this.colixEdge = 0;
this.visibilityFlags = 0;
this.colix = 23;
this.modelIndex = -2147483648;
this.offset = null;
this.scale = 1;
this.pointScale = 0;
this.elemNos = null;
Clazz.instantialize (this, arguments);
}, J.shapespecial, "Polyhedron");
Clazz.makeConstructor (c$, 
function () {
});
Clazz.defineMethod (c$, "set", 
function (id, modelIndex, atomOrPt, points, nPoints, vertexCount, triangles, triangleCount, faces, normals, bsFlat, collapsed, distanceRef, pointScale) {
this.pointScale = pointScale;
this.distanceRef = distanceRef;
if (id == null) {
this.centralAtom = atomOrPt;
this.modelIndex = this.centralAtom.mi;
} else {
this.id = id;
this.center = atomOrPt;
this.modelIndex = modelIndex;
}this.nVertices = vertexCount;
this.vertices =  new Array (nPoints + 1);
this.normals =  new Array (triangleCount);
this.faces = faces;
this.bsFlat = bsFlat;
this.triangles = JU.AU.newInt2 (triangleCount);
for (var i = nPoints + 1; --i >= 0; ) this.vertices[i] = points[i];

for (var i = triangleCount; --i >= 0; ) this.normals[i] = JU.V3.newV (normals[i]);

for (var i = triangleCount; --i >= 0; ) this.triangles[i] = triangles[i];

this.collapsed = collapsed;
return this;
}, "~S,~N,JU.P3,~A,~N,~N,~A,~N,~A,~A,JU.BS,~B,~N,~N");
Clazz.defineMethod (c$, "setInfo", 
function (info, at) {
try {
this.collapsed = info.containsKey ("collapsed");
this.id = (info.containsKey ("id") ? info.get ("id").asString () : null);
if (this.id == null) {
this.centralAtom = at[info.get ("atomIndex").intValue];
this.modelIndex = this.centralAtom.mi;
} else {
this.center = JU.P3.newP (JS.SV.ptValue (info.get ("center")));
this.modelIndex = info.get ("modelIndex").intValue;
this.colix = JU.C.getColixS (info.get ("color").asString ());
this.colixEdge = JU.C.getColixS (info.get ("colorEdge").asString ());
if (info.containsKey ("offset")) this.offset = JU.P3.newP (JS.SV.ptValue (info.get ("offset")));
if (info.containsKey ("scale")) this.scale = JS.SV.fValue (info.get ("scale"));
}var lst = info.get ("vertices").getList ();
var vc = info.get ("vertexCount");
if (vc == null) {
this.nVertices = lst.size ();
this.vertices =  new Array (this.nVertices + 1);
this.vertices[this.nVertices] = JS.SV.ptValue (info.get ("ptRef"));
} else {
this.nVertices = vc.intValue;
this.vertices =  new Array (lst.size ());
vc = info.get ("r");
if (vc != null) this.distanceRef = vc.asFloat ();
}for (var i = lst.size (); --i >= 0; ) this.vertices[i] = JS.SV.ptValue (lst.get (i));

lst = info.get ("elemNos").getList ();
for (var i = this.nVertices; --i >= 0; ) {
var n = lst.get (i).intValue;
if (n > 0) {
var p =  new JU.Point3fi ();
p.setT (this.vertices[i]);
p.sD = n;
this.vertices[i] = p;
}}
if (info.containsKey ("pointScale")) this.pointScale = Math.max (0, JS.SV.fValue (info.get ("pointScale")));
var faces = info.get ("faces");
var o = info.get ("triangles");
if (o == null) {
o = faces;
} else {
this.faces = this.toInt2 (faces);
}this.triangles = this.toInt2 (o);
this.normals =  new Array (this.triangles.length);
var vAB =  new JU.V3 ();
for (var i = this.triangles.length; --i >= 0; ) {
this.normals[i] =  new JU.V3 ();
var a = this.triangles[i];
JU.Measure.getNormalThroughPoints (this.vertices[a[0]], this.vertices[a[1]], this.vertices[a[2]], this.normals[i], vAB);
}
this.bsFlat = JS.SV.getBitSet (info.get ("bsFlat"), false);
} catch (e) {
if (Clazz.exceptionOf (e, Exception)) {
return null;
} else {
throw e;
}
}
return this;
}, "java.util.Map,~A");
Clazz.defineMethod (c$, "toInt2", 
 function (o) {
var lst = o.getList ();
var ai = JU.AU.newInt2 (lst.size ());
for (var i = ai.length; --i >= 0; ) {
var lst2 = lst.get (i).getList ();
var a = ai[i] =  Clazz.newIntArray (lst2.size (), 0);
for (var j = a.length; --j >= 0; ) a[j] = lst2.get (j).intValue;

}
return ai;
}, "JS.SV");
Clazz.defineMethod (c$, "getInfo", 
function (vwr, isState) {
if (!isState && this.info != null && !JU.Logger.debugging) return this.info;
var info =  new java.util.Hashtable ();
info.put ("vertexCount", Integer.$valueOf (this.nVertices));
var nv = (isState ? this.vertices.length : this.nVertices);
var pts =  new Array (nv);
for (var i = 0; i < nv; i++) pts[i] = JU.P3.newP (this.vertices[i]);

info.put ("vertices", pts);
info.put ("elemNos", this.getElemNos ());
if (this.id == null) {
info.put ("atomIndex", Integer.$valueOf (this.centralAtom.i));
} else {
info.put ("id", this.id);
info.put ("center", JU.P3.newP (this.center));
info.put ("color", JU.C.getHexCode (this.colix));
info.put ("colorEdge", JU.C.getHexCode (this.colixEdge == 0 ? this.colix : this.colixEdge));
if (this.offset != null) info.put ("offset", this.offset);
if (this.scale != 1) info.put ("scale", Float.$valueOf (this.scale));
}if (this.id != null || !isState) info.put ("modelIndex", Integer.$valueOf (this.modelIndex));
if (!isState) {
this.info = info;
if (this.id == null) {
info.put ("center", JU.P3.newP (this.centralAtom));
info.put ("modelNumber", Integer.$valueOf (this.centralAtom.getModelNumber ()));
info.put ("atomNumber", Integer.$valueOf (this.centralAtom.getAtomNumber ()));
info.put ("atomName", this.centralAtom.getInfo ());
info.put ("element", this.centralAtom.getElementSymbol ());
var energy = vwr.ms.getInfo (this.centralAtom.mi, "Energy");
if (energy != null) info.put ("energy", energy);
}info.put ("triangleCount", Integer.$valueOf (this.triangles.length));
info.put ("volume", this.getVolume ());
var names =  new Array (this.nVertices);
var indices =  Clazz.newIntArray (this.nVertices, 0);
for (var i = this.nVertices; --i >= 0; ) {
var pt = this.vertices[i];
var isNode = Clazz.instanceOf (pt, JU.Node);
names[i] = (isNode ? (pt).getAtomName () : Clazz.instanceOf (pt, JU.Point3fi) ? JU.Elements.elementSymbolFromNumber ((pt).sD) : "");
indices[i] = (isNode ? (pt).getIndex () : -1);
}
info.put ("atomNames", names);
info.put ("vertexIndices", indices);
if (this.faces != null) info.put ("faceCount", Integer.$valueOf (this.faces.length));
if (this.smarts != null) info.put ("smarts", this.smarts);
if (this.smiles != null) info.put ("smiles", this.smiles);
if (this.polySmiles != null) info.put ("polySmiles", this.polySmiles);
if (this.pointGroup != null) info.put ("pointGroup", this.pointGroup.getPointGroupName ());
if (this.pointGroupFamily != null) info.put ("pointGroupFamily", this.pointGroupFamily.getPointGroupName ());
}if (this.pointScale > 0) info.put ("pointScale", Float.$valueOf (this.pointScale));
if (this.faces != null) info.put ("faces", this.faces);
if (isState || JU.Logger.debugging) {
info.put ("bsFlat", this.bsFlat);
if (this.collapsed) info.put ("collapsed", Boolean.$valueOf (this.collapsed));
if (this.distanceRef != 0) info.put ("r", Float.$valueOf (this.distanceRef));
var n =  new Array (this.normals.length);
for (var i = n.length; --i >= 0; ) n[i] = JU.P3.newP (this.normals[i]);

info.put ("normals", n);
info.put ("triangles", JU.AU.arrayCopyII (this.triangles, this.triangles.length));
}return info;
}, "JV.Viewer,~B");
Clazz.defineMethod (c$, "getElemNos", 
function () {
if (this.elemNos == null) {
this.elemNos =  Clazz.newIntArray (this.nVertices, 0);
for (var i = 0; i < this.nVertices; i++) {
var pt = this.vertices[i];
this.elemNos[i] = (Clazz.instanceOf (pt, JU.Node) ? (pt).getElementNumber () : Clazz.instanceOf (pt, JU.Point3fi) ? (pt).sD : -2);
}
}return this.elemNos;
});
Clazz.defineMethod (c$, "getSymmetry", 
function (vwr, withPointGroup) {
if (this.id == null) {
if (this.smarts == null) {
this.info = null;
var sm = vwr.getSmilesMatcher ();
try {
var details = (this.distanceRef <= 0 ? null : "r=" + this.distanceRef);
this.smarts = sm.polyhedronToSmiles (this.centralAtom, this.faces, this.nVertices, null, 8192, null);
this.smiles = sm.polyhedronToSmiles (this.centralAtom, this.faces, this.nVertices, this.vertices, 1, null);
this.polySmiles = sm.polyhedronToSmiles (this.centralAtom, this.faces, this.nVertices, this.vertices, 196609, details);
} catch (e) {
if (Clazz.exceptionOf (e, Exception)) {
} else {
throw e;
}
}
}}if (!withPointGroup) return null;
if (this.pointGroup == null) {
var pts =  new Array (this.nVertices);
for (var i = pts.length; --i >= 0; ) pts[i] = this.vertices[i];

this.pointGroup = vwr.ms.getSymTemp (true).setPointGroup (null, null, pts, null, false, vwr.getFloat (570425382), vwr.getFloat (570425384), true);
for (var i = pts.length; --i >= 0; ) pts[i] = JU.P3.newP (this.vertices[i]);

this.pointGroupFamily = vwr.ms.getSymTemp (true).setPointGroup (null, null, pts, null, false, vwr.getFloat (570425382), vwr.getFloat (570425384), true);
}return (this.center == null ? this.centralAtom : this.center) + "    \t" + this.pointGroup.getPointGroupName () + "\t" + this.pointGroupFamily.getPointGroupName ();
}, "JV.Viewer,~B");
Clazz.defineMethod (c$, "getVolume", 
 function () {
if (this.volume != null) return this.volume;
var vAB =  new JU.V3 ();
var vAC =  new JU.V3 ();
var vTemp =  new JU.V3 ();
var v = 0;
if (this.bsFlat.cardinality () < this.triangles.length) for (var i = this.triangles.length; --i >= 0; ) {
var face = this.triangles[i];
v += this.triangleVolume (face[0], face[1], face[2], vAB, vAC, vTemp);
}
return Float.$valueOf (v / 6);
});
Clazz.defineMethod (c$, "triangleVolume", 
 function (i, j, k, vAB, vAC, vTemp) {
vAB.setT (this.vertices[i]);
vAC.setT (this.vertices[j]);
vTemp.cross (vAB, vAC);
vAC.setT (this.vertices[k]);
return vAC.dot (vTemp);
}, "~N,~N,~N,JU.V3,JU.V3,JU.V3");
Clazz.defineMethod (c$, "getState", 
function (vwr) {
var ident = (this.id == null ? "({" + this.centralAtom.i + "})" : "ID " + JU.Escape.e (this.id));
return "  polyhedron @{" + JU.Escape.e (this.getInfo (vwr, true)) + "} " + (this.isFullyLit ? " fullyLit" : "") + ";" + (this.visible ? "" : "polyhedra " + ident + " off;") + "\n";
}, "JV.Viewer");
Clazz.defineMethod (c$, "move", 
function (mat, bsMoved) {
this.info = null;
for (var i = 0; i < this.nVertices; i++) {
var p = this.vertices[i];
if (Clazz.instanceOf (p, JM.Atom)) {
if (bsMoved.get ((p).i)) continue;
p = this.vertices[i] = JU.P3.newP (p);
}mat.rotTrans (p);
}
for (var i = this.normals.length; --i >= 0; ) mat.rotate (this.normals[i]);

this.normixes = null;
}, "JU.M4,JU.BS");
Clazz.defineMethod (c$, "getNormixes", 
function () {
if (this.normixes == null) {
this.normixes =  Clazz.newShortArray (this.normals.length, 0);
var bsTemp =  new JU.BS ();
for (var i = this.normals.length; --i >= 0; ) this.normixes[i] = (this.bsFlat.get (i) ? JU.Normix.get2SidedNormix (this.normals[i], bsTemp) : JU.Normix.getNormixV (this.normals[i], bsTemp));

}return this.normixes;
});
Clazz.defineMethod (c$, "setOffset", 
function (value) {
if (this.center == null) return;
var v = JU.P3.newP (value);
if (this.offset != null) value.sub (this.offset);
this.offset = v;
this.center.add (value);
for (var i = this.vertices.length; --i >= 0; ) this.vertices[i].add (value);

}, "JU.P3");
});
