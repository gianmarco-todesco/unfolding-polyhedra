
from geometry import Point3, Point2


class UnfoldedFace:
    """Unfolded face: 2d-points and some info about the original polyhedron."""

    class Vertex:
        def __init__(self, p, index):
            self.p = p
            self.index = index
            
    def __init__(self, face_index, face_color, vertices):
        self.face_index = face_index
        self.face_color = face_color
        self.vertices = vertices
        self.m = len(vertices)
        

        
def unfold_first_face(ph, face_index, offset=Point2(0,0)):
    face = ph.faces[face_index]
    pts = [ph.vertices[i] for i in face.indices]
    face_center = ph.get_face_center(face_index)
    e0 = (pts[0] - face_center).normalize()
    e1 = pts[1] - face_center
    e1 = (e1 - e0*Point3.dot(e0,e1)).normalize()
    vv = []
    for i in range(face.m):
        p = pts[i] - face_center
        p2 = Point2(Point3.dot(p,e0), Point3.dot(p,e1)) + offset
        vv.append(UnfoldedFace.Vertex(p2, face.indices[i]))

    return UnfoldedFace(face_index, face.color, vv)

    
def unfold_face(ph, face_index, edge_index, p0, p1):
    face = ph.faces[face_index]
    ii = [face.indices[(edge_index+i)%face.m] for i in range(face.m)]
    pts = [ph.vertices[i] for i in ii]
    q0 = pts[0]
    e0 = (pts[1] - q0).normalize()
    e1 = pts[2] - q0
    e1 = (e1 - e0*Point3.dot(e0,e1)).normalize()
    
    f0 = (p1-p0).normalize()
    f1 = f0.rotate90()

    vv = []
    for i in range(face.m):
        p = pts[i] - q0
        u,v = Point3.dot(p,e0), Point3.dot(p,e1)
        p2 = p0 + f0*u + f1*v
        vv.append(UnfoldedFace.Vertex(p2, ii[i]))

    return UnfoldedFace(face_index, face.color, vv)
        
    
def unfold(ph, first_face=0, center = Point2(0,0)):
    if ph.faces == []:
        return []
    assert 0<=first_face<len(ph.faces)

    result = []
    visited = set()
    visited.add(first_face)
        
    unfolded_face = unfold_first_face(ph, first_face, center)
    result.append(unfolded_face)

    todo = [unfolded_face]
    while todo != []:
        unfolded_face = todo.pop()
        vv = unfolded_face.vertices
        m = unfolded_face.m
        for i in range(m):
            
            a,b = vv[i].index, vv[(i+1)%m].index
            fi,ei = ph.get_face_by_edge(b,a)
            if fi == None or fi in visited: continue
            face_fi = ph.faces[fi]
            assert face_fi.indices[ei]==b
            assert face_fi.indices[(ei+1)%face_fi.m]==a
            
            visited.add(fi)
            
            p0 = vv[i].p
            p1 = vv[(i+1)%m].p
            new_unfolded_face = unfold_face(ph, fi,ei, p1, p0)
            assert new_unfolded_face.vertices[0].index == b
            assert new_unfolded_face.vertices[1].index == a
            
            result.append(new_unfolded_face)
            todo.append(new_unfolded_face)
        

    return result

    