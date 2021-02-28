
class Point3:
    """A point in three-dimensional space."""
    def __init__(self, x=0, y=0, z=0):
        self.x=x
        self.y=y
        self.z=z

    def __repr__(self):
        return f"({self.x},{self.y},{self.z})"

    def __str__(self):
        return self.__repr__()

    def __add__(self, other):
        return Point3(self.x+other.x, self.y+other.y, self.z+other.z)

    def __sub__(self, other):
        return Point3(self.x-other.x, self.y-other.y, self.z-other.z)

    def __mul__(self, k):
        return Point3(k*self.x, k*self.y, k*self.z)

    def __neg__(self):
        return Point3(-self.x, -self.y, -self.z)

    def length(self):
        return (self.x**2+self.y**2+self.z**2)**0.5

    def normalize(self):
        s = 1/self.length()
        return Point3(s*self.x, s*self.y, s*self.z)

    def tuple(self):
        return (self.x, self.y, self.z)
    
    @staticmethod
    def dot(a,b):
        return a.x*b.x + a.y*b.y + a.z*b.z

    @staticmethod
    def cross(a,b):
        return Point3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x)
    
    @staticmethod
    def lerp(a,b,t):
        it = 1.0-t
        return Point3(it*a.x+t*b.x, it*a.y+t*b.y, it*a.z+t*b.z)
        
class Point2:
    """A point in two-dimensional space."""
    def __init__(self, x=0, y=0):
        self.x=x
        self.y=y

    def __repr__(self):
        return f"({self.x},{self.y})"

    def __str__(self):
        return self.__repr__()

    def __add__(self, other):
        return Point2(self.x+other.x, self.y+other.y)

    def __sub__(self, other):
        return Point2(self.x-other.x, self.y-other.y)

    def __mul__(self, k):
        return Point2(k*self.x, k*self.y)

    def length(self):
        return (self.x**2+self.y**2)**0.5

    def normalize(self):
        s = 1/self.length()
        return Point2(s*self.x, s*self.y)

    def __neg__(self):
        return Point2(-self.x, -self.y)

    def rotate90(self):
        return Point2(-self.y, self.x)

    def rotate270(self):
        return Point2(self.y, -self.x)

    def tuple(self):
        return (self.x, self.y)
    
    @staticmethod
    def dot(a,b):
        return a.x*b.x + a.y*b.y

    @staticmethod
    def lerp(a,b,t):
        it = 1.0-t
        return Point2(it*a.x+t*b.x, it*a.y+t*b.y)
    

class Plane:
    """A plane three-dimensional space; set of (x,y,z): a*x+b*y+c*z+d==0."""
    def __init__(self, a,b,c,d):
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    def comp(self, p):
        return self.a*p.x+self.b*p.y+self.c*p.z+self.d
        
    def __neg__(self):
        return Plane(-self.a, -self.b, -self.c, -self.d) 


class Polyhedron:
    """A polyhedron."""
    
    class Face:
        """Polyhedron face: a list of m indices and a color."""
        def __init__(self, indices, color=1):
            self.color = color
            self.indices = indices
            self.m = len(indices)
        def __repl__(self):
            ii = ",".join([str(i) for i in self.indices])
            return f"face({ii})"
        def __str__(self):
            return self.__repl__()
            
    def __init__(self, vertices, faces):
        self.vertices = vertices
        self.faces = faces

    def get_face_center(self, face_index):
        assert 0<=face_index<len(self.faces)
        face = self.faces[face_index]
        sc = 1/face.m
        p = sum([self.vertices[i] for i in face.indices], Point3(0,0,0))*sc
        return p

    def get_face_by_edge(self, a, b): # return (face-index, edge-index)
        for i,face in enumerate(self.faces):
            indices = face.indices
            for j in range(face.m):
                c,d = indices[j], indices[(j+1)%face.m]
                if (c,d) == (a,b):
                    return (i,j)
        return None,None

    def assert_integrity(self):
        edges = set()
        for face in self.faces:
            m = face.m
            assert len(face.indices) == m
            assert m>=3
            assert len(set(face.indices)) == m
            for i in range(m):
                a, b = face.indices[i], face.indices[(i+1)%m]
                edge = a,b
                assert edge not in edges, "Polyhedron with duplicated edges or misoriented faces"
                edges.add(edge)
        for a,b in edges:
            # assert (b,a) in edges, "Polyhedron with holes"
            pass

        for face in self.faces:
            pts = [self.vertices[j] for j in face.indices]
            fc = sum(pts,Point3(0,0,0)) * (1.0/len(pts))
            e0 = (pts[0]-fc).normalize()
            e1 = pts[1]-fc
            e1 = (e1-e0*Point3.dot(e0,e1)).normalize()
            e2 = Point3.cross(e0,e1).normalize()
            assert Point3.dot(Point3(0,0,0) - fc, e2) > 0.0
            for i,p in enumerate(self.vertices):
                w = Point3.dot(p-fc, e2)
                assert w >= -1.0e-10, "Bad face orientation"


# create a cube 
def create_tetrahedron():
    u = 1
    vertices = [Point3(x,y,z) for (x,y,z) in [
        ( u, u, u),(-u,-u, u),(-u, u,-u),( u,-u,-u)]]
    faces = [Polyhedron.Face(ii) for ii in [
        (0,1,2),(1,0,3),(2,1,3),(0,2,3)]]

    ph = Polyhedron(vertices, faces)
    ph.assert_integrity()
    return ph

def create_cube():
    u = 1
    vertices = [Point3(x,y,z) for (x,y,z) in [
        (-u,-u,-u),( u,-u,-u),(-u, u,-u),(u, u,-u),
        (-u,-u, u),( u,-u, u),(-u, u, u),(u, u, u)]]
    faces = [Polyhedron.Face(ii) for ii in [
        (0,1,3,2),(4,6,7,5),(1,0,4,5),
        (3,1,5,7),(2,3,7,6),(0,2,6,4)]]

    ph = Polyhedron(vertices, faces)
    ph.assert_integrity()
    return ph

def create_octahedron():
    u = 1
    vertices = [Point3(x,y,z) for (x,y,z) in [
        ( 0, u, 0),( u, 0, 0),(-u, 0, 0),
        ( 0, 0, u),( 0, 0,-u),( 0,-u, 0)]]
    faces = [Polyhedron.Face(ii) for ii in [
        (0,1,3),(0,3,2),(0,2,4),(0,4,1),
        (5,3,1),(5,2,3),(5,4,2),(5,1,4)]]

    ph = Polyhedron(vertices, faces)
    ph.assert_integrity()
    return ph

def create_icosahedron():
    u = 1
    f = (-1+5**0.5)/2
    vertices = [Point3(x,y,z) for (x,y,z) in [
        (-f, u, 0),( f, u, 0),(-f,-u, 0),( f,-u, 0),
        (-1, 0, f),(-1, 0,-f),( 1, 0, f),( 1, 0,-f),
        ( 0,-f, 1),( 0, f, 1),( 0,-f,-1),( 0, f,-1)]]
    faces = [Polyhedron.Face(ii) for ii in [
        (9,0,1),(9,1,6),(9,6,8),(9,8,4),(9,4,0),
        (8,6,3),(8,3,2),(8,2,4),(2,3,10),(0,11,1),
        (1,11,7),(6,1,7),(3,6,7),(10,3,7),(11,10,7),
        (5,11,0),(5,0,4),(5,4,2),(5,2,10),(5,10,11)
        ]]

    ph = Polyhedron(vertices, faces)
    ph.assert_integrity()
    return ph



class Clipper:
    """The class cuts polyhedron across a plane."""


    def clip(self, ph, plane):
        """Main method. Returns the (possibly empty) polyhedron part below the plane."""
        self.ph = ph
        self.plane = plane
        self.ws = [plane.comp(p) for p in ph.vertices]
        epsilon = 1.0e-8
        self.ws = [w if abs(w)>epsilon else 0.0 for w in self.ws]
        self.vertices = []
        self.faces = []

        self.vtb = {} # old-vertex => new-vertex
        self.etb = {} # old-edge(a,b) => new-vertex

        self.section_face_tb = {} # 

        for i,face in enumerate(ph.faces):
            self._process_face(i,face)

        self._add_section_face()
        result = Polyhedron(self.vertices, self.faces)
        return result


    def _vp(self, i):
        if i in self.vtb: return self.vtb[i]
        j = len(self.vertices)
        self.vertices.append(self.ph.vertices[i])
        self.vtb[i] = j
        return j

    def _ep(self, a,b):
        if a>b: a,b=b,a
        edge = (a,b)
        wa = self.ws[a]
        wb = self.ws[b]        
        assert wa * wb < 0.0
        if edge in self.etb: return self.etb[edge]
        pa = self.ph.vertices[a]
        pb = self.ph.vertices[b]
        t = -wa/(wb-wa)
        j = len(self.vertices)
        self.vertices.append(Point3.lerp(pa, pb, t))
        self.etb[edge] = j
        return j

    def _process_zero_face(self, i, face):
        new_face = [self._vp(j) for j in face.indices]
        self.faces.append(Polyhedron.Face(new_face, face.color))


    def _process_face(self, i, face):
        m = face.m
        ws = self.ws
        face_ws = [ws[j] for j in face.indices]
        zm,inm,outm = 0,0,0
        for w in  [ws[j] for j in face.indices]:
            if w>0.0: outm += 1
            elif w<0.0: inm += 1
            else: zm += 0
        if zm >= 3:
            # face belong to cutting plane => pass it
            assert inm == outm == 0
            self._process_zero_face(i, face)
            return

        new_face = []
        new_a, new_b = None, None
        for j in range(m):
            a,b = face.indices[j], face.indices[(j+1)%m]

            if ws[a]<0:
                
                if ws[b]<0:
                    j = self._vp(b)
                    new_face.append(j)
                elif ws[b]>0:
                    j = self._ep(a, b)
                    new_face.append(j)
                    assert new_b == None
                    new_b = j
                else: # ws[b]==0
                    j = self._vp(b)
                    new_face.append(j)
                    assert new_b == None
                    new_b = j
                    
            elif ws[a]>0:                
                
                if ws[b]<0:
                    j = self._ep(a, b) 
                    new_face.append(j)
                    jb = self._vp(b)
                    new_face.append(jb)
                    assert new_a == None
                    new_a = j
                elif ws[b]==0:
                    jb = self._vp(b)
                    new_face.append(jb)
                    assert new_a == None
                    new_a = jb
            
            else: # ws[a] == 0

                if ws[b]<0:
                    j = self._vp(b) 
                    new_face.append(j)
                elif ws[b] == 0:
                    # segment a,b belongs to the cutting plane
                    if inm > 0:
                        j = self._vp(b) 
                        new_face.append(j)
                        
                        # assert new_b == None
                        new_a = j

        if new_face != [] and len(new_face)>=3:
            self.faces.append(Polyhedron.Face(new_face, face.color))
            # print("  new face : ", ", ".join([str(self.vertices[j]) for j in new_face]))

        # assert (new_a is None) == (new_b is None)
        # print("  new edge : ", new_a, new_b)
        if new_a is not None and new_b is not None and new_a != new_b:
            assert new_a not in self.section_face_tb
            self.section_face_tb[new_a] = new_b

    def _add_section_face(self):
        # print("section face: ", self.section_face_tb)
        if len(self.section_face_tb) >= 3:
            section_face = [next(iter(self.section_face_tb))]
            while True:
                j = self.section_face_tb.get(section_face[-1],None)
                if j == None: break
                if j == section_face[0]: break
                section_face.append(j)
            if j == section_face[0]:
                self.faces.append(Polyhedron.Face(section_face, 2))
        

def clip(ph, plane):
    clipper = Clipper()
    return clipper.clip(ph,plane)
    
