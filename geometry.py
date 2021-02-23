
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

    def __neg__(self, k):
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
        return Point(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x)
    
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

    def __neg__(self, k):
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


# create a cube 
def create_cube(r=1):
    vertices = [Point3(x,y,z) for (x,y,z) in [
        (-r,-r,-r),( r,-r,-r),(-r,r,-r),(r,r,-r),
        (-r,-r,r),( r,-r,r),(-r,r,r),(r,r,r)]]
    faces = [Polyhedron.Face(ii) for ii in [
        (0,1,3,2),(4,6,7,5),(1,0,4,5),
        (3,1,5,7),(2,3,7,6),(0,2,6,4)]]

    return Polyhedron(vertices, faces)



        
def clip(ph, plane):
    ws = [plane.comp(p) for p in ph.vertices]
    vertices = []
    faces = []

    vtb = {}
    etb = {}

    section_face_tb = {}
    
    def vp(i):
        if i in vtb: return vtb[i]
        j = len(vertices)
        vertices.append(ph.vertices[i])
        vtb[i] = j
        return j

    def ep(a,b):
        if a>b: a,b=b,a
        if (a,b) in etb: return etb[(a,b)]
        assert ws[a]*ws[b]<0.0
        pa = ph.vertices[a]
        pb = ph.vertices[b]
        t = -ws[a]/(ws[b]-ws[a])
        j = len(vertices)
        vertices.append(Point3.lerp(pa, pb, t))
        etb[(a,b)] = j
        return j
    
    for i,face in enumerate(ph.faces):
        m = face.m
        new_face = []
        new_a, new_b = None, None
        for j in range(m):
            a,b = face.indices[j], face.indices[(j+1)%m]
            if ws[a]<0:
                if ws[b]<0:
                    new_face.append(vp(b))
                else:
                    j = ep(a, b)
                    assert new_b == None
                    new_b = j
                    new_face.append(j)
            else:
                if ws[b]<0:
                    j = ep(a, b)
                    assert new_a == None
                    new_a = j
                    new_face.append(j)
                    new_face.append(vp(b))
                else:
                    pass
        assert new_face==[] or len(new_face)>=3
        if new_face != []:
            faces.append(Polyhedron.Face(new_face, face.color))
        assert (new_a is None) == (new_b is None)
        if new_a is not None:
            assert new_a not in section_face_tb
            section_face_tb[new_a] = new_b

    if len(section_face_tb) >= 3:
        section_face = [next(iter(section_face_tb))]
        while True:
            j = section_face_tb.get(section_face[-1],None)
            if j == None: break
            if j == section_face[0]: break
            section_face.append(j)
        if j == section_face[0]:
            faces.append(Polyhedron.Face(section_face, 2))
        
    result = Polyhedron(vertices, faces)
    return result

    
    
