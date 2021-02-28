from geometry import *
from unfolding import unfold
import svg


print("tet")
tetrahedron = create_tetrahedron()
print("cube")
cube = create_cube()
print("octahedron")
octahedron = create_octahedron()
print("icosahedron")
icosahedron = create_icosahedron()

svg.write_unfolded_faces('tetrahedron.svg', unfold(tetrahedron))
svg.write_unfolded_faces('cube.svg', unfold(cube))
svg.write_unfolded_faces('octahedron.svg', unfold(octahedron))
svg.write_unfolded_faces('icosahedron.svg', unfold(icosahedron))

ph = clip(icosahedron, Plane(0,1,0.05,0))

svg.write_unfolded_faces('icosahedron_1.svg', unfold(ph))

# plane = Plane(1,1,1,0)

#svg.write_unfolded_faces('test1a.svg', unfold(clip(cube,plane)))
#svg.write_unfolded_faces('test1b.svg', unfold(clip(cube,-plane)))




