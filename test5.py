from geometry import create_cube, clip, Plane
from unfolding import unfold
import svg



cube = create_cube()
cube.assert_integrity()
plane = Plane(1,1,0,0)

ph = clip(cube,plane)
ph.assert_integrity()

svg.write_unfolded_faces('test5.svg', unfold(ph))




