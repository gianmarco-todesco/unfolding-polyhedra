from geometry import create_cube, clip, Plane
from unfolding import unfold
import svg


cube = create_cube()
plane1 = Plane(0,0,1, 1/3)
plane2 = Plane(0,0,1,-1/3)

ph12 = clip(cube,plane2)
ph1 = clip(ph12,plane1)
ph2 = clip(ph12,-plane1)
ph3 = clip(cube,-plane2)


svg.write_unfolded_faces('test2a.svg', unfold(ph1))
svg.write_unfolded_faces('test2b.svg', unfold(ph2))
svg.write_unfolded_faces('test2c.svg', unfold(ph3))




