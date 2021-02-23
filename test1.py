from geometry import create_cube, clip, Plane
from unfolding import unfold
import svg



cube = create_cube()
plane = Plane(1,1,1,0)

svg.write_unfolded_faces('test1a.svg', unfold(clip(cube,plane)))
svg.write_unfolded_faces('test1b.svg', unfold(clip(cube,-plane)))




