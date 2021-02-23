from geometry import Point3, create_cube, clip, Plane
from unfolding import unfold
import svg
import math

m = 5

cube = create_cube()

rest = cube

pieces = []

for i in range(m):
    theta = math.pi/6 * (1-2*i/(m-1))

    u = Point3(math.sin(theta),0,math.cos(theta))
    plane = Plane(u.x,u.y,u.z, -2*u.x)

    pieces.append(clip(rest, -plane))
    rest = clip(rest, plane)

pieces.append(rest)

for i,piece in enumerate(pieces):
    svg.write_unfolded_faces(f"test3-{i+1}.svg", unfold(piece))




