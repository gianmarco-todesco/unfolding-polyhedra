from geometry import Point3, create_cube, clip, Plane
from unfolding import unfold
import svg
import math

m = 5

cube = create_cube()

rest = cube

pieces = []

phi = math.pi/4
cs_phi = math.cos(phi)
sn_phi = math.sin(phi)

for i in range(m):
    theta = math.pi/6 * (1-2*i/(m-1))

    u = Point3(math.sin(theta),0,math.cos(theta))
    u.x,u.y = cs_phi*u.x, sn_phi*u.x
    plane = Plane(u.x,u.y,u.z, -2*(u.x+u.y))

    pieces.append(clip(rest, -plane))
    rest = clip(rest, plane)

pieces.append(rest)

for i,piece in enumerate(pieces):
    svg.write_unfolded_faces(f"test4-{i+1}.svg", unfold(piece))




