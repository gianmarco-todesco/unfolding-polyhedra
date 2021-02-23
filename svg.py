import svgwrite
from geometry import *


def write_unfolded_faces(filename, unfolded_faces):
    unit = 100
    mrg = 50
    inf = 10**9
    xmin,ymin,xmax,ymax=inf,inf,-inf,-inf
    for face in unfolded_faces:
        for v in face.vertices:
            x = v.p.x * unit
            y = v.p.y * unit
            if x < xmin: xmin = x
            elif x > xmax: xmax = x
            if y < ymin: ymin = y
            elif y > ymax: ymax = y
    lx = xmax-xmin+2*mrg
    ly = ymax-ymin+2*mrg
    offset = Point2(mrg-xmin, mrg-ymin)

    dwg = svgwrite.Drawing(filename, width=lx, height=ly)
    for face in unfolded_faces:
        pts = [(v.p*unit + offset).tuple() for v in face.vertices]
        color = "rgb(200,230,200)" if face.face_color==1 else "rgb(200,100,200)"
        dwg.add(dwg.polygon(
                points=pts,
                fill=color))

    for index, face in enumerate(unfolded_faces):
        pts = [(v.p*unit + offset).tuple() for v in face.vertices]
        pts.append(pts[0])
        if index>0:
            pts=pts[1:]
        dwg.add(dwg.polyline(
                points=pts,
                fill="none",
                stroke='black',
                stroke_width=1,
                stroke_linejoin="round",
                stroke_linecap="round")) 
    dwg.save()
    print(f"written {filename}")

