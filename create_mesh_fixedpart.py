import numpy as np
import pygmsh

inflow_marker = 1
outflow_marker = 2
wall_marker = 3
obstacle_marker = 4
# resolution
resolution = 0.2 #1 # 0.005 #0.1

# geometric properties
L = 10 #2.5 #20            # length of channel
H = 6 #0.4 #6           # heigth of channel
c = [5, 3, 0]  #[0.2, 0.2, 0] #[10, 3, 0]  # position of object
r = 0.5 #0.05 #0.5 # radius of object


def create_stokes_mesh(res):
    """
    Creates a mesh containing a circular obstacle.
    Inlet marker (left): 1
    Oulet marker (right): 2
    Wall marker (top/bottom): 3
    Obstacle marker: 4
    """
    try:
        import gmsh
        import meshio

    except ImportError:
        print("meshio and/or gmsh not installed. Requires the non-python libraries:\n",
              "- libglu1\n - libxcursor-dev\n - libxinerama1\n And Python libraries:\n"
              " - h5py",
              " (pip3 install --no-cache-dir --no-binary=h5py h5py)\n",
              "- gmsh \n - meshio")
        exit(1)

    # Initialize empty geometry using the build in kernel in GMSH
    geometry = pygmsh.geo.Geometry()
    # Fetch model we would like to add data to
    model = geometry.__enter__()
    # Add circle
    pc = model.add_point(c)
    sin = 0.5  # sin(30°)
    cos = np.sqrt(3) / 2  # cos(30°)
    pc0 = model.add_point(c)
    pc1 = model.add_point((c[0] - r, c[1], 0), mesh_size=0.1 * resolution)
    pc2 = model.add_point((c[0] + cos * r, c[1] + sin * r, 0), mesh_size=0.1 * resolution)
    pc3 = model.add_point((c[0] + cos * r, c[1] - sin * r, 0), mesh_size=0.1 * resolution)
    circle1 = model.add_circle_arc(pc2, pc0, pc1)
    circle2 = model.add_circle_arc(pc1, pc0, pc3)
    circle3 = model.add_circle_arc(pc3, pc0, pc2)
    circle = model.add_curve_loop([circle1, circle2, circle3])

    # Add points with finer resolution on left side
    points = [model.add_point((0, 0, 0), mesh_size=resolution),
              model.add_point((L, 0, 0), mesh_size=resolution),  # 5*resolution
              model.add_point((L, H, 0), mesh_size=resolution),  # 5*resolution
              model.add_point((0, H, 0), mesh_size=resolution)]

    # Add lines between all points creating the rectangle
    channel_lines = [model.add_line(points[i], points[i + 1])
                     for i in range(-1, len(points) - 1)]

    # Create a line loop and plane surface for meshing
    channel_loop = model.add_curve_loop(channel_lines)
    plane_surface = model.add_plane_surface(
        channel_loop, holes=[circle])

    # Call gmsh kernel before add physical entities
    model.synchronize()

    volume_marker = 6
    model.add_physical([channel_lines[0]], "inflow")  # mark inflow boundary with 1
    model.add_physical([channel_lines[2]], "outflow")  # mark outflow boundary with 2
    model.add_physical([channel_lines[1], channel_lines[3], circle3], "walls")  # mark walls with 3
    model.add_physical([circle1, circle2], "obstacle")  # mark obstacle with 4
    model.add_physical([plane_surface], "Volume")

    geometry.generate_mesh(dim=2)
    import gmsh
    gmsh.write("mesh.msh")
    gmsh.clear()
    geometry.__exit__()

    # Read in and convert mesh to msh
    msh = meshio.read("mesh.msh")

    line_cells = []
    for cell in msh.cells:
        if cell.type == "triangle":
            triangle_cells = cell.data
        elif cell.type == "line":
            if len(line_cells) == 0:
                line_cells = cell.data
            else:
                line_cells = np.vstack([line_cells, cell.data])

    line_data = []
    for key in msh.cell_data_dict["gmsh:physical"].keys():
        if key == "line":
            if len(line_data) == 0:
                line_data = msh.cell_data_dict["gmsh:physical"][key]
            else:
                line_data = np.vstack([line_data, msh.cell_data_dict["gmsh:physical"][key]])
        elif key == "triangle":
            triangle_data = msh.cell_data_dict["gmsh:physical"][key]

    triangle_mesh = meshio.Mesh(points=msh.points[:, :2], cells={"triangle": triangle_cells},
                                cell_data={"name_to_read": [triangle_data]})

    line_mesh = meshio.Mesh(points=msh.points[:, :2], cells=[("line", line_cells)],
                            cell_data={"name_to_read": [line_data]})
    meshio.write("mesh.xdmf", triangle_mesh)
    meshio.write("mf.xdmf", line_mesh)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--res", default=resolution, type=np.float64, dest="res",
                        help="Resolution near circular obstacle")
    args = parser.parse_args()
    create_stokes_mesh(args.res)
