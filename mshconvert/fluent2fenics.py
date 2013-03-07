#!/usr/bin python
import sys
import re
from sets import Set
from copy import copy
from numpy import zeros, array, sign, cross, dot, ones, arctan, sin, cos, pi, mod, sqrt
from numpy.linalg import norm
from pylab import find
from operator import add
from scipy.optimize import fsolve
from dolfin import File, MeshEditor, Mesh, Cell, facets

print "Converting from ANSYS Fluent format (.msh) to FEniCS format"

# Use regular expressions to identify sections and tokens found in a fluent file
re_dimline  = re.compile(r"\(2\s(\d)\)")
re_comment  = re.compile(r"\(0\s.*")
re_zone0    = re.compile(r"\(10\s\(0\s(\w+)\s(\w+)\s(\d+)\s(\d+)\)\)")
re_zone     = re.compile(r"\(10\s\((\w+)\s(\w+)\s(\w+)\s(\d+)\s(\d)\)(\(|)")
re_face0    = re.compile(r"\(13(\s*)\(0\s+(\w+)\s+(\w+)\s+(0|0 0)\)\)")
re_face     = re.compile(r"\(13(\s*)\((\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\)(\s*)(\(|)")
re_periodic = re.compile(r"\(18.*\((\w+)\s+(\w+)\s+(\w+)\s+(\w+)\).*\(")
re_pfaces   = re.compile(r"((^\s)|)(\w+)(\s*)(\w+)")
re_cells0   = re.compile(r"\(12(\s*)\(0(\s+)(\w+)(\s+)(\w+)(\s+)(0|0 0)\)\)")
re_cells    = re.compile(r"\(12.*\((\w+)\s+(\w+)\s+(\w+)\s+(\d+)\s+(\d+)\)\)")
re_cells2   = re.compile(r"\(12(\s*)\((\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\)(\s*)(\(|)")
re_zones    = re.compile(r"\((45|39)\s+\((\d+)\s+(\S+)\s+(\S+).*\)\((.*|[0-9]+[\.]*[0-9]*)\)\)")
re_parant   = re.compile(r"(^\s*\)(\s*)|^\s*\)\)(\s*)|^\s*\(\s*)")

# The fluent mesh (the .msh file) is basically stored as a list of nodes, and then a 
# list of faces for each zone of the mesh, the interior and the boundaries.

# Declare som maps that will be built when reading in the lists of nodes and faces:
cell_map = {}               # Maps cell id with nodes
cell_face_map = {}          # Maps cell id with faces
face_cell_map = {}          # Maps face id with two cells
face_list = []              # List of faces [[id, 2-4 nodes, 2 connecting cells and type]]
face_map = {}               # For each cell a dictionary with key=face and val=local face number
nodes = None                # Will be numpy array of nodes

# Information about connectivity and boundaries
boundary_cells = {}         # List of cells attached to a boundary facet. Key is zone id

# Some global values
num_cells = {}              # Total number of cells in different zones
zones = {}                  # zone information
zone_number_of_faces = {}   # number of faces for each zone

def read_periodic(ifile, periodic_dx):
    """Scan past periodic section"""
    while 1:
        line = ifile.readline()
        a = re.search(re_pfaces, line)
        if a:
            continue
        break

def read_zone_nodes(dim, Nmin, Nmax, ifile):
    """Scan lines for nodes and return in an array."""
    line = ifile.readline()
    readline = False
    if re.search(re_parant, line): # check for initial paranthesis
        readline = True
        #dummy = lines.pop(0)
    global nodes 
    nodes = zeros((dim, Nmax - Nmin + 1))
    for i in range(Nmin, Nmax + 1):
        if readline:
            line = ifile.readline()
        readline = True
        nodes[:, i - Nmin] = [eval(x) for x in line.split()]
    
def read_faces(zone_id, Nmin, Nmax, bc_type, face, ifile):
    """Read all faces and create cell_face_map + some boundary maps."""    
    line = ifile.readline()
    readline = False
    if re.search(re_parant, line): # check for initial paranthesis
        readline = True

    ls = []
    for i in range(Nmin, Nmax + 1):
        if readline:
            line = ifile.readline()
        readline = True
        ln = line.split()
        if face == 0:
            nd = int(ln[0]) # Number of nodes
            nds = [int(x, 16) for x in ln[1:(nd + 1)]]
            cells = [int(x, 16) for x in ln[(nd + 1):]]
        else:
            nd = face
            nds = [int(x, 16) for x in ln[:nd]]
            cells = [int(x, 16) for x in ln[nd:]]
                    
        if min(cells) == 0: # A boundary zone
            if zone_id in boundary_cells:
                boundary_cells[zone_id][max(cells)] = array(nds)
            else:
                boundary_cells[zone_id] = {max(cells): array(nds)}

        for c in cells:
            if c > 0:    
                if not c in cell_map:
                    cell_map[c] = copy(nds)                    
                else:
                    cell_map[c] = list(Set(cell_map[c] + nds))
                                    
def scan_fluent_mesh(ifile):  
    """Scan fluent mesh and generate numerous maps."""
    # Warning! Not yet tested for multiple interior zones
    dim = 0
    one = 0
    num_faces = 0
    while 1:
        line = ifile.readline()
        if len(line) == 0:
            print 'Finished reading file\n'
            break

        if dim == 0: # Dimension usually comes first
            a = re.search(re_dimline, line)
            if a: 
                print 'Reading dimensions\n'
                dim = int(a.group(1))
                print 'Mesh is ' + str(dim) + 'D\n'
                continue
        
        if one == 0: # The total number of nodes
            a = re.search(re_zone0, line)
            if a:
                print 'Reading zone info\n'
                one, num_vertices, dummy1, dummy2 = int(a.group(1)), \
                     int(a.group(2), 16), int(a.group(3), 16), int(a.group(4))
                continue
            
        a = re.search(re_zone, line) # Nodes
        if a:
            zone_id, first_id, last_id, dummy1, dummy2 = int(a.group(1), 16), \
                int(a.group(2), 16), int(a.group(3), 16), int(a.group(4)), \
                int(a.group(5))
            print 'Reading ', last_id - first_id + 1,' nodes in zone ', zone_id + 1, '\n'
            read_zone_nodes(dim, first_id, last_id, ifile)
            continue
            
        a = re.search(re_zones,line) # Zone info
        if a:
            print 'Reading zone ', line
            dummy, zone_id, zone_type, zone_name, radius =  \
                       int(a.group(1)), int(a.group(2)),  a.group(3), \
                       a.group(4), a.group(5)
            zones[zone_id] = [zone_type, zone_name, radius]
            continue
        
        a = re.search(re_cells0, line) # Get total number of cells/elements
        if a:
            print 'Reading cell info ', line
            first_id, tot_num_cells = int(a.group(3),16), int(a.group(5), 16)
            continue

        a = re.search(re_cells,line) # Get the cell info.
        if a:
            zone_id, first_id, last_id, bc_type, element_type = \
                int(a.group(1)), int(a.group(2), 16), int(a.group(3), 16), \
                int(a.group(4), 16), int(a.group(5), 16)
            print 'Reading ', last_id - first_id + 1,' cells in zone ', zone_id, '\n'
            if last_id == 0:
                raise TypeError("Zero elements!")
            num_cells[zone_id] = [first_id, last_id, bc_type, element_type]
            continue

        a = re.search(re_cells2,line) # Get the cell info.
        if a:
            raise TypeError("Wrong cell type. Can only handle one single cell type")

        a = re.search(re_face0, line)
        if a:
            print 'Reading total number of faces\n', line
            num_faces = int(a.group(3),16)
            continue
            
        a = re.search(re_face, line)
        if a:
            print 'Reading faces ', line
            zone_id, first_id, last_id, bc_type, face_type = \
                 int(a.group(2), 16), int(a.group(3), 16), int(a.group(4), 16), \
                 int(a.group(5), 16), int(a.group(6), 16)
            read_faces(zone_id, first_id, last_id, bc_type, face_type, ifile)
            zone_number_of_faces[zone_id] = last_id - first_id + 1
            continue
        
        a = re.search(re_periodic, line)
        if a:
            print 'Reading periodic connectivity\n', line
            read_periodic(ifile, periodic_dx)
            continue
        
        #print 'Line = ',line
        if any([re.search(st, line) for st in (re_parant, re_comment)]) or \
                                                             not line.strip():
            continue
                
        # Should not make it here
        print 'Line = ',line
        raise IOError('Something went wrong reading fluent mesh.')
    

def write_fenics_file(dim, ofilename):
    ofile  = File(ofilename + '.xml')
    mesh = Mesh()
    editor = MeshEditor()
    editor.open(mesh, dim, dim)
    editor.init_vertices(nodes.shape[1])
    editor.init_cells(len(cell_map))
    
    for i in range(nodes.shape[1]):
        if dim == 2:
            editor.add_vertex(i, nodes[0, i], nodes[1, i])
        else:
            editor.add_vertex(i, nodes[0, i], nodes[1, i], nodes[2, i])
            
    for i in range(1, len(cell_map)+1):
        if dim == 2:
            editor.add_cell(i-1, cell_map[i][0]-1, cell_map[i][1]-1, cell_map[i][2]-1)
        else:
            editor.add_cell(i-1, cell_map[i][0]-1, cell_map[i][1]-1, cell_map[i][2]-1, cell_map[i][3]-1)
    
    mesh.order()
    mvc = mesh.domains().markers(dim-1)
    for zone, cells in boundary_cells.iteritems():
        for cell, nds in cells.iteritems():
            dolfin_cell = Cell(mesh, cell-1)
            nodes_of_cell = dolfin_cell.entities(0)
            #print cell
            #print nodes_of_cell
            nodes_of_face = nds - 1
            #print nodes_of_face
            for jj, ff in enumerate(facets(dolfin_cell)):
                facet_nodes = ff.entities(0)
                #print facet_nodes
                if all(map(lambda x: x in nodes_of_face, facet_nodes)):
                    local_index = jj
                    break
            mvc.set_value(cell-1, local_index, zone)
        
    ofile << mesh        
    from dolfin import plot
    plot(mesh, interactive=True)
    print 'Finished writing FEniCS mesh\n'
    
def convert(fluentmesh):
    """Converts a fluent mesh to a mesh format that can be used by FEniCS. 
         
         fluentmesh = fluent mesh (*.msh file)         
    """
    ofilename = fluentmesh[:-4]
    ifile  = open(fluentmesh, "r")
    scan_fluent_mesh(ifile)
    write_fenics_file(nodes.shape[0], ofilename)
    ifile.close()
